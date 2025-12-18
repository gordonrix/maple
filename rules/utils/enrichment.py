"""
This script calculates enrichment scores from counts of barcodes or genotypes.

It takes as input a CSV file containing the counts of each entity for each sample, and outputs
a CSV file containing the enrichment scores for each sample.

The enrichment score is calculated as the slope of the log-linear regression line fitted to the
counts of each entity across the samples, normalized to the total counts (or reference entity
counts) for each timepoint.

Two modes of operation:
- 'demux': Calculate enrichment for barcodes or references (from demux-stats.csv)
- 'genotype': Calculate enrichment for genotypes (from mutation_data/{tag}/{barcode}/ directories)

See Rubin, ... Fowler. Genome Biol 2017. (https://doi.org/10.1186/s13059-017-1272-5)
for derivations of enrichment score and standard error calculations.
"""

import pandas as pd
import numpy as np
import holoviews as hv
from holoviews.operation import histogram
import hvplot.pandas
hv.extension('bokeh')
import scipy
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import statsmodels.api as sm
from collections import Counter
import argparse
import itertools
import os
from concurrent.futures import ProcessPoolExecutor


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def filter_by_quantile(df, column, quantile, keep_upper=True):
    """Filter DataFrame by quantile threshold on a column."""
    threshold = df[column].quantile(quantile)
    if keep_upper:
        return df[df[column] >= threshold]
    else:
        return df[df[column] <= threshold]


# =============================================================================
# REGRESSION FUNCTIONS (shared by all modes)
# =============================================================================

def apply_regression(row, x, timepoints, weighted=False):
    """
    Apply linear regression to a single entity across timepoints.

    Used as pandas apply function for linear regression using statsmodels.
    x values and weights are taken from each individual row.

    Args:
        row (pandas.Series): row of a pandas DataFrame, must contain columns corresponding
            to the timepoints and, if weighted, the weights for each timepoint
        x (np.array): array of normalized x values (range 0-1)
        timepoints (list): list of non-normalized timepoints, used to find column names
        weighted (bool): whether to use weights for the regression

    Returns:
        tuple(float, float): (slope, standard_error) of the regression line
    """
    norm_cols = [f"{tp}_normalized" for tp in timepoints]

    y = row[norm_cols].astype(float).to_numpy()
    not_nan = ~np.isnan(y)
    not_nan_count = np.count_nonzero(not_nan)
    not_nan_idx = np.where(not_nan)[0]

    if not_nan_count < 2:
        slope, se = np.nan, np.nan
    # if there's only two timepoints, then return slope as normalized log ratio of final / first timepoint counts
    #   if more than 2, use weighted linear regression
    elif not_nan_count == 2:
        count_cols = [f"{tp}_count" for tp in timepoints]
        reference_cols = [f"{tp}_reference" for tp in timepoints]
        y = y[not_nan_idx]
        slope = y[1] - y[0]
        counts = row[count_cols].astype(int).to_numpy()[not_nan_idx]
        reference_counts = row[reference_cols].astype(int).to_numpy()[not_nan_idx]
        se = np.sqrt( (1/ (reference_counts[0]+0.5) ) + (1/ (reference_counts[1]+0.5) ) + (1/ (counts[0]+0.5) ) + (1/ (counts[1]+0.5) ) )
    else:
        x = sm.add_constant(x)
        w = 1
        if weighted:
            weight_cols = [f"{tp}_weight" for tp in timepoints]
            w = row[weight_cols].astype(float).to_numpy()
        results = sm.WLS(y, x, weights=w, missing='drop').fit()
        slope = results.params[1]  # Slope coefficient
        se = results.bse[1]  # SE of slope (not intercept)

    return slope, se


def regression_worker(df_chunk, x, timepoints, weighted=False):
    """
    Worker function for parallelizing the regression apply function.

    Args:
        df_chunk (pd.DataFrame): chunk of dataframe to process
        x (np.array): array of normalized x values (range 0-1)
        timepoints (list): list of non-normalized timepoints
        weighted (bool): whether to use weights for regression

    Returns:
        pd.DataFrame: DataFrame with enrichment_score and standard_error columns
    """
    slopes, ses = zip(*df_chunk.apply(apply_regression, args=[x, timepoints, weighted], axis=1))
    return pd.DataFrame({"enrichment_score": slopes, "standard_error": ses}, index=df_chunk.index)


def regression_parallelize(n_threads, df, x, timepoints, weighted=False):
    """
    Parallelize the regression apply function and add the results to new columns in the dataframe.

    Args:
        n_threads (int): number of threads to use
        df (pandas.DataFrame): dataframe containing the data to be regressed
        x (np.array): array of normalized x values (range 0-1)
        timepoints (list): list of non-normalized timepoints, used to find column names
        weighted (bool): whether to use weights for the regression

    Returns:
        pd.DataFrame: original dataframe with two additional columns for the slope and standard error
                     of the regression for each row
    """
    if len(df) < 2*n_threads:
        n_threads = 1

    # Split DataFrame into chunks manually to avoid FutureWarning from np.array_split
    chunk_size = len(df) // n_threads + (1 if len(df) % n_threads else 0)
    df_chunks = [df.iloc[i:i + chunk_size] for i in range(0, len(df), chunk_size)]
    with ProcessPoolExecutor(max_workers=n_threads) as executor:
        results = executor.map(regression_worker, df_chunks, itertools.repeat(x),
                             itertools.repeat(timepoints), itertools.repeat(weighted))

    # concatenate results and add as new columns to the dataframe
    results = pd.concat(list(results))
    df = pd.concat([df, results], axis=1)

    return df


# =============================================================================
# DATA NORMALIZATION (shared by all modes)
# =============================================================================

def normalize_counts(counts_pivot, entity_columns, reference_entity=''):
    """
    Apply log-normalization relative to reference entity and prepare data for regression.

    Takes counts in wide format with metadata and timepoint columns, normalizes by sample/replicate,
    and prepares data for enrichment score calculation via linear regression.

    Args:
        counts_pivot (pd.DataFrame): Wide-format DataFrame with columns:
            - group, replicate, sample_label: Sample metadata
            - entity_id: Composite entity identifier
            - entity_columns: Columns defining entity identity (reference_name, mutations, etc.)
            - Numeric timepoint columns (e.g., 0, 24, 48): Count values
        entity_columns (list): Column names that define entity identity
        reference_entity (str): Entity to normalize against. Options:
            - 'all': Normalize to total counts across all entities
            - '': Auto-select a stable, abundant entity
            - Specific entity_id: Use this entity as reference

    Returns:
        pd.DataFrame: Long-format with one row per (group, replicate, sample_label, entity) containing:
            - group, replicate, sample_label: Sample metadata
            - entity_id and entity columns
            - {tp}_count columns: Raw counts at each timepoint
            - {tp}_normalized columns: Log-normalized frequencies
            - {tp}_weight columns: Regression weights
            - {tp}_reference columns: Reference entity counts (for SE calculation)
            - enrichment_score, standard_error: Added by regression (not present yet)
    """
    # Identify timepoint columns by numeric dtype
    timepoint_cols = [col for col in counts_pivot.columns if pd.api.types.is_numeric_dtype(counts_pivot[col])]

    # Sort timepoint columns numerically for proper time ordering
    timepoint_cols = sorted(timepoint_cols)

    # Auto-select reference entity if needed
    # Must find an entity that meets criteria across ALL samples
    if reference_entity != 'all' and not reference_entity:
        # For each sample, find entities present at all timepoints
        sample_candidates = []
        for sample_label_val in counts_pivot['sample_label'].unique():
            sample_data = counts_pivot[counts_pivot['sample_label'] == sample_label_val]
            entity_counts = sample_data.groupby('entity_id')[timepoint_cols].sum()
            # Entities with counts > 0 at all timepoints in this sample
            present_all_tp = entity_counts[(entity_counts > 0).all(axis=1)].index.tolist()
            sample_candidates.append(set(present_all_tp))

        # Intersection: entities present at all timepoints in ALL samples
        if sample_candidates:
            common_entities = set.intersection(*sample_candidates)

            if common_entities:
                # Calculate mean counts across all samples for filtering
                common_entity_data = counts_pivot[counts_pivot['entity_id'].isin(common_entities)]
                entity_means = common_entity_data.groupby('entity_id')[timepoint_cols].mean()

                t0 = timepoint_cols[0]
                tLast = timepoint_cols[-1]

                # Apply filtering criteria
                candidates = entity_means.copy()
                candidates = filter_by_quantile(candidates, t0, 0.25, keep_upper=True)
                candidates = filter_by_quantile(candidates, tLast, 0.25, keep_upper=False)

                if len(candidates) > 0:
                    reference_entity = candidates[tLast].idxmax()
                    print(f"[NOTICE] enrichment.py: Auto-selected enrichment reference {reference_entity} (well-represented in all samples, not among top enriched).\n")
                else:
                    # Fallback: most abundant at tLast among common entities
                    reference_entity = entity_means[tLast].idxmax()
                    print(f"[NOTICE] enrichment.py: Filtering criteria too strict. Using most abundant common entity at tLast, {reference_entity}.\n")
            else:
                reference_entity = 'all'
                print("[WARNING] enrichment.py: No entities present at all timepoints in all samples. Using 'all' as reference.\n")
        else:
            reference_entity = 'all'

    if reference_entity == 'all':
        print("[NOTICE] enrichment.py: Using all entities as the enrichment reference (per-sample totals).\n")
    elif reference_entity and reference_entity in counts_pivot['entity_id'].values:
        print(f"[NOTICE] enrichment.py: Using {reference_entity} as the enrichment reference.\n")
    elif reference_entity:
        print(f"[WARNING] enrichment.py: Specified reference {reference_entity} not found. Using 'all' as reference.\n")
        reference_entity = 'all'

    # Process each sample separately
    sample_dfs = []
    for sample_label_val in counts_pivot['sample_label'].unique():
        sample_df = counts_pivot[counts_pivot['sample_label'] == sample_label_val].copy()

        # Extract counts as numpy array (rows=entities, cols=timepoints)
        counts = sample_df[timepoint_cols].to_numpy().astype(float)

        # Calculate reference counts for THIS sample
        if reference_entity == 'all':
            # Sum across all entities in this sample
            reference_counts = counts.sum(axis=0)
        else:
            # Get counts for specific reference entity in this sample
            ref_row = sample_df[sample_df['entity_id'] == reference_entity]
            if len(ref_row) > 0:
                reference_counts = ref_row[timepoint_cols].to_numpy().astype(float)[0]
            else:
                # Fallback to 'all' if reference not found in this sample
                print(f"[WARNING] enrichment.py: Reference {reference_entity} not found in sample {sample_label_val}, using 'all' for this sample.\n")
                reference_counts = counts.sum(axis=0)

        # Calculate normalized values: log((count+0.5) / (reference+0.5))
        # Set counts to nan after first zero (entity dropped out)
        counts_pruned = counts.copy()
        mask = (np.cumsum(counts_pruned == 0, axis=1) > 1)
        counts_pruned[mask] = np.nan
        normalized = np.log((counts_pruned + 0.5) / (reference_counts + 0.5))

        # Calculate weights: reference_proportion / (1/(count+0.5) + 1/(reference+0.5))
        reference_proportion = reference_counts / reference_counts.sum()
        weights = reference_proportion / ((1 / (counts + 0.5)) + (1 / (reference_counts + 0.5)))

        # Add new columns to dataframe
        for i, tp in enumerate(timepoint_cols):
            sample_df[f'{tp}_count'] = counts[:, i].astype(int)
            sample_df[f'{tp}_normalized'] = normalized[:, i]
            sample_df[f'{tp}_weight'] = weights[:, i]
            sample_df[f'{tp}_reference'] = reference_counts[i]

        # Drop original timepoint columns (now have _count versions)
        sample_df = sample_df.drop(columns=timepoint_cols)

        sample_dfs.append(sample_df)

    # Concatenate all samples
    enrichment_df = pd.concat(sample_dfs, ignore_index=True)

    # Reorder columns: metadata, entity info, then timepoint data grouped by type
    # (all counts, all normalized, all weights, all reference)
    entity_columns_list = entity_columns if isinstance(entity_columns, list) else [entity_columns]
    metadata_cols = ['group', 'replicate', 'sample_label', 'entity_id'] + entity_columns_list

    count_cols = [f'{tp}_count' for tp in timepoint_cols]
    normalized_cols = [f'{tp}_normalized' for tp in timepoint_cols]
    weight_cols = [f'{tp}_weight' for tp in timepoint_cols]
    reference_cols = [f'{tp}_reference' for tp in timepoint_cols]

    enrichment_df = enrichment_df[metadata_cols + count_cols + normalized_cols + weight_cols + reference_cols]

    # Move reference entity to top of dataframe if specified (not 'all')
    if reference_entity != 'all' and reference_entity:
        ref_rows = enrichment_df[enrichment_df['entity_id'] == reference_entity]
        other_rows = enrichment_df[enrichment_df['entity_id'] != reference_entity]
        enrichment_df = pd.concat([ref_rows, other_rows], ignore_index=True)

    return enrichment_df


# =============================================================================
# MODE-SPECIFIC DATA LOADING FUNCTIONS
# =============================================================================

def load_demux_data(demux_stats_path, timepoints_df, screen_failures=True):
    """
    Load and pivot barcode counts from demux-stats.csv for demux enrichment mode.

    Reads demux-stats.csv to get entity counts by tag_bcGroup, then uses timepoints DataFrame
    to map tag_bcGroups to (group, replicate, sample_label, timepoint) metadata.

    Args:
        demux_stats_path (str): Path to demux-stats.csv file
        timepoints_df (pd.DataFrame): Must contain exactly three metadata columns (group, replicate,
            sample_label) followed by numeric timepoint columns. Timepoint columns contain tag_bcGroup
            identifiers (tag_barcode format). Each row represents one sample replicate across all timepoints.
        screen_failures (bool): Whether to exclude reads with 'fail' in label_barcodes column

    Returns:
        tuple: (counts_pivot, entity_columns)
            - counts_pivot (pd.DataFrame): Wide-format counts with columns [group, replicate,
              sample_label, entity_id, entity_columns..., timepoint1, timepoint2, ...]. Each row
              represents one (group, replicate, sample_label, entity) combination with counts at
              each timepoint.
            - entity_columns (list): Column names defining entity identity ['reference_name', 'label_barcodes']
    """
    all_demux_df = pd.read_csv(demux_stats_path, index_col=False)

    # Check that required columns exist
    required_cols = ['tag', 'reference_name', 'output_file_barcodes', 'label_barcodes', 'count']
    missing_cols = [col for col in required_cols if col not in all_demux_df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns in demux-stats.csv: {missing_cols}")

    # Remove failures if requested
    if screen_failures:
        n_before = len(all_demux_df)
        all_demux_df = all_demux_df[~all_demux_df['label_barcodes'].str.contains('fail', na=False)].copy()
        n_after = len(all_demux_df)
        n_removed = n_before - n_after
        if n_removed > 0:
            print(f"[NOTICE] enrichment.py: Removed {n_removed} rows with 'fail' in label_barcodes\n")

    # Create tag_bcGroup identifier
    all_demux_df['tag_bcGroup'] = all_demux_df['tag'] + '_' + all_demux_df['output_file_barcodes']

    # Build mapping DataFrame from tag_bcGroup to metadata (allows multiple samples per tag_bcGroup)
    metadata_cols = ['group', 'replicate', 'sample_label']
    timepoint_cols = [col for col in timepoints_df.columns if col not in metadata_cols]

    mapping_rows = []
    for _, row in timepoints_df.iterrows():
        for tp_col in timepoint_cols:
            tag_bcgroup = row[tp_col]
            if pd.notna(tag_bcgroup) and str(tag_bcgroup).strip() != '':
                mapping_rows.append({
                    'tag_bcGroup': tag_bcgroup,
                    'group': row['group'],
                    'replicate': row['replicate'],  # Already str from timepoints_df
                    'sample_label': row['sample_label'],
                    'timepoint': float(tp_col)  # Numeric for easy identification later
                })

    mapping_df = pd.DataFrame(mapping_rows)

    # Entity identity columns
    entity_columns = ['reference_name', 'label_barcodes']

    # Aggregate by entity + tag_bcGroup
    counts_aggregated = all_demux_df.groupby(entity_columns + ['tag_bcGroup'])['count'].sum().reset_index()

    # Merge with metadata mapping (this will duplicate rows for shared tag_bcGroups)
    counts_with_metadata = counts_aggregated.merge(mapping_df, on='tag_bcGroup', how='inner')

    # Aggregate by entity + metadata + timepoint (in case of duplicates)
    groupby_cols = entity_columns + ['group', 'replicate', 'sample_label', 'timepoint']
    counts_aggregated = counts_with_metadata.groupby(groupby_cols)['count'].sum().reset_index()

    # Pivot to get (group, replicate, sample_label, entity) × timepoints
    index_cols = entity_columns + ['group', 'replicate', 'sample_label']
    counts_pivot = counts_aggregated.pivot(index=index_cols, columns='timepoint', values='count')
    counts_pivot = counts_pivot.fillna(0).reset_index()

    # Create entity_id
    counts_pivot['entity_id'] = counts_pivot['reference_name'] + '|' + counts_pivot['label_barcodes']

    # Reorder columns: metadata, entity_id, entity columns, timepoint columns
    timepoint_value_cols = [col for col in counts_pivot.columns
                            if col not in ['group', 'replicate', 'sample_label', 'entity_id'] + entity_columns]
    counts_pivot = counts_pivot[['group', 'replicate', 'sample_label', 'entity_id'] + entity_columns + timepoint_value_cols]

    print(f"[NOTICE] enrichment.py: Loaded {len(counts_pivot)} unique (group, replicate, sample_label, entity) combinations\n")
    print(f"[NOTICE] enrichment.py: Timepoints: {timepoint_value_cols}\n")

    return counts_pivot, entity_columns


def load_genotype_data(timepoints_df):
    """
    Load and aggregate genotype counts from multiple CSV files for genotype enrichment mode.

    Reads the timepoints DataFrame to identify which genotype files belong to which samples
    and timepoints, then aggregates counts by genotype identity (combination of mutation columns).

    Args:
        timepoints_df (pd.DataFrame): Must contain exactly three metadata columns (group, replicate,
            sample_label) followed by numeric timepoint columns. Timepoint columns contain file paths
            to genotype CSV files. Each row represents one sample replicate across all timepoints.

    Returns:
        tuple: (counts_pivot, entity_columns)
            - counts_pivot (pd.DataFrame): Wide-format counts with columns [group, replicate,
              sample_label, entity_id, entity_columns..., timepoint1, timepoint2, ...]. Each row
              represents one (group, replicate, sample_label, entity) combination with counts at
              each timepoint.
            - entity_columns (list): Column names defining genotype identity (reference_name,
              NT_substitutions, NT_insertions, NT_deletions, and optionally AA columns if present
              in all input files).
    """
    # Entity columns for genotypes (NT columns always required, AA columns optional)
    nt_entity_columns = ['reference_name', 'NT_substitutions', 'NT_insertions', 'NT_deletions']
    aa_entity_columns = ['AA_substitutions_nonsynonymous', 'AA_substitutions_synonymous']

    # Collect all genotype dataframes
    all_genotype_dfs = []

    # Identify metadata vs timepoint columns
    metadata_cols = ['group', 'replicate', 'sample_label']
    timepoint_cols = [col for col in timepoints_df.columns if col not in metadata_cols]

    # Track which files have AA columns
    files_with_aa = 0
    files_without_aa = 0

    # Iterate through timepoints DataFrame to collect all files
    for _, row in timepoints_df.iterrows():
        group = row['group']
        replicate = row['replicate']
        sample_label = row['sample_label']

        for timepoint_col in timepoint_cols:
            file_path = row[timepoint_col]

            # Skip empty cells
            if pd.isna(file_path) or str(file_path).strip() == '':
                continue

            # Read genotype file
            if os.path.exists(file_path):
                geno_df = pd.read_csv(file_path)

                # Check that required NT columns exist
                required_cols = nt_entity_columns + ['count']
                missing_cols = [col for col in required_cols if col not in geno_df.columns]
                if missing_cols:
                    print(f"[WARNING] enrichment.py: Missing required columns {missing_cols} in {file_path}, skipping\n")
                    continue

                # Check if AA columns are present
                has_aa = all(col in geno_df.columns for col in aa_entity_columns)
                if has_aa:
                    files_with_aa += 1
                else:
                    files_without_aa += 1

                # Add metadata columns to genotype dataframe
                geno_df = geno_df.copy()
                geno_df['group'] = group
                geno_df['replicate'] = replicate  # Already str from timepoints_df
                geno_df['sample_label'] = sample_label
                geno_df['timepoint'] = float(timepoint_col)  # Numeric for easy identification later

                # Store dataframe with metadata
                all_genotype_dfs.append((geno_df, has_aa))
            else:
                print(f"[WARNING] enrichment.py: Genotype file not found: {file_path}\n")

    if not all_genotype_dfs:
        raise ValueError(f"No genotype files could be loaded from timepoints DataFrame")

    # Determine whether to include AA columns (only if ALL files have them)
    if files_without_aa > 0:
        if files_with_aa > 0:
            print(f"[NOTICE] enrichment.py: {files_without_aa} file(s) missing AA columns. Excluding AA columns from analysis.\n")
        entity_columns = nt_entity_columns
    else:
        print(f"[NOTICE] enrichment.py: All genotype files contain AA substitution columns. Including in analysis.\n")
        entity_columns = nt_entity_columns + aa_entity_columns

    # Process genotype dataframes with final column set
    processed_dfs = []
    for geno_df, has_aa in all_genotype_dfs:
        # Select entity columns + metadata + count
        cols_to_keep = entity_columns + ['group', 'replicate', 'sample_label', 'timepoint', 'count']
        geno_df = geno_df[cols_to_keep].copy()
        processed_dfs.append(geno_df)

    # Concatenate all genotype dataframes
    all_genotypes = pd.concat(processed_dfs, ignore_index=True)

    # Replace empty strings with a placeholder TEMPORARILY to avoid NaN issues in groupby/pivot
    # We'll convert back to empty strings after pivoting
    placeholder = '__EMPTY__'
    for col in entity_columns:
        all_genotypes[col] = all_genotypes[col].fillna('').astype(str)
        # Replace empty string with placeholder temporarily
        all_genotypes[col] = all_genotypes[col].replace('', placeholder)

    # Aggregate by entity columns + metadata + timepoint
    groupby_cols = entity_columns + ['group', 'replicate', 'sample_label', 'timepoint']
    counts_aggregated = all_genotypes.groupby(groupby_cols)['count'].sum().reset_index()

    # Pivot to get entities + metadata as rows, timepoints as columns
    index_cols = entity_columns + ['group', 'replicate', 'sample_label']
    counts_pivot = counts_aggregated.pivot(index=index_cols, columns='timepoint', values='count')
    counts_pivot = counts_pivot.fillna(0).reset_index()

    # Convert placeholder back to empty strings in entity columns
    for col in entity_columns:
        counts_pivot[col] = counts_pivot[col].replace(placeholder, '')

    # Create entity_id as a helper column for operations
    # Format: "ref|NT_subs|NT_ins|NT_dels|AA_nonsyn|AA_syn"
    # Empty values will just appear as empty between pipes (e.g., "ref|A1G|||M1V|")
    id_components = []
    for col in entity_columns:
        id_components.append(counts_pivot[col].fillna('').astype(str))
    counts_pivot['entity_id'] = id_components[0]
    for component in id_components[1:]:
        counts_pivot['entity_id'] = counts_pivot['entity_id'] + '|' + component

    # Reorder columns: metadata, entity_id, entity columns, then timepoint columns
    timepoint_value_cols = [col for col in counts_pivot.columns
                            if col not in ['group', 'replicate', 'sample_label', 'entity_id'] + entity_columns]
    counts_pivot = counts_pivot[['group', 'replicate', 'sample_label', 'entity_id'] + entity_columns + timepoint_value_cols]

    print(f"[NOTICE] enrichment.py: Loaded {len(counts_pivot)} unique (group, replicate, sample_label, genotype) combinations\n")
    print(f"[NOTICE] enrichment.py: Timepoints: {timepoint_value_cols}\n")

    return counts_pivot, entity_columns


# =============================================================================
# MAIN ENRICHMENT CALCULATION
# =============================================================================

def calculate_enrichment(mode='demux', output_csv='', **kwargs):
    """
    Main entry point for enrichment calculation with mode-based dispatch.

    Args:
        mode (str): 'demux' or 'genotype'
        output_csv (str): Path to save enrichment scores CSV
        **kwargs: Mode-specific and common parameters

    Mode-specific parameters:
        For 'demux' mode (barcode or reference enrichment from demux-stats.csv):
            - demux_stats (str): path to demux-stats.csv
            - screen_failures (bool): exclude reads with 'fail' in label_barcodes (default: True)

            Note: Handles both barcode enrichment (label_barcodes = barcode sequences)
                  and reference-only enrichment (label_barcodes = 'none')

        For 'genotype' mode (genotype enrichment from mutation_data directory):
            - timepoints_path (str): path to timepoints CSV (used to identify genotype files)
            - timepoint_sample (str, optional): specific timepoint sample to process

            Note: Genotype files are auto-discovered from mutation_data/{tag}/{barcode}/ directories.
                  Entity columns include NT mutations and optionally AA mutations (if present in all files)

    Common parameters:
        - timepoints_path (str): path to timepoints CSV file
        - reference_entity (str): specific entity to use as reference ('', 'all', or specific entity name)
        - n_threads (int): number of parallel workers (default: 1)

    Returns:
        pd.DataFrame: enrichment scores with columns:
            - sample info (sample, replicate)
            - entity identifier (reference_name|label_barcodes for demux mode, genotype for genotype mode)
            - counts per timepoint
            - normalized values per timepoint
            - weights per timepoint
            - enrichment_score
            - standard_error
    """
    # Load timepoints CSV
    # Keep replicate/group/sample_label as str so only timepoint columns are numeric (for easy identification)
    timepoints_df = pd.read_csv(kwargs['timepoints_path'],
                                dtype={'replicate': str, 'group': str, 'sample_label': str})

    # Filter to specific timepoint sample if provided
    if kwargs.get('timepoint_sample'):
        timepoint_sample = kwargs['timepoint_sample']
        if timepoint_sample not in timepoints_df['sample_label'].values:
            raise ValueError(f"Timepoint sample '{timepoint_sample}' not found in timepoints CSV")
        timepoints_df = timepoints_df[timepoints_df['sample_label'] == timepoint_sample].copy()
        print(f"[NOTICE] enrichment.py: Processing single timepoint sample: {timepoint_sample}\n")

    # Validate sample_label uniqueness
    if not timepoints_df['sample_label'].is_unique:
        duplicate_labels = timepoints_df[timepoints_df['sample_label'].duplicated()]['sample_label'].tolist()
        raise ValueError(f"sample_label values must be unique. Found duplicates: {duplicate_labels}")

    # Timepoints are all columns after group, replicate, sample_label (starting at index 3)
    timepoints = timepoints_df.columns.to_list()[3:]

    # Load data based on mode
    if mode == 'demux':
        counts_pivot, entity_columns = load_demux_data(
            kwargs['demux_stats'],
            timepoints_df,
            screen_failures=kwargs.get('screen_failures', True)
        )
    elif mode == 'genotype':
        counts_pivot, entity_columns = load_genotype_data(timepoints_df)
    else:
        raise ValueError(f"Invalid mode: {mode}. Must be 'demux' or 'genotype'")

    # Normalize counts
    enrichment_df = normalize_counts(
        counts_pivot,
        entity_columns,
        reference_entity=kwargs.get('reference_entity', '')
    )

    # Check if there are any entities to analyze
    if len(enrichment_df) == 0:
        print(f"[WARNING] enrichment.py: No entities remain after filtering (no entities with counts > 0 at first timepoint). Cannot calculate enrichment scores.\n")
        # Return empty dataframe with expected columns
        expected_cols = list(enrichment_df.columns) + ['enrichment_score', 'standard_error']
        empty_df = pd.DataFrame(columns=expected_cols)
        if output_csv:
            empty_df.to_csv(output_csv, index=False)
        return empty_df

    # Calculate enrichment scores via weighted linear regression
    # timepoints are already floats from column names
    tp_vals = [float(tp) for tp in timepoints]
    x_normalized = [x / max(tp_vals) for x in tp_vals]
    enrichment_df = regression_parallelize(kwargs.get('n_threads', 1), enrichment_df, x_normalized, tp_vals, weighted=True)

    # Filter out entities with ≤1 non-zero count across timepoints or NA enrichment score
    count_cols = [f'{tp}_count' for tp in tp_vals]
    initial_count = len(enrichment_df)

    # Count non-zero values across timepoint count columns for each row
    non_zero_counts = (enrichment_df[count_cols] > 0).sum(axis=1)

    enrichment_df = enrichment_df[
        (non_zero_counts > 1) &
        (enrichment_df['enrichment_score'].notna())
    ]

    filtered_count = initial_count - len(enrichment_df)
    if filtered_count > 0:
        print(f"[NOTICE] enrichment.py: Filtered out {filtered_count} entities with ≤1 non-zero count across timepoints or NA enrichment score.\n", file=sys.stderr)

    # Remove reference columns (used for calculation but not needed in output)
    ref_cols = [f"{tp}_reference" for tp in tp_vals]
    enrichment_df = enrichment_df.drop(columns=ref_cols)

    # Round very small values to zero (avoid floating point precision errors like 7e-16)
    enrichment_df['enrichment_score'] = enrichment_df['enrichment_score'].apply(lambda x: 0.0 if abs(x) < 1e-10 else x)
    enrichment_df['standard_error'] = enrichment_df['standard_error'].apply(lambda x: 0.0 if abs(x) < 1e-10 else x)

    # Save output
    if output_csv:
        enrichment_df.to_csv(output_csv, index=False)

    return enrichment_df


# =============================================================================
# FILTERING AND PLOTTING UTILITIES (not directly used by this script)
# =============================================================================

def enrichment_mean_filter(enrichment_df, SE_filter=0, t0_filter=0, score_filter=False,
                           filtered_csv='', mean_csv='', include_mean_normalized=False):
    """
    Filter enrichment scores by standard error then pivot and calculate average across replicates.

    Args:
        enrichment_df (pd.DataFrame): DataFrame with enrichment scores for each sample
        SE_filter (float): proportion of standard error to filter out.
            rows with standard error in the top SE_filter proportion will be filtered out
        t0_filter (float): proportion of t0 counts to filter out.
            rows with t0 counts in the bottom t0_filter proportion will be filtered out
        score_filter (float or False): minimum enrichment score to keep (or False for no filter)
        filtered_csv (str): path to save filtered enrichment scores
        include_mean_normalized (bool): if True, also average normalized frequency columns across replicates (default: False)
        mean_csv (str): path to save mean enrichment scores

    Returns:
        tuple(pd.DataFrame, pd.DataFrame): (filtered enrichment_df, mean_enrichment_df)
    """
    filter_list = []
    if 0 < SE_filter < 1:  # filter by SE
        filter_list.append(('standard_error', SE_filter, False))
    if 0 < t0_filter < 1:  # filter by t0
        # Find first count column (ends with '_count')
        first_count = [col for col in enrichment_df.columns if col.endswith('_count')][0]
        filter_list.append((first_count, t0_filter, True))

    # Get column names: sample_label, replicate, entity_id are first 3
    sample_label = enrichment_df.columns[0]
    enrichment_entity = 'entity_id'
    if filter_list:
        # Calculate per-group quantiles and create filter mask
        mask = pd.Series(True, index=enrichment_df.index)
        for column, quantile, take_upper in filter_list:
            if not take_upper:
                quantile = 1 - quantile
            # Calculate threshold for each group
            thresholds = enrichment_df.groupby([sample_label, 'replicate'])[column].transform(lambda x: x.quantile(quantile))
            # Apply filter
            if take_upper:
                mask &= (enrichment_df[column] >= thresholds)
            else:
                mask &= (enrichment_df[column] <= thresholds)
        enrichment_df = enrichment_df[mask].reset_index(drop=True)

    if score_filter:  # filter by enrichment score
        enrichment_df = enrichment_df[enrichment_df['enrichment_score'] >= score_filter]

    # Pivot enrichment scores to get replicates as columns and calculate mean and SD
    mean_enrichment_df = enrichment_df.pivot(index=[sample_label, enrichment_entity],
                                            columns='replicate', values='enrichment_score')
    mean_enrichment_df.columns = [f'replicate_{rep}' for rep in mean_enrichment_df.columns]
    count_col = mean_enrichment_df.count(axis=1)
    mean_enrichment_df['mean_enrichment_score'] = mean_enrichment_df.mean(axis=1)
    mean_enrichment_df['std_enrichment_score'] = mean_enrichment_df.std(axis=1, ddof=1)  # Sample std (N-1)
    mean_enrichment_df['valid_replicates'] = count_col
    mean_enrichment_df = mean_enrichment_df.reset_index()

    # Optionally average normalized frequency columns (e.g., 0_normalized, 24_normalized, etc.)
    if include_mean_normalized:
        normalized_cols = [col for col in enrichment_df.columns if col.endswith('_normalized')]
        for norm_col in normalized_cols:
            pivot_norm = enrichment_df.pivot(index=[sample_label, enrichment_entity],
                                            columns='replicate', values=norm_col)
            mean_norm = pivot_norm.mean(axis=1)
            mean_enrichment_df = mean_enrichment_df.set_index([sample_label, enrichment_entity])
            mean_enrichment_df[f'mean_{norm_col}'] = mean_norm
            mean_enrichment_df = mean_enrichment_df.reset_index()

    if filtered_csv:
        enrichment_df.to_csv(filtered_csv, index=False)
    if mean_csv:
        mean_enrichment_df.to_csv(mean_csv, index=False)

    return enrichment_df, mean_enrichment_df


def plot_enrichment(enrichment_df, plots_out):
    """
    Plot enrichment scores for each sample and save to file.

    Args:
        enrichment_df (pd.DataFrame): DataFrame with enrichment scores for each sample
        plots_out (str): path to save plots
    """
    sample_label, _, enrichment_entity = enrichment_df.columns[:3]
    samples = list(enrichment_df[sample_label].unique())
    samples.sort()
    plots = {}

    # Use same min/max values for fit lines and label location
    min_value = enrichment_df['enrichment_score'].quantile(0.01)
    max_value = enrichment_df['enrichment_score'].quantile(0.99)

    for sample in samples:
        sample_df = enrichment_df[enrichment_df[sample_label] == sample]
        replicates = list(sample_df['replicate'].unique())
        replicates.sort()
        sample_plots = []

        # If there are no replicates, just plot the distribution of scores
        if len(replicates) == 1:
            sample_plots.append(sample_df.hvplot.hist(y='enrichment_score', ylabel='sequences count',
                                                     bins=20, title=f'{sample}', color='grey')
                              .opts(fontsize={'title':16,'labels':14,'xticks':10,'yticks':10},
                                   width=500, height=500))

        else:
            for replicate1, replicate2 in itertools.combinations(replicates, 2):
                # Filter for the two replicates, pivot to get both replicates in the same row
                pivot_df = sample_df[sample_df['replicate'].isin([replicate1, replicate2])].pivot(
                    index=enrichment_entity, columns='replicate', values='enrichment_score').dropna()

                # Plot points
                kdims = ["replicate " + str(replicate1), "replicate " + str(replicate2)]
                points = hv.Points((pivot_df[replicate1], pivot_df[replicate2]), kdims=kdims)
                points.opts(title=f'{sample}, Replicates {replicate1} and {replicate2}',
                          width=400, height=400,
                          fontsize={'title':16,'labels':14,'xticks':10,'yticks':10},
                          size=2, color='black', alpha=0.1)

                # Use statsmodels OLS to calculate R² between the two replicates
                y = pivot_df[replicate2].to_numpy()
                x = sm.add_constant(pivot_df[replicate1].to_numpy())
                fit = sm.OLS(y,x).fit()
                fit_r2 = fit.rsquared
                intercept = fit.params[0]
                slope = fit.params[1]

                # Use hv.Curve to show the fit line
                trendline = hv.Curve(([min_value, max_value],
                                     [min_value * slope + intercept, max_value * slope + intercept])).opts(
                    xlabel=kdims[0], ylabel=kdims[1], color="lightgrey", line_dash='dashed')

                # Also calculate the R² for fit to x=y line
                r2 = r2_score(pivot_df[replicate1], pivot_df[replicate2])
                r2_label = hv.Text(min_value, max_value,
                                  f' slope 1, R² = {r2:.2f}\n slope {slope:.2f}, R² = {fit_r2:.2f}',
                                  halign='left', valign='top').opts(color='black')

                # Distributions
                x_hist, y_hist = (histogram(points, dimension=k).opts(
                    fontsize={'title':16,'labels':14,'xticks':10,'yticks':10}) for k in kdims)

                sample_plots.append((trendline * points * r2_label) <<
                                  y_hist.opts(width=125, color='grey', tools=['hover']) <<
                                  x_hist.opts(height=125, color='grey', tools=['hover']))

        plots[sample] = hv.Layout(sample_plots).cols(1)

    final_layout = hv.Layout(plots.values()).cols(len(samples))
    hv.save(final_layout, plots_out, backend='bokeh')


# =============================================================================
# MAIN ENTRY POINT (Snakemake or standalone)
# =============================================================================

if __name__ == '__main__':

    # Check if running from Snakemake
    if 'snakemake' in globals():
        # Running from Snakemake
        mode = snakemake.params.get('mode', 'demux')

        args_dict = {
            'mode': mode,
            'timepoints_path': snakemake.input.timepoints,
            'reference_entity': snakemake.params.get('reference_entity', ''),
            'n_threads': snakemake.threads,
            'output_csv': snakemake.output.scores
        }

        # Mode-specific parameters
        if mode == 'demux':
            args_dict['demux_stats'] = snakemake.input.CSV
            args_dict['screen_failures'] = snakemake.params.get('screen_failures', True)
        elif mode == 'genotype':
            args_dict['timepoint_sample'] = snakemake.params.timepoint_sample
    else:
        # Running as standalone script
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument('--demux_stats',
                          help='Path to demux-stats.csv file (required for demux mode)')
        parser.add_argument('--timepoints', required=True,
                          help='Path to timepoints.csv file')
        parser.add_argument('--output', required=True,
                          help='Output CSV file for enrichment scores')
        parser.add_argument('--reference_entity', default='',
                          help='Reference entity to use (empty = auto-select most abundant)')
        parser.add_argument('--screen_failures', action='store_true', default=True,
                          help='Exclude reads with "fail" in label_barcodes (default: True)')
        parser.add_argument('--no_screen_failures', dest='screen_failures', action='store_false',
                          help='Include reads with "fail" in label_barcodes')
        parser.add_argument('--threads', type=int, default=1,
                          help='Number of threads for parallel processing')
        parser.add_argument('--mode', default='demux', choices=['demux', 'genotype'],
                          help='Enrichment mode (default: demux)')
        parser.add_argument('--timepoint_sample',
                          help='Specific timepoint sample to process (for genotype mode)')
        args = parser.parse_args()

        # Validate arguments
        if args.mode == 'demux' and not args.demux_stats:
            parser.error('--demux_stats is required for demux mode')

        args_dict = {
            'mode': args.mode,
            'demux_stats': args.demux_stats,
            'timepoints_path': args.timepoints,
            'reference_entity': args.reference_entity,
            'screen_failures': args.screen_failures,
            'n_threads': args.threads,
            'output_csv': args.output
        }

        # Add optional genotype-specific parameters
        if args.timepoint_sample:
            args_dict['timepoint_sample'] = args.timepoint_sample

    # Run enrichment calculation
    enrichment_df = calculate_enrichment(**args_dict)