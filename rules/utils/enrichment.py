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

    y = row[norm_cols].astype(float).to_numpy().transpose()
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
        counts = row[count_cols].astype(int).to_numpy().transpose()[not_nan_idx]
        reference_counts = row[reference_cols].astype(int).to_numpy().transpose()[not_nan_idx]
        se = np.sqrt( (1/ (reference_counts[0]+0.5) ) + (1/ (reference_counts[1]+0.5) ) + (1/ (counts[0]+0.5) ) + (1/ (counts[1]+0.5) ) )
    else:
        x = sm.add_constant(x)
        w = 1
        if weighted:
            weight_cols = [f"{tp}_weight" for tp in timepoints]
            w = row[weight_cols].astype(float).to_numpy().transpose()
        results = sm.WLS(y, x, weights=w, missing='drop').fit()
        slope = results.params[1]
        se = results.bse[0]

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
    df_chunks = np.array_split(df, n_threads)
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

def normalize_counts(counts_pivot, timepoints_df, entity_columns, reference_entity=''):
    """
    Apply log-normalization relative to reference entity and prepare data for regression.

    This function takes a pivot table of counts (entities × samples) and normalizes each
    sample/replicate combination according to the timepoints mapping.

    Args:
        counts_pivot (pd.DataFrame): DataFrame with entity_id + entity_columns + sample columns
        timepoints_df (pd.DataFrame): Timepoints mapping with sample/replicate info
        entity_columns (list): List of column names that define entity identity (e.g., ['reference_name', 'label_barcodes'])
        reference_entity (str): Specific entity_id to use as reference, or '' for auto-selection, or 'all' for all entities

    Returns:
        pd.DataFrame: DataFrame with rows per (sample, replicate, entity) and columns for:
                     - sample info
                     - entity_id
                     - entity columns (reference_name, label_barcodes, etc.)
                     - counts per timepoint
                     - normalized values per timepoint
                     - weights per timepoint
                     - reference counts per timepoint (for SE calculation)
    """
    # Determine reference entity if not specified
    # Get sample label to identify which columns are metadata vs sample counts
    sample_label = timepoints_df.columns[0]

    # Get all unique sample names from time points
    timepoint_cols = [col for col in timepoints_df.columns if col != sample_label]
    all_sample_names = set()
    for col in timepoint_cols:
        all_sample_names.update(timepoints_df[col].unique())

    # Sample columns are those that appear in the timepoints data (excluding entity_id and entity_columns)
    sample_cols = [col for col in counts_pivot.columns if col in all_sample_names]

    # Use entity_id as index for operations, get just numeric sample columns
    counts_pivot_numeric = counts_pivot.set_index('entity_id')[sample_cols]
    counts_pivot_no_na = counts_pivot_numeric.dropna()

    # Determine the reference entity to use
    if reference_entity != 'all':
        if reference_entity and reference_entity in counts_pivot_no_na.index:
            print(f"[NOTICE] enrichment.py: Provided enrichment reference {reference_entity} found in all samples. Using this entity as the enrichment reference.\n")
        else:
            if reference_entity:
                print(f"[WARNING] enrichment.py: Provided enrichment reference {reference_entity} not found in all samples.\n")

            # Auto-select reference using filtering criteria
            if counts_pivot_no_na.shape[0] > 0:
                # Identify t0 and tLast samples
                timepoint_values = [int(col) for col in timepoint_cols]
                t0 = str(min(timepoint_values))
                tLast = str(max(timepoint_values))

                # Get samples for each timepoint
                t0_samples = timepoints_df[t0].unique()
                tLast_samples = timepoints_df[tLast].unique()

                # Calculate mean counts at t0 and tLast for each entity
                candidate_entities = counts_pivot_no_na.copy()
                candidate_entities['t0_mean'] = candidate_entities[t0_samples].mean(axis=1)
                candidate_entities['tLast_mean'] = candidate_entities[tLast_samples].mean(axis=1)

                # Filter: keep entities with t0_mean >= 25th percentile (exclude bottom 75%)
                candidate_entities = filter_by_quantile(candidate_entities, 't0_mean', 0.25, keep_upper=True)

                # Filter: keep entities with tLast_mean <= 25th percentile (exclude top 75%)
                candidate_entities = filter_by_quantile(candidate_entities, 'tLast_mean', 0.25, keep_upper=False)

                # Select most abundant at tLast from filtered candidates
                if len(candidate_entities) > 0:
                    reference_entity = candidate_entities['tLast_mean'].idxmax()
                    print(f"[NOTICE] enrichment.py: Auto-selected enrichment reference {reference_entity} (well-represented throughout, not among top enriched).\n")
                else:
                    # Fallback 1: Most abundant at tLast among entities present in all timepoints
                    if len(counts_pivot_no_na) > 0:
                        reference_entity = counts_pivot_no_na[tLast_samples].mean(axis=1).idxmax()
                        print(f"[NOTICE] enrichment.py: Filtering criteria too strict. Using most abundant entity at tLast among those present in all timepoints, {reference_entity}, as the enrichment reference.\n")
                    else:
                        # Fallback 2: Use 'all' if no entities present in all timepoints
                        reference_entity = 'all'
            else:
                reference_entity = 'all'

    if reference_entity == 'all':
        print("[NOTICE] enrichment.py: Using all entities as the enrichment reference.\n")

    # Add replicate column to timepoints_df
    replicate_col = []
    sample_count = Counter()
    for _, row in timepoints_df.iterrows():
        sample_count[row[sample_label]] += 1
        replicate_col.append(sample_count[row[sample_label]])
    timepoints_df = timepoints_df.copy()
    timepoints_df['replicate'] = replicate_col
    timepoints_df = timepoints_df.set_index([sample_label, 'replicate'])

    # Loop through each sample/replicate and normalize
    sample_df_list = []
    # After set_index, columns only contain timepoints (sample_label is now in index)
    tp_cols = timepoints_df.columns.to_list()

    cols_dict = {tp_type:[f'{x}_{tp_type}' for x in tp_cols] for tp_type in ['count', 'normalized', 'weight', 'reference']}

    for idx_num, (index, row) in enumerate(timepoints_df.iterrows()):
        # Get the tag_bcGroup names for each timepoint (exclude sample label column)
        slice_cols = row[tp_cols].to_list()

        # Check that all required columns exist
        for col in slice_cols:
            if col not in counts_pivot.columns:
                print(f'[WARNING] enrichment.py: No entities above the threshold were identified for sample `{col}`, setting counts for this sample to 0\n')
                counts_pivot = counts_pivot.assign(**{col:0})

        # Select entity_id + entity_columns + sample counts for this timepoint
        sample_df = counts_pivot.loc[:,['entity_id'] + entity_columns + slice_cols].copy()

        # Remove rows with first timepoint counts of 0
        sample_df = sample_df[sample_df[slice_cols[0]] > 0]

        # Calculate reference counts
        if reference_entity == 'all':
            # Get the total counts for all entities
            reference_counts = sample_df[slice_cols].sum(axis=0).to_numpy()
        elif reference_entity:
            # Get the counts for the specific reference entity
            reference_counts = sample_df[sample_df['entity_id'] == reference_entity][slice_cols].sum(axis=0).to_numpy()
        else:
            # Get the total counts for only entities that appear in all timepoints
            reference_counts = sample_df[(sample_df[slice_cols] > 0).all(axis=1)][slice_cols].sum(axis=0).to_numpy()

        reference_proportion = reference_counts / reference_counts.sum()
        counts = sample_df[slice_cols].to_numpy().astype(float)

        # Set counts following 0 counts to nan then calculate normalized log transformed y values
        counts_pruned = counts.copy()
        mask = (np.cumsum(counts_pruned == 0, axis=1) > 1)
        counts_pruned[mask] = np.nan
        y_normalized = np.log( (counts_pruned+0.5) / (reference_counts+0.5) )

        # Calculate weights as reference_proportion / ( 1/(counts)+0.5 + 1/(reference_counts)+0.5) )
        weights = np.repeat(reference_proportion[np.newaxis,], counts.shape[0], axis=0) / ( ( 1/ (counts + 0.5) ) + np.repeat(( 1/ (reference_counts + 0.5) )[np.newaxis,], counts.shape[0], axis=0) )

        # Add normalized values, weights, reference count to the sample_df
        sample_df[cols_dict['normalized']] = y_normalized
        sample_df[cols_dict['weight']] = weights
        sample_df[cols_dict['reference']] = np.repeat(reference_counts[np.newaxis,], counts.shape[0], axis=0)

        # Add sample and replicate info, rename count columns
        sample_cols = {sample_label:index[0], 'replicate':index[1]}
        rename_dict = {old:new for old, new in zip(slice_cols, cols_dict['count'])}
        sample_df = sample_df.assign(**sample_cols).rename(columns=rename_dict)
        sample_df = sample_df.astype({count_col:'int64' for count_col in cols_dict['count']})
        # Reorder columns: sample info, entity_id, entity_columns, then count/normalized/weight/reference data
        sample_df = sample_df[[sample_label, 'replicate', 'entity_id'] + entity_columns + cols_dict['count'] + cols_dict['normalized'] + cols_dict['weight'] + cols_dict['reference']]
        sample_df_list.append(sample_df)

    enrichment_df = pd.concat(sample_df_list)

    # Bring the reference entity to the top of the dataframe if specified
    if reference_entity != 'all':
        ref_df = enrichment_df[enrichment_df['entity_id'] == reference_entity]
        no_ref_df = enrichment_df[enrichment_df['entity_id'] != reference_entity]
        enrichment_df = pd.concat([ref_df, no_ref_df])

    return enrichment_df


# =============================================================================
# MODE-SPECIFIC DATA LOADING FUNCTIONS
# =============================================================================

def load_demux_data(demux_stats_path, screen_failures=True):
    """
    Load and pivot barcode counts from demux-stats.csv for demux enrichment mode.

    Groups by (reference_name, label_barcodes) to create entity identities.

    Structure:
    - output_file_barcodes: partition group (defines the sample)
    - label_barcodes: entity to track (barcode sequence or group name)
    - Numerator: count of each label_barcodes
    - Denominator: total count within each tag/output_file_barcodes

    Args:
        demux_stats_path (str): Path to demux-stats.csv file
        screen_failures (bool): Whether to exclude reads with 'fail' in label_barcodes column

    Returns:
        tuple: (counts_pivot, entity_columns)
            - counts_pivot: DataFrame with entity_id (for operations) + entity columns (for data) + sample columns
            - entity_columns: list of column names that define entity identity ['reference_name', 'label_barcodes']
    """
    all_demux_df = pd.read_csv(demux_stats_path, index_col=False)

    # Check that required columns exist
    required_cols = ['tag', 'reference_name', 'output_file_barcodes', 'label_barcodes', 'count']
    missing_cols = [col for col in required_cols if col not in all_demux_df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns in demux-stats.csv: {missing_cols}")

    # Remove failures if requested
    # label_barcodes can be 'fail' or contain 'fail-' or '-fail' if multiple barcodes
    if screen_failures:
        n_before = len(all_demux_df)
        all_demux_df = all_demux_df[~all_demux_df['label_barcodes'].str.contains('fail', na=False)].copy()
        n_after = len(all_demux_df)
        n_removed = n_before - n_after
        if n_removed > 0:
            print(f"[NOTICE] enrichment.py: Removed {n_removed} rows with 'fail' in label_barcodes\n")

    # Create tag_bcGroup identifier for samples (partition groups within tags)
    all_demux_df['tag_bcGroup'] = all_demux_df['tag'] + '_' + all_demux_df['output_file_barcodes']

    # Entity identity columns
    entity_columns = ['reference_name', 'label_barcodes']

    # Aggregate by entity columns + sample
    counts_aggregated = all_demux_df.groupby(entity_columns + ['tag_bcGroup'])['count'].sum().reset_index()

    # Pivot to get entities as rows, samples as columns
    counts_pivot = counts_aggregated.pivot(index=entity_columns, columns='tag_bcGroup', values='count')
    counts_pivot = counts_pivot.fillna(0).reset_index()

    # Create entity_id as a helper column for operations (but keep original columns)
    # Format: "ref|label"
    counts_pivot['entity_id'] = counts_pivot['reference_name'] + '|' + counts_pivot['label_barcodes']

    # Reorder columns: entity_id, entity columns, then samples
    sample_cols = [col for col in counts_pivot.columns if col not in entity_columns + ['entity_id']]
    counts_pivot = counts_pivot[['entity_id'] + entity_columns + sample_cols]

    print(f"[NOTICE] enrichment.py: Loaded {len(counts_pivot)} unique entities (reference × label_barcodes combinations)\n")
    print(f"[NOTICE] enrichment.py: Found {len(sample_cols)} samples: {sample_cols[:5]}{'...' if len(sample_cols) > 5 else ''}\n")

    return counts_pivot, entity_columns


def load_genotype_data(timepoints_df, timepoints_path):
    """
    Load and aggregate genotype counts from multiple CSV files for genotype enrichment mode.

    Uses timepoints DataFrame to identify which genotype files correspond to which samples,
    then aggregates counts by genotype (combination of mutation columns).

    Args:
        timepoints_df (pd.DataFrame): Timepoints DataFrame with columns:
            - First column: sample names
            - Other columns: timepoint values, with cells containing 'tag_barcode' identifiers
        timepoints_path (str): Path to timepoints CSV file, used to locate mutation_data directory

    Returns:
        tuple: (counts_pivot, entity_columns)
            - counts_pivot: DataFrame with entity_id (for operations) + entity columns (for data) + sample columns
            - entity_columns: list of column names that define genotype identity:
                ['reference_name', 'NT_substitutions', 'NT_insertions', 'NT_deletions',
                 'AA_substitutions_nonsynonymous', 'AA_substitutions_synonymous']

    Expected directory structure:
        mutation_data/{tag}/{barcode}/{tag}_{barcode}_genotypes.csv
    """
    # Entity columns for genotypes (NT columns always required, AA columns optional)
    nt_entity_columns = ['reference_name', 'NT_substitutions', 'NT_insertions', 'NT_deletions']
    aa_entity_columns = ['AA_substitutions_nonsynonymous', 'AA_substitutions_synonymous']

    # Collect all genotype dataframes
    all_genotype_dfs = []

    # Get base directory (mutation_data should be in same directory as timepoints file)
    base_dir = os.path.dirname(timepoints_path)
    mutation_data_dir = os.path.join(base_dir, 'mutation_data')

    if not os.path.exists(mutation_data_dir):
        raise ValueError(f"mutation_data directory not found at {mutation_data_dir}")

    # Iterate through timepoints to find all genotype files
    sample_col = timepoints_df.columns[0]
    timepoint_cols = timepoints_df.columns[1:]

    # Collect unique tag_barcode identifiers from timepoints
    tag_barcodes = set()
    for _, row in timepoints_df.iterrows():
        for timepoint_col in timepoint_cols:
            tag_barcodes.add(row[timepoint_col])

    # Track which files have AA columns
    files_with_aa = 0
    files_without_aa = 0

    # Parse tag_barcode format: "tag_barcode"
    # Split on underscore to get tag and barcode
    for tag_barcode in tag_barcodes:
        parts = tag_barcode.split('_')
        if len(parts) != 2:
            print(f"[WARNING] enrichment.py: Could not parse tag_barcode '{tag_barcode}', expected format 'tag_barcode'\n")
            continue

        tag, barcode = parts

        # Construct file path: mutation_data/{tag}/{barcode}/{tag}_{barcode}_genotypes.csv
        genotype_file = os.path.join(mutation_data_dir, tag, barcode, f"{tag}_{barcode}_genotypes.csv")

        if os.path.exists(genotype_file):
            # Read genotype file
            geno_df = pd.read_csv(genotype_file)

            # Check that required NT columns exist
            required_cols = nt_entity_columns + ['count']
            missing_cols = [col for col in required_cols if col not in geno_df.columns]
            if missing_cols:
                print(f"[WARNING] enrichment.py: Missing required columns {missing_cols} in {genotype_file}, skipping\n")
                continue

            # Check if AA columns are present
            has_aa = all(col in geno_df.columns for col in aa_entity_columns)
            if has_aa:
                files_with_aa += 1
            else:
                files_without_aa += 1

            # Store dataframe with AA status
            all_genotype_dfs.append((geno_df, tag_barcode, has_aa))
        else:
            print(f"[WARNING] enrichment.py: Genotype file not found: {genotype_file}\n")

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
    for geno_df, tag_barcode, has_aa in all_genotype_dfs:
        # Select entity columns + count
        geno_df = geno_df[entity_columns + ['count']].copy()

        # Add tag_barcode identifier
        geno_df['tag_barcode'] = tag_barcode

        processed_dfs.append(geno_df)

    if not processed_dfs:
        raise ValueError(f"No genotype files found in {mutation_data_dir}")

    # Concatenate all genotype dataframes
    all_genotypes = pd.concat(processed_dfs, ignore_index=True)

    # Replace empty strings with a placeholder TEMPORARILY to avoid NaN issues in groupby/pivot
    # We'll convert back to empty strings after pivoting
    placeholder = '__EMPTY__'
    for col in entity_columns:
        all_genotypes[col] = all_genotypes[col].fillna('').astype(str)
        # Replace empty string with placeholder temporarily
        all_genotypes[col] = all_genotypes[col].replace('', placeholder)

    # Aggregate by entity columns + tag_barcode
    counts_aggregated = all_genotypes.groupby(entity_columns + ['tag_barcode'])['count'].sum().reset_index()

    # Pivot to get entities as rows, samples as columns
    counts_pivot = counts_aggregated.pivot(index=entity_columns, columns='tag_barcode', values='count')
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

    # Reorder columns: entity_id, entity columns, then samples
    sample_cols = [col for col in counts_pivot.columns if col not in entity_columns + ['entity_id']]
    counts_pivot = counts_pivot[['entity_id'] + entity_columns + sample_cols]

    print(f"[NOTICE] enrichment.py: Loaded {len(counts_pivot)} unique genotypes\n")
    print(f"[NOTICE] enrichment.py: Found {len(sample_cols)} samples: {sample_cols[:5]}{'...' if len(sample_cols) > 5 else ''}\n")

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
    # Load timepoints CSV once
    timepoints_df = pd.read_csv(kwargs['timepoints_path']).astype(str)

    # Filter to specific timepoint sample if provided
    if kwargs.get('timepoint_sample'):
        timepoint_sample = kwargs['timepoint_sample']
        if timepoint_sample not in timepoints_df.iloc[:, 0].values:
            raise ValueError(f"Timepoint sample '{timepoint_sample}' not found in timepoints CSV")
        timepoints_df = timepoints_df[timepoints_df.iloc[:, 0] == timepoint_sample].copy()
        print(f"[NOTICE] enrichment.py: Processing single timepoint sample: {timepoint_sample}\n")

    # Timepoints are all columns except the first (sample label)
    timepoints = timepoints_df.columns.to_list()[1:]

    # Load data based on mode
    if mode == 'demux':
        counts_pivot, entity_columns = load_demux_data(
            kwargs['demux_stats'],
            screen_failures=kwargs.get('screen_failures', True)
        )
    elif mode == 'genotype':
        counts_pivot, entity_columns = load_genotype_data(
            timepoints_df,
            kwargs['timepoints_path']
        )
    else:
        raise ValueError(f"Invalid mode: {mode}. Must be 'demux' or 'genotype'")

    # Normalize counts
    enrichment_df = normalize_counts(
        counts_pivot,
        timepoints_df,
        entity_columns,
        reference_entity=kwargs.get('reference_entity', '')
    )

    # Calculate enrichment scores via weighted linear regression
    tp_vals = [int(x) for x in timepoints]
    x_normalized = [float(x)/float(max(tp_vals)) for x in tp_vals]
    enrichment_df = regression_parallelize(kwargs.get('n_threads', 1), enrichment_df, x_normalized, timepoints, weighted=True)

    # Remove reference columns (used for calculation but not needed in output)
    ref_cols = [f"{tp}_reference" for tp in timepoints]
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