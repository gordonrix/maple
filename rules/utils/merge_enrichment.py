"""
part of the maple pipeline, written by Gordon Rix

combines genotypes with enrichment scores (genotype-centric merge)
- Start with genotypes (multiple sample files concatenated and deduplicated)
- Add barcode information via seq_IDs
- Merge enrichment scores, keeping all genotypes even without scores
"""

import pandas as pd

if __name__=='__main__':
    if type(snakemake.input.enrichment) == str:
        input_list = [snakemake.input.enrichment]
    else:
        input_list = snakemake.input.enrichment

    # Step 1: Load seq_IDs to get the mapping between genotype_ID, reference_name, and label_barcodes
    seq_ids = pd.read_csv(snakemake.input.seq_ids, index_col=False)
    # Keep only genotype_ID, reference_name, and label_barcodes, drop duplicates
    genotype_barcode_map = seq_ids[['genotype_ID', 'reference_name', 'label_barcodes']].drop_duplicates()

    # Step 2: Load sample genotypes, concatenate, and deduplicate
    sample_files = snakemake.input.sample_genotypes
    if isinstance(sample_files, str):
        sample_files = [sample_files]
    genotypes_list = [pd.read_csv(f, index_col=False) for f in sample_files]
    genotypes = pd.concat(genotypes_list, ignore_index=True)

    # Deduplicate on genotype columns
    dedupe_cols = snakemake.params.dedupe_columns
    genotypes = genotypes.drop_duplicates(subset=dedupe_cols)

    # Record original dtypes to restore after merge (prevents int->float conversion from NaN)
    original_dtypes = genotypes.dtypes.to_dict()

    genotypes_with_barcodes = pd.merge(genotype_barcode_map, genotypes, on=['genotype_ID', 'reference_name'], how='left')

    # Step 3: Load and process enrichment scores
    for mean_csv in input_list:
        tag = mean_csv.replace('enrichment/','').replace('_enrichment-scores-mean.csv','')
        mean_enrichment = pd.read_csv(mean_csv, index_col=False)
        sample_label, entity_id = list(mean_enrichment.columns)[:2]
        mean_enrichment = mean_enrichment.pivot(index=entity_id, columns=sample_label, values='mean_enrichment_score')
        if snakemake.params.filter_missing_replicates:
            mean_enrichment = mean_enrichment.dropna(how='any')
        mean_enrichment.columns = [f'mean_enrichment_score_{tag}_' + str(sample) for sample in mean_enrichment.columns]
        mean_enrichment.reset_index(inplace=True)

        # Parse entity_id (format: 'reference|barcode') into reference_name and label_barcodes
        mean_enrichment['reference_name'] = mean_enrichment[entity_id].str.split('|').str[0]
        mean_enrichment['label_barcodes'] = mean_enrichment[entity_id].str.split('|').str[1]
        mean_enrichment = mean_enrichment.drop(columns=[entity_id])

        # Step 4: Merge enrichment with genotypes (genotype-centric: keep all genotypes)
        # Merge on both reference_name and label_barcodes to ensure correct matching
        genotypes_with_barcodes = pd.merge(genotypes_with_barcodes, mean_enrichment,
                                          on=['reference_name', 'label_barcodes'], how='left')

    # Restore original integer dtypes for columns that were converted to float due to NaN
    for col, dtype in original_dtypes.items():
        if col in genotypes_with_barcodes.columns:
            if pd.api.types.is_integer_dtype(dtype):
                # Use nullable integer type to preserve NaN while keeping int dtype
                genotypes_with_barcodes[col] = genotypes_with_barcodes[col].astype('Int64')

    genotypes_with_barcodes.to_csv(snakemake.output.genotypes_enrichment, index=False)