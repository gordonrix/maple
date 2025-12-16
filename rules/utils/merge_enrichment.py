"""
part of the maple pipeline, written by Gordon Rix

combines genotypes with enrichment scores (enrichment-centric merge)
- Start with enrichment scores
- Add barcode information via seq_IDs
- Add genotype details via genotypes file
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

    # Step 2: Load genotypes and merge with barcode mapping
    genotypes = pd.read_csv(snakemake.input.genotypes, index_col=False)
    genotypes_with_barcodes = pd.merge(genotype_barcode_map, genotypes, on='genotype_ID', how='left')

    # Step 3: Load and process enrichment scores
    for mean_csv in input_list:
        tag = mean_csv.replace('enrichment/','').replace('_enrichment-scores-mean.csv','')
        mean_enrichment = pd.read_csv(mean_csv, index_col=False)
        sample_label, barcode = list(mean_enrichment.columns)[:2]
        mean_enrichment = mean_enrichment.pivot(index=barcode, columns=sample_label, values='mean_enrichment_score')
        if snakemake.params.filter_missing_replicates:
            mean_enrichment = mean_enrichment.dropna(how='any')
        mean_enrichment.columns = [f'mean_enrichment_score_{tag}_' + str(sample) for sample in mean_enrichment.columns]
        mean_enrichment.reset_index(inplace=True)
        # Barcode column is already 'label_barcodes' from enrichment

        # Step 4: Merge enrichment with genotypes (enrichment-centric: keep all enrichment entities)
        # Merge on both reference_name and label_barcodes to ensure correct matching
        genotypes_with_barcodes = pd.merge(mean_enrichment, genotypes_with_barcodes,
                                          on=['reference_name', 'label_barcodes'], how='left')

    genotypes_with_barcodes.to_csv(snakemake.output.genotypes_enrichment, index=False)