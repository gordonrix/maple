"""
part of the maple pipeline, written by Gordon Rix

combines enrichment scores with a genotypes dataframe
"""

import pandas as pd

if __name__=='__main__':
    if type(snakemake.input.enrichment) == str:
        input_list = [snakemake.input.enrichment]
    else:
        input_list = snakemake.input.enrichment
        
    genotypes = pd.read_csv(snakemake.input.genotypes, index_col=False)
    # pivot mean enrichment scores from each tag and relabel them so that they can be merged with the genotypes dataframe sequentially
    for mean_csv in input_list:
        tag = mean_csv.replace('enrichment/','').replace('_enrichment-scores-mean.csv','')
        mean_enrichment = pd.read_csv(mean_csv, index_col=False)
        sample_label, barcode = list(mean_enrichment.columns)[:2]
        mean_enrichment = mean_enrichment.pivot(index=barcode, columns=sample_label, values='mean_enrichment_score')
        if snakemake.params.filter_missing_replicates:
            mean_enrichment = mean_enrichment.dropna(how='any')
        mean_enrichment.columns = [f'mean_enrichment_score_{tag}_' + str(sample) for sample in mean_enrichment.columns]
        mean_enrichment.reset_index(inplace=True)
        # rename barcode column to match genotypes barcode column
        mean_enrichment.rename(columns={barcode: 'barcode(s)'}, inplace=True)

        # merge genotypes and mean_enrichment
        genotypes = pd.merge(genotypes, mean_enrichment, on='barcode(s)', how='left')
    genotypes.to_csv(snakemake.output.genotypes_enrichment, index=False)