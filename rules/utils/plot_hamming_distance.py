#
#  DESCRIPTION   : Script for maple pipeline. Generates bar plots for hamming distance
#                       distributions
#
#  AUTHOR(S)     : Gordon Rix
#

import pandas as pd
import numpy as np
import holoviews as hv
import hvplot.pandas
hv.extension('bokeh')

def main():
    plots = [plot_HD_dist(csv) for csv in snakemake.input]
    plots = hv.Layout(plots).cols(1)
    hvplot.save(plots, snakemake.output.hamDistPlot)

def plot_HD_dist(CSV):
    tag, barcodes, NTorAA = CSV.split('/')[-1].split('_')
    NTorAA = NTorAA.strip('-hamming-distance-distribution.csv')
    HDdistDF = pd.read_csv(CSV)
    totalPairs = int(HDdistDF['number_of_sequence_pairs_with_n_hamming_distance'].sum())

    hvDF = hv.Dataset(HDdistDF, kdims=['n'], vdims=['proportion_of_sequence_pairs_with_n_hamming_distance','number_of_sequence_pairs_with_n_hamming_distance'])
    return hv.Histogram(hvDF).opts(hv.opts.Histogram(tools=['hover'], color='grey', width=800, height=600,
                                    title=f"hamming distance distribution: {tag}, {barcodes}",
                                    ylabel=f'proportion of sequence pairs with n hamming distance\n(total pairs: {totalPairs})',
                                    fontsize={'title':16,'labels':14,'xticks':10,'yticks':10}))

if __name__ == '__main__':
    main()