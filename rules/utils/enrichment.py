"""
This script calculates enrichment scores from counts of barcodes. It takes as input a CSV file containing the counts of each barcode for each sample, and outputs a CSV file containing the enrichment scores for each sample.

The enrichment score is calculated as the slope of the log-linear regression line fitted to the counts of each barcode across the samples,
normalized to the total barcode counts for each timepoint. The script uses the `scikit-learn` library to perform the linear regression.

The input CSV file should have the following format:
XXX

The output CSV file will have one row per barcode, 2 columns for each timepoint (total count and normalized count), and one column for the enrichment score.

Usage: python enrichment.py input_file.csv output_file.csv
"""

import pandas as pd
import numpy as np
import holoviews as hv
from holoviews.operation import histogram
hv.extension('bokeh')
import scipy
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from collections import Counter
import argparse
import itertools

def calculate_enrichment(demux_stats, timepoints, barcode_info, barcode_groups, scores_csv='', mean_scores_csv='', screen_no_group=True):
    """
    Calculate enrichment scores from counts of barcodes.

    Args:
        demux_stats (str): Path to the demux_stats.csv file output by `demux.py`.
        timepoints (str): Path to the timepoints.csv file output by `demux.py`.
        barcode_info (dict): Dictionary of barcode information from config file
        barcode_groups (dict): Dictionary of barcode groups from config file
        output_csv (str): Path to the output CSV file.
        screen_no_group (bool): Whether to screen out barcodes that were not assigned to a barcode group. Default: True.

    Returns:
        enrichment_df (pandas.DataFrame): DataFrame containing the enrichment scores for each barcode/sample combination.
    """

    demux_stats_DF = pd.read_csv(demux_stats,index_col=False)

    barcode_types = list(barcode_info.keys())
    grouped_barcode_types = list(barcode_groups[list(barcode_groups.keys())[0]].keys())

    ungrouped_barcode_types = [bc_type for bc_type in barcode_types if bc_type not in grouped_barcode_types]
    if len(ungrouped_barcode_types) != 1:
        raise ValueError('Enrichment currently requires exactly one barcode type not used for barcode groups.')
    enrichment_bc = '_'.join(ungrouped_barcode_types)

    ### QC ###
    # remove any barcode_type columns that are 'fail'
    drop_idx = []
    for bc_type in barcode_types:
        drop_idx.extend(demux_stats_DF[(demux_stats_DF[bc_type] == 'fail')].index.tolist())
    if screen_no_group:
        # remove any sequences that were not assigned a barcode group
        demux_stats_DF['named_by_group'] = demux_stats_DF.apply(lambda row:
            False if row['output_file_barcodes'] not in barcode_groups.keys() else True, axis=1)
        drop_idx.extend(demux_stats_DF[(demux_stats_DF['named_by_group'] == False)].index.tolist())
    demux_stats_DF.drop(drop_idx, inplace=True)
    ### END QC ###

    # add a demuxed count column for all enrichment barcodes that did not fail demux
    demux_stats_DF['enrichment_bc_no_fail_count'] = demux_stats_DF.groupby('output_file_barcodes')['barcodes_count'].transform('sum')

    demux_stats_DF['tag_bcGroup'] = demux_stats_DF.apply(lambda row:
        row['tag'] + '_' + row['output_file_barcodes'], axis=1)
    
    # pivot the counts to give one row per enrichment barcode, with columns for counts of that enrichment barcode
    #   per each barcode group, then convert to fraction of counts for sample/timepoint
    counts_pivot = demux_stats_DF.pivot(index=enrichment_bc, columns='tag_bcGroup',values='barcodes_count').fillna(0)
    counts_pivot = counts_pivot.reset_index().rename(columns={'index':enrichment_bc})
    # sample_total_counts = counts_pivot.sum(axis=0)
    # # get relative counts for each barcode in each sample to use for weights
    # sample_counts_relative = sample_total_counts / sample_total_counts.mean()
    # count_fractions_pivot = counts_pivot.divide(sample_total_counts, axis=1)
    # count_fractions_pivot = count_fractions_pivot.reset_index().rename(columns={'index':enrichment_bc})

    # Grab the timepoints DF and use this to convert tag_bcGroup column names to timepoint names
    timepoints_DF = pd.read_csv(timepoints, header=1).astype(str)

    # get time unit and sample label
    top_row = [x for x in pd.read_csv(timepoints).columns if 'Unnamed: ' not in x]
    if len(top_row) > 0:
        time_unit = top_row[0]
    else:
        time_unit = 'unit time'    # default label for time unit if none is provided in the first row of timepoints CSV
    sample_label = timepoints_DF.columns[0]

    # add a column to the timepoints DF that gives the replicate number for each sample
    replicateCol = []
    sample_count = Counter()
    for _, row in timepoints_DF.iterrows():
        sample_count[row[sample_label]] += 1
        replicateCol.append(sample_count[row[sample_label]])
    timepoints_DF['replicate'] = replicateCol
    timepoints_DF.set_index([sample_label, 'replicate'], inplace=True)

    # loop through each sample of the timepoints DF and use the columns to slice the enrichment DF
    #   such that we get a DF with (n= (samples per timepoint * replicates) ) rows per enrichment barcode, and one column per timepoint
    #   with the count fraction for that barcode as values, then concatenate to get a single dataframe
    sample_df_list = []
    tp_cols = timepoints_DF.columns.to_list()
    tp_vals = [int(x) for x in tp_cols]
    for index, row in timepoints_DF.iterrows():
        slice_cols = row.to_list()
        for col in slice_cols:
            if col not in counts_pivot.columns:
                print(f'[WARNING] No genotype barcodes above the threshold were identified for tag_barcodeGroup `{col}`, setting counts for this sample to 0')
                counts_pivot = counts_pivot.assign(**{col:0})
        sample_df = counts_pivot[[enrichment_bc] + slice_cols]
        # remove rows with 0 at first timepoint
        # sample_df = sample_df[sample_df[slice_cols[0]] > 0]
        sample_df = sample_df[(sample_df[slice_cols] > 0).all(axis=1)]
        # filter out rows in the bottom 25% of counts at the first timepoint
        sample_df = sample_df[sample_df[slice_cols[0]] > sample_df[slice_cols[0]].quantile(0.25)]
        # get the total counts for only barcodes that appear in all timepoints
        survivors_sum = sample_df[(sample_df[slice_cols] > 0).all(axis=1)][slice_cols].sum(axis=0).to_numpy()
        # ln(count+0.5/survivors_sum+0.5)
        y_normalized = np.log( (sample_df[slice_cols].to_numpy()+0.5) / (survivors_sum+0.5) )

        # time is normalized to the maximum timepoint
        x_normalized = np.array([float(x) for x in tp_cols])
        x_normalized = (x_normalized / max(x_normalized)).reshape(-1,1)

        # perform weighted linear regression for each barcode / sample combination
        fit = LinearRegression().fit(x_normalized, y_normalized.transpose(), sample_weight=survivors_sum)
        
        # add normalized values and fit values to the sample_df
        new_cols = [f'{x}_normalized' for x in tp_cols]
        sample_df[new_cols] = y_normalized
        sample_df['enrichment_score'] = fit.coef_

        sample_cols = {sample_label:index[0], 'replicate':index[1]}
        rename_dict = {old:new for old, new in zip(slice_cols, timepoints_DF.columns.to_list())}
        sample_df = sample_df.assign(**sample_cols).rename(columns=rename_dict)
        sample_df = sample_df[[sample_label, 'replicate', enrichment_bc] + tp_cols + new_cols + ['enrichment_score']]
        sample_df_list.append(sample_df)

    enrichment_df = pd.concat(sample_df_list)

    # # calculate the enrichment score for each barcode / sample combination using linear regression on log transformed data
    # x_time = np.array([float(x) for x in tp_cols]).reshape(-1,1)
    # proportions = enrichment_df[tp_cols].to_numpy()

    # # add a psuedocount to avoid log(0) errors
    # y_proportions_relative =  proportions + 10**-10

    # # get proportions relative to first timepoint, take log of this
    # y_proportions_relative = np.log2( y_proportions_relative / y_proportions_relative[:,0].reshape(-1,1) )

    # # linear regression on log transformed data
    # #   0 values for any timepoint should exclude all subsequent timepoints from the regression for that sample
    # #   to do this without looping through each sample, we can apply the regression for 2, 3, ..., n timepoints separately,
    # #   but parallelized across all samples. We then assign the fits to samples depending on when the first 0 count timepoint appears
    # n_samples, n_timepoints = y_proportions_relative.shape
    # slopes = np.full((n_samples, n_timepoints), np.nan)
    # for n in range(1, n_timepoints):
    #     fit = LinearRegression(fit_intercept=False).fit(x_time[:n+1], y_proportions_relative[:, :n+1].transpose())
    #     slopes[:, n] = fit.coef_.reshape(-1)
    
    # # assign slopes based on where the first 0 appears. If the first 0 appears at timepoint 2, then the slope is the fit for timepoints 1 and 2, etc
    # slope_idx = np.argmax(proportions == 0, axis=1)

    # # For those samples where a zero never occurs, we'll set their index
    # #   so that they always take the last calculated slope
    # slope_idx[np.all(proportions != 0, axis=1)] = n_timepoints - 1
    # slopes_selected = slopes[np.arange(n_samples), slope_idx]

    # enrichment_df[[col + '_relative_enrichment' for col in tp_cols]] = y_proportions_relative
    # enrichment_df['enrichment_score'] = slopes_selected

    if scores_csv:
        enrichment_df.to_csv(scores_csv, index=False)
    if mean_scores_csv:
        mean_enrichment_df = enrichment_df.pivot(index=[sample_label, enrichment_bc], columns='replicate', values='enrichment_score').reset_index()
        mean_enrichment_df['mean_enrichment_score'] = mean_enrichment_df[mean_enrichment_df.columns[2:]].mean(axis=1)
        mean_enrichment_df.to_csv(mean_scores_csv, index=False)
    return enrichment_df

def plot_enrichment(enrichment_df, plots_out):
    """
    Plot enrichment scores for each sample and save to file

    Args:
        enrichment_df (pd.DataFrame): DataFrame with enrichment scores for each sample
        plots_out (str): path to save plots
    """


    sample_label, _, enrichment_bc = enrichment_df.columns[:3]
    samples = enrichment_df[sample_label].unique()
    plots = {}

    for sample in samples:
        sample_df = enrichment_df[enrichment_df[sample_label] == sample]
        replicates = sample_df['replicate'].unique()
        sample_plots = []

        for replicate1, replicate2 in itertools.combinations(replicates, 2):

            # filter for the two replicates, pivot to get both replicates in the same row, remove rows with NaNs
            pivot_df = sample_df[sample_df['replicate'].isin([replicate1, replicate2])].pivot(index=enrichment_bc, columns='replicate', values='enrichment_score').dropna()

            # plot points
            kdims = ["replicate " + str(replicate1), "replicate " + str(replicate2)]
            points = hv.Points((pivot_df[replicate1], pivot_df[replicate2]), kdims=kdims)
            points.opts(title=f'{sample}, Replicates {replicate1} and {replicate2}', width=400, height=400,
                        fontsize={'title':16,'labels':14,'xticks':10,'yticks':10}, size=2, color='black', alpha=0.1)

            # add x=y line and R² fit
            min_value = np.min([pivot_df[replicate1].min(), pivot_df[replicate2].min()])
            max_value = np.max([pivot_df[replicate1].max(), pivot_df[replicate2].max()])
            top_quartile_y = pivot_df[replicate2].quantile(0.75)
            trendline = hv.Curve(([min_value, max_value], [min_value, max_value])).opts(
                xlabel=kdims[0], ylabel=kdims[1], color="lightgrey", line_dash='dashed')
            r2 = r2_score(pivot_df[replicate1], pivot_df[replicate2])
            r2_label = hv.Text(max_value, top_quartile_y, f'R² = {r2:.2f}', halign='right', valign='top').opts(color='black')

            # distributions
            x_hist, y_hist = (histogram(points, dimension=k).opts(fontsize={'title':16,'labels':14,'xticks':10,'yticks':10}) for k in kdims)

            sample_plots.append((trendline * points * r2_label) << y_hist.opts(width=125, color='grey', tools=['hover']) << x_hist.opts(height=125, color='grey', tools=['hover']))
        
        plots[sample] = hv.Layout(sample_plots).cols(1)
    
    final_layout = hv.Layout(plots.values()).cols(len(samples))

    hv.save(final_layout, plots_out, backend='bokeh')



def main():

    enrichment_df = calculate_enrichment(snakemake.input.CSV, snakemake.input.timepoints, snakemake.params.barcodeInfo, snakemake.params.barcodeGroups, scores_csv=snakemake.output.scores, mean_scores_csv=snakemake.output.mean, screen_no_group=snakemake.params.screen_no_group)
    plot_enrichment(enrichment_df, snakemake.output.plots)

if __name__ == '__main__':
    main()