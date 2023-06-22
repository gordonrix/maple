"""
This script calculates enrichment scores from counts of barcodes. It takes as input a CSV file containing the counts of each barcode for each sample, and outputs a CSV file containing the enrichment scores for each sample.
It can also be used to get mean enrichment scores and to plot the enrichment scores in some capacity

The enrichment score is calculated as the slope of the log-linear regression line fitted to the counts of each barcode across the samples,
normalized to the total barcode counts for each timepoint.

The input CSV file should have the following format:
XXX

See Rubin, ... Fowler. Genome Biol 2017. (https://doi.org/10.1186/s13059-017-1272-5) for derivations of enrichment score and standard error calculations. 
"""

import pandas as pd
import numpy as np
import holoviews as hv
from holoviews.operation import histogram
hv.extension('bokeh')
import scipy
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import statsmodels.api as sm
from collections import Counter
import argparse
import itertools
from  concurrent.futures import ProcessPoolExecutor

def print_(string):
    """
    print function that writes to a file
    """
    with open("enrichment_log.txt", "a") as f:
        f.write(string + "\n")

def apply_regression(row, x, timepoints, weighted=False):
    """
    pandas apply function for linear regression using statsmodels.
    x values and weights are taken from each individual row

    Args:
        row (pandas.Series): row of a pandas DataFrame, must contain columns corresponding
            to the timepoints and, if weighted, the weights for each timepoint
        x (np.array): array of normalized x values (range 0-1)
        timepoints (list): list of non-normalized timepoints, used to find column names
        weighted (bool): whether to use weights for the regression
    
    Returns:
        slope (float): slope of the regression line
        se (float): standard error of the slope
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
    worker function for parallelizing the regression apply function
    """
    slopes, ses = zip(*df_chunk.apply(apply_regression, args=[x, timepoints, weighted], axis=1))
    return pd.DataFrame({"enrichment_score": slopes, "standard_error": ses}, index=df_chunk.index)

def regression_parallelize(n_threads, df, x, timepoints, weighted=False):
    """
    parallelize the regression apply function and add the results to new columns in the dataframe

    Args:
        n_threads (int): number of threads to use
        df (pandas.DataFrame): dataframe containing the data to be regressed
        x (np.array): array of normalized x values (range 0-1)
        timepoints (list): list of non-normalized timepoints, used to find column names
        weighted (bool): whether to use weights for the regression

    Returns:
        df (pandas.DataFrame): original dataframe with two additional columns for the slope and standard error of the regression for each row
    """
    df_chunks = np.array_split(df, n_threads)
    with ProcessPoolExecutor(max_workers=n_threads) as executor:
        results = executor.map(regression_worker, df_chunks, itertools.repeat(x), itertools.repeat(timepoints), itertools.repeat(weighted))

    # concatenate results and add as new columns to the dataframe
    results = pd.concat(list(results))
    df = pd.concat([df, results], axis=1)

    return df


def calculate_enrichment(demux_stats, n_threads, timepoints, barcode_info, barcode_groups, reference_bc='', scores_csv='', screen_no_group=True):
    """
    Calculate enrichment scores from counts of barcodes.

    Args:
        demux_stats (str): Path to the demux_stats.csv file output by `demux.py`.
        n_threads (int): Number of threads to use for linear regression parallelization.
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

    # add a demuxed count column for all enrichment barcodes that did not fail demux
    demux_stats_DF['enrichment_bc_no_fail_count'] = demux_stats_DF.groupby('output_file_barcodes')['barcodes_count'].transform('sum')

    demux_stats_DF['tag_bcGroup'] = demux_stats_DF.apply(lambda row:
        row['tag'] + '_' + row['output_file_barcodes'], axis=1)
    

    # pivot the counts to give one row per enrichment barcode, with columns for counts of that enrichment barcode
    #   per each barcode group, then convert to fraction of counts for sample/timepoint
    counts_pivot = demux_stats_DF.pivot(index=enrichment_bc, columns='tag_bcGroup',values='barcodes_count')
    sample_sums = counts_pivot.sum(axis=0)
    counts_pivot_no_na = counts_pivot.dropna()
    counts_pivot = counts_pivot.fillna(0)
    counts_pivot = counts_pivot.reset_index().rename(columns={'index':enrichment_bc})

    # identify the barcode that should be used as a reference
    if reference_bc:
        if reference_bc in counts_pivot_no_na.index:
            print(f"[NOTICE] enrichment.py: Provided enrichment reference_bc {reference_bc} found in all samples. Using this barcode as the enrichment reference.\n")
        else:
            print(f"[WARNING] enrichment.py: Provided enrichment reference_bc {reference_bc} not found in all samples.\n")
            reference_bc = ''
    if reference_bc == '':
        # if pivot no na isnt empty, use the most abundant barcode as the reference
        if counts_pivot_no_na.shape[0] > 0:
            reference_bc = counts_pivot_no_na.divide( counts_pivot_no_na.mean(axis=0) ).sum(axis=1).sort_values(ascending=False).idxmax()
            print(f"[NOTICE] enrichment.py: Using the most abundant barcode, {reference_bc}, as the enrichment reference.\n")
        else:
            print("[NOTICE] enrichment.py: No barcodes found in all samples. Using the total count of all barcodes found in each set of timepoints as the reference.\n")

    # Grab the timepoints DF and use this to convert tag_bcGroup column names to timepoint names
    timepoints_DF = pd.read_csv(timepoints, header=1).astype(str)

    # # get time unit and sample label
    # top_row = [x for x in pd.read_csv(timepoints).columns if 'Unnamed: ' not in x]
    # if len(top_row) > 0:
    #     time_unit = top_row[0]
    # else:
    #     time_unit = 'unit time'    # default label for time unit if none is provided in the first row of timepoints CSV
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
    #   such that we get a DF with (n= (samples per timepoint * replicates) ) rows per enrichment barcode,
    #   and three columns per timepoint (count, normalized and log transformed count, and weights)
    sample_df_list = []
    tp_cols = timepoints_DF.columns.to_list()
    cols_dict = {tp_type:[f'{x}_{tp_type}' for x in tp_cols] for tp_type in ['count', 'normalized', 'weight', 'reference']}
    for index, row in timepoints_DF.iterrows():
        slice_cols = row.to_list()
        for col in slice_cols:
            if col not in counts_pivot.columns:
                print(f'[WARNING] enrichment.py: No genotype barcodes above the threshold were identified for tag_barcodeGroup `{col}`, setting counts for this sample to 0\n')
                counts_pivot = counts_pivot.assign(**{col:0})
        sample_df = counts_pivot.loc[:,[enrichment_bc] + slice_cols]
        
        if reference_bc:
            # get the total counts for the reference barcode
            reference_counts = sample_df[sample_df[enrichment_bc] == reference_bc][slice_cols].sum(axis=0).to_numpy()
        else:
            # get the total counts for only barcodes that appear in all timepoints, used as a reference
            reference_counts = sample_df[(sample_df[slice_cols] > 0).all(axis=1)][slice_cols].sum(axis=0).to_numpy()
        counts = sample_df[slice_cols].to_numpy()


        # set counts following 0 counts to nan then calculate normalized log transformed y values
        counts_pruned = counts.copy()
        mask = (np.cumsum(counts_pruned == 0, axis=1) > 1)
        counts_pruned[mask] = np.nan
        y_normalized = np.log( (counts_pruned+0.5) / (reference_counts+0.5) )

        # calculate weights as 1/ 1/(count proportion)+0.5 + 1/(reference count proportion)+0.5)
        #   weight from survivor sum needs to be broadcasted to the same shape as counts
        count_proportion = counts / np.repeat( sample_sums[slice_cols].to_numpy() [np.newaxis,], counts.shape[0], axis=0)
        reference_proportion = reference_counts / sample_sums[slice_cols].to_numpy()
        weights = 1 / ( ( 1/ (count_proportion + 0.5) ) + np.repeat(( 1/ (reference_proportion + 0.5) )[np.newaxis,], counts.shape[0], axis=0) )
        
        # add normalized values, weights, reference count to the sample_df. reference will be removed before output
        sample_df[cols_dict['normalized']] = y_normalized
        sample_df[cols_dict['weight']] = weights
        sample_df[cols_dict['reference']] = np.repeat(reference_counts[np.newaxis,], counts.shape[0], axis=0)

        sample_cols = {sample_label:index[0], 'replicate':index[1]}
        rename_dict = {old:new for old, new in zip(slice_cols, cols_dict['count'])}
        sample_df = sample_df.assign(**sample_cols).rename(columns=rename_dict)
        sample_df = sample_df.astype({count_col:'int64' for count_col in cols_dict['count']})
        sample_df = sample_df[[sample_label, 'replicate', enrichment_bc] + cols_dict['count'] + cols_dict['normalized'] + cols_dict['weight'] + cols_dict['reference']]
        sample_df_list.append(sample_df)

    enrichment_df = pd.concat(sample_df_list)

    # time is normalized to the maximum timepoint
    tp_vals = [int(x) for x in tp_cols]
    x_normalized = [float(x)/float(max(tp_vals)) for x in tp_vals]

    enrichment_df = regression_parallelize(n_threads, enrichment_df, x_normalized, tp_cols, weighted=True)
    enrichment_df = enrichment_df.drop(columns=cols_dict['reference'])

    if reference_bc: # bring the reference to the top of the dataframe
        ref_df = enrichment_df[enrichment_df[enrichment_bc] == reference_bc]
        no_ref_df = enrichment_df[enrichment_df[enrichment_bc] != reference_bc]
        enrichment_df = pd.concat([ref_df, no_ref_df])


    if scores_csv:
        enrichment_df.to_csv(scores_csv, index=False)

    return enrichment_df

def filter_by_proportions(group, column_proportion_upper):
    """
    Filter a groupby dataframe by one or more columns' quantiles so that each group and each column is filtered to the same extent.
        e.g. if column_proportion_upper = [('count', 0.1, True), ('standard_error', 0.15, False)], then rows with 'count' values in
        the bottom 10% and 'standard_error' values in the top 15% will be removed, but only a total of 15% of rows will be removed
        from each group and each column.

    Args:
        df (pd.DataFrame): DataFrame to filter
        column_proportion_upper list(tuple(str, float, bool): list of tuples of column names (str), proportion (float),
            and whether to keep the upper (True) or lower (False) fraction of columns (bool)
        take_upper (bool): if True, filter out the bottom proportion, if False, filter out the top proportion
    
    Returns:
        pd.DataFrame: filtered DataFrame
    """
    # need to convert the proportion to a specific threshold for each column
    column_threshold_upper = []
    for column, quantile, take_upper in column_proportion_upper:
        if not take_upper: # convert to opposite quantile if take_upper is False
            quantile = 1-quantile
        threshold = group[column].quantile(quantile)
        column_threshold_upper.append((column, threshold, take_upper))
    for column, threshold, take_upper in column_threshold_upper:
        if take_upper:
            group = group[group[column] >= threshold]
        else:
            group = group[group[column] <= threshold]
    return group

def enrichment_mean_filter(enrichment_df, SE_filter=0, t0_filter=0, score_filter=False, filtered_csv='', mean_csv=''):
    """
    Filter enrichment scores by standard error then pivot and calculate average across replicates, then save to csv

    Args:
        enrichment_df (pd.DataFrame): DataFrame with enrichment scores for each sample
        SE_filter (float): proportion of standard error to filter out.
            rows with standard error in the top SE_filter proportion will be filtered out
        t0_filter (float): proportion of t0 counts to filter out.
            rows with t0 counts in the bottom t0_filter proportion will be filtered out
        csv_out (str): path to save csv    
    
    Returns:
        tuple(pd.DataFrame, pd.DataFrame): tuple of (filtered enrichment_df, mean_enrichment_df)
    """
    filter_list = []
    if 0 < SE_filter < 1: # filter by SE
        filter_list.append(('standard_error', SE_filter, False))
    if 0 < t0_filter < 1: # filter by t0
        first_count = enrichment_df.columns[3]
        filter_list.append((first_count, t0_filter, True))
    
    sample_label, _, enrichment_bc = enrichment_df.columns[:3]
    if filter_list:
        enrichment_df = enrichment_df.groupby([sample_label, 'replicate']).apply(filter_by_proportions, filter_list).reset_index(drop=True)

    if score_filter: # filter by enrichment score
        enrichment_df = enrichment_df[enrichment_df['enrichment_score'] >= score_filter]

    mean_enrichment_df = enrichment_df.pivot(index=[sample_label, enrichment_bc], columns='replicate', values='enrichment_score')
    mean_enrichment_df.columns = [f'replicate_{rep}' for rep in mean_enrichment_df.columns]
    count_col = mean_enrichment_df.count(axis=1)
    mean_enrichment_df['mean_enrichment_score'] = mean_enrichment_df.mean(axis=1)
    mean_enrichment_df['valid_replicates'] = count_col
    mean_enrichment_df.reset_index(inplace=True)

    if filtered_csv:
        enrichment_df.to_csv(filtered_csv, index=False)
    if mean_csv:
        mean_enrichment_df.to_csv(mean_csv, index=False)

    return enrichment_df, mean_enrichment_df

def plot_enrichment(enrichment_df, plots_out):
    """
    Plot enrichment scores for each sample and save to file

    Args:
        enrichment_df (pd.DataFrame): DataFrame with enrichment scores for each sample
        plots_out (str): path to save plots
    """

    sample_label, _, enrichment_bc = enrichment_df.columns[:3]
    samples = list(enrichment_df[sample_label].unique())
    samples.sort()
    plots = {}
    # use same min/max values for fit lines and label location
    min_value = enrichment_df['enrichment_score'].quantile(0.01)
    max_value = enrichment_df['enrichment_score'].quantile(0.99)

    for sample in samples:
        sample_df = enrichment_df[enrichment_df[sample_label] == sample]
        replicates = list(sample_df['replicate'].unique())
        replicates.sort()
        sample_plots = []

        for replicate1, replicate2 in itertools.combinations(replicates, 2):

            # filter for the two replicates, pivot to get both replicates in the same row, remove rows with NaNs
            pivot_df = sample_df[sample_df['replicate'].isin([replicate1, replicate2])].pivot(index=enrichment_bc, columns='replicate', values='enrichment_score').dropna()

            # plot points
            kdims = ["replicate " + str(replicate1), "replicate " + str(replicate2)]
            points = hv.Points((pivot_df[replicate1], pivot_df[replicate2]), kdims=kdims)
            points.opts(title=f'{sample}, Replicates {replicate1} and {replicate2}', width=400, height=400,
                        fontsize={'title':16,'labels':14,'xticks':10,'yticks':10}, size=2, color='black', alpha=0.1)

            # use statsmodels OLS to calculate R² between the two replicates
            y = pivot_df[replicate2].to_numpy()
            x = sm.add_constant(pivot_df[replicate1].to_numpy())
            fit = sm.OLS(y,x).fit()
            fit_r2 = fit.rsquared
            intercept = fit.params[0]
            slope = fit.params[1]

            # use hv.Curve to show the fit line
            trendline = hv.Curve(([min_value, max_value], [min_value * slope + intercept, max_value * slope + intercept])).opts(
                xlabel=kdims[0], ylabel=kdims[1], color="lightgrey", line_dash='dashed')

            # also calculate the R² for fit to x=y line
            r2 = r2_score(pivot_df[replicate1], pivot_df[replicate2])
            r2_label = hv.Text(min_value, max_value, f' slope 1, R² = {r2:.2f}\n slope {slope:.2f}, R² = {fit_r2:.2f}', halign='left', valign='top').opts(color='black')

            # distributions
            x_hist, y_hist = (histogram(points, dimension=k).opts(fontsize={'title':16,'labels':14,'xticks':10,'yticks':10}) for k in kdims)

            sample_plots.append((trendline * points * r2_label) << y_hist.opts(width=125, color='grey', tools=['hover']) << x_hist.opts(height=125, color='grey', tools=['hover']))
        
        plots[sample] = hv.Layout(sample_plots).cols(1)
    
    final_layout = hv.Layout(plots.values()).cols(len(samples))

    hv.save(final_layout, plots_out, backend='bokeh')


if __name__ == '__main__':
    enrichment_df = calculate_enrichment(snakemake.input.CSV, snakemake.threads, snakemake.input.timepoints, snakemake.params.barcodeInfo, snakemake.params.barcodeGroups, reference_bc=snakemake.params.reference_bc, scores_csv=snakemake.output.scores, screen_no_group=snakemake.params.screen_no_group)