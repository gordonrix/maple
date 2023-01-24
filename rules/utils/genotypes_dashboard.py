import pathlib
import pandas as pd
import panel as pn
import holoviews as hv
import numpy as np
import argparse

pn.extension('tabulator', template='material')

import colorcet as cc
import hvplot.pandas # noqa

from holoviews.streams import Selection1D, BoundsXY
from holoviews.operation.datashader import datashade, rasterize, dynspread
from bokeh.models import HoverTool
from bokeh import palettes

from common import dist_to_DF
from dimension_reduction_genotypes import seq_array_from_genotypes, SequenceEncoder
from plot_distribution import plot_dist

from Bio import SeqIO
from Bio.Seq import Seq

from timeit import default_timer as now

## define functions that will be used by pipeline

def add_MOI_column(ds, onehotArrayDict, muts, NTorAA, max_groups):
    """
    adds a column to an hv.Dataset that tracks the NT or AA mutations of interest as both a comma separated list
        and as a count of the number of mutations of interest contained within each sequence
        
    parameters:
        ds (hv.Dataset):        genotypes.csv file converted into an hv.Dataset object
        onehotArrayDict (dict): dictionary containing onehot encoded sequences and related information
                                    onehot array is shape (N,L,C) where N = len(genotypes), L = len(refSeq),
                                    and C = len(characters) (e.g. 21 for amino acids+stop)
        muts (str):             comma separated list of mutations of interest
                                    can be a full mutation, or just the position, or just a
                                    mutated AA/NT. For example, "A100G, 150, G_, _C" would mark
                                    any sequences that have the specific A100G mutation or a
                                    mutation to position 150, or a mutation to any G, or any mutation
                                    that produces a C
       NTorAA (str):            "NT" or "AA"
       max_groups (int):        maximum number of unique groups two include, decided based on frequency
                                    the least common groups will be combined under the name 'other',
                                    and sequences without any of the listed mutations will be given the name 'none'
    returns:
        hv.Dataset with one additional column
    """

    df = ds.data
    
    idx_all = np.arange(onehotArrayDict[NTorAA]['array'].shape[0])
    selected_idx = np.array(ds.data.index)
    arrIdx_in_selected = np.where(np.isin(idx_all, selected_idx))[0]
    onehot_selected = onehotArrayDict[NTorAA]['array'][arrIdx_in_selected,:,:]
    
    
    # onehot_selected = onehotArrayDict[NTorAA]['array']
    
    refSeq = onehotArrayDict[NTorAA]['refSeq']
    all_L_idx = [x for x in range(0,len(refSeq))]
    charactersDict = {'NT':'ATGC-', 'AA':'ACDEFGHIKLMNPQRSTVWY*'}
    characters = charactersDict[NTorAA]
    all_C_idx = [i for i in range(0,len(characters))]
    
    muts_list = muts.replace(' ','').split(',')
    numbers = [str(x) for x in range(0,10)]
    
    # build a 2D array of shape (L,C) of onehot encoded mutations of interes
    #    to be broadcast and multiplied by the 3D onehot array of shape (N,L,C)
    #    to find genotypes that have one or more of the mutations in muts
    all_muts_2D_arr = np.zeros((len(refSeq), len(charactersDict[NTorAA])))
    
    for mut_input in muts_list:
        mut = mut_input
        position_specific, ID_specific = True, True
        
        if mut[0] not in numbers:
            if mut[0] == '_':
                position_specific = False
            else:
                wt_input = mut[0].upper()
            mut = mut[1:]
        if mut[-1] not in numbers:
            if mut[-1] == '_':
                position_specific = False
            else:
                C_idx = charactersDict[NTorAA].find(mut[-1])
            mut = mut[:-1]
        else:
            ID_specific = False
            C_idx = all_C_idx
            
        if position_specific:
            L_idx = int(mut)-1
            wt = refSeq[L_idx].upper()
            if wt != wt_input:
                print(f'[ERROR] User provided wt nucleotide/residue for {NTorAA} mutation {mut_input} does not match the wt nucleotide/residue {wt}. Please double check indexing.')
        else:
            L_idx = all_L_idx
            
        all_muts_2D_arr[L_idx, C_idx] = 1
        
    # zero out all non-mutations by subtracting out the onehot encoded wild type sequence
    #    then converting any -1s to 0
    wt_onehot = SequenceEncoder(refSeq, characters).onehot
    all_muts_2D_arr = np.subtract( all_muts_2D_arr, wt_onehot )
    all_muts_2D_arr[all_muts_2D_arr<0] = 0
        
    # perform matrix multiplication between the 2D onehot array representing all mutations 
    #    of interest and the 3D onehot array representing all sequences resulting 3D array
    #    will have ones at positions where a sequence contained a mutation of interest, 0s elsewhere
    all_muts_2D_arr = (all_muts_2D_arr==1) # mask non-1 values as False so these indices will be 0 in the result
    all_muts_3D_arr = np.einsum('NLC,LC->NLC', onehot_selected, all_muts_2D_arr)
    
    # find all unique combinations of mutations of interest, sort by count, then loop through the most common
    #    ones, assign strings to them, then assign these strings to the correct row of the dataset
    unique_MOI_combo_arrs, idx_original, counts = np.unique(all_muts_3D_arr, axis=0, return_counts=True, return_inverse=True)
    sorted_idx = np.argsort(-counts)
    MOI_combo_strings_sorted = []
    none_found = False # flag for doing a search for the MOI_combo without any mutations of interest
    
    for num_combos, unique_idx in enumerate(sorted_idx):
        
        # if no more groups need to be made, fill remaining values with 'other'
        if num_combos > (max_groups-1): # convert to 1-indexed
            remaining_genotypes = len(sorted_idx) - len(MOI_combo_strings_sorted)
            MOI_combo_strings_sorted.extend(['other']*(remaining_genotypes))
            break

        combo_arr = unique_MOI_combo_arrs[unique_idx,:,:]
        mut_indices = np.argwhere(combo_arr == 1)
        if len(mut_indices) > 0:
            muts = []
            for L_idx, C_idx in mut_indices:
                wt = refSeq[L_idx]
                mut = characters[C_idx]
                muts.append(wt + str(L_idx+1) + mut)
            mut_string = ', '.join(muts)
        else:
            max_groups += 1
            none_found = True
            mut_string = 'none'

        MOI_combo_strings_sorted.append(mut_string)
        
    MOI_combo_strings_sorted = np.array(MOI_combo_strings_sorted)
    MOI_combo_strings_unsorted = MOI_combo_strings_sorted[sorted_idx] # convert to original order of unique arrays with reverse indexing

    # search for unique array without any MOIs if it hasn't already been found
    if not none_found:
        none_arr = np.zeros( (len(refSeq), len(AAs)) )
        none_arr_matches = np.any(np.all( unique_MOI_combo_arrs == none_arr, axis=(1,2)), axis=0)
        none_idx = np.argmax(none_arr_matches) if none_arr_matches.any() else None
        if none_idx!=None:
            MOI_combo_strings_unsorted[none_idx] = 'none'

    MOI_combo_strings_original = MOI_combo_strings_unsorted[idx_original] # convert to original order of all arrays/genotypes with reverse indexing
    
    df[f'{NTorAA}_muts_of_interest'] = MOI_combo_strings_original
    return hv.Dataset(df)

def aggregate_mutations(onehotArr, arr_idx, selected_idx, NTorAA, refSeq):
    """
    Generate a tidy format pd.DataFrame of counts of all mutations within the provided onehot-encoded
        sequences among the selected indices
    
    Parameters:
        onehotArr (np.array):   3D array of onehot encoded sequences, output by seq_array_from_genotypes
        arrIdx (np.array):      1d array of indices of onehot encoded sequences for matching to selection indices
        selectedIdx (np.array): 1d array of indices of a selection
        NTorAA (str):           'NT' or 'AA', type of analysis
        refSeq (str):           NT or AA sequence to be used for reference to know the WT NT or AA at each position
    
    Returns:
        pd.DataFrame: A tidy dataframe that tabulates counts for all possible mutations
        int:          total number of sequences analyzed
    """

    idx_all = np.arange(onehotArrayDict[NTorAA]['array'].shape[0])
    # selected_idx = np.array(ds.data.index)
    arrIdx_in_selected = np.where(np.isin(idx_all, selected_idx))[0]
    onehot_selected = onehotArrayDict[NTorAA]['array'][arrIdx_in_selected,:,:]
    
    # arrIdx_in_selected = np.where(np.isin(arr_idx, selected_idx))[0]
    # onehot_selected = onehotArr[arrIdx_in_selected,:,:]
    total_seqs = len(onehot_selected)
    
    aggregatedMuts = np.sum(onehot_selected, axis=0, dtype=np.int32)

    charactersDict = {'NT':'ATGC-', 'AA':'ACDEFGHIKLMNPQRSTVWY*'}
    rows = []

    chars = charactersDict[NTorAA]
    # loop through all possible mutations and generate a dataframe of all, including those with count of 0
    for index, count in np.ndenumerate(aggregatedMuts):
        posi, charIdx = index
        wt = refSeq[posi].upper()
        posi += 1
        mut = chars[charIdx]
        if wt != mut:
            rows.append([wt, posi, mut, count])
            
    df = pd.DataFrame(rows, columns=['wt', 'position', 'mutation', 'total_count'])
    df['proportion_of_seqs'] = df['total_count']/total_seqs

    return df, total_seqs

def conspicuous_mutations(df, num_positions, colormap, most_common=True):
    """
    produces a bar plot of the most or least frequent mutations
    
    parameters:
        df (pd.DataFrame):   dataframe of aggregated mutations output by aggregate_mutations
        num_positions (int): number of mutations to include in the bar plot output
        colormap (dict):     AA/NT letter : color hex code key:value pairs to use for the plot
        most_common (bool):  if True/False, output the most/False commonly mutated positions
        
    returns:
        hv.Bars object showing the topN most frequently observed mutations in the aggregated
            mutations dataframe
    """
    
    df = df.sort_values(['total_count','position'], ascending=[(not most_common),True])
    df_grouped = df.groupby('position', as_index=False).sum().sort_values(['total_count','position'],ascending=[not most_common, True])
    positions = df_grouped['position'].iloc[:num_positions]
    df = df[df['position'].isin(positions)]
    df = df.sort_values('position', ascending=True)
    df['position'] = df['wt'] + df['position'].astype(str)
    plot = hv.Bars(df, kdims=['position','mutation'], vdims=['proportion_of_seqs', 'total_count']).opts(
                    show_legend=False, height=500, xlabel='position',
                    xrotation=40, stacked=True, cmap=colormap, tools=['hover'])
    return plot

def plot_mutations_aggregated(ds, onehotArrayDict, idx, num_positions, NTorAA, plot_type):
    """
    produces plots that describe aggregated mutations. plot is chosen based on the text in the plot_type input
    """
    
    NTorAA = NTorAA[-3:-1]
    if plot_type.endswith('frequent'):
        if plot_type.startswith('most'):
            most_common=True
        elif plot_type.startswith('least'):
            most_common=False
        selected_idx = np.array(ds.data.index)
        total_seqs = len(selected_idx)
        df, total_seqs = aggregate_mutations(
                                 onehotArrayDict[NTorAA]['array'],
                                 arr_idx=idx,
                                 selected_idx=selected_idx,
                                 NTorAA=NTorAA,
                                 refSeq=onehotArrayDict[NTorAA]['refSeq'])
        plot = conspicuous_mutations(df, num_positions, onehotArrayDict[NTorAA]['colormap'],
                                     most_common=most_common).opts(ylabel=f"frequency (n={total_seqs})")
    
    return plot




parser = argparse.ArgumentParser()

parser.add_argument('--genotypes', type=str, help='Name of the csv file that describes genotypes')
parser.add_argument('--reference', type=str, help='Name of the fasta file that contains alignment, nucleotide, and ORF sequenes')
args = parser.parse_args()

# load in data
data = pd.read_csv(pathlib.Path(args.genotypes))
embedding_options = [col.split('_')[0] for col in list(data.columns) if col.endswith('_PaCMAP1')]
do_AA_analysis = 'AA' in embedding_options

# remove genotypes that don't have all 2D embedding options
for embedding in embedding_options:
    data = data[data[embedding+'_PaCMAP1'].notna()]
data = data.reset_index().drop('index', axis='columns')

## Create two 3D arrays of onehot-encoded sequences, one for NT and one for AA. These will be used for vectorized
#   queries of mutations. Index corresponds to the index of 'data'. Also adding colormaps and reference sequences

# nested dictionary to hold the arrays for both NT and AA analysis
onehotArrayDict = {'NT':{'colormap':{'A':palettes.Greens[3][1], #take the middle color from the 3 length color list
                                    'T':palettes.Reds[3][1],
                                    'G':'#000000',           #black
                                    'C':palettes.Blues[3][1],
                                    '-':'#d3d3d3'}}}          #grey

if do_AA_analysis:
    AAs_by_group = [['K','R','H'],              # positive, red
                    ['D','E'],                  # negative, green
                    ['F','Y','W'],              # aromatic, purple
                    ['A','V','I','L'],          # small hydrophobic, blue
                    ['C','M','S','T','N','Q']]  # sulfurous and polar uncharged, yellow->orange

    colorDictList= [palettes.Reds,
                    palettes.Greens,
                    palettes.Purples,
                    palettes.Blues,
                    palettes.YlOrBr]

    amino_acid_colormap = {}
    for AAs,colorDict in zip(AAs_by_group,colorDictList):
        colorList = colorDict[len(AAs)+1][::-1]
        for i,AA in enumerate(AAs):
            amino_acid_colormap[AA] = colorList[i+1]
    amino_acid_colormap.update({'P':'#FA11F2','G':'#FEFBEA','*':'#000000','-':'#d3d3d3'}) # pink and cream for proline and glycine, black for stop, grey for gap
    
    onehotArrayDict.update({'AA':{'colormap':amino_acid_colormap}})

for seqIndex, seqType in enumerate(onehotArrayDict.keys(), start=1):
    refSeq = str(list(SeqIO.parse(args.reference, 'fasta'))[seqIndex].seq).upper()
    if seqType == 'AA':
        refSeq = Seq.translate(refSeq)
    onehot3DArray, genotypesDF = seq_array_from_genotypes(refSeq, data, seqType, onehot=True)
    onehotArrayDict[seqType].update( {'refSeq':refSeq, 'array':onehot3DArray} )



## Convert into a hv.Dataset, then use some widgets to downsample data and add a column for the size of points
downsample_slider = pn.widgets.IntSlider(name='downsample',
                                        start=1000,
                                        end=len(data),
                                        step=1000,
                                        value=10000)

size_column_select = pn.widgets.Select(name='point size column', options=['NT_substitutions_count', 'AA_substitutions_nonsynonymous_count', 'count'])
size_range_slider = pn.widgets.IntRangeSlider(name='point size range',
                                    start=1,
                                    end=100,
                                    step=1,
                                    value=(10,60))

def initialize_ds(ds, downsample, size_column, size_range):
    df = ds.data
    
    if downsample < len(df):
        df = df.sample(n=downsample, random_state=0)
    
    # add a size column scaled to user provided point size values
    minSize, maxSize = size_range
    sizeCol = df[size_column]
    maxSizeCol = sizeCol.max()
    minSizeCol = sizeCol.min()
    if minSizeCol==maxSizeCol:
        df['size'] = minSize
    else:
        df['size'] = ( (sizeCol-minSizeCol) / (maxSizeCol-minSizeCol) ) * (maxSize-minSize) + minSize
    
    return hv.Dataset(df)

downsampled = hv.Dataset(data).apply(initialize_ds,
                                downsample=downsample_slider,
                                size_column=size_column_select,
                                size_range=size_range_slider)



## add columns for coloring by nucleotide or amino acid mutations of interest

NT_muts_text = pn.widgets.TextInput(name='NT muts of interest', placeholder='Input comma separated list of muts e.g. "A100T, 150, G_, _C"')
AA_muts_text = pn.widgets.TextInput(name='AA muts of interest', placeholder='Input comma separated list of muts e.g. "A100T, 150, G_, _C"', disabled=not do_AA_analysis)
max_mut_combos_slider = pn.widgets.IntSlider(name='maximum number of mutation of interest combinations',
                                        start=1,
                                        end=100,
                                        step=1,
                                        value=10)

def add_muts_of_interest_columns(ds, onehotArrayDict, NT_muts, AA_muts, max_groups, do_AA_analysis):
    if NT_muts != '':
        ds = add_MOI_column(ds, onehotArrayDict, NT_muts, 'NT', max_groups)
    if AA_muts != '' and do_AA_analysis:
        ds = add_MOI_column(ds, onehotArrayDict, AA_muts, 'AA', max_groups)
    return ds

dataset = downsampled.apply(add_muts_of_interest_columns,
                                onehotArrayDict=onehotArrayDict,
                                NT_muts=NT_muts_text,
                                AA_muts=AA_muts_text,
                                max_groups=max_mut_combos_slider,
                                do_AA_analysis=do_AA_analysis)



## define some widgets for filtering and apply filters

filter_by_select = pn.widgets.Select(name='substitutions count filter', options=['NT_substitutions_count', 'AA_substitutions_nonsynonymous_count'], value='NT_substitutions_count')

filter_range_slider = pn.widgets.IntRangeSlider(name='NT mutations range',
                                    start=round(data['NT_substitutions_count'].min()),
                                    end=round(data['NT_substitutions_count'].max()),
                                    step=1,
                                    value=( round(data['NT_substitutions_count'].min()),round(data['NT_substitutions_count'].max()) ))

# update the range of the filter based on which column is selected
def range_widget_callback(IntRangeSlider, event):
    column = event.new
    NTorAA = column[:2]
    minVal = round(data[column].min())
    maxVal = round(data[column].max())
    IntRangeSlider.name  = f"{NTorAA} substitutions range"
    IntRangeSlider.start = minVal
    IntRangeSlider.end   = maxVal
    IntRangeSlider.value = (minVal, maxVal)
    
filter_by_select.link(filter_range_slider, callbacks={'value':range_widget_callback})

count_range_slider = pn.widgets.IntRangeSlider(name='genotype count range',
                                    start=round(data['count'].min()),
                                    end=round(data['count'].max()),
                                    step=1,
                                    value=( round(data['count'].min()),round(data['count'].max()) ))

# multi select widget that depends on the kind of genotypes CSV file columns
for choice in ['barcode_group','timepoint']:
    if choice in data.columns:
        choice_col = choice
group_choice = pn.widgets.MultiChoice( name=f'select {choice_col}', value=list(data[choice_col].unique()),
                                        options=list(data[choice_col].unique()) )

def filter_ds(ds, filter_by, filter_range, count_range, group):
    df = ds.data
    df = df[
        (df['count'] >= count_range[0]) &
        (df['count'] <= count_range[1])]

    df = df[
        (df[filter_by] >= filter_range[0]) &
        (df[filter_by] <= filter_range[1])]
    
    df = df[df[choice_col].isin(group)]

    return hv.Dataset(df)

filtered_ds = dataset.apply(filter_ds,
                            filter_by=filter_by_select,
                            filter_range=filter_range_slider,
                            count_range=count_range_slider,
                            group=group_choice)



## define widgets for the static points plot, then make a points plot that will be used for tracking selections, but
#   will just show grey empty circles which will remain when a subset are selected. selected points will then
#   be used to further filter the data via a selection stream, and this selected data will be used by other plots

color_options = ['NT_substitutions_count', 'AA_substitutions_nonsynonymous_count', 'count', 'NT_muts_of_interest', choice_col]
if do_AA_analysis:
    color_options.append('AA_muts_of_interest')
embedding_select = pn.widgets.Select(name='sequence embedding', options=embedding_options, value=embedding_options[0])

tools = ['box_select', 'lasso_select']

def points(ds, embedding, tools):
    dim1, dim2 = f"{embedding}_PaCMAP1", f"{embedding}_PaCMAP2"
    return ds.data.hvplot(kind='points', x=dim1, y=dim2, size='size', hover_cols=color_options, xticks=[100], yticks=[100]).opts(
        xlabel=dim1, ylabel=dim2, height=600, width=800, tools=tools)

static_points = dataset.apply(points, embedding=embedding_select, tools=tools)
index_stream = Selection1D(source=static_points)

# filter dataset by selection if any points are selected
def select_ds(ds, idx_selected):
    df = ds.data
    df2 = df[df.index.isin(idx_selected)]
    if len(df2)==0:
        return hv.Dataset(df)
    else:
        return hv.Dataset(df2)

selected_ds = filtered_ds.apply(select_ds, idx_selected=index_stream.param.index)



## define widgets for the dynamic points plot, then use pn.bind to combine this with the static points plot

color_by_select = pn.widgets.Select(name='color by', options=color_options)

cmaps  = ['kbc', 'fire', 'bgy', 'bgyw', 'bmy', 'gray', 'rainbow4']
cmaps_r = [c+'_r' for c in cmaps]
cmaps = cmaps_r+cmaps
cmap_selector = pn.widgets.Select(name='colormap', options=cmaps)

selected_alpha_slider = pn.widgets.FloatSlider(name='selected point opacity',
                                    start=0.01,
                                    end=1.00,
                                    step=0.01,
                                    value=1.0 )

unselected_alpha_slider = pn.widgets.FloatSlider(name='unselected point opacity',
                                    start=0.00,
                                    end=1.00,
                                    step=0.01,
                                    value=0.1 )

# I have to use panel.bind for the dynamic scatter plot because I want to be able to switch between numerical and categorical color labelling, but there's a bug when using apply to do that. see https://github.com/holoviz/holoviews/issues/5591
def points_bind(ds, embedding, filtered_ds, points_unfiltered, colorby, cmap, selected_alpha, unselected_alpha):
    points_unfiltered.opts(fill_alpha=0, line_alpha=unselected_alpha, nonselection_line_alpha=unselected_alpha, color='grey')
    points_filtered = filtered_ds.apply(points, embedding=embedding_select, tools=[]).apply.opts(alpha=selected_alpha, nonselection_alpha=unselected_alpha, color=colorby, cmap=cmap)
    
    return points_unfiltered*points_filtered

dynamic_points = pn.bind(points_bind, ds=dataset, embedding=embedding_select,
                        points_unfiltered=static_points,
                        filtered_ds=filtered_ds,
                        colorby=color_by_select,
                        cmap=cmap_selector,
                        unselected_alpha=unselected_alpha_slider,
                        selected_alpha=selected_alpha_slider)



## histogram for a column chosen by a widget. selection will be highlighted by combining a low opacity histogram
#   of the pre-filter dataset with a high opacity histogram of the selected dataset

hist_col_selector = pn.widgets.Select(name='histogram column', options=['count', 'NT_substitutions_count', 'AA_substitutions_nonsynonymous_count'], value='NT_substitutions_count')

def histogram(ds, histCol):
    # calculate bin range based on maximum value, 1:1 bin:mutations ratio
    max_x_val = ds.data[histCol].max()
    bins = np.linspace(-0.5, max_x_val-0.5, max_x_val+1)
    return ds.data.hvplot(y=histCol, kind='hist', width=500,height=500, bins=bins, xlabel=histCol, ylabel='total sequences', color='grey')

dynamic_hist = selected_ds.apply(histogram, histCol=hist_col_selector).opts(alpha=1)
static_hist = dataset.apply(histogram, histCol=hist_col_selector).opts(alpha=0.1)



## plot that describes mutations aggregated over the selected genotypes
# Current options:
#   - bar plot of most frequent mutations
#   - bar plot of least frequent mutations
agg_muts_width_slider = pn.widgets.IntSlider(name='plot width',
                                start=160,
                                end=1200,
                                step=10,
                                value=590 )

num_positions_slider = pn.widgets.IntSlider(name='number of positions',
                                    start=5,
                                    end=100,
                                    step=1,
                                    value=10 )

NTorAA_radio_options = ['nucleotide (NT)']
if do_AA_analysis:
    NTorAA_radio_options.append('amino acid (AA)')
NTorAA_radio = pn.widgets.RadioButtonGroup(
    name='bar plot mutation type', options=NTorAA_radio_options, value=NTorAA_radio_options[-1], disabled=(not do_AA_analysis), button_type='default')
agg_plot_type_selector = pn.widgets.Select(name='aggregated mutations plot type', options=['most frequent', 'least frequent'])
aggregated_muts_plot = selected_ds.apply(plot_mutations_aggregated,
                                        onehotArrayDict=onehotArrayDict,
                                        idx=data.index.to_numpy(),
                                        num_positions=num_positions_slider,
                                        NTorAA=NTorAA_radio,
                                        plot_type=agg_plot_type_selector)
aggregated_muts_panel = pn.panel(aggregated_muts_plot, width=agg_muts_width_slider.value)
def mod_width_callback(muts_panel, event):
    muts_panel.width=event.new
agg_muts_width_slider.link(aggregated_muts_panel, callbacks={'value':mod_width_callback})

snakemake_dir = pathlib.Path(__file__).parent.parent.parent

## build the layout
layout = pn.Column(
pn.Row(
    dynamic_points,
    pn.Column(
        downsample_slider,
        pn.Row(selected_alpha_slider,unselected_alpha_slider),
        pn.Row(size_column_select,size_range_slider),
        pn.Row(pn.Column(embedding_select,count_range_slider),group_choice),
        pn.Row(filter_by_select,filter_range_slider),
        pn.Row(color_by_select,cmap_selector),
        pn.Row(NT_muts_text,AA_muts_text),
        max_mut_combos_slider),
    ),
pn.Row(
    pn.panel(pathlib.Path(snakemake_dir/'images'/'dashboard_legend.png'),width=200,align='start'),
    pn.Column(
        aggregated_muts_panel, NTorAA_radio,
        agg_plot_type_selector,
        num_positions_slider,
        agg_muts_width_slider,
        ),
    pn.Column(
        dynamic_hist*static_hist,
        hist_col_selector)
    ))

layout.servable(title=f'maple dashboard, file: {args.genotypes.split("/")[-1]}')