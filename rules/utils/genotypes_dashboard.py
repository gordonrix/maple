#
#  DESCRIPTION   : Custom python script used in the maple snakemake pipeline.
#                   Runs an interactive dashboard for visualizing genotypes.
#
#  RESTRICTIONS  : none
#
#  AUTHOR(S)     : Gordon Rix
#

import pathlib
import pandas as pd
import panel as pn
import holoviews as hv
import numpy as np
import datashader as ds
import colorcet as cc
import argparse

pn.extension('tabulator', template='material')

import hvplot.pandas # noqa

from holoviews.streams import Selection1D
from holoviews.selection import link_selections
from holoviews.operation.datashader import datashade, rasterize, dynspread
from datashader.colors import Sets1to3
from bokeh.models import HoverTool
from bokeh.models import LinearColorMapper, ColorBar, Label
from bokeh.plotting import figure
import matplotlib
from colorcet import palette

from common import conspicuous_mutations, colormaps, export_svg_plots, cmap_dict
from SequenceAnalyzer import SequenceAnalyzer
from plot_distribution import plot_dist

from Bio import SeqIO
from Bio.Seq import Seq

from timeit import default_timer as now
import datetime

parser = argparse.ArgumentParser()

parser.add_argument('--genotypes', type=str, help='Name of the csv file that describes genotypes')
parser.add_argument('--reference', type=str, help='Name of the fasta file that contains alignment, nucleotide, and ORF sequenes')
parser.add_argument('--exclude_indels', action='store_true', help='''Whether to use sequences with indels or not.
                        Default behavior is to include sequences that have indels, though insertions will not influence the 
                        location of a sequence in the 2D plot''')
args = parser.parse_args()

# Use SequenceAnalyzer to handle the data
all_data = SequenceAnalyzer(reference_fasta=args.reference, genotypesCSV=args.genotypes, exclude_indels=args.exclude_indels)
do_AA_analysis = all_data.do_AA_analysis
if 'NT_muts_of_interest' not in all_data.genotypes.columns:
    all_data.genotypes['NT_muts_of_interest'] = 'none'
if do_AA_analysis:
    if 'AA_muts_of_interest' not in all_data.genotypes.columns:
        all_data.genotypes['AA_muts_of_interest'] = 'none'

## Convert into a hv.Dataset, then use some widgets to downsample data and add a column for the size of points

ds_checkbox = pn.widgets.Checkbox(name='datashade (unticking may impede performance)', value=True)

size_column_select = pn.widgets.Select(name='point size column', options=['NT_substitutions_count', 'AA_substitutions_nonsynonymous_count', 'count'], disabled=ds_checkbox.value)
size_range_slider = pn.widgets.IntRangeSlider(name='point size range',
                                    start=1,
                                    end=100,
                                    step=1,
                                    value=(10,60),
                                    disabled=ds_checkbox.value)

# size of points can be controlled only when datashading is off
def ds_widget_callback(widget, event):
    do_ds = event.new
    widget.disabled  = do_ds

ds_checkbox.link(size_column_select, callbacks={'value':ds_widget_callback})
ds_checkbox.link(size_range_slider, callbacks={'value':ds_widget_callback})

downsample_slider = pn.widgets.IntSlider(name='downsample',
                                        start=1000,
                                        end=len(all_data.genotypes),
                                        step=1000,
                                        value=len(all_data.genotypes) )

num_barcodes_slider = pn.widgets.IntSlider(name='maximum number of barcodes',
                                            start = 1,
                                            end = 100,
                                            step = 1,
                                            value = 10,
                                            disabled = not 'barcode(s)' in all_data.genotypes.columns)

def initialize_ds(dataset, all_data, downsample, size_column, size_range, num_barcodes, NT_muts, AA_muts, max_groups):
    idx = all_data.downsample(downsample)['df'].index
    df = dataset.data
    
    # add a size column scaled to user provided point size values
    minSize, maxSize = size_range
    sizeCol = df[size_column]
    maxSizeCol = sizeCol.max()
    minSizeCol = sizeCol.min()
    if minSizeCol==maxSizeCol:
        df.loc[idx, 'size'] = minSize
    else:
        df.loc[idx, 'size'] = ( (sizeCol-minSizeCol) / (maxSizeCol-minSizeCol) ) * (maxSize-minSize) + minSize

    # add a new column for the subset of the most common barcodes
    if not num_barcodes_slider.disabled:
        barcode_counts = df['barcode(s)'].value_counts()
        if len(barcode_counts) > num_barcodes:
            top_barcodes = barcode_counts[:num_barcodes].index
            df.loc[idx, 'barcode(s)_subset'] = df['barcode(s)'].where(df['barcode(s)'].isin(top_barcodes), 'other')
        else:
            df.loc[idx, 'barcode(s)_subset'] = df['barcode(s)']

    # add new columns for mutations of interest
    df.loc[df.index,'NT_muts_of_interest'] = all_data.get_mutations_of_interest('NT', NT_muts, max_groups, idx=df.index)
    df.loc[df.index,'AA_muts_of_interest'] = all_data.get_mutations_of_interest('AA', AA_muts, max_groups, idx=df.index)
    return hv.Dataset(df)


## add columns for coloring by nucleotide or amino acid mutations of interest

NT_muts_text = pn.widgets.TextInput(name='NT muts of interest', placeholder='Input comma separated list of muts e.g. "A100T, 150, G_, _C"')
AA_muts_text = pn.widgets.TextInput(name='AA muts of interest', placeholder='Input comma separated list of muts e.g. "A100T, 150, G_, _C"', disabled=not do_AA_analysis)
max_mut_combos_slider = pn.widgets.IntSlider(name='maximum number of mutation of interest combinations',
                                        start=1,
                                        end=100,
                                        step=1,
                                        value=10)

downsampled_MOI = hv.Dataset(all_data.genotypes).apply(initialize_ds,
                                all_data=all_data,
                                downsample=downsample_slider,
                                size_column=size_column_select,
                                size_range=size_range_slider,
                                num_barcodes=num_barcodes_slider,
                                NT_muts=NT_muts_text,
                                AA_muts=AA_muts_text,
                                max_groups=max_mut_combos_slider)


## define some widgets for filtering and apply filters

# filter by any numeric columns (muts of interest are classified as numeric sometimes, not sure why)
numerical_columns = [x for x in list(all_data.genotypes.select_dtypes(include=['number']).columns) if 'muts_of_interest' not in x]
filter_by_select = pn.widgets.Select(name='filter column', options=numerical_columns, value='NT_substitutions_count')

filter_range_slider = pn.widgets.RangeSlider(name='NT mutations range',
                                    start=round(all_data.genotypes['NT_substitutions_count'].min()),
                                    end=round(all_data.genotypes['NT_substitutions_count'].max()),
                                    step=1,
                                    value=( round(all_data.genotypes['NT_substitutions_count'].min()),round(all_data.genotypes['NT_substitutions_count'].max()) ))

# update the range of the filter based on which column is selected
def range_widget_callback(RangeSlider, event):
    column = event.new
    col_dtype = all_data.genotypes[column].dtype
    if np.issubdtype(col_dtype, np.integer):
        step = 1
    elif np.issubdtype(col_dtype, np.floating):
        step = 0.1
    minVal = all_data.genotypes[column].min()
    maxVal = all_data.genotypes[column].max()
    RangeSlider.name  = column + ' range'
    RangeSlider.start = minVal
    RangeSlider.end   = maxVal
    RangeSlider.value = (minVal, maxVal)
    RangeSlider.step  = step
    
filter_by_select.link(filter_range_slider, callbacks={'value':range_widget_callback})

count_range_slider = pn.widgets.IntRangeSlider(name='genotype count range',
                                    start=round(all_data.genotypes['count'].min()),
                                    end=round(all_data.genotypes['count'].max()),
                                    step=1,
                                    value=( round(all_data.genotypes['count'].min()),round(all_data.genotypes['count'].max()) ))

# multi select widget that depends on the kind of genotypes CSV file columns
for choice in ['barcode_group','timepoint']:
    if choice in all_data.genotypes.columns:
        choice_col = choice
group_choice = pn.widgets.MultiChoice( name=f'select {choice_col}', value=list(all_data.genotypes[choice_col].unique()),
                                        options=list(all_data.genotypes[choice_col].unique()) )

def filter_dataset(dataset, filter_by, filter_range, count_range, group):
    df = dataset.data
    df = df[
        (df['count'] >= count_range[0]) &
        (df['count'] <= count_range[1])]

    df = df[
        (df[filter_by] >= filter_range[0]) &
        (df[filter_by] <= filter_range[1])]
    
    df = df[df[choice_col].isin(group)]

    return hv.Dataset(df)

filtered = downsampled_MOI.apply(filter_dataset,
                            filter_by=filter_by_select,
                            filter_range=filter_range_slider,
                            count_range=count_range_slider,
                            group=group_choice)


## define widgets for the static points plot, then make a points plot that will be used for tracking selections, but
#   will just show grey empty circles which will remain when a subset are selected. selected points will then
#   be used to further filter the data via a selection stream, and this selected data will be used by other plots

color_options = ['NT_substitutions_count', 'NT_muts_of_interest', choice_col] + [col for col in numerical_columns if col != 'NT_substitutions_count']
if do_AA_analysis:
    color_options.append('AA_muts_of_interest')

if 'barcode(s)' in all_data.genotypes.columns:
    color_options.append('barcode(s)_subset')

embedding_options = [col.split('_')[0] for col in list(all_data.genotypes.columns) if col.endswith('_PaCMAP1')]
points_x_select = pn.widgets.Select(name='column for x axis', options=numerical_columns, value=embedding_options[0] + '_PaCMAP1')
points_y_select = pn.widgets.Select(name='column for y axis', options=numerical_columns, value=embedding_options[0] + '_PaCMAP2')

#TODO: add toggle switch for hover
# hover = HoverTool(tooltips=[('index','@index'),('count','@count'),('NT mutations count','@NT_substitutions_count'),('AA mutations count','@AA_substitutions_nonsynonymous_count'),
#                             ('NT mutations','@NT_substitutions'),('AA mutations','@AA_substitutions_nonsynonymous')])
tools = ['box_select', 'lasso_select']

def points(dataset, x_column, y_column, tools, link=None):
    # dim1, dim2 = f"{embedding}_PaCMAP1", f"{embedding}_PaCMAP2"
    # ticks outside range to hide them
    plot = dataset.data.hvplot(kind='points', x=x_column, y=y_column, hover_cols=color_options ).opts(
        xlabel=x_column, ylabel=y_column, tools=tools, height=600, width=770)
    
    return plot

downsampled_points = downsampled_MOI.apply(points, x_column=points_x_select, y_column=points_y_select, tools=tools)

ls = link_selections.instance()

# filter dataset by selection if any points are selected and clear selection if no points are selected
def select_dataset(dataset, selection_expr):
    idx = dataset.select(selection_expr).data.index
    if len(idx) == 0:
        print('clearing selection')
        ls.selection_expr = None
    # stores the most recently selected indices in the all_data SequenceAnalyzer object and retrieve the indices, which will be the same as before unless nothing is selected
    idx = all_data.select(idx=idx)['df'].index
    return hv.Dataset(dataset.data.loc[idx])

selected = filtered.apply(select_dataset, selection_expr=ls.param.selection_expr)

# make a hv.Table of the selected points, default to randomly sampling 50 sequences
sample_size_slider = pn.widgets.IntSlider(name='table sample size', start=5, step=1,
                                        end=500,
                                        value=50)

## define widgets for the filtered points plot, then use pn.bind to combine this with the unfiltered points plot

color_by_select = pn.widgets.Select(name='color by', options=color_options, value=color_options[0])
colormap_dict = cmap_dict()
cmap_selector = pn.widgets.Select(name='colormap', options=list(colormap_dict.keys()), value='kbc_r')

# manual cmap for datashading
def get_color_key(categories, cmap):
    """ Returns a dictionary mapping supplied categories to colors from a supplied colormap
    """
    all_colors = colormap_dict[cmap]
    num_categories = len(categories)

    # Create an array of evenly spaced indices for selecting colors from the colormap
    indices = np.linspace(0, len(all_colors) - 1, num_categories, dtype=int)
    subset_colors = [all_colors[i] for i in indices]

    cmap_dict = {cat: color for cat, color in zip(categories, subset_colors)}
    return cmap_dict

def float_to_gray(value):
    """ converts a float value between 0 and 1 to a hex color between white and grey
    """
    grey = (128, 128, 128)  # RGB values for #808080 (grey)
    white = (255, 255, 255)  # RGB values for #ffffff (white)

    def interpolate_channel(channel1, channel2, value):
        return int(channel1 + (channel2 - channel1) * value)

    red = interpolate_channel(white[0], grey[0], value)
    green = interpolate_channel(white[1], grey[1], value)
    blue = interpolate_channel(white[2], grey[2], value)

    return matplotlib.colors.rgb2hex((red / 255, green / 255, blue / 255))


# alpha for selected points, disabled if datashader is on
selected_alpha_slider = pn.widgets.FloatSlider(name='selected point opacity',
                                    start=0.01,
                                    end=1.00,
                                    step=0.01,
                                    value=1.0,
                                    disabled=ds_checkbox.value )
ds_checkbox.link(selected_alpha_slider, callbacks={'value':ds_widget_callback})

# alpha for unselected points, used to determine shade of grey if datashader is on
unselected_alpha_slider = pn.widgets.FloatSlider(name='unselected point opacity',
                                    start=0.00,
                                    end=1.00,
                                    step=0.01,
                                    value=0.1 )


# I have to use panel.bind for the dynamic scatter plot because I want to be able to switch between numerical and categorical color labelling, but there's a bug when using apply to do that. see https://github.com/holoviz/holoviews/issues/5591
def points_bind(link, x_col, y_col, filtered_dataset, points_unfiltered, colorby, cmap, selected_alpha, unselected_alpha, do_datashade, NT_muts, AA_muts, max_muts):
    # widgets that aren't used in this function are just there to trigger the callback

    if do_datashade:
        points_unfiltered = datashade(points_unfiltered, aggregator=ds.reductions.any(), cmap=[float_to_gray(unselected_alpha)])
        points_unfiltered = dynspread(points_unfiltered, threshold=1, max_px=2).opts(height=600, width=600)

        clims = None
        cnorm = 'eq_hist'
        min_alpha = 40
        color_key = None
        # decide how to aggregate based on if the colorby column is categorical or numerical
        is_not_numeric = not pd.api.types.is_numeric_dtype(all_data.genotypes[colorby].dtype)
        if is_not_numeric: # highest per category count within a pixel
            agg = ds.by(colorby)
            color_key = get_color_key(all_data.get_selection()['df'][colorby].unique(), cmap)
            min_alpha = 255 # prevents coloring with non-discrete colors
        elif colorby == 'count':
            agg = ds.reductions.sum(colorby)
        else:          # average value within a pixel
            agg = ds.reductions.mean(colorby)
            clims = all_data.genotypes[colorby].min(), all_data.genotypes[colorby].max()
            cnorm = 'linear'

        points_filtered = filtered_dataset.apply(points, x_column=x_col, y_column=y_col, tools=[])
        points_filtered = datashade(points_filtered, aggregator=agg, cmap=cmap, color_key=color_key, cnorm=cnorm, clims=clims, min_alpha=min_alpha)
        if '_PaCMAP' in x_col:
            points_filtered = points_filtered.opts(xticks=[100])
        if '_PaCMAP' in y_col:
            points_filtered = points_filtered.opts(yticks=[100])
        points_filtered = dynspread(points_filtered, threshold=1, max_px=2).opts(height=600, width=600)
        plot = points_unfiltered * points_filtered

    else:
        points_unfiltered.opts(fill_alpha=0, line_alpha=unselected_alpha, nonselection_line_alpha=unselected_alpha, color='grey')
        points_filtered = filtered_dataset.apply(points, x_column=x_col, y_column=y_col, tools=[]).opts(colorbar=True).apply.opts(alpha=selected_alpha, nonselection_alpha=unselected_alpha, color=colorby, cmap=cmap, colorbar_opts={'title':colorby})
        if '_PaCMAP' in x_col:
            points_filtered = points_filtered.opts(xticks=[100])
        if '_PaCMAP' in y_col:
            points_filtered = points_filtered.opts(yticks=[100])
        plot = points_unfiltered * points_filtered
    
    if unselected_alpha == 0: # if unselected points are invisible, don't plot them to save compute. to maintain selections, link_selections is applied to the filtered points
        plot = link(points_filtered)
    else:
        link(points_unfiltered)
    return plot

dynamic_points = pn.bind(points_bind,
                                link=ls,
                                x_col=points_x_select,
                                y_col=points_y_select,
                                points_unfiltered=downsampled_points,
                                filtered_dataset=selected,
                                colorby=color_by_select,
                                cmap=cmap_selector,
                                unselected_alpha=unselected_alpha_slider,
                                selected_alpha=selected_alpha_slider,
                                do_datashade=ds_checkbox,
                                NT_muts=NT_muts_text,
                                AA_muts=AA_muts_text,
                                max_muts=max_mut_combos_slider)


# datashader requires some help to make the colorbar for numerical data
def colorbar_bind(colorby, cmap, do_datashade):

    color_bar = False
    if ( colorby in all_data.get_selection()['df'].columns) & do_datashade:
        is_not_numeric = not pd.api.types.is_numeric_dtype(all_data.genotypes[colorby].dtype)
        if is_not_numeric: # no colorbar for categorical data
            color_bar = False
        elif colorby == 'count': # unclear how to set up colorbar for count data
            color_bar = False
        else: 
            color_mapper = LinearColorMapper(palette=colormap_dict[cmap], low=all_data.genotypes[colorby].min(), high=all_data.genotypes[colorby].max())
            color_bar = ColorBar(color_mapper=color_mapper)

            # Wrap the color_bar in a Bokeh Figure so that Panel can render it
            p = figure(width=100, height=550, toolbar_location=None, min_border=0)
            p.add_layout(color_bar, 'right')
            p.axis.visible = False
            p.xgrid.grid_line_color = None
            p.ygrid.grid_line_color = None
            p.outline_line_color = None

            label = Label(
                text=colorby,
                x=10,
                y=200,
                x_units="screen",
                y_units="screen",
                angle=0.5 * np.pi,
                text_align="center",
                text_baseline="middle",
            )
            p.add_layout(label)

            color_bar = pn.panel(p)
    
    if not color_bar:
        color_bar = pn.Spacer(width=0)

    return color_bar

colorbar = pn.bind(colorbar_bind, colorby=color_by_select, cmap=cmap_selector, do_datashade=ds_checkbox)



## histogram for a column chosen by a widget. selection will be highlighted by combining a low opacity histogram
#   of the pre-filter dataset with a high opacity histogram of the selected dataset

hist_col_selector = pn.widgets.Select(name='histogram column', options=numerical_columns, value='NT_substitutions_count')
hist_bins_slider = pn.widgets.IntSlider(name='histogram bins', start=1, step=1,
                                        end=round(all_data.genotypes['NT_substitutions_count'].max()),
                                        value=round(all_data.genotypes['NT_substitutions_count'].max()))

# update the range of bins based on the maximum value of the column
def bins_widget_callback(IntSlider, event):
    column = event.new
    if all_data.genotypes[column].dtype == 'int64':
        maxVal = round(all_data.genotypes[column].max())
        IntSlider.end = maxVal
        IntSlider.value = maxVal
    elif all_data.genotypes[column].dtype == 'float64':
        IntSlider.end = 30
        IntSlider.value = 10
    
hist_col_selector.link(hist_bins_slider, callbacks={'value':bins_widget_callback})

def histogram(dataset, histCol, num_bins):
    # calculate bin range based on maximum value, allowing for bins=max value for integer columns
    if all_data.genotypes[histCol].dtype == 'int64':
        max_x_val = all_data.genotypes[histCol].max()
        min_x_val = 0
    elif all_data.genotypes[histCol].dtype == 'float64':
        max_x_val = all_data.genotypes[histCol].max()
        min_x_val = all_data.genotypes[histCol].min()
    bins = np.linspace(min_x_val-0.5, max_x_val-0.5, num_bins+1)
    return dataset.data.hvplot(y=histCol, kind='hist', width=500,height=500, bins=bins, xlabel=histCol, ylabel='total sequences', color='grey')

dynamic_hist = selected.apply(histogram, histCol=hist_col_selector, num_bins=hist_bins_slider).opts(alpha=1)
static_hist = downsampled_MOI.apply(histogram, histCol=hist_col_selector, num_bins=hist_bins_slider).opts(alpha=0.1)



## plot that describes mutations aggregated over the selected genotypes
# Current options:
#   - bar plot of most frequent mutations
#   - bar plot of least frequent mutations
#   - heatmap of most frequent mutations
#   - heatmap of least frequent mutations
agg_muts_width_slider = pn.widgets.IntSlider(name='plot width',
                                start=160,
                                end=1600,
                                step=10,
                                value=500 )

max_positions = all_data.integer_matrix['NT'].shape[1]
NTorAA_radio_options = ['nucleotide (NT)']
if do_AA_analysis:
    NTorAA_radio_options.append('amino acid (AA)')
    max_positions = all_data.integer_matrix['AA'].shape[1]
NTorAA_radio = pn.widgets.RadioButtonGroup(
    name='bar plot mutation type', options=NTorAA_radio_options, value=NTorAA_radio_options[-1], disabled=(not do_AA_analysis), button_type='default')
agg_plot_type_selector = pn.widgets.Select(name='aggregated mutations plot type', options=['most frequent, heatmap', 'least frequent, heatmap', 'most frequent, bars', 'least frequent, bars'])
num_positions_slider = pn.widgets.IntSlider(name='number of positions',
                                    start=5,
                                    end=max_positions,
                                    step=1,
                                    value=10 )
# link the number of positions slider to the NTorAA radio button
# update the range of bins based on the maximum value of the column
def num_positions_callback(IntSlider, event):
    value = event.new
    NTorAA = value[-3:-1]
    max_positions = all_data.integer_matrix[NTorAA].shape[1]
    IntSlider.end = max_positions
    
NTorAA_radio.link(num_positions_slider, callbacks={'value':num_positions_callback})

def plot_mutations_aggregated(filtered_dataset, all_data, num_positions, NTorAA, plot_type, cmap):
    """
    produces plots that describe aggregated mutations. plot is chosen based on the text in the plot_type input. link is just used to trigger the callback when a selection is made
    """

    # apply is faster but can't deal with both heatmap and bars, so we wrap within a pn.bind function
    def apply_plot(dataset, all_data, num_positions, NTorAA, most_common, heatmap, cmap):
        selected_idx = dataset.data.index
        # selected_idx = all_data.selected_idx
        total_seqs = all_data.get_count(selected_idx)
        df = all_data.aggregate_mutations(
                                    NTorAA=NTorAA,
                                    idx=selected_idx)

        plot = conspicuous_mutations(df, total_seqs, num_positions=num_positions,
                                    colormap=cmap, most_common=most_common, heatmap=heatmap)
        if heatmap: # normally the plot would have n in the clabel but this doesn't get updated after a selection
            plot.opts(title=str(f'n={total_seqs}'), clabel='frequency', fontsize={'title':10})
        return plot
    
    NTorAA = NTorAA[-3:-1] # get NT or AA from the text in the radio button
    if plot_type.endswith('heatmap'):
        heatmap = True
    elif plot_type.endswith('bars'):
        heatmap = False
        cmap = colormaps[NTorAA]
    if plot_type.startswith('most'):
        most_common=True
    elif plot_type.startswith('least'):
        most_common=False
    
    plot = filtered_dataset.apply(apply_plot, all_data=all_data, num_positions=num_positions, NTorAA=NTorAA, most_common=most_common, heatmap=heatmap, cmap=cmap)

    return plot

aggregated_muts_plot = pn.bind(plot_mutations_aggregated,
                                        filtered_dataset=selected,
                                        all_data=all_data,
                                        num_positions=num_positions_slider,
                                        NTorAA=NTorAA_radio,
                                        plot_type=agg_plot_type_selector,
                                        cmap=cmap_selector)

# Create a panel with a slider to control the width of the plot
aggregated_muts_panel = pn.panel(aggregated_muts_plot, width=agg_muts_width_slider.value)
def mod_width_callback(muts_panel, event):
    muts_panel.width=event.new
agg_muts_width_slider.link(aggregated_muts_panel, callbacks={'value':mod_width_callback})

# make a text input widget for naming exported files
selection_name_text = pn.widgets.TextInput(name='selection name', value='unnamed')

# make a button that reruns the aggregate_mutations function and exports the resulting dataframe to a csv file
export_agg_muts_button = pn.widgets.Button(name='export aggregated mutations to .CSV', button_type='primary')
def export_agg_muts(event):
    # get date/time for filename to prevent overwriting
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

    # get the most recently selected indices
    selected_idx = all_data.selected_idx

    # get the aggregated mutations dataframe
    NTorAA = NTorAA_radio.value[-3:-1]
    df = all_data.aggregate_mutations(NTorAA=NTorAA, idx=selected_idx)

    # export to csv
    df.to_csv(f'dashboard/dashboard-{selection_name_text.value}-{timestamp}_mutations-aggregated.csv', index=False)
export_agg_muts_button.on_click(export_agg_muts)

# make a table to display the selected genotypes
def tabulate(dataset, sample_size):
    df = dataset.data.drop(['barcode(s)_subset','size'], axis=1)
    if sample_size < len(df):
        df = df.sample(n=sample_size, random_state=0)
    if df.empty:
        table = hv.Table(pd.DataFrame())
    else:
        table = hv.Table(df)
    table.opts(width=1000, height=450)
    return table

selected_table = selected.apply(tabulate, sample_size=sample_size_slider)

# add the legend for the bar plots only if the bar plot is being displayed
def get_bar_legend(agg_plot_type):
    if agg_plot_type.endswith('heatmap'):
        return pn.Spacer(width=0)
    elif agg_plot_type.endswith('bars'):
        snakemake_dir = pathlib.Path(__file__).parent.parent.parent
        return pn.panel(pathlib.Path(snakemake_dir/'images'/'AA_NT_colormap.png'),width=200,align='start')
bars_legend = pn.bind(get_bar_legend, agg_plot_type=agg_plot_type_selector)

# Add a button that exports all three plots to SVG format
SVG_button = pn.widgets.Button(name='export plots to .SVG', button_type='primary')
def export_svg(event):
    # get date/time for filename
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

    points = dynamic_points()

    # get current plots
    export_svg_plots([points, aggregated_muts_plot().opts(width=agg_muts_width_slider.value), dynamic_hist.last*static_hist.last], f'dashboard/dashboard-plot_{selection_name_text.value}_{timestamp}.html', ['points', 'mutations-frequencies', 'histogram']) # "html" gets clipped and replaced with "SVG"
    print('plots exported as SVGs to dashboard folder')
SVG_button.on_click(export_svg)

# add button that exports the selected genotypes to a csv file
genotypes_CSV_button = pn.widgets.Button(name='export selected genotypes to .csv', button_type='primary')
def export_csv(event):
    # get date/time for filename
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    filename = f'dashboard/dashboard-{selection_name_text.value}-{timestamp}_genotypes.csv'
    selected_genotypes = all_data.get_selection()['df'].drop(['barcode(s)_subset','size'], axis=1)
    pathlib.Path(filename).parent.absolute().mkdir(parents=True, exist_ok=True)
    selected_genotypes.to_csv(filename)
    print(f'selected genotypes exported to {filename}')
genotypes_CSV_button.on_click(export_csv)

# add buttons that write sequences corresponding to selected indices to a fasta file
NT_fasta_button = pn.widgets.Button(name='export selected NT sequences to .fasta', button_type='primary')
def export_NT_fasta(event):
    # get date/time for filename
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    filename = f'dashboard/dashboard-{selection_name_text.value}-{timestamp}_NT-seqs.fasta'
    all_data.write_fasta(filename, 'NT', idx=all_data.selected_idx)
    print(f'wrote {len(all_data.selected_idx)} NT sequences to {filename}')
NT_fasta_button.on_click(export_NT_fasta)

AA_fasta_button = pn.widgets.Button(name='export selected AA sequences to .fasta', button_type='primary', disabled=(not do_AA_analysis))
def export_AA_fasta(event):
    # get date/time for filename
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    filename = f'dashboard/dashboard-{selection_name_text.value}-{timestamp}_AA-seqs.fasta'
    all_data.write_fasta(filename, 'AA', idx=all_data.selected_idx)
    print(f'wrote {len(all_data.selected_idx)} AA sequences to {filename}')
AA_fasta_button.on_click(export_AA_fasta)

# add buttons that aggregates mutations in selected sequences and uses them to produce a consensus sequence, which is then written to a fasta file
NT_consensus_button = pn.widgets.Button(name='export NT consensus of selected sequences to .fasta', button_type='primary')
def export_NT_consensus(event):
    # get date/time for filename
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    filename = 'dashboard/dashboard-NT-consensus.fasta'
    all_data.get_consensus('NT', idx=all_data.selected_idx, write_to=filename, append=True, name=f'{selection_name_text.value}_{timestamp}')
    print(f'appended 1 NT consensus sequence to {filename}')
NT_consensus_button.on_click(export_NT_consensus)

AA_consensus_button = pn.widgets.Button(name='export AA consensus of selected sequences to .fasta', button_type='primary', disabled=(not do_AA_analysis))
def export_AA_consensus(event):
    # get date/time for filename
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    filename = 'dashboard/dashboard-AA-consensus.fasta'
    all_data.get_consensus('AA', idx=all_data.selected_idx, write_to=filename, append=True, name=f'{selection_name_text.value}_{timestamp}')
    print(f'appended 1 AA consensus sequence to {filename}')
AA_consensus_button.on_click(export_AA_consensus)

# make a button to reset all the widgets to their default values
widgets = [ds_checkbox, size_column_select, size_range_slider, downsample_slider, NT_muts_text, AA_muts_text, max_mut_combos_slider, filter_by_select, filter_range_slider, count_range_slider, group_choice, points_x_select, points_y_select, selected_alpha_slider, unselected_alpha_slider, color_by_select, cmap_selector, NTorAA_radio, agg_plot_type_selector, sample_size_slider, agg_muts_width_slider]
widgets_default_values = {widget:widget.value for widget in widgets}
reset_button = pn.widgets.Button(name='reset all widgets', button_type='primary')
def reset(event):
    for widget, value in widgets_default_values.items():
        widget.value = value
reset_button.on_click(reset)

## build the layout
layout = pn.Column(
    pn.Row(
    pn.Column(
        pn.Row(downsample_slider, ds_checkbox),
        pn.Row(selected_alpha_slider,unselected_alpha_slider),
        pn.Row(size_column_select,size_range_slider),
        pn.Row(pn.Column(points_x_select,points_y_select,count_range_slider),group_choice),
        pn.Row(filter_by_select,filter_range_slider),
        pn.Row(color_by_select,cmap_selector),
        pn.Row(NT_muts_text,AA_muts_text),
        pn.Row(SVG_button, max_mut_combos_slider),
        pn.Row(selection_name_text, num_barcodes_slider),
        pn.Row(reset_button)),
    pn.Row(dynamic_points, colorbar),
    pn.Column(sample_size_slider,
              selected_table,
              genotypes_CSV_button,
              pn.Row(NT_fasta_button, AA_fasta_button),
              pn.Row(NT_consensus_button, AA_consensus_button))
    ),
pn.Row(
    pn.Column(
        dynamic_hist*static_hist,
        pn.Row(hist_col_selector, hist_bins_slider)),
    pn.Column(
        aggregated_muts_panel, NTorAA_radio,
        agg_plot_type_selector,
        num_positions_slider,
        agg_muts_width_slider,
        export_agg_muts_button),
    bars_legend
    ))

layout.servable(title=f'maple dashboard, file: {args.genotypes.split("/")[-1]}')