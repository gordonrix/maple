#
#  DESCRIPTION   : Script for maple pipeline. Declares functions and classes that are used
#                   in several scripts throughout the pipeline
#
#  AUTHOR(S)     : Gordon Rix
#

import pandas as pd
import numpy as np
import pathlib
import holoviews as hv

from bokeh.io import export_svgs
from selenium import webdriver as wd
from selenium.webdriver.chrome.service import Service
from bokeh import palettes
from colorcet import palette
from natsort import natsorted

# define colormaps
colormaps = {'NT':{'A':palettes.Greens[3][1], #take the middle color from the 3 length color list
                                'T':palettes.Reds[3][1],
                                'G':'#000000',           #black
                                'C':palettes.Blues[3][1],
                                '-':'#d3d3d3'}}           #grey

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

colormaps.update({'AA':amino_acid_colormap})

def sort_barcodes(bcList, bcGroupsDict):
    """
    bcList:         list of strings, barcode groups
    bcGroupsList:   dict, barcode groups dictionary for a specific tag from the config file,
                        or an empty dictionary

    returns the same list of barcode sorted first according to user provided barcode groups,
        then by natural sorting (see natsort package)
    """
    defined = [bc for bc in bcGroupsDict.keys() if bc in bcList]
    undefined = natsorted([bc for bc in bcList if bc not in bcGroupsDict.keys()])
    return defined+undefined

def cmap_dict():
    """
    basically just a wrapper for the colorcet palette function that returns a dictionary
    (allowing for key retrieval) with a limited number of color palettes that includes reversed color palettes
    """

    cmaps_fwd = ['kbc', 'fire', 'bgy', 'bgyw', 'bmy', 'gray', 'rainbow4']

    # need a dict to map color names to lists of hex color values
    colormap_dict = {}

    # add reversed colormaps to cmap dict
    for cmap in cmaps_fwd:
        cmap_rvs = cmap+'_r'
        colors = palette[cmap]
        colormap_dict[cmap] = colors
        colormap_dict[cmap_rvs] = colors[::-1]

    return colormap_dict

def get_colors(labels, cmap, background=False):
    """
    retrieve a list of colors for each of the samples in the input list

    args:
        labels (list(str)):    list of labels, used to determine background sample index
        cmap (str):           string that describes the colormap to use
        background (str):     string for the background sample if it is being used,
                                will assign grey to the same position in the list as in the labels, if it is present
                                must only appear once in the list of labels

    returns:
        colors:         a list of hex color values of length = len(labels)
    """
    labels = labels.copy() # don't modify the input list

    if background:
        if background in labels: # set aside background sample while we determine colors
            background_idx = labels.index(background)
            labels.remove(background)
        else:
            background = False

    colormap = cmap_dict()[cmap]
    steps = len(labels)-1 # number of steps between colors
    colors = []
    for i in range(0,len(labels)):
        if steps > 0:
            color_idx = int(round(i * (len(colormap) - 1) / steps))
        else:
            color_idx = len(colormap)//2 # use a color in the middle of the colormap
        colors.append(colormap[color_idx])
        
    # insert grey into list at position of the background sample
    if background:
        colors.insert(background_idx, 'grey')

    return colors


def dist_to_DF(dist, x, y):
    """given a np.bincount output, i.e. a distribution of values in which each value is
    the number of observations corresponding to the value's position in the distribution
    (e.g. [0,3,1] is a distribution for 4 values in which 3 of the values were '1'),
    this will calculate the distribution as a proportion and as a cumulative proportion
    and produce a DataFrame from these three numpy arrays that describe the distribution
    
    dist:       np.array of shape (maximumValue)
    x, y:       strings that describes the x and y variable"""

    maxVal = dist.shape[0]
    dist = dist.reshape(1,maxVal)

    values = np.arange(maxVal).reshape(1,maxVal)
    proportion = np.divide(dist, dist.sum())
    cumsum = np.cumsum(proportion).reshape(1,maxVal)

    df = pd.DataFrame( np.concatenate((values, dist, proportion, cumsum), axis=0).T, columns = [x, f"total {y}", f"proportion of {y}", f"cumulative proportion of {y}"] )
    df[df.columns[[0,1]]] = df[df.columns[[0,1]]].apply(lambda x: x.astype(int)) # convert x value and total counts columns to int

    return df

def export_svg_plots(plots, file_name, labels=[], export=True):
    """
    attempts to export individual bokeh plots from a list of holoviews plots
    with the provided file name and the index of the plot in the list.
    If it fails, it will print a warning and return without exporting any plots or raising an error.
    
    plots:      list of holoviews bokeh plots
    file_name   file name being used to save the plots. must end in '.html'
    labels      list of strings to be appended to the end of the file name for 
                    each of the plots, which are exported as separate files
    export:     bool or string. if True, export all plots. if False, export none.
                    if string, export only plots that contain the string in the file name or label
    """

    # unless labels are provided, name plots just using their order in the plot list
    if not labels:
        labels = [str(i) for i in range(0,len(plots))]
    file_name_base = file_name[:-5] # remove '.html' from file name

    pathlib.Path(file_name_base).parent.absolute().mkdir(parents=True, exist_ok=True)

    # Unless export is just True, only export plots that contain the export string in the file name or label
    if type(export) == bool:
        if not export:
            return
    elif type(export) == str:
        if export not in file_name:
            labels = [label for label in labels if export in label]
            plots = [plots[i] for i,label in enumerate(labels) if export in label]

    if plots:

        try:

            options = wd.ChromeOptions()
            options.add_argument('--headless')
            options.add_argument('--disable-gpu')
            options.add_argument("--no-sandbox")
            options.add_argument("--window-size=2000x2000")
            options.add_argument('--disable-dev-shm-usage')
            webdriver = wd.Chrome(service=Service(), options=options)
            
            for plot, label in zip(plots, labels):
                fName = f'{file_name_base}_{label}.svg'
                p = hv.render(plot,backend='bokeh')
                p.output_backend='svg'
                export_svgs(p, 
                    filename=fName,
                    webdriver=webdriver)
                
        except:
            print(f"""\n[ERROR] SVG export for {file_name} failed. Pipeline continuing but the SVG version of this plot was not generated. 
                    This usually doesn't indicate a consistent issue. Try again by deleting the associated .html file and rerunning snakemake.\n""")
    return
        
def conspicuous_mutations(df, total_seqs, num_positions=None, colormap='kbc_r', most_common=True, heatmap=False):
    """
    produces a bar plot of the most or least frequent mutations
    
    parameters:
        df (pd.DataFrame):   dataframe of aggregated mutations output by aggregate_mutations
        total_seqs (int):    total number of sequences in the sample, just used for labelling the plot
        num_positions (int): number of mutations to include in the bar plot output
        colormap (dict):     AA/NT letter : color hex code key:value pairs to use for the stacked bars plot output
                                or a name of the colormap to use for the heatmap output
        most_common (bool):  if True/False, output the most/least commonly mutated positions
        
    returns:
        hv.Bars or hv.heatmap object showing the topN most frequently observed mutations
            in the aggregated mutations dataframe
    """
    if num_positions is None:
        num_positions = len(df['position'].unique())

    df = df.sort_values(['total_count','position'], ascending=[(not most_common),True])
    df_grouped = df.groupby('position', as_index=False).sum().sort_values(['total_count','position'],ascending=[not most_common, True])
    positions = df_grouped['position'].iloc[:num_positions]
    df = df[df['position'].isin(positions)]
    df = df.sort_values('position', ascending=True)
    df['WT_position'] = df['wt'] + df['position'].astype(str)

    if heatmap:
        AAs = list('AILPVFWYNQSTCMDEHKRG*')
        order = AAs + [m for m in df['mutation'].unique().tolist() if m not in AAs]
        df['mutation'] = pd.Categorical(df['mutation'], categories=order, ordered=True)
        df = df.sort_values(['position','mutation'], ascending=[True,False])

        plot = hv.HeatMap(df, kdims=['WT_position','mutation'], vdims=['proportion_of_seqs', 'total_count']).opts(
                    colorbar=True, clabel=f"frequency (n={total_seqs})", ylabel="mutation")
    else:
        plot = hv.Bars(df, kdims=['WT_position','mutation'], vdims=['proportion_of_seqs', 'total_count']).opts(
                    show_legend=False, ylabel=f"frequency (n={total_seqs})", stacked=True)
                    
    plot = plot.opts(height=500, width=1000, xrotation=90, tools=['hover'],  cmap=colormap, xlabel='position', fontsize={'title':16,'labels':14,'xticks':10,'yticks':10})
    return plot

def dashboard_input(wildcards, config):
    """
    returns the input files for the dashboard based on the sample name provided by the user
    or the first run tag in the config file if no sample name is provided

    args:
        wildcards (snakemake object):   snakemake wildcards object, not used but required for unpacking in rule input
        config (dict):                  dictionary of config file contents

    returns:
        inputDict (dict):               dictionary of input files for the dashboard, genotypes and refFasta
    """
    sample = config.get('dashboard_input', False)
    # use the first run tag as the sample for the dashboard if no sample is provided by user
    if not sample:
        sample = config['runs'].values()[0]
    if sample in config['runs']:
        genotypes = f'mutation_data/{sample}/{sample}_genotypes-reduced-dimensions.csv'
        if 'enrichment' in config['runs'][sample]:
            genotypes = genotypes[:-4] + '-enrichment.csv'
        inputDict = {'genotypes': genotypes, 'refFasta': config['runs'][sample]['reference']}
    elif sample in config.get('timepointsInfo', ''):
        inputDict = {'genotypes': f'mutation_data/timepoints/{sample}_merged-timepoint_genotypes-reduced-dimensions.csv',
                    'refFasta': config['timepointsInfo'][sample]['reference']}
    elif '_' in sample and (sample.split('_')[0] in config['runs']):
        print(f'[NOTICE] dashboard_input {sample} contains an underscore, so tag_barcode input is assumed. This may cause an error if this is not a valid tag_barcode combination.')
        tag, barcodes = sample.split('_')
        inputDict = {'genotypes': f'mutation_data/{tag}/{barcodes}/{tag}_{barcodes}_genotypes-reduced-dimensions.csv',
                    'refFasta': config['runs'][tag]['reference']}
    else:
        print(f'[NOTICE] dashboard_input {sample} does not conform to a valid tag, timepoint name, or tag_barcode combination. Dashboard input will not be generated.')
        return None
    return inputDict