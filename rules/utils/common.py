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
from bokeh.models import HoverTool
from selenium import webdriver as wd
from selenium.webdriver.chrome.service import Service
from bokeh import palettes
from colorcet import palette
import re
import os
import glob
import xml.etree.ElementTree as ET
from natsort import natsorted
import sys

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

WILDCARD_CONSTRAINTS = {
    "tag": "[^\/_]+",
    "barcodes": "[^\/_]+",
    "NTorAA": "[^\/_-]+",
    "timepointsGroup": "[^\/_]+",
}

def get_demuxed_barcodes(tag, bcGroupsDict):
    """
    tag:            string, run tag
    bcGroupsDict:   dictionary, barcodeGroups dict from config file for a specific tag,
                        if it exists, otherwise an empty dict if no sorting is needed

    grabs all barcodes for a given tag that were properly demultiplexed
            if that tag was demultiplexed, then sorts according to the
            barcode groups in the config file"""

    if config['do_demux'][tag]:
        checkpoint_demux_output = checkpoints.demultiplex.get(tag=tag).output[0]
        checkpoint_demux_prefix = checkpoint_demux_output.split('demultiplex')[0]
        checkpoint_demux_files = checkpoint_demux_prefix.replace('.','') + '{BCs}.bam'
        BCs = glob_wildcards(checkpoint_demux_files).BCs
        out = sort_barcodes(BCs, bcGroupsDict)
    else:
        out = ['all']
    return out

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

def export_svg_plots(plots, file_name, labels=[], export=True, dimensions=None):
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
    dimensions: tuple of ints (width, height), dimensions to use for the plots. if None, use default dimensions
    """

    # unless labels are provided, name plots just using their order in the plot list
    if not labels:
        labels = [str(i) for i in range(0,len(plots))]
    out_dir = file_name[:-5]+'_svg' # remove '.html' from file name and add '_svg'

    # Unless export is just True, only export plots that contain the export string in the file name or label
    if type(export) == bool:
        if not export:
            return
    elif type(export) == str:
        if export not in file_name:
            labels = [label for label in labels if export in label]
            plots = [plots[i] for i,label in enumerate(labels) if export in label]

    # Function to normalize SVG units to standardize the SVG output across different programs
    # Ensures that svg outputs look the same as html rendering
    def normalize_svg_units(fName):
        tree = ET.parse(fName)
        root = tree.getroot()
        NS = {'svg': 'http://www.w3.org/2000/svg'}

        # make width/height explicit px
        w = root.attrib.get('width')
        h = root.attrib.get('height')
        if w and not w.endswith('px'):
            root.set('width',  w + 'px')
        if h and not h.endswith('px'):
            root.set('height', h + 'px')

        # esure viewBox matches pixel dimensions
        if 'viewBox' not in root.attrib and w and h:
            wi = float(w.rstrip('px'))
            hi = float(h.rstrip('px'))
            root.set('viewBox', f"0 0 {int(wi)} {int(hi)}")

        # optionally convert pt→px (1 pt ≈ 1.333px at 96 ppi)
        for text in root.findall('.//svg:text', NS) + root.findall('.//svg:tspan', NS):
            fs = text.attrib.get('font-size') or text.attrib.get('style','')
            # handle style="font-size:16pt;…"
            if 'pt' in fs:
                import re
                def repl(m):
                    pt = float(m.group(1))
                    px = pt * 96.0/72.0
                    return f"{px:.1f}px"
                newstyle = re.sub(r'(\d+(\.\d+)?)pt', repl, fs)
                if 'font-size' in text.attrib:
                    text.set('font-size', newstyle)
                else:
                    text.set('style', newstyle)

        tree.write(fName, encoding='utf-8', xml_declaration=True)


    if plots:

        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        try:

            options = wd.ChromeOptions()
            options.add_argument('--headless')
            options.add_argument('--disable-gpu')
            options.add_argument("--no-sandbox")
            options.add_argument("--window-size=2000x2000")
            options.add_argument('--disable-dev-shm-usage')
            webdriver = wd.Chrome(service=Service(), options=options)
            
            for plot, label in zip(plots, labels):
                fName = os.path.join(out_dir, f'{label}.svg')
                plot = plot.opts(show_legend=False) # legends break svg export with bokeh
                p = hv.render(plot,backend='bokeh')
                if dimensions:
                    p.width = dimensions[0]
                    p.height = dimensions[1]
                p.output_backend='svg'
                export_svgs(p, 
                    filename=fName,
                    webdriver=webdriver)
                
                # after export_svgs(…) has written fName:
                tree = ET.parse(fName)
                root = tree.getroot()

                # handle SVG namespace
                ns = {'svg': 'http://www.w3.org/2000/svg'}

                # bokeh svg export draws text twice in different ways, so remove stroke attributes only on text or tspan elements
                for tag in ('text', 'tspan'):
                    for elem in root.findall(f'.//svg:{tag}', ns):
                        for attr in ['stroke', 'stroke-width', 'stroke-linecap', 'stroke-linejoin']:
                            if attr in elem.attrib:
                                del elem.attrib[attr]

                # write back, preserving the XML prolog
                tree.write(fName, encoding='utf-8', xml_declaration=True)

                normalize_svg_units(fName)
                
        except:
            print(f"""\n[ERROR] SVG export for {file_name} failed. Pipeline continuing but the SVG version of this plot was not generated. 
                    This usually doesn't indicate a consistent issue. Try again by deleting the associated .html file and rerunning snakemake.\n""")
    return
        
def parse_position_ranges(range_str):
    """
    Parse position range string like '20-30' or '20-30,50-60' into a list of positions

    Args:
        range_str (str): Range string like '20-30' or '20-30,50-60'

    Returns:
        list: List of positions to include
    """
    positions = []
    if not range_str.strip():
        return positions

    ranges = range_str.split(',')
    for range_part in ranges:
        range_part = range_part.strip()
        if '-' in range_part:
            try:
                start, end = range_part.split('-')
                start, end = int(start.strip()), int(end.strip())
                positions.extend(range(start, end + 1))
            except ValueError:
                print(f"[WARNING] Invalid range '{range_part}' in position_ranges '{range_str}'. Skipping this range.")
                continue  # Skip invalid ranges
        else:
            try:
                positions.append(int(range_part.strip()))
            except ValueError:
                print(f"[WARNING] Invalid position '{range_part}' in position_ranges '{range_str}'. Skipping this position.")
                continue  # Skip invalid single positions

    return sorted(list(set(positions)))  # Remove duplicates and sort

def conspicuous_mutations(df, total_seqs, num_positions=None, colormap='kbc_r', most_common=True, heatmap=False, fontsize={'title':18,'labels':18,'xticks':16,'yticks':16, 'legend': 16}, position_ranges=None, axis_type='categorical'):
    """
    produces a plot of the most or least frequent mutations

    parameters:
        df (pd.DataFrame):   dataframe of aggregated mutations output by aggregate_mutations
        total_seqs (int):    total number of sequences in the sample, just used for labelling the plot
        num_positions (int): number of mutations to include in the bar plot output
        colormap (dict):     AA/NT letter : color hex code key:value pairs to use for the stacked bars plot output
                                or a name of the colormap to use for the heatmap output
        most_common (bool):  if True/False, output the most/least commonly mutated positions
        position_ranges (str): optional position range string like '20-30' or '20-30,50-60' to filter specific positions
        axis_type (str): 'numerical' or 'categorical' - determines x-axis type

    returns:
        hv.Bars or hv.heatmap object showing the topN most frequently observed mutations
            in the aggregated mutations dataframe
    """
    # Always do the standard data preparation (sorting and grouping)
    df = df.sort_values(['total_count','position'], ascending=[(not most_common),True])
    df_grouped = df.groupby('position', as_index=False).sum().sort_values(['total_count','position'],ascending=[not most_common, True])

    # Determine which positions to keep based on filtering mode
    if position_ranges:
        # Filter by specific position ranges
        positions_to_keep = parse_position_ranges(position_ranges)
        positions = df_grouped[df_grouped['position'].isin(positions_to_keep)]['position'] if positions_to_keep else df_grouped['position']
    elif num_positions is None:
        # Show all positions
        positions = df_grouped['position']
    else:
        # Show top N most/least common positions
        positions = df_grouped['position'].iloc[:num_positions]

    # Apply the position filter
    df = df[df['position'].isin(positions)]

    # Set axis type based on user preference
    if axis_type == 'numerical':
        xrotation = 0
        kdims = ['position','mutation'] # Use numerical positions
        vdims = ['proportion_of_seqs', 'total_count', 'WT_position']
    else:  # categorical
        xrotation = 90
        kdims = ['WT_position','mutation'] # Use WT_position labels
        vdims = ['proportion_of_seqs', 'total_count']

    df = df.sort_values('position', ascending=True)
    df['WT_position'] = df['wt'] + df['position'].astype(str)

    df = df.drop_duplicates(subset=['position', 'mutation', 'wt'])

    # Determine if this is NT or AA data based on wildtype alphabet size
    is_AA = len(df['wt'].unique()) > 5

    # Filter out ambiguous characters and deletions from plot
    if is_AA:
        df = df[~df['mutation'].isin(['X', '-'])]
    else:
        df = df[~df['mutation'].isin(['N', '-'])]

    df.loc[df['mutation'] == df['wt'], 'proportion_of_seqs'] = -.000001 # set the proportion of the WT to -.000001 to make it white in the plot

    if heatmap:
        if is_AA:
            y_values = list('AILMPVFWYNQSTDEHKRGC*')
        else:
            y_values = list('ATGC')

        # Get all unique values from both wt and mutation columns, avoiding duplicates
        all_values = list(set(df['wt'].unique().tolist() + df['mutation'].unique().tolist()))
        order = y_values + [m for m in all_values if m not in y_values]
        df['mutation'] = pd.Categorical(df['mutation'], categories=order, ordered=True)
        df = df.sort_values(['position','mutation'], ascending=[True,False])

        # Ensure unique index by dropping duplicates before creating HeatMap
        df = df.drop_duplicates(subset=kdims)

        plot = hv.HeatMap(df, kdims=kdims, vdims=vdims).opts(clipping_colors={'min':'white'}, clim=(0,None),
                    colorbar=True, clabel=f"frequency (n={total_seqs})", ylabel="mutation")
    else:
        plot = hv.Bars(df, kdims=kdims, vdims=vdims).opts(
                    show_legend=False, ylabel=f"frequency (n={total_seqs})", stacked=True, ylim=(-0.01, None))

    plot = plot.opts(height=500, width=1000, xrotation=xrotation, tools=['hover'],
                     cmap=colormap, xlabel='position', fontsize=fontsize)
    return plot

def str_to_bool(value):
    """
    Convert a string to a boolean if possible, raise error if not.
    """
    if isinstance(value, str):
        if value.lower() == 'true':
            return True
        elif value.lower() == 'false':
            return False
    raise ValueError(f"Invalid boolean string: {value}")


def load_csv_as_dict(csv_path, required=[], lists=[], tag=False, defaults={}):
    """
    Load the given csv file as a nested dictionary in which the first column is used as the outer key 
    and columns are used as the inner key, with values as the inner values.
    If a tag is given then the tag column will be used to filter the dataframe, then removed, 
    and the first of the remaining columns will be used as the outer key.

    params:
        csv_path: str, path to the csv file
        required: list, optional, list of required column values to check for in each row in the csv file
        lists: list, optional, list of column names that should be converted to lists using '__' as the delimiter
        tag: str, optional, the tag to filter the dataframe by. If this is used but no tag column is present,
            no filtering will be done. This is to allow for the same file to be used by multiple tags
        defaults: dict, optional, default values to use if a column is not present or if a value is missing
    
    returns:
        dict, the csv file as a nested dictionary
    """
    if os.path.isfile(csv_path):
        # Read CSV with object dtype to avoid pandas type inference issues
        df = pd.read_csv(csv_path, index_col=False, dtype=str)
        columns = df.columns.tolist()
        
        if ('tag' in df.columns) and tag:
            df = df.loc[df['tag']==tag]
            df = df.drop(columns='tag')
            
        if any([c.startswith('Unnamed: ') for c in df.columns]):
            print(f"[WARNING] Column name beginning with 'Unnamed: ' detected in csv {csv_path}. This usually results from erroneous whitespace characters.\n")
            
        for req in required:
            if req not in columns:
                print(f"[WARNING] Required value `{req}` not found in csv `{csv_path}`. Cannot use this file as a dictionary.\n")
                return {}
            elif df[req].isnull().values.any():
                print(f"[WARNING] Some rows are missing value for required key `{req}` in csv `{csv_path}`. Cannot use this file as a dictionary.\n")
                return {}
                
        # Apply defaults before converting to dict to avoid dtype issues
        for default_col, default_val in defaults.items():
            if default_col not in columns:
                df[default_col] = str(default_val)
            else:
                # Fill NaN/null values with the default, converting to string
                df[default_col] = df[default_col].fillna(str(default_val))

        # issue a warning if the first column is not unique
        if not df[df.columns[0]].is_unique:
            duplicate_values = df[df.duplicated(subset=[df.columns[0]], keep=False)][df.columns[0]].to_list()
            dup_string = '\n'.join(duplicate_values)
            print(f"""[WARNING] The first column `{df.columns[0]}` in csv `{csv_path}` is not unique. This will cause an error.\n
                  Duplicate values:\n{dup_string}\n
                  Please ensure that the first column contains unique values.\n""")
                
        df = df.set_index(df.columns[0])
        csv_dict = df.to_dict('index')
        
        # Clean up the dictionary and convert types appropriately
        cleaned_dict = {}
        for key in csv_dict:
            row_dict = {}
            for k, v in csv_dict[key].items():
                # Skip truly null values (but keep string representations like 'False')
                if pd.isnull(v) or v == 'nan':
                    continue
                    
                # Convert string representations back to appropriate types
                if v == 'True':
                    row_dict[k] = True
                elif v == 'False':
                    row_dict[k] = False
                elif v == 'all':
                    row_dict[k] = 'all'
                else:
                    # Try to convert to numeric if possible, otherwise keep as string
                    try:
                        # Try integer first
                        if '.' not in str(v):
                            row_dict[k] = int(v)
                        else:
                            row_dict[k] = float(v)
                    except (ValueError, TypeError):
                        # Keep as string if conversion fails
                        row_dict[k] = v
                        
            cleaned_dict[key] = row_dict
        
        # Convert list columns to lists
        for col in lists:
            for key in cleaned_dict:
                if col in cleaned_dict[key]:
                    cleaned_dict[key][col] = str(cleaned_dict[key][col]).split('__')

        return cleaned_dict
    else:
        print(f"[WARNING] Provided string `{csv_path}` is not a csv file. Please provide a YAML-formatted dictionary or a csv file.\n")
        return {}


def validate_bc(bc_info_dict, bc, metadata_dir, errors_list=None):
    """
    Validate and convert the barcode info dictionary values for a given barcode type.
    
    This function ensures that string representations from CSV files are converted to 
    appropriate Python types and performs a subset of data validations.
    
    Parameters:
        bc_info_dict (dict): Dictionary containing barcode information for all barcode types
        bc (str): The specific barcode type name to validate
        metadata_dir (str): Path to the metadata directory, used to append to fasta file paths
        errors_list (list): Optional list to append errors to. If None, errors are printed.
    
    Returns:
        tuple: (validated_dict, errors) - The bc_info_dict with validated values and list of errors
    
    Validates the following fields:
        - fasta: appends the provided metadata folder and ensures the file exists
        - reverse_complement: Must be boolean (converts string "True"/"False" to bool)
        - label_only: Must be boolean (converts string "True"/"False" to bool)
        - generate: Must be False, "all", or an integer (converts string representations)
        - hamming_distance: Must be a non-negative integer (converts string to int)
    """
    # Track whether we need to print errors
    print_errors = errors_list is None
    
    # Create a local list to collect errors
    local_errors = []

    bc_data = bc_info_dict.get(bc, {})

    if not bc_data:
        error_msg = f"[ERROR] Barcode `{bc}` not found in barcode_info dictionary.\n"
        local_errors.append(error_msg)
        if print_errors:
            print(error_msg, file=sys.stderr)
        else:
            errors_list.append(error_msg)
        return bc_info_dict, (None if print_errors else errors_list)
    
    # Handle fasta field
    if 'fasta' in bc_data and isinstance(bc_data['fasta'], str):
        bc_fasta = bc_data['fasta']
        if not bc_fasta.startswith(metadata_dir):
            bc_fasta = os.path.join(metadata_dir, bc_fasta)
            bc_data['fasta'] = bc_fasta 
        if not os.path.isfile(bc_fasta):
            local_errors.append(f"[ERROR] Fasta file `{bc_fasta}` for barcode `{bc}` does not exist. Please provide a valid fasta file path.\n")
    
    # Handle reverse_complement field
    if 'reverse_complement' in bc_data and isinstance(bc_data['reverse_complement'], str):
        try:
            bc_data['reverse_complement'] = str_to_bool(bc_data['reverse_complement'])
        except ValueError:
            local_errors.append(f"[ERROR] Invalid boolean string for `reverse_complement` in barcode_info for barcode `{bc}`: {bc_data['reverse_complement']}. Please provide a boolean string (True/False).\n")
    
    # Handle label_only field
    if 'label_only' in bc_data and isinstance(bc_data['label_only'], str):
        try:
            bc_data['label_only'] = str_to_bool(bc_data['label_only'])
        except ValueError:
            local_errors.append(f"[ERROR] Invalid boolean string for `label_only` in barcode_info for barcode `{bc}`: {bc_data['label_only']}. Please provide a boolean string (True/False).\n")
    
    # Handle generate field - can be False, 'all', or an integer
    if 'generate' in bc_data:
        gen_value = bc_data['generate']
        
        if gen_value is False or gen_value == 'all':
            pass
        elif isinstance(gen_value, str):
            if gen_value.lower() == 'false':
                bc_data['generate'] = False
            elif gen_value.lower() == 'all':
                bc_data['generate'] = 'all'
            else:
                try:
                    bc_data['generate'] = int(gen_value)
                except ValueError:
                    local_errors.append(f"[ERROR] Invalid value for `generate` in barcode_info for barcode `{bc}`. Must be False, 'all', or an integer (got '{gen_value}').\n")
        elif isinstance(gen_value, (int, float)):
            bc_data['generate'] = int(gen_value)
        else:
            errors_list.append(f"[ERROR] Invalid type for `generate` in barcode_info for barcode `{bc}`. Must be False, 'all', or an integer (got type {type(gen_value).__name__}).\n")
    
    # Handle hamming_distance field
    if 'hamming_distance' in bc_data:
        ham_value = bc_data['hamming_distance']
        error_msg = None
        try:
            # Try to convert to int regardless of type
            if isinstance(ham_value, (int, float, str)):
                converted_value = int(ham_value)
                if converted_value < 0:
                    error_msg = f"[ERROR] Invalid value for `hamming_distance` in barcode_info for barcode `{bc}`. Must be a non-negative integer (got {ham_value}).\n"
                else:
                    bc_data['hamming_distance'] = converted_value
            else:
                error_msg = f"[ERROR] Invalid type for `hamming_distance` in barcode_info for barcode `{bc}`. Must be a non-negative integer (got type {type(ham_value).__name__}).\n"
        except (ValueError, OverflowError):
            error_msg = f"[ERROR] Invalid value for `hamming_distance` in barcode_info for barcode `{bc}`. Must be a non-negative integer (got '{ham_value}').\n"
        if error_msg: local_errors.append(error_msg)

    # Handle errors based on whether errors_list was provided
    if print_errors:
        # Print all errors to stderr and return None for errors
        for error in local_errors:
            print(error, file=sys.stderr)
        return bc_info_dict, None
    else:
        # Append all errors to the provided list
        errors_list.extend(local_errors)
        return bc_info_dict, errors_list

def dashboard_input(wildcards, config):
    """
    returns the input files for the dashboard based on the sample name provided by the user
    or the first run tag in the config file if no sample name is provided

    args:
        wildcards (snakemake object):   snakemake wildcards object, not used but required for unpacking in rule input
        config (dict):                  dictionary of config file contents

    returns:
        inputDict (dict):               dictionary of input files for the dashboard, genotypes and refFasta, and optionally pdbFile
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
        # Add PDB file if it exists in the tags, otherwise empty string
        inputDict['pdb'] = config['runs'][sample].get('pdb', '')
    elif sample in config.get('timepointsInfo', ''):
        inputDict = {'genotypes': f'mutation_data/timepoints/{sample}_merged-timepoint_genotypes-reduced-dimensions.csv',
                    'refFasta': config['timepointsInfo'][sample]['reference']}
        tag = config['timepointsInfo'][sample]['tag']
        if 'enrichment' in config['runs'][tag]:
            inputDict['genotypes'] = inputDict['genotypes'][:-4] + '-enrichment.csv'
        # Add PDB file if it exists for the first tag in the timepoint, otherwise empty string
        inputDict['pdb'] = config['runs'][tag].get('pdb', '')
    elif '_' in sample and (sample.split('_')[0] in config['runs']):
        print(f'[NOTICE] dashboard_input {sample} contains an underscore, so tag_barcode input is assumed. This may cause an error if this is not a valid tag_barcode combination.')
        tag, barcodes = sample.split('_')
        inputDict = {'genotypes': f'mutation_data/{tag}/{barcodes}/{tag}_{barcodes}_genotypes-reduced-dimensions.csv',
                    'refFasta': config['runs'][tag]['reference']}
        # Add PDB file if it exists in the tags, otherwise empty string
        inputDict['pdb'] = config['runs'][tag].get('pdb', '')
    else:
        # print(f'[NOTICE] dashboard_input {sample} does not conform to a valid tag, timepoint name, or tag_barcode combination. Dashboard input will not be generated.')
        return None
    return inputDict

def _is_fastq_file(filename):
    """Check if a filename is a FASTQ file."""
    return filename.endswith(('.fastq', '.fastq.gz', '.fq.gz'))

def _find_fastq_files_in_path(base_path, extensions=('*.fastq.gz', '*.fq.gz', '*.fastq')):
    """Find all FASTQ files in a given path."""
    files = []
    for ext in extensions:
        pattern = os.path.join(base_path, ext)
        files.extend(glob.glob(pattern, recursive=True))
    return files

def _search_for_specific_file(root_folder, filename):
    """Search for a specific FASTQ file in the root folder structure."""
    pattern = os.path.join(root_folder, '**', filename)
    matches = glob.glob(pattern, recursive=True)

    if len(matches) == 0:
        raise FileNotFoundError(f"File '{filename}' not found in '{root_folder}'")
    elif len(matches) > 1:
        raise ValueError(f"Multiple files found matching '{filename}': {matches}")

    return matches

def _validate_folder_structure(root_folder, folder, subfolders):
    """Validate and report on folder structure."""
    warnings = []

    # Check if folder exists
    folder_pattern = os.path.join(root_folder, '**', folder)
    if not glob.glob(folder_pattern, recursive=True):
        warnings.append(f"Folder '{folder}' not found in root folder '{root_folder}'")
        return warnings, False

    # Check for each subfolder
    for sub in subfolders:
        sub_pattern = os.path.join(root_folder, '**', folder, sub)
        if not glob.glob(sub_pattern, recursive=True):
            warnings.append(f"Subfolder '{sub}' not found under folder '{folder}' in '{root_folder}'")

    return warnings, True

def retrieve_fastqs(root_folder, list_of_folders_or_fqs, subfolder_string, select=''):
    """
    Search for a uniquely named folder within the root folder, and retrieve the file names
    ending in '.fastq.gz' for all files within a list of subfolders.

    Parameters:
        root_folder (str):          The path to the root folder to search in.
        list_of_folders_or_fqs:     Either a list of folders to search within the root directory,
                                   or a list of fastq filenames to find directly in the root folder.
                                   Can also be a single string (folder or fastq filename).
        subfolder_string (str):     A string of comma separated subfolder names to search for within
                                   the uniquely named folder. At least one must be present, but all
                                   don't need to be.
        select (str):              An optional string to select a specific file name within the
                                   subfolders. If not specified, all files ending in '.fastq.gz'
                                   will be retrieved.

    Returns:
        List[str]: A list of the full file paths for all matching FASTQ files.
    """
    # Early validation
    if not os.path.isdir(root_folder):
        print(f"[WARNING] Provided root folder '{root_folder}' does not exist.")
        return []

    # Convert single string to list for uniform handling
    if isinstance(list_of_folders_or_fqs, str):
        list_of_folders_or_fqs = [list_of_folders_or_fqs]

    # Check if all items are fastq files
    if all(_is_fastq_file(item) for item in list_of_folders_or_fqs):
        # Handle list of fastq files
        found_files = []
        for fq_file in list_of_folders_or_fqs:
            pattern = os.path.join(root_folder, '**', fq_file)
            matches = glob.glob(pattern, recursive=True)

            if len(matches) == 0:
                raise FileNotFoundError(f"File '{fq_file}' not found in '{root_folder}'")
            elif len(matches) > 1:
                raise ValueError(f"Multiple files found matching '{fq_file}': {matches}")

            found_files.extend(matches)

        # Ensure we found exactly the same number of files as requested
        if len(found_files) != len(list_of_folders_or_fqs):
            raise ValueError(f"Expected {len(list_of_folders_or_fqs)} files but found {len(found_files)}")

        return found_files

    # If not fastq files, treat as folder names
    folder_list = list_of_folders_or_fqs

    # Parse subfolders
    subfolders = [s.strip() for s in subfolder_string.split(',')]
    file_paths = []
    select_matches = 0
    all_warnings = []

    # Process each folder
    for folder in folder_list:
        # Validate folder structure
        warnings, folder_exists = _validate_folder_structure(root_folder, folder, subfolders)
        for warning in warnings:
            if "Subfolder" in warning:
                print(f"[NOTICE] {warning}")
            else:
                print(f"[WARNING] {warning}")

        if not folder_exists:
            continue

        # Collect FASTQ files from all subfolders
        folder_hits = []
        for sub in subfolders:
            # Build search patterns
            search_paths = glob.glob(os.path.join(root_folder, '**', folder, sub), recursive=True)

            for search_path in search_paths:
                files = _find_fastq_files_in_path(search_path)

                # Apply select filter if specified
                if select:
                    files = [f for f in files if os.path.basename(f) == select]
                    select_matches += len(files)

                folder_hits.extend(files)

        # Report findings for this folder
        if folder_hits:
            if len(folder_hits) > 1 and select:
                print(f"[WARNING] More than one file found that matches select '{select}':")
                for f in folder_hits:
                    print(f"  {f}")
                print(f"Using all {len(folder_hits)} matching files.")
            file_paths.extend(folder_hits)
        else:
            print(f"[WARNING] No .fastq.gz files found within subfolders {subfolders} "
                  f"within folder '{folder}' in root folder '{root_folder}'.")

    # Final reporting
    if select and select_matches == 0:
        print(f"[ERROR] No fastq file matching '{select}' was found.")
    if not file_paths:
        print('[NOTICE] Did not find any sequences to import.')

    return file_paths
