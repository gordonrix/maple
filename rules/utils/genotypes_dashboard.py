#!/usr/bin/env python3
"""
Maple Genotypes Dashboard - 2025 Rebuild
Modern reactive architecture using Panel 1.7+ and HoloViews
"""

import argparse
import datetime
import os
import pathlib
import pickle
import traceback

import datashader as ds
import holoviews as hv
import hvplot.pandas
import numpy as np
import pandas as pd
import panel as pn
import param
import spatialpandas as spd
from bokeh.models import LinearColorMapper, ColorBar, Label
from bokeh.plotting import figure
from holoviews.operation.datashader import datashade as hv_datashade, rasterize, dynspread
from holoviews.selection import link_selections
from holoviews.streams import Selection1D
from holoviews.util.transform import dim
from natsort import natsorted, index_natsorted

# Import existing utilities
from SequenceAnalyzer import SequenceAnalyzer
from common import cmap_dict, conspicuous_mutations, colormaps, export_svg_plots, load_references_from_csv
from structure_heatmap import create_structure_pane

# Configure Panel
pn.extension('tabulator')

# Plot styling constants (matching original dashboard font sizes)
PLOT_FONT_SIZES = {
    'title': 18,
    'labels': 17,
    'xticks': 13,  # Reduced as requested
    'yticks': 13,  # Reduced as requested
    'legend': 13,  # Reduced to 13 as requested
    'clabel': 13
}

def create_histogram_plot(data, hist_column, bins, alpha):
    """Helper function to create histogram plots with consistent styling"""
    return data.hvplot.hist(
        y=hist_column,
        bins=bins,
        width=500,
        height=400,
        xlabel=hist_column,
        ylabel='Total sequences',
        color='grey',
        alpha=alpha,
        fontsize=PLOT_FONT_SIZES
    )


def process_selection(selection_expr, dashboard_instance):
    """
    Process selection for all dashboard components.
    Handles both datashaded (HoloViews dim) and non-datashaded (list) selections.
    Updates analyzer.selected_idx for mutation analysis and table components.
    Returns indices for histogram component.
    """
    # Clear cache whenever selection changes (idx parameter changes)
    dashboard_instance.data_manager.clear_agg_muts_cache()

    # Check type first to avoid weird comparison issues with dim objects
    # Datashaded: HoloViews dim expression
    if str(type(selection_expr)) == "<class 'holoviews.util.transform.dim'>":
        # Get current data for conversion
        current_data = getattr(dashboard_instance, 'current_data', dashboard_instance.data_manager.genotypes)
        dataset = hv.Dataset(current_data)
        selected_df = dataset.select(selection_expr).data

        if len(selected_df) > 0:
            indices = selected_df.index.tolist()
            # Update analyzer for mutation analysis and table
            selected_data = dashboard_instance.data_manager.analyzer.select(idx=indices)
            dashboard_instance.data_manager.analyzer.selected_idx = selected_data['df'].index.tolist()
            # Return indices for histogram
            return indices
        else:
            # Empty selection
            dashboard_instance.data_manager.analyzer.selected_idx = None
            return None

    # Non-datashaded: list of indices
    elif isinstance(selection_expr, list) and len(selection_expr) > 0:
        # Update analyzer for mutation analysis and table
        selected_data = dashboard_instance.data_manager.analyzer.select(idx=selection_expr)
        dashboard_instance.data_manager.analyzer.selected_idx = selected_data['df'].index.tolist()
        # Return indices for histogram
        return selection_expr

    # No selection cases (None or empty list)
    else:
        # Use current filtered data indices instead of None (all data)
        current_data = getattr(dashboard_instance, 'current_data', None)
        if current_data is not None:
            # Set to filtered data indices
            dashboard_instance.data_manager.analyzer.selected_idx = current_data.index.tolist()
            return current_data.index.tolist()
        else:
            # Fallback to None if no current_data
            dashboard_instance.data_manager.analyzer.selected_idx = None
            return None


class DataManager:
    """Handles data loading and basic processing"""

    def __init__(self, genotypes_file, reference_file, exclude_indels=False, downsample_size=None):
        self.reference_file = reference_file
        self.exclude_indels = exclude_indels

        # Load reference sequences from CSV
        self.references_dict, errors, _, _ = load_references_from_csv(reference_file)
        if errors:
            raise ValueError(f"Error loading references: {errors}")

        self.current_reference = list(self.references_dict.keys())[0]

        # Load only genotypes DataFrame (lightweight)
        genotypes = pd.read_csv(genotypes_file, low_memory=False)
        if exclude_indels:
            genotypes = genotypes[genotypes['NT_insertions'].isna()]
            genotypes = genotypes[genotypes['NT_deletions'].isna()]
        genotypes = genotypes.reset_index().drop('index', axis='columns')

        self.do_AA_analysis = 'AA_substitutions_nonsynonymous' in genotypes.columns

        if 'NT_muts_of_interest' not in genotypes.columns:
            genotypes['NT_muts_of_interest'] = 'none'
        if self.do_AA_analysis and 'AA_muts_of_interest' not in genotypes.columns:
            genotypes['AA_muts_of_interest'] = 'none'

        self.genotypes = genotypes
        self.numerical_columns = genotypes.select_dtypes(include=[np.number]).columns.tolist()
        all_embeddings = [col for col in self.numerical_columns if 'PaCMAP' in col or 'UMAP' in col or 'tSNE' in col]
        self.embedding_columns = sorted(all_embeddings, key=lambda x: (0 if x.startswith('AA_') else 1, x))

        # Downsample before creating SequenceAnalyzer
        if downsample_size and downsample_size < len(genotypes):
            np.random.seed(0)
            counts = genotypes['count'].values if 'count' in genotypes.columns else np.ones(len(genotypes))
            probabilities = counts / counts.sum()
            sampled_idx = np.sort(np.random.choice(np.arange(len(genotypes)), size=downsample_size, replace=False, p=probabilities))
            sampled_genotypes = genotypes.iloc[sampled_idx]
        else:
            sampled_genotypes = genotypes

        # Create SequenceAnalyzer with downsampled data only
        ref_data = self.references_dict[self.current_reference]
        reference_sequences = [
            ref_data['alignment_seq'],
            ref_data['NT_seq'],
            ref_data['coding_seq']
        ]

        self.analyzer = SequenceAnalyzer(
            reference_sequences=reference_sequences,
            genotypesCSV=sampled_genotypes,
            exclude_indels=False
        )

        # Cache for aggregate_mutations to avoid duplicate computation
        self._agg_muts_cache = {}

    def get_structure_file(self):
        """Get structure file for current reference"""
        ref_data = self.references_dict.get(self.current_reference, {})
        structure_file = ref_data.get('structure_file', '')
        return structure_file if structure_file else None

    def get_aggregated_mutations(self, NTorAA, idx=None, unique_only=False):
        """
        Get aggregated mutations with caching to avoid duplicate computation.

        Parameters:
            NTorAA (str): 'NT' or 'AA'
            idx (list): Indices to aggregate, or None for all
            unique_only (bool): Whether to count unique sequences only

        Returns:
            pd.DataFrame: Aggregated mutations dataframe
        """
        # Create cache key from parameters
        idx_key = tuple(sorted(idx)) if idx is not None and len(idx) > 0 else None
        cache_key = (NTorAA, idx_key, unique_only)

        # Return cached result if available
        if cache_key in self._agg_muts_cache:
            return self._agg_muts_cache[cache_key]

        # Compute and cache result
        result = self.analyzer.aggregate_mutations(
            NTorAA=NTorAA,
            idx=idx,
            unique_only=unique_only
        )
        self._agg_muts_cache[cache_key] = result

        return result

    def clear_agg_muts_cache(self):
        """Clear the aggregated mutations cache"""
        self._agg_muts_cache = {}



class ScatterPlotComponent(pn.viewable.Viewer):
    """Interactive scatter plot with datashading support"""

    # Parameters for the component
    x_column = param.String(default='NT_PaCMAP1')
    y_column = param.String(default='NT_PaCMAP2')
    color_by = param.String(default='NT_substitutions_count')
    colormap = param.String(default='kbc_r')
    datashade = param.Boolean(default=True)
    point_size = param.Integer(default=5, bounds=(1, 20))
    alpha = param.Number(default=0.8, bounds=(0, 1))

    def __init__(self, data_manager, filtered_data=None, **params):
        super().__init__(**params)
        self.data_manager = data_manager
        self.filtered_data = filtered_data or data_manager.genotypes

        # Set up default columns if they exist
        if data_manager.embedding_columns:
            self.x_column = data_manager.embedding_columns[0]
            if len(data_manager.embedding_columns) > 1:
                self.y_column = data_manager.embedding_columns[1]

    @param.depends('x_column', 'y_column', 'color_by', 'datashade', 'point_size', 'alpha')
    def plot(self):
        """Create the main scatter plot using hvplot - reactive method"""

        # Get current data (reactive or static)
        if hasattr(self.filtered_data, 'rx'):
            # Reactive data
            current_data = self.filtered_data.rx.value
        else:
            # Static data
            current_data = self.filtered_data

        if self.datashade:
            # Use hvplot with datashader
            return current_data.hvplot.scatter(
                x=self.x_column,
                y=self.y_column,
                c=self.color_by,
                datashade=True,
                width=800,
                height=600,
                fontsize=PLOT_FONT_SIZES
            )
        else:
            # Regular scatter plot
            return current_data.hvplot.scatter(
                x=self.x_column,
                y=self.y_column,
                c=self.color_by,
                size=self.point_size,
                alpha=self.alpha,
                width=800,
                height=600,
                fontsize=PLOT_FONT_SIZES
            )

    def __panel__(self):
        """Return the panel object"""
        return pn.pane.HoloViews(
            self.plot,
            linked_axes=False,
            sizing_mode='stretch_width'
        )


class HistogramComponent(pn.viewable.Viewer):
    """Interactive histogram with selection overlay"""

    # Parameters for the component
    hist_column = param.String(default='NT_substitutions_count')
    num_bins = param.Integer(default=20, bounds=(1, 100))

    def __init__(self, data_manager, **params):
        super().__init__(**params)
        self.data_manager = data_manager

        # Set default column - prioritize AA_substitutions_nonsynonymous_count, then NT_substitutions_count
        if 'AA_substitutions_nonsynonymous_count' in data_manager.numerical_columns:
            self.hist_column = 'AA_substitutions_nonsynonymous_count'
        elif 'NT_substitutions_count' in data_manager.numerical_columns:
            self.hist_column = 'NT_substitutions_count'
        elif data_manager.numerical_columns:
            self.hist_column = data_manager.numerical_columns[0]


class TableComponent(pn.viewable.Viewer):
    """Interactive data table for displaying selected genotypes"""

    # Parameters
    sample_size = param.Integer(default=50, bounds=(5, 500), doc="Number of genotypes to display in table")

    def __init__(self, data_manager, **params):
        super().__init__(**params)
        self.data_manager = data_manager

    @param.depends('sample_size')
    def table(self):
        """Create data table - reactive method"""
        # Get selected data from the analyzer
        selected_idx = getattr(self.data_manager.analyzer, 'selected_idx', None)

        if selected_idx is not None and len(selected_idx) > 0:
            # Use selected data from analyzer's downsampled genotypes
            df = self.data_manager.analyzer.genotypes.loc[selected_idx].copy()
        else:
            # Use all analyzer data if no selection (already downsampled)
            df = self.data_manager.analyzer.genotypes.copy()

        # Drop unnecessary columns like original dashboard
        drop_cols = ['size', 'index']  # Hide index column from table
        # Check if barcode column should be dropped (if num_barcodes control exists and is disabled)
        if 'barcode(s)_subset' in df.columns:
            drop_cols.append('barcode(s)_subset')

        df = df.drop(columns=[col for col in drop_cols if col in df.columns], errors='ignore')

        # Sort by timepoint then genotype_ID if timepoint exists, otherwise just genotype_ID
        if 'timepoint' in df.columns and 'genotype_ID' in df.columns:
            # Sort by timepoint first, then by genotype_ID within each timepoint
            df = df.sort_values(['timepoint', 'genotype_ID'], key=lambda x: np.argsort(index_natsorted(x)))
        elif 'genotype_ID' in df.columns:
            # Sort by genotype_ID only if no timepoint column
            df = df.iloc[index_natsorted(df['genotype_ID'])]

        # Take first N rows if larger than sample size
        if self.sample_size < len(df) and len(df) > 0:
            df = df.head(self.sample_size)

        # Return empty table if no data
        if df.empty:
            return pn.widgets.Tabulator(pd.DataFrame(), width=400, height=600)

        # Reset index to hide original index and create tabulator widget
        # Width adjusted: 662 + 46 = 708
        return pn.widgets.Tabulator(
            df.reset_index(drop=True),
            width=708,
            height=600,  # Same height as scatter plot
            pagination='remote',
            page_size=20,
            sizing_mode='fixed',
            show_index=False  # Don't show the index column
        )

    def __panel__(self):
        """Return the panel object"""
        return self.table


class MutationAnalysisComponent(pn.viewable.Viewer):
    """Component for displaying aggregated mutation analysis"""

    # Parameters
    most_common = param.Boolean(default=True, doc="Show most common (True) or least common (False) mutations")
    num_positions = param.Integer(default=20, bounds=(1, None), doc="Number of mutation positions to display")
    colormap = param.String(default='kbc_r', doc="Colormap for heatmap")
    NTorAA = param.String(default='NT', doc="Nucleotide (NT) or Amino Acid (AA) analysis")
    plot_type = param.String(default='heatmap', doc="Plot type: heatmap or bars")
    position_mode = param.String(default='most common', doc="How to select positions: most common, least common, or range")
    position_range = param.String(default='', doc="Position range(s) to display, e.g., '20-30' or '20-30,50-60'")
    axis_type = param.String(default='categorical', doc="Axis type: numerical or categorical")
    plot_width = param.Integer(default=500, bounds=(300, 1200), doc="Width of the mutation analysis plot")

    def __init__(self, data_manager, **params):
        super().__init__(**params)
        self.data_manager = data_manager

        # Set default NTorAA to AA if available, otherwise NT
        if 'NTorAA' not in params:
            has_aa = (hasattr(data_manager.analyzer, 'integer_matrix') and
                      'AA' in data_manager.analyzer.integer_matrix)
            self.NTorAA = 'AA' if has_aa else 'NT'

        # Set max positions based on sequence length
        self._update_max_positions()

    def _update_max_positions(self):
        """Update max positions based on current NT/AA selection"""
        if hasattr(self.data_manager.analyzer, 'integer_matrix') and self.NTorAA in self.data_manager.analyzer.integer_matrix:
            max_positions = self.data_manager.analyzer.integer_matrix[self.NTorAA].shape[1]
            self.param.num_positions.bounds = (1, max_positions)

    @param.depends('most_common', 'num_positions', 'colormap', 'NTorAA', 'plot_type', 'position_mode', 'position_range', 'axis_type', 'plot_width')
    def plot(self):
        """Create aggregated mutations heatmap - reactive method"""

        # Update max positions when NTorAA changes
        self._update_max_positions()

        # Use slider width (dynamic width calculation removed since slider is always used)
        actual_width = self.plot_width

        # Use selected data if available, otherwise use all data
        # Get selected indices from the data manager's analyzer (updated by dashboard)
        selected_idx = getattr(self.data_manager.analyzer, 'selected_idx', None)

        # Get aggregated mutations using cached method
        try:
            agg_df = self.data_manager.get_aggregated_mutations(
                NTorAA=self.NTorAA,  # Use parameter instead of hardcoded
                idx=selected_idx,    # Use selected sequences or None (all) if no selection
                unique_only=False
            )

            total_seqs = self.data_manager.analyzer.get_count(idx=selected_idx)

            # Create the plot using conspicuous_mutations with global font sizes
            is_heatmap = (self.plot_type == 'heatmap')

            # Use different colormaps for bars vs heatmaps
            if is_heatmap:
                colormap_to_use = self.colormap  # Use the selected colormap for heatmaps
            else:
                # Use the original dashboard's colormaps for bars
                colormap_to_use = colormaps[self.NTorAA]

            # Set parameters based on position_mode
            if self.position_mode == 'range':
                num_positions_arg = None
                most_common_arg = True  # Default for range mode
                position_ranges_arg = self.position_range if self.position_range.strip() else None
            else:
                num_positions_arg = self.num_positions
                most_common_arg = (self.position_mode == 'most common')
                position_ranges_arg = None

            plot = conspicuous_mutations(
                df=agg_df,
                total_seqs=total_seqs,
                num_positions=num_positions_arg,
                colormap=colormap_to_use,
                most_common=most_common_arg,
                heatmap=is_heatmap,
                fontsize=PLOT_FONT_SIZES,  # Use global font sizes
                position_ranges=position_ranges_arg,
                axis_type=self.axis_type
            )

            return plot.opts(
                height=428,
                width=actual_width
                # Don't override fontsize - already set in conspicuous_mutations call
            )

        except Exception as e:
            # Return error message if something goes wrong
            return hv.Text(0.5, 0.5, f"Error: {str(e)}").opts(
                width=600,
                height=400,
                xlim=(0, 1),
                ylim=(0, 1)
            )

    def __panel__(self):
        """Return the panel object"""
        return pn.pane.HoloViews(
            self.plot,
            linked_axes=False,
            sizing_mode='stretch_width'
        )


class StructureHeatmapComponent(pn.viewable.Viewer):
    """Component for displaying 3D protein structure viewer with mutation frequency coloring"""

    # Parameters
    colormap = param.String(default='kbc_r', doc="Colormap for structure heatmap")
    label_step = param.Selector(default=50, objects=[None, 10, 50, 100], doc="Label every Nth residue (None = no labels)")
    frequency_floor = param.Number(default=0.0, bounds=(0, 1), doc="Minimum frequency threshold for coloring residues")

    def __init__(self, data_manager, **params):
        super().__init__(**params)
        self.data_manager = data_manager

    def plot(self, position_range=''):
        """Create structure viewer - method that takes position_range as argument"""

        # Get structure file from data_manager based on current reference
        structure_file = self.data_manager.get_structure_file()

        # Return empty spacer if no structure file
        if not structure_file:
            return pn.Spacer(width=0)

        # Use selected data if available, otherwise use all data
        selected_idx = getattr(self.data_manager.analyzer, 'selected_idx', None)

        # Get aggregated mutations using cached method
        try:
            agg_df = self.data_manager.get_aggregated_mutations(
                NTorAA='AA',  # Structure viewer always uses AA
                idx=selected_idx,
                unique_only=False
            )

            # Create structure viewer pane
            structure_pane = create_structure_pane(
                structure_file=structure_file,
                agg_df=agg_df,
                colormap=self.colormap,
                width=500,
                height=428,
                label_step=self.label_step,
                min_identity=0.95,
                frequency_floor=self.frequency_floor,
                position_range=position_range
            )

            return structure_pane

        except Exception as e:
            # Return error message if something goes wrong
            return pn.pane.Markdown(f"**Error creating structure viewer:**\n\n{str(e)}", width=800, height=600)

    def __panel__(self):
        """Return the panel object"""
        return self.plot


class ControlsPanel(pn.viewable.Viewer):
    """Control widgets for the dashboard"""

    def __init__(self, data_manager, scatter_component, histogram_component, mutation_analysis_component, table_component, downsample_slider, structure_component=None, dataset_selector=None, reference_selector=None, **params):
        super().__init__(**params)
        self.data_manager = data_manager
        self.scatter = scatter_component
        self.histogram = histogram_component
        self.mutation_analysis = mutation_analysis_component
        self.table = table_component
        self.downsample_slider = downsample_slider
        self.structure = structure_component
        self.dataset_selector = dataset_selector
        self.reference_selector = reference_selector

        # Create control widgets
        self._create_widgets()

    def _create_widgets(self):
        """Create all control widgets"""

        # Data controls (downsample_slider passed from main dashboard)

        self.datashade_toggle = pn.widgets.Checkbox(
            name='Enable Datashading',
            value=True,
            width=200
        )

        # Visualization controls
        self.x_select = pn.widgets.Select(
            name='X Column',
            options=self.data_manager.numerical_columns,
            value=self.scatter.x_column,
            width=200
        )

        self.y_select = pn.widgets.Select(
            name='Y Column',
            options=self.data_manager.numerical_columns,
            value=self.scatter.y_column,
            width=200
        )

        # Build color options including numerical columns and mutations of interest
        color_options = self.data_manager.numerical_columns.copy()
        color_options.append('NT_muts_of_interest')
        if self.data_manager.do_AA_analysis:
            color_options.append('AA_muts_of_interest')

        self.color_select = pn.widgets.Select(
            name='Color By',
            options=color_options,
            value=self.scatter.color_by,
            width=200
        )

        # Create colormap selector like original dashboard
        colormap_dict = cmap_dict()
        self.colormap_select = pn.widgets.Select(
            name='Colormap',
            options=list(colormap_dict.keys()),
            value=self.scatter.colormap,
            width=200
        )

        self.size_slider = pn.widgets.IntSlider(
            name='Point Size',
            start=1,
            end=20,
            value=self.scatter.point_size,
            width=200
        )

        self.alpha_slider = pn.widgets.FloatSlider(
            name='Opacity',
            start=0.1,
            end=1.0,
            step=0.1,
            value=self.scatter.alpha,
            width=200
        )

        # Histogram controls
        self.hist_column_select = pn.widgets.Select(
            name='Histogram Column',
            options=self.data_manager.numerical_columns,
            value=self.histogram.hist_column,
            width=200
        )

        self.hist_bins_slider = pn.widgets.IntSlider(
            name='Number of Bins',
            start=1,
            end=50,
            value=self.histogram.num_bins,
            width=200
        )

        # Selection opacity controls
        self.selected_alpha_slider = pn.widgets.FloatSlider(
            name='Selected Opacity',
            start=0.0,
            end=1.0,
            step=0.1,
            value=1.0,
            width=200
        )

        self.unselected_alpha_slider = pn.widgets.FloatSlider(
            name='Unselected Opacity',
            start=0.0,
            end=1.0,
            step=0.1,
            value=0.2,
            width=200
        )

        # Mutation analysis controls
        self.freq_select = pn.widgets.Select(
            name='Mutation Positions',
            options=['Most Common', 'Least Common', 'Range'],
            value='Most Common',
            width=200
        )

        # Check if AA analysis is available
        has_aa = (hasattr(self.data_manager.analyzer, 'integer_matrix') and
                  'AA' in self.data_manager.analyzer.integer_matrix)

        NTorAA_options = ['NT']
        if has_aa:
            NTorAA_options.append('AA')

        # Use the mutation analysis component's default (which is already set correctly)
        self.NTorAA_toggle = pn.widgets.RadioButtonGroup(
            name='Analysis Type',
            options=NTorAA_options,
            value=self.mutation_analysis.NTorAA,  # Get default from component
            button_type='default',
            width=200
        )

        self.plot_type_toggle = pn.widgets.RadioButtonGroup(
            name='Plot Type',
            options=['heatmap', 'bars'],
            value=self.mutation_analysis.plot_type,
            button_type='default',
            width=200
        )

        self.axis_type_toggle = pn.widgets.RadioButtonGroup(
            name='Axis Type',
            options=['numerical', 'categorical'],
            value=self.mutation_analysis.axis_type,
            button_type='default',
            width=200
        )

        # Get max positions for the current NT/AA selection
        current_NTorAA = self.NTorAA_toggle.value
        max_pos = self.data_manager.analyzer.integer_matrix[current_NTorAA].shape[1] if has_aa or current_NTorAA == 'NT' else 50

        self.num_positions_text = pn.widgets.TextInput(
            name='Number of positions to display',
            placeholder=f'from 5 to {max_pos}',
            value=str(self.mutation_analysis.num_positions),
            width=200
        )

        self.range_input = pn.widgets.TextInput(
            name='Range(s)',
            placeholder='e.g., 20-30 or 20-30,50-60',
            value='',
            width=200
        )

        # Link visibility of widgets to the frequency select and update axis options
        def update_widget_visibility(event):
            is_range = event.new == 'Range'
            self.num_positions_text.visible = not is_range
            self.range_input.visible = is_range
            # Reset range input when switching away from Range mode
            if not is_range:
                self.range_input.value = ''
            # Force sync the component parameter after changing axis options
            self._update_axis_type_options()
            self.mutation_analysis.param.trigger('axis_type')

        self.freq_select.param.watch(update_widget_visibility, 'value')

        # Watch range input to update axis options when range changes
        self.range_input.param.watch(lambda event: self._update_axis_type_options(), 'value')

        # Set initial visibility
        is_range = self.freq_select.value == 'Range'
        self.num_positions_text.visible = not is_range
        self.range_input.visible = is_range

        self.plot_width_slider = pn.widgets.IntSlider(
            name='Plot Width',
            start=300,
            end=1200,
            value=self.mutation_analysis.plot_width,
            width=200
        )

        # Table controls
        self.sample_size_slider = pn.widgets.IntSlider(
            name='Table Sample Size',
            start=5,
            end=500,
            step=1,
            value=self.table.sample_size,
            width=200
        )

        # Mutations of interest controls
        self.NT_muts_text = pn.widgets.TextInput(
            name='NT mutations of interest',
            placeholder='e.g., A100T, 150, G_, _C',
            value='',
            width=200
        )

        self.AA_muts_text = pn.widgets.TextInput(
            name='AA mutations of interest',
            placeholder='e.g., A100T, 150, G_, _C',
            value='',
            width=200,
            disabled=not self.data_manager.do_AA_analysis
        )

        self.max_mut_combos_slider = pn.widgets.IntSlider(
            name='Max mutation combinations',
            start=1,
            end=20,
            step=1,
            value=10,
            width=200
        )

        # MultiChoice widgets for filtering by mutations
        self.NT_muts_multichoice = pn.widgets.MultiChoice(
            name='Select NT mutations to filter',
            options=[],  # Will be populated dynamically
            value=[],
            width=200
        )

        self.AA_muts_multichoice = pn.widgets.MultiChoice(
            name='Select AA mutations to filter',
            options=[],  # Will be populated dynamically
            value=[],
            width=200,
            disabled=not self.data_manager.do_AA_analysis
        )

        # Export controls
        self.selection_name_text = pn.widgets.TextInput(
            name='Selection name',
            value='unnamed',
            width=200
        )

        self.export_genotypes_button = pn.widgets.Button(
            name='Genotypes CSV',
            button_type='primary',
            width=200
        )

        self.export_NT_fasta_button = pn.widgets.Button(
            name='NT FASTA',
            button_type='primary',
            width=95
        )

        self.export_AA_fasta_button = pn.widgets.Button(
            name='AA FASTA',
            button_type='primary',
            width=95,
            disabled=not self.data_manager.do_AA_analysis
        )

        self.export_NT_consensus_button = pn.widgets.Button(
            name='NT consensus',
            button_type='primary',
            width=95
        )

        self.export_AA_consensus_button = pn.widgets.Button(
            name='AA consensus',
            button_type='primary',
            width=95,
            disabled=not self.data_manager.do_AA_analysis
        )

        self.export_agg_muts_button = pn.widgets.Button(
            name='Aggregated mutations CSV',
            button_type='primary',
            width=200
        )

        self.export_plots_button = pn.widgets.Button(
            name='Plots SVG',
            button_type='primary',
            width=200
        )

        # Data filter controls - two independent numerical column filters
        numerical_columns = [col for col in self.data_manager.numerical_columns if 'muts_of_interest' not in col]

        # Filter 1
        self.filter1_column_select = pn.widgets.Select(
            name='Filter 1 Column',
            options=numerical_columns,
            value=numerical_columns[0] if numerical_columns else None,
            width=200
        )

        self.filter1_range_slider = pn.widgets.RangeSlider(
            name='Filter 1 Range',
            start=0,
            end=100,
            value=(0, 100),
            width=200
        )

        # Filter 2
        self.filter2_column_select = pn.widgets.Select(
            name='Filter 2 Column',
            options=numerical_columns,
            value=numerical_columns[1] if len(numerical_columns) > 1 else (numerical_columns[0] if numerical_columns else None),
            width=200
        )

        self.filter2_range_slider = pn.widgets.RangeSlider(
            name='Filter 2 Range',
            start=0,
            end=100,
            value=(0, 100),
            width=200
        )

        # Structure viewer controls (only create if structure component exists)
        if self.structure is not None:
            self.structure_label_step_select = pn.widgets.Select(
                name='Residue Label Step',
                options=[None, 10, 50, 100],
                value=self.structure.label_step,
                width=200
            )

            self.structure_frequency_floor_slider = pn.widgets.FloatSlider(
                name='Frequency Floor',
                start=0.0,
                end=1.0,
                step=0.05,
                value=self.structure.frequency_floor,
                width=200
            )

        # Initialize range sliders with actual data ranges
        self._update_filter_ranges()

        # Link widgets to component parameters
        self._link_widgets()

    def _link_widgets(self):
        """Link control widgets to component parameters"""

        # Link scatter plot controls
        self.x_select.link(self.scatter, value='x_column')
        self.y_select.link(self.scatter, value='y_column')
        self.color_select.link(self.scatter, value='color_by')
        self.colormap_select.link(self.scatter, value='colormap')
        self.datashade_toggle.link(self.scatter, value='datashade')
        self.size_slider.link(self.scatter, value='point_size')
        # Opacity is handled directly in plot creation function via selected_alpha/unselected_alpha parameters

        # Link histogram controls
        self.hist_column_select.link(self.histogram, value='hist_column')
        self.hist_bins_slider.link(self.histogram, value='num_bins')

        # Link mutation analysis controls
        # Convert 'Most Common'/'Least Common'/'Range' to position_mode
        def freq_callback(target, event):
            if event.new == 'Most Common':
                target.position_mode = 'most common'
            elif event.new == 'Least Common':
                target.position_mode = 'least common'
            else:  # Range
                target.position_mode = 'range'
        self.freq_select.param.watch(lambda event: freq_callback(self.mutation_analysis, event), 'value')

        # Link NTorAA toggle
        self.NTorAA_toggle.link(self.mutation_analysis, value='NTorAA')

        # Link plot type toggle
        self.plot_type_toggle.link(self.mutation_analysis, value='plot_type')

        # Link axis type toggle
        self.axis_type_toggle.link(self.mutation_analysis, value='axis_type')

        # Handle text input for positions with validation
        def positions_callback(target, event):
            try:
                value = int(event.new)
                bounds = target.param.num_positions.bounds
                if bounds[0] <= value <= bounds[1]:
                    target.num_positions = value
                else:
                    # Reset to valid range if out of bounds
                    self.num_positions_text.value = str(target.num_positions)
            except ValueError:
                # Reset to current value if not a valid integer
                self.num_positions_text.value = str(target.num_positions)
        self.num_positions_text.param.watch(lambda event: positions_callback(self.mutation_analysis, event), 'value')

        # Link range input directly
        self.range_input.link(self.mutation_analysis, value='position_range')

        # Initial axis type options update
        self._update_axis_type_options()

    def _update_axis_type_options(self):
        """Update axis type options based on current position selection mode"""
        if self.freq_select.value == 'Range':
            range_str = self.range_input.value.strip()
            if not range_str or ',' not in range_str:
                # Empty range or single range - default to numerical, allow both
                self.axis_type_toggle.value = 'numerical'
                self.axis_type_toggle.disabled = False
            else:
                # Multiple ranges - force categorical and disable
                self.axis_type_toggle.value = 'categorical'
                self.axis_type_toggle.disabled = True
        else:
            # Most/Least Common modes - force categorical and disable
            self.axis_type_toggle.value = 'categorical'
            self.axis_type_toggle.disabled = True

        # Update text input prompt when NTorAA changes
        def update_positions_prompt(event):
            new_NTorAA = event.new
            if hasattr(self.data_manager.analyzer, 'integer_matrix') and new_NTorAA in self.data_manager.analyzer.integer_matrix:
                max_pos = self.data_manager.analyzer.integer_matrix[new_NTorAA].shape[1]
                self.num_positions_text.placeholder = f'from 5 to {max_pos}'
        self.NTorAA_toggle.param.watch(update_positions_prompt, 'value')

        # Link plot width slider directly
        self.plot_width_slider.link(self.mutation_analysis, value='plot_width')
        # Set flag to use slider width instead of dynamic calculation
        self.mutation_analysis._use_slider_width = True

        self.colormap_select.link(self.mutation_analysis, value='colormap')  # Share colormap with scatter plot

        # Link table controls
        self.sample_size_slider.link(self.table, value='sample_size')

        # Link structure viewer controls (if structure component exists)
        if self.structure is not None:
            self.structure_label_step_select.link(self.structure, value='label_step')
            self.structure_frequency_floor_slider.link(self.structure, value='frequency_floor')
            self.colormap_select.link(self.structure, value='colormap')  # Share colormap with scatter and mutation plots

        # Set up filter column callbacks to update range sliders
        def update_filter1_range(event):
            if event.new and event.new in self.data_manager.genotypes.columns:
                col_data = self.data_manager.genotypes[event.new]
                self.filter1_range_slider.start = np.floor(col_data.min())
                self.filter1_range_slider.end = np.ceil(col_data.max())
                self.filter1_range_slider.value = (np.floor(col_data.min()), np.ceil(col_data.max()))
                self.filter1_range_slider.name = f'{event.new} range'

        def update_filter2_range(event):
            if event.new and event.new in self.data_manager.genotypes.columns:
                col_data = self.data_manager.genotypes[event.new]
                self.filter2_range_slider.start = np.floor(col_data.min())
                self.filter2_range_slider.end = np.ceil(col_data.max())
                self.filter2_range_slider.value = (np.floor(col_data.min()), np.ceil(col_data.max()))
                self.filter2_range_slider.name = f'{event.new} range'

        self.filter1_column_select.param.watch(update_filter1_range, 'value')
        self.filter2_column_select.param.watch(update_filter2_range, 'value')

        # Link export button callbacks (plot export callback will be set by Dashboard)
        self.export_genotypes_button.on_click(self._export_genotypes)
        self.export_NT_fasta_button.on_click(self._export_NT_fasta)
        self.export_AA_fasta_button.on_click(self._export_AA_fasta)
        self.export_NT_consensus_button.on_click(self._export_NT_consensus)
        self.export_AA_consensus_button.on_click(self._export_AA_consensus)
        self.export_agg_muts_button.on_click(self._export_agg_muts)


    def _export_genotypes(self, event):
        """Export selected genotypes to CSV"""
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        filename = f'dashboard/dashboard-{self.selection_name_text.value}-{timestamp}_genotypes.csv'

        # Create dashboard directory
        pathlib.Path(filename).parent.absolute().mkdir(parents=True, exist_ok=True)

        # Get selected data from analyzer
        selected_idx = self.data_manager.analyzer.selected_idx
        if selected_idx is not None and len(selected_idx) > 0:
            selected_genotypes = self.data_manager.genotypes.loc[selected_idx]
        else:
            selected_genotypes = self.data_manager.genotypes

        # Drop temporary columns
        drop_cols = ['size', 'index']  # Exclude index column from export
        if 'barcode(s)_subset' in selected_genotypes.columns:
            drop_cols.append('barcode(s)_subset')
        selected_genotypes = selected_genotypes.drop(columns=[col for col in drop_cols if col in selected_genotypes.columns], errors='ignore')

        # Export
        selected_genotypes.to_csv(filename)
        print(f'Exported {len(selected_genotypes)} genotypes to {filename}')

    def _export_NT_fasta(self, event):
        """Export NT sequences to FASTA"""
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        filename = f'dashboard/dashboard-{self.selection_name_text.value}-{timestamp}_NT-seqs.fasta'

        selected_idx = self.data_manager.analyzer.selected_idx
        self.data_manager.analyzer.write_fasta(filename, 'NT', idx=selected_idx)
        count = len(selected_idx) if selected_idx is not None else len(self.data_manager.genotypes)
        print(f'Exported {count} NT sequences to {filename}')

    def _export_AA_fasta(self, event):
        """Export AA sequences to FASTA"""
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        filename = f'dashboard/dashboard-{self.selection_name_text.value}-{timestamp}_AA-seqs.fasta'

        selected_idx = self.data_manager.analyzer.selected_idx
        self.data_manager.analyzer.write_fasta(filename, 'AA', idx=selected_idx)
        count = len(selected_idx) if selected_idx is not None else len(self.data_manager.genotypes)
        print(f'Exported {count} AA sequences to {filename}')

    def _export_NT_consensus(self, event):
        """Export NT consensus sequence to FASTA"""
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        filename = 'dashboard/dashboard-NT-consensus.fasta'

        selected_idx = self.data_manager.analyzer.selected_idx
        self.data_manager.analyzer.get_consensus('NT', idx=selected_idx, write_to=filename, append=True,
                                                   name=f'{self.selection_name_text.value}_{timestamp}')
        print(f'Appended NT consensus sequence to {filename}')

    def _export_AA_consensus(self, event):
        """Export AA consensus sequence to FASTA"""
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        filename = 'dashboard/dashboard-AA-consensus.fasta'

        selected_idx = self.data_manager.analyzer.selected_idx
        self.data_manager.analyzer.get_consensus('AA', idx=selected_idx, write_to=filename, append=True,
                                                   name=f'{self.selection_name_text.value}_{timestamp}')
        print(f'Appended AA consensus sequence to {filename}')

    def _export_agg_muts(self, event):
        """Export aggregated mutations to CSV"""
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        filename = f'dashboard/dashboard-{self.selection_name_text.value}-{timestamp}_mutations-aggregated.csv'

        # Create dashboard directory
        pathlib.Path(filename).parent.absolute().mkdir(parents=True, exist_ok=True)

        # Get selected data and aggregate using cached method
        selected_idx = self.data_manager.analyzer.selected_idx
        NTorAA = self.NTorAA_toggle.value  # Use current NT/AA toggle setting
        df = self.data_manager.get_aggregated_mutations(NTorAA=NTorAA, idx=selected_idx)

        # Export
        df.to_csv(filename, index=False)
        print(f'Exported aggregated {NTorAA} mutations to {filename}')


    def _update_filter_ranges(self):
        """Initialize filter range sliders with actual data ranges"""
        # Update filter 1
        if self.filter1_column_select.value and self.filter1_column_select.value in self.data_manager.genotypes.columns:
            col_data = self.data_manager.genotypes[self.filter1_column_select.value]
            self.filter1_range_slider.start = np.floor(col_data.min())
            self.filter1_range_slider.end = np.ceil(col_data.max())
            self.filter1_range_slider.value = (np.floor(col_data.min()), np.ceil(col_data.max()))
            self.filter1_range_slider.name = f'{self.filter1_column_select.value} range'

        # Update filter 2
        if self.filter2_column_select.value and self.filter2_column_select.value in self.data_manager.genotypes.columns:
            col_data = self.data_manager.genotypes[self.filter2_column_select.value]
            self.filter2_range_slider.start = np.floor(col_data.min())
            self.filter2_range_slider.end = np.ceil(col_data.max())
            self.filter2_range_slider.value = (np.floor(col_data.min()), np.ceil(col_data.max()))
            self.filter2_range_slider.name = f'{self.filter2_column_select.value} range'

    def __panel__(self):
        """Return the controls panel layout"""

        # Build data controls with dataset and reference selectors at top
        data_controls_widgets = [
            pn.Row(self.dataset_selector, self.reference_selector),
            pn.Row(self.downsample_slider, self.datashade_toggle),
            pn.Row(self.filter1_column_select, self.filter1_range_slider),
            pn.Row(self.filter2_column_select, self.filter2_range_slider)
        ]

        data_controls = pn.Card(
            pn.Column(*data_controls_widgets),
            title="Data Controls",
            collapsed=False
        )

        viz_controls = pn.Card(
            pn.Column(
                pn.Row(self.x_select, self.y_select),
                pn.Row(self.color_select, self.colormap_select),
                self.size_slider
            ),
            title="Plot Controls",
            collapsed=False
        )

        muts_of_interest_controls = pn.Card(
            pn.Column(
                pn.Row(self.NT_muts_text, self.AA_muts_text),
                self.max_mut_combos_slider,
                pn.pane.Markdown("**Filter by mutations:**"),
                pn.Row(self.NT_muts_multichoice, self.AA_muts_multichoice)
            ),
            title="Mutations of Interest Controls",
            collapsed=True
        )

        hist_controls = pn.Card(
            pn.Column(
                pn.Row(self.hist_column_select, self.hist_bins_slider),
                pn.Row(self.selected_alpha_slider, self.unselected_alpha_slider)
            ),
            title="Histogram Controls",
            collapsed=False
        )

        mutation_analysis_controls = pn.Card(
            pn.Column(
                pn.Row(self.NTorAA_toggle, self.plot_type_toggle),  # Radio buttons on top
                pn.Row(self.freq_select, self.axis_type_toggle),    # Position mode and axis type
                pn.Row(self.num_positions_text, self.range_input, self.plot_width_slider)  # Number positions or range and plot width
            ),
            title="Aggregated Mutations Controls",
            collapsed=False
        )

        table_controls = pn.Card(
            self.sample_size_slider,
            title="Data Table Controls",
            collapsed=False
        )

        export_controls = pn.Card(
            pn.Column(
                pn.Row(self.selection_name_text, self.export_plots_button),
                pn.Row(self.export_genotypes_button, self.export_agg_muts_button),
                pn.Row(self.export_NT_fasta_button, self.export_AA_fasta_button,
                       self.export_NT_consensus_button, self.export_AA_consensus_button)
            ),
            title="Export Controls",
            collapsed=True
        )

        # Structure viewer controls (only show if structure component exists)
        control_columns = [
            data_controls,
            viz_controls,
            muts_of_interest_controls,
            hist_controls,
            mutation_analysis_controls,
            table_controls
        ]

        if self.structure is not None:
            structure_controls = pn.Card(
                pn.Row(
                    self.structure_label_step_select,
                    self.structure_frequency_floor_slider
                ),
                title="Structure Viewer Controls",
                collapsed=False
            )
            control_columns.append(structure_controls)

        control_columns.append(export_controls)

        return pn.Column(
            *control_columns,
            width=500,  # Increased width as requested
            sizing_mode='stretch_height'
        )


class MapleDashboard:
    """Main dashboard class orchestrating all components"""

    def __init__(self, datasets_config, initial_downsample=100000):
        """
        Initialize dashboard with unified datasets_config.

        datasets_config: list of dicts with keys: name, genotypes_file, reference_file, exclude_indels
        Note: structure_file is now obtained from the reference CSV file
        """
        self.datasets_config = datasets_config
        dataset_names = [d['name'] for d in datasets_config]

        # Create dataset selector (disabled if only one option)
        self.dataset_selector = pn.widgets.Select(
            name='Dataset',
            options=dataset_names,
            value=dataset_names[0],
            width=200,
            disabled=(len(dataset_names) == 1)
        )

        # Load first dataset
        first_dataset = datasets_config[0]
        self.data_manager = DataManager(
            first_dataset['genotypes_file'],
            first_dataset['reference_file'],
            first_dataset['exclude_indels'],
            downsample_size=initial_downsample
        )
        self.genotypes_filename = pathlib.Path(first_dataset['genotypes_file']).name

        # Create reference selector (populated from loaded data)
        reference_names = self._get_reference_names()
        self.reference_selector = pn.widgets.Select(
            name='Reference',
            options=reference_names,
            value=reference_names[0],
            width=200,
            disabled=(len(reference_names) <= 1)
        )

        # Cache for saved selection
        self.saved_selection = None

        # Create downsample slider with command line argument support
        max_size = len(self.data_manager.genotypes)
        initial_value = min(initial_downsample, max_size)

        self.downsample_slider = pn.widgets.IntSlider(
            name='Downsample',
            start=1000,
            end=max_size,
            step=1000,
            value=initial_value,
            width=200
        )

        # Create components first
        self.scatter_component = ScatterPlotComponent(self.data_manager)
        self.histogram_component = HistogramComponent(self.data_manager)
        self.mutation_analysis_component = MutationAnalysisComponent(self.data_manager)
        self.table_component = TableComponent(self.data_manager)
        self.structure_component = StructureHeatmapComponent(self.data_manager)
        self.controls_panel = ControlsPanel(
            self.data_manager,
            self.scatter_component,
            self.histogram_component,
            self.mutation_analysis_component,
            self.table_component,
            self.downsample_slider,
            structure_component=self.structure_component,
            dataset_selector=self.dataset_selector,
            reference_selector=self.reference_selector
        )

        # Reactive data pipeline with filtering (already downsampled at load time)
        self.filtered_data = pn.rx(self._get_filtered_data)(
            self.reference_selector.param.value,
            self.controls_panel.filter1_column_select.param.value,
            self.controls_panel.filter1_range_slider.param.value,
            self.controls_panel.filter2_column_select.param.value,
            self.controls_panel.filter2_range_slider.param.value,
            self.controls_panel.NT_muts_text.param.value,
            self.controls_panel.AA_muts_text.param.value,
            self.controls_panel.max_mut_combos_slider.param.value,
            self.controls_panel.NT_muts_multichoice.param.value,
            self.controls_panel.AA_muts_multichoice.param.value
        )

        # Create shared selection state like original dashboard
        self.ls = link_selections.instance()
        self.selected_indices = []  # Shared selection state for all plots to use

        # Create selection control buttons
        self.clear_selection_button = pn.widgets.Button(
            name='Clear Selection',
            button_type='primary',
            width=133
        )
        self.clear_selection_button.on_click(self._clear_selection)

        self.save_selection_button = pn.widgets.Button(
            name='Save Selection',
            button_type='primary',
            width=133
        )
        self.save_selection_button.on_click(self._save_selection)

        self.load_selection_button = pn.widgets.Button(
            name='Load Selection',
            button_type='primary',
            width=133
        )
        self.load_selection_button.on_click(self._load_selection)

        # Create reactive plots with shared selection (pass link as parameter like original)
        self.reactive_plot = pn.bind(self._create_scatter_plot_with_selection,
                                   link=self.ls,  # Pass link_selections as parameter
                                   data=self.filtered_data,
                                   x_col=self.scatter_component.param.x_column,
                                   y_col=self.scatter_component.param.y_column,
                                   color_by=self.scatter_component.param.color_by,
                                   colormap=self.scatter_component.param.colormap,
                                   datashade=self.scatter_component.param.datashade,
                                   point_size=self.scatter_component.param.point_size,
                                   alpha=self.scatter_component.param.alpha,
                                   # Make sure opacity changes trigger plot updates
                                   selected_alpha=self.controls_panel.selected_alpha_slider.param.value,
                                   unselected_alpha=self.controls_panel.unselected_alpha_slider.param.value)

        # Create reactive mutation analysis that updates with selections and filtered data
        self.reactive_mutation_analysis = pn.bind(self._create_mutation_analysis_with_selection,
                                                data=self.filtered_data,
                                                selection_expr=self.ls.param.selection_expr,
                                                mutation_component=self.mutation_analysis_component)

        # Create conditional legend for bar charts
        self.reactive_legend = pn.bind(self._create_mutation_legend,
                                     plot_type=self.mutation_analysis_component.param.plot_type)

        # Create reactive histogram with selection overlay
        self.reactive_histogram = pn.bind(self._create_histogram_with_selection,
                                        data=self.filtered_data,
                                        hist_column=self.histogram_component.param.hist_column,
                                        num_bins=self.histogram_component.param.num_bins,
                                        selected_alpha=self.controls_panel.selected_alpha_slider.param.value,
                                        unselected_alpha=self.controls_panel.unselected_alpha_slider.param.value,
                                        selection_expr=self.ls.param.selection_expr)

        # Create colorbar/legend for datashaded plots
        self.reactive_colorbar = pn.bind(self._create_color_legend,
                                       data=self.filtered_data,
                                       selection_expr=self.ls.param.selection_expr,
                                       color_by=self.scatter_component.param.color_by,
                                       colormap=self.scatter_component.param.colormap,
                                       datashade=self.scatter_component.param.datashade,
                                       NT_muts=self.controls_panel.NT_muts_text.param.value,
                                       AA_muts=self.controls_panel.AA_muts_text.param.value,
                                       max_groups=self.controls_panel.max_mut_combos_slider.param.value)

        # Create reactive table that updates with selections and filtered data
        self.reactive_table = pn.bind(self._create_table_with_selection,
                                    data=self.filtered_data,
                                    selection_expr=self.ls.param.selection_expr,
                                    table_component=self.table_component)

        # Create reactive structure viewer that updates with selections and filtered data (only if PDB file provided)
        # Note: Structure file availability is checked dynamically in StructureHeatmapComponent
        self.reactive_structure = pn.bind(self._create_structure_with_selection,
                                        data=self.filtered_data,
                                        selection_expr=self.ls.param.selection_expr,
                                        structure_component=self.structure_component,
                                        position_range=self.mutation_analysis_component.param.position_range,
                                        colormap=self.structure_component.param.colormap,
                                        label_step=self.structure_component.param.label_step,
                                        frequency_floor=self.structure_component.param.frequency_floor)

        # Link export plots button callback (button is in ControlsPanel, logic is here)
        self.controls_panel.export_plots_button.on_click(self._export_plots)

        # Link mutations of interest callbacks to update MultiChoice options
        self.controls_panel.NT_muts_text.param.watch(self._update_NT_muts_options, 'value')
        self.controls_panel.AA_muts_text.param.watch(self._update_AA_muts_options, 'value')
        self.controls_panel.max_mut_combos_slider.param.watch(self._update_NT_muts_options, 'value')
        self.controls_panel.max_mut_combos_slider.param.watch(self._update_AA_muts_options, 'value')

        # Link dataset selector callback (only active if more than one dataset)
        if len(datasets_config) > 1:
            self.dataset_selector.param.watch(self._switch_dataset, 'value')

        # Link reference selector callback to same method (reference switch is same as dataset switch)
        self.reference_selector.param.watch(self._switch_dataset, 'value')

        # Create the dashboard template
        self._create_template()

    def _get_reference_names(self):
        """Get list of reference names sorted by abundance from current data_manager"""
        if 'reference_name' in self.data_manager.genotypes.columns:
            ref_abundance = self.data_manager.genotypes.groupby('reference_name')['count'].sum().sort_values(ascending=False)
            return ref_abundance.index.tolist()
        else:
            return ['N/A']

    def _get_filtered_data(self, reference_name, filter1_col, filter1_range, filter2_col, filter2_range, NT_muts, AA_muts, max_groups, NT_muts_filter, AA_muts_filter):
        """Get filtered data (already downsampled at load time)"""
        # Clear cache when filters change (data subset changes)
        self.data_manager.clear_agg_muts_cache()

        # Start with analyzer genotypes (already downsampled)
        df = self.data_manager.analyzer.genotypes.copy()

        # Apply reference filter if reference_name column exists and not N/A
        if 'reference_name' in df.columns and reference_name != 'N/A':
            df = df[df['reference_name'] == reference_name]

        # Apply filter 1 if column and range are valid
        if filter1_col and filter1_col in df.columns and filter1_range:
            df = df[
                (df[filter1_col] >= filter1_range[0]) &
                (df[filter1_col] <= filter1_range[1])
            ]

        # Apply filter 2 if column and range are valid
        if filter2_col and filter2_col in df.columns and filter2_range:
            df = df[
                (df[filter2_col] >= filter2_range[0]) &
                (df[filter2_col] <= filter2_range[1])
            ]

        # Add mutations of interest columns
        df.loc[df.index, 'NT_muts_of_interest'] = self.data_manager.analyzer.get_mutations_of_interest('NT', NT_muts, max_groups, idx=df.index)
        if self.data_manager.do_AA_analysis:
            df.loc[df.index, 'AA_muts_of_interest'] = self.data_manager.analyzer.get_mutations_of_interest('AA', AA_muts, max_groups, idx=df.index)

        # Apply mutations of interest filter
        if NT_muts_filter and len(NT_muts_filter) > 0:
            df = df[df['NT_muts_of_interest'].isin(NT_muts_filter)]

        if self.data_manager.do_AA_analysis and AA_muts_filter and len(AA_muts_filter) > 0:
            df = df[df['AA_muts_of_interest'].isin(AA_muts_filter)]

        # Add index as a column for programmatic selection expressions
        df['index'] = df.index

        return df

    def _get_categorical_color_mapping(self, data, color_by, colormap):
        """
        Helper method to generate categorical color mapping.
        Used by both scatter plot and colorbar/legend creation.

        Returns: (categories, color_key) tuple
        """
        categories = sorted(data[color_by].unique())

        # Get colormap colors
        colormap_dict = cmap_dict()
        cmap_colors = colormap_dict[colormap]

        # Create color mapping
        num_categories = len(categories)
        indices = np.linspace(0, len(cmap_colors) - 1, num_categories, dtype=int)
        subset_colors = [cmap_colors[i] for i in indices]
        color_key = {cat: color for cat, color in zip(categories, subset_colors)}

        return categories, color_key

    def _clear_selection(self, event):
        """Clear all selections manually and reset state"""
        # Set to empty list instead of None to properly clear visual selection
        self.ls.selection_expr = []
        self.selected_indices = []
        # Clear the analyzer's selected_idx to show all data in mutation analysis
        self.data_manager.analyzer.selected_idx = None
        # Clear cache when selection changes
        self.data_manager.clear_agg_muts_cache()
        print("Selection manually cleared")

    def _save_selection(self, event):
        """Save current selection to cache and pickle file"""
        selection_expr = self.ls.selection_expr

        if selection_expr is None or selection_expr == []:
            print("No selection to save")
            return

        # Cache the selection expression
        self.saved_selection = selection_expr

        # Save to pickle file
        pathlib.Path('dashboard').mkdir(parents=True, exist_ok=True)
        save_data = {
            'filename': self.genotypes_filename,
            'selection_expr': selection_expr
        }

        with open('dashboard/selection.pkl', 'wb') as f:
            pickle.dump(save_data, f)

        # Get count for user feedback
        selected_idx = self.data_manager.analyzer.selected_idx
        count = len(selected_idx) if selected_idx else 0
        print(f"Saved selection with {count} sequences")

    def _load_selection(self, event):
        """Load selection from cache or pickle file"""
        # Try cache first
        if self.saved_selection is not None:
            print("Loading selection from cache")
            self._apply_saved_selection(self.saved_selection)
            return

        # Try pickle file if no cache
        selection_file = pathlib.Path('dashboard/selection.pkl')
        if not selection_file.exists():
            print("No saved selection found (no cache or file)")
            return

        # Load and validate pickle
        try:
            with open(selection_file, 'rb') as f:
                save_data = pickle.load(f)

            # Validate filename
            saved_filename = save_data.get('filename')
            if saved_filename != self.genotypes_filename:
                print(f"WARNING: Saved selection is from different file!")
                print(f"  Current: {self.genotypes_filename}")
                print(f"  Saved:   {saved_filename}")
                print("Selection NOT loaded")
                return

            # Load selection expression
            selection_expr = save_data.get('selection_expr')
            if selection_expr is not None:
                print("Loading selection from file")
                self._apply_saved_selection(selection_expr)
            else:
                print("Invalid selection file (no selection_expr)")

        except Exception as e:
            print(f"Error loading selection file: {e}")

    def _apply_saved_selection(self, selection_expr):
        """Apply a saved selection expression to the dashboard"""
        # Update link_selections with the saved expression
        self.ls.selection_expr = selection_expr
        print("Selection applied")

    def _switch_dataset(self, event):
        """Switch to a different dataset or reference"""
        # If triggered by reference selector, use current dataset. Otherwise use new dataset.
        if event.obj == self.reference_selector:
            dataset_name = self.dataset_selector.value
        else:
            dataset_name = event.new
        dataset_config = next(d for d in self.datasets_config if d['name'] == dataset_name)

        # Clear selection state BEFORE switching data manager
        self.selected_indices = []
        if self.data_manager.analyzer.selected_idx is not None:
            self.data_manager.analyzer.selected_idx = None

        # Clear cache from old data manager before switching
        self.data_manager.clear_agg_muts_cache()

        # Clear selection expr to prevent reactive updates during switch
        self.ls.selection_expr = []

        # Create new data manager with new reference (old data_manager will be garbage collected)
        current_downsample = self.downsample_slider.value
        self.data_manager = DataManager(
            dataset_config['genotypes_file'],
            dataset_config['reference_file'],
            dataset_config['exclude_indels'],
            downsample_size=current_downsample
        )
        self.genotypes_filename = pathlib.Path(dataset_config['genotypes_file']).name

        # Update data_manager reference in all components
        self.scatter_component.data_manager = self.data_manager
        self.histogram_component.data_manager = self.data_manager
        self.mutation_analysis_component.data_manager = self.data_manager
        self.table_component.data_manager = self.data_manager
        self.structure_component.data_manager = self.data_manager
        self.controls_panel.data_manager = self.data_manager

        # Update downsample slider range
        max_size = len(self.data_manager.genotypes)
        self.downsample_slider.end = max_size
        self.downsample_slider.value = min(current_downsample, max_size)

        # Update reference selector for new dataset
        reference_names = self._get_reference_names()
        self.reference_selector.options = reference_names
        self.reference_selector.value = reference_names[0]
        self.reference_selector.disabled = (len(reference_names) <= 1)

        # Update filter ranges for new dataset (do this last to trigger reactive update with fresh data)
        self.controls_panel._update_filter_ranges()

    def _update_NT_muts_options(self, event):
        """Update NT mutations MultiChoice options based on current filtered data"""
        current_data = self.filtered_data.rx.value
        if 'NT_muts_of_interest' in current_data.columns:
            unique_muts = current_data['NT_muts_of_interest'].unique()
            options = [m for m in unique_muts if m != 'none']
            self.controls_panel.NT_muts_multichoice.options = sorted(options)

    def _update_AA_muts_options(self, event):
        """Update AA mutations MultiChoice options based on current filtered data"""
        current_data = self.filtered_data.rx.value
        if self.data_manager.do_AA_analysis and 'AA_muts_of_interest' in current_data.columns:
            unique_muts = current_data['AA_muts_of_interest'].unique()
            options = [m for m in unique_muts if m != 'none']
            self.controls_panel.AA_muts_multichoice.options = sorted(options)


    def _create_mutation_analysis_with_selection(self, data, selection_expr, mutation_component):
        """Update mutation analysis component based on current selection and filtered data"""
        # Store current data for process_selection to use
        self.current_data = data
        # Process selection using centralized helper (updates analyzer.selected_idx)
        process_selection(selection_expr, self)

        # Return the mutation component's plot method (not called) so Panel can track dependencies
        return mutation_component.plot

    def _create_mutation_legend(self, plot_type):
        """Create legend image for bar charts, spacer for heatmaps"""
        if plot_type == 'bars':
            # Show legend image for bar charts
            legend_path = os.path.join(os.path.dirname(__file__), '..', '..', 'images', 'AA_NT_colormap.png')
            if os.path.exists(legend_path):
                return pn.pane.PNG(legend_path, width=150, height=400)
            else:
                return pn.pane.Markdown(f"**Legend image not found**\n\nExpected location:\n`maple/images/AA_NT_colormap.png`", width=150, height=400)
        else:
            # Show spacer for heatmaps (they don't need legends)
            return pn.Spacer(width=0)

    def _create_table_with_selection(self, data, selection_expr, table_component):
        """Update table component based on current selection and filtered data"""
        # Store current data for process_selection to use
        self.current_data = data
        # Process selection using centralized helper (updates analyzer.selected_idx)
        process_selection(selection_expr, self)

        # Return the table component's table (which will now use the updated selected_idx)
        return table_component.table

    def _create_structure_with_selection(self, data, selection_expr, structure_component, position_range, colormap, label_step, frequency_floor):
        """Update structure viewer component based on current selection and filtered data"""
        # Store current data for process_selection to use
        self.current_data = data
        # Process selection using centralized helper (updates analyzer.selected_idx)
        process_selection(selection_expr, self)

        # Call plot with position_range argument
        # colormap, label_step, frequency_floor are in the signature to trigger updates when they change
        return structure_component.plot(position_range=position_range)


    def _create_color_legend(self, data, selection_expr, color_by, colormap, datashade, NT_muts=None, AA_muts=None, max_groups=None):
        """Create colorbar for numerical data or legend for categorical data"""
        # NT_muts, AA_muts, max_groups are just to trigger updates when they change
        if not datashade:
            return pn.Spacer(width=0)  # No colorbar for non-datashaded plots

        # Process selection to get the data we should use for the colorbar
        # Store current data for process_selection to use
        self.current_data = data
        selected_indices = process_selection(selection_expr, self)

        # Determine which data to use: selected data if there's a selection, otherwise filtered data
        if selected_indices is not None and len(selected_indices) > 0:
            # Use selected data (subset of filtered data)
            display_data = data.loc[data.index.isin(selected_indices)]
        else:
            # Use all filtered data
            display_data = data

        # Check if column is categorical
        is_categorical = not pd.api.types.is_numeric_dtype(self.data_manager.genotypes[color_by].dtype)

        if is_categorical:
            # Create HoloViews-based legend with colored circles (exportable to SVG)
            if display_data is None or color_by not in display_data.columns:
                return pn.Spacer(width=0)

            # Use helper method to get color mapping (same logic as scatter plot)
            categories, color_key = self._get_categorical_color_mapping(display_data, color_by, colormap)

            # Create HoloViews elements for each category
            elements = []
            for i, (cat, color) in enumerate(color_key.items()):
                # Create point (circle) with no border
                point = hv.Points([(0, i)]).opts(
                    color=color,
                    size=10,
                    line_color=None,  # Remove border
                    line_width=0
                )
                elements.append(point)

                # Create text label with correct font size
                label = hv.Text(0.3, i, str(cat)).opts(
                    text_font_size=f'{PLOT_FONT_SIZES["legend"]}pt',  # Match colorbar label size (13pt)
                    text_align='left',
                    text_baseline='middle'
                )
                elements.append(label)

            # Calculate height: 33px per line of text (excluding title)
            legend_height = max(33 * len(categories), 100)  # Minimum 100px

            # Combine all elements into overlay
            legend = hv.Overlay(elements).opts(
                width=200,
                height=legend_height,
                title=color_by,
                xlim=(-0.2, 2.5),
                ylim=(-0.5, len(categories) - 0.5),
                xaxis=None,
                yaxis=None,
                show_frame=False,
                toolbar=None,  # Disable pan tool
                fontsize={'title': PLOT_FONT_SIZES['clabel']}  # Match colorbar title size (13pt)
            )

            return legend

        # Create colorbar using Bokeh for all numerical columns (including count)
        # Use the color range from the display data (filtered/selected)
        color_min = display_data[color_by].min()
        color_max = display_data[color_by].max()

        # Use the same colormap as the plot like original dashboard
        colormap_dict = cmap_dict()
        cmap_colors = colormap_dict[colormap]
        color_mapper = LinearColorMapper(palette=cmap_colors, low=color_min, high=color_max)

        color_bar = ColorBar(
            color_mapper=color_mapper,
            label_standoff=8,
            major_label_text_font_size=f"{PLOT_FONT_SIZES['yticks']}px"
        )

        # Wrap in Bokeh figure like original
        p = figure(width=100, height=550, toolbar_location=None, min_border=0)
        p.add_layout(color_bar, 'right')
        p.axis.visible = False
        p.xgrid.grid_line_color = None
        p.ygrid.grid_line_color = None
        p.outline_line_color = None
        p.background_fill_color = None
        p.border_fill_color = None

        # Add label for column name with consistent font size (moved 10px left)
        label = Label(
            text=color_by,
            x=10,  # 20 - 10 = 10 pixels from left
            y=20,
            x_units="screen",
            y_units="screen",
            angle=90,
            angle_units="deg",
            text_align="left",
            text_baseline="middle",
            text_font_size=f"{PLOT_FONT_SIZES['clabel']}pt"
        )
        p.add_layout(label)

        return pn.panel(p)

    def _update_selection_state(self, event):
        """Update shared selection state when any plot selection changes"""
        self.selected_indices = event.new if event.new else []

    def _sync_selection_to_ls(self, event):
        """Sync non-datashaded selection to ls.selection_expr"""
        if event.new:
            # Convert index selection to selection_expr format
            self.ls.selection_expr = event.new
        else:
            # Clear selection properly to refresh histogram
            self.ls.selection_expr = None
            print("Selection cleared")

    def _create_scatter_plot_with_selection(self, link, data, x_col, y_col, color_by, colormap, datashade, point_size, alpha, selected_alpha, unselected_alpha):
        """Create scatter plot with selection support - different approaches for datashaded vs non-datashaded"""

        if datashade:
            # For datashaded: use link_selections approach
            tools = ['box_select', 'lasso_select', 'reset']

            # Check if categorical using full dataset (not filtered)
            is_categorical = not pd.api.types.is_numeric_dtype(self.data_manager.genotypes[color_by].dtype)

            if is_categorical:
                # For categorical: no color range
                color_range = None
            else:
                # For numerical: use filtered data min/max to update with filter changes
                color_range = (data[color_by].min(), data[color_by].max())

            # Store this range for the colorbar
            self.current_color_range = color_range

            # Create basic scatter plot first (without datashading)
            scatter = data.hvplot.scatter(
                x=x_col,
                y=y_col,
                c=color_by,
                width=800,
                height=600,
                tools=tools
            )

            # Apply datashader with appropriate settings
            # Get colormap colors from cmap_dict like original dashboard
            colormap_dict = cmap_dict()
            cmap_colors = colormap_dict[colormap]

            if is_categorical:
                # For categorical: use ds.by aggregator with color_key
                categories, color_key = self._get_categorical_color_mapping(data, color_by, colormap)

                scatter = hv_datashade(scatter,
                                  aggregator=ds.by(color_by),
                                  color_key=color_key,
                                  cnorm='eq_hist',
                                  min_alpha=255)  # Prevents blending for categorical
            else:
                # For numerical: use mean aggregator with fixed clims
                scatter = hv_datashade(scatter,
                                  aggregator=ds.reductions.mean(color_by),
                                  cmap=cmap_colors,
                                  cnorm='linear',
                                  clims=color_range)

            # Apply dynspread with point_size controlling spread amount
            # Convert point_size (1-20) to reasonable max_px range (1-6)
            max_px = max(1, min(6, int(point_size / 3)))
            scatter = dynspread(scatter, threshold=1, max_px=max_px)

            # Apply sizing options that were lost when we manually applied datashade
            scatter = scatter.opts(width=730, height=600, fontsize=PLOT_FONT_SIZES)  # Reduced datashaded plot by 70

            # Apply link_selections inside function like original
            plot = link(scatter)
            self.current_data = data
            return plot

        else:
            # For non-datashaded: also use link_selections approach for consistency
            tools = ['box_select', 'lasso_select']  # Remove reset from tools like original

            # Check if categorical
            is_categorical = not pd.api.types.is_numeric_dtype(self.data_manager.genotypes[color_by].dtype)

            # Get colormap colors from cmap_dict using the widget parameter
            colormap_dict = cmap_dict()
            cmap_colors = colormap_dict[colormap]

            if is_categorical:
                # For categorical: try using points instead of scatter
                categories, color_map = self._get_categorical_color_mapping(data, color_by, colormap)

                scatter = data.hvplot(
                    kind='points',
                    x=x_col,
                    y=y_col,
                    c=color_by,
                    s=point_size,
                    width=838,
                    height=600,
                    tools=tools,
                    cmap=color_map,  # Pass the category -> color mapping
                    legend='right',
                    fontsize=PLOT_FONT_SIZES
                ).opts(
                    legend_opts={
                        'title': color_by,
                        'title_text_font_size': f'{PLOT_FONT_SIZES["clabel"]}pt'  # Match colorbar title size (13pt)
                    }
                )
            else:
                # For numerical: use c= parameter for continuous color mapping
                scatter = data.hvplot(
                    kind='points',
                    x=x_col,
                    y=y_col,
                    c=color_by,
                    s=point_size,
                    width=838,  # Increased by 38px as requested
                    height=600,
                    tools=tools,
                    cmap=cmap_colors,  # Use selected colormap from widget
                    colorbar=True,  # Add colorbar
                    fontsize=PLOT_FONT_SIZES  # Add standardized font sizes
                ).opts(
                    colorbar_opts={
                        'title': color_by,
                        'title_text_font_size': f'{PLOT_FONT_SIZES["legend"]}px',
                        'major_label_text_font_size': f'{PLOT_FONT_SIZES["legend"]}px'
                    }
                )
            # NOTE: Opacity controls don't work for hvplot scatter plots

            # Apply link_selections inside function like original
            plot = link(scatter)
            self.current_data = data
            return plot

    def _create_histogram_with_selection(self, data, hist_column, num_bins, selected_alpha, unselected_alpha, selection_expr):
        """Create histogram with selection overlay using modern patterns"""
        # Calculate bin range based on column type
        if data[hist_column].dtype == 'int64':
            max_x_val = data[hist_column].max()
            min_x_val = 0
        else:  # float64
            max_x_val = data[hist_column].max()
            min_x_val = data[hist_column].min()

        bins = np.linspace(min_x_val - 0.5, max_x_val + 0.5, num_bins + 1)

        # Static histogram - opacity depends on whether there's a selection
        hist_alpha = unselected_alpha if (selection_expr is not None and selection_expr != []) else selected_alpha

        static_hist = create_histogram_plot(data, hist_column, bins, hist_alpha)

        # Create selected data histogram using centralized helper
        indices = process_selection(selection_expr, self)
        if indices is not None:
            # Use loc instead of iloc since indices are absolute (not positional within filtered data)
            # Only select indices that exist in the current filtered data
            valid_indices = [idx for idx in indices if idx in data.index]
            if len(valid_indices) > 0:
                selected_data = data.loc[valid_indices]
                selected_hist = create_histogram_plot(selected_data, hist_column, bins, selected_alpha)
                return static_hist * selected_hist

        # Return just static histogram if no selection
        return static_hist

    def _export_plots(self, event):
        """Export current plots to SVG"""
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        selection_name = self.controls_panel.selection_name_text.value
        filename = f'dashboard/{selection_name}_{timestamp}_plots.html'

        scatter_plot = self.reactive_plot()
        histogram_plot = self.reactive_histogram()
        mutation_plot_method = self.reactive_mutation_analysis()
        mutation_plot = mutation_plot_method()
        colorbar = self.reactive_colorbar()

        export_svg_plots(
            [scatter_plot, mutation_plot, histogram_plot, colorbar],
            filename,
            ['scatter', 'mutations', 'histogram', 'colorbar']
        )
        print(f'Exported plots to {filename[:-5]}_svg/ directory')

    def _create_template(self):
        """Create the main dashboard template"""

        # Quick start panel in collapsible card
        quick_start_panel = pn.Card(
            pn.pane.Markdown("""
        - Use **datashading** for large datasets (>10k points)
        - Select **X/Y columns** to change the view
        - **Color by** different mutation counts
        - Use **box select** or **lasso select** to highlight points
        - Refreshing the page or clicking the header resets the dashboard
        - Note that the datashaded plot might look weird if your screen has an unusual aspect ratio
        """),
            title="Quick Start",
            collapsed=True
        )

        # Main plot panel with table to the right and colorbar (like original dashboard layout)
        plot_panel = pn.Row(
            self.reactive_plot,
            self.reactive_colorbar,
            self.reactive_table,
            sizing_mode='stretch_width'
        )


        # Create Selection Controls card
        selection_controls_card = pn.Card(
            pn.Row(
                self.save_selection_button,
                self.load_selection_button,
                self.clear_selection_button
            ),
            title="Selection Controls",
            collapsed=False
        )

        # Sidebar items (dataset selector is now in controls panel)
        sidebar_items = [quick_start_panel, self.controls_panel, selection_controls_card]

        # Create template with wider sidebar
        self.template = pn.template.FastListTemplate(
            title="",  # Remove main title
            sidebar=sidebar_items,
            main=[
                pn.Column(
                    plot_panel,
                    pn.Row(
                        self.reactive_histogram,
                        self.reactive_mutation_analysis,
                        self.reactive_legend,
                        self.reactive_structure,
                        sizing_mode='stretch_width'
                    ),
                    sizing_mode='stretch_both'
                )
            ],
            header_background='#2596be',
            accent_base_color='#2596be',
            sidebar_width=450  # Set sidebar width directly on template
        )

    def show(self):
        """Display the dashboard"""
        self.template.show(autoreload=True)
        return self.template

    def servable(self, title=None):
        """Make the dashboard servable"""
        if title:
            self.template.title = title
        return self.template.servable()


# Safe argument parsing and dashboard creation
try:
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Maple Genotypes Dashboard - 2025 Rebuild')
    parser.add_argument('--genotypes', type=str,
                       help='Path to genotypes CSV file (single dataset mode)')
    parser.add_argument('--reference', type=str,
                       help='Path to reference FASTA file (single dataset mode)')
    parser.add_argument('--exclude_indels', action='store_true',
                       help='Exclude sequences with indels (single dataset mode)')
    parser.add_argument('--downsample', type=int, default=100000,
                       help='Initial downsample value (default: 100000)')
    parser.add_argument('--structure_file', type=str, default='',
                       help='Path to structure file (PDB/CIF) for structure viewer (optional, single dataset mode)')
    parser.add_argument('--datasets', type=str,
                       help='Path to datasets CSV file (multi-dataset mode)')

    args, unknown = parser.parse_known_args()

    print("Starting Maple Dashboard 2025...")

    # Build datasets_config from either CSV or command line args
    if args.datasets:
        print(f"Multi-dataset mode: {args.datasets}")
        # Load datasets CSV
        datasets_df = pd.read_csv(args.datasets)

        # Validate required columns
        required_cols = ['genotypes_csv', 'reference_csv']
        missing_cols = [col for col in required_cols if col not in datasets_df.columns]
        if missing_cols:
            raise ValueError(f"Datasets CSV missing required columns: {missing_cols}")

        # Create dataset configurations
        datasets_config = []
        for _, row in datasets_df.iterrows():
            # Extract dataset name from genotypes filename
            genotypes_path = row['genotypes_csv']
            dataset_name = pathlib.Path(genotypes_path).name.replace('_genotypes.csv', '')

            # Handle optional columns
            exclude_indels = False
            if 'exclude_indels' in row and pd.notna(row['exclude_indels']):
                exclude_indels = str(row['exclude_indels']).lower() in ['true', '1', 'yes']

            datasets_config.append({
                'name': dataset_name,
                'genotypes_file': genotypes_path,
                'reference_file': row['reference_csv'],
                'exclude_indels': exclude_indels
            })

    else:
        # Single dataset mode - convert to datasets_config format
        if not args.genotypes or not args.reference:
            raise ValueError("Single dataset mode requires --genotypes and --reference arguments")

        print(f"Single dataset mode")
        print(f"Genotypes: {args.genotypes}")
        print(f"Reference: {args.reference}")

        # Build single-item datasets_config
        # Note: structure_file now comes from reference CSV, not command line
        datasets_config = [{
            'name': pathlib.Path(args.genotypes).name.replace('_genotypes.csv', ''),
            'genotypes_file': args.genotypes,
            'reference_file': args.reference,
            'exclude_indels': args.exclude_indels
        }]

    # Create dashboard with unified datasets_config
    dashboard = MapleDashboard(
        datasets_config=datasets_config,
        initial_downsample=args.downsample
    )

    layout = dashboard.servable(title='Maple Dashboard')

except Exception as e:
    print(f"Error creating dashboard: {e}")
    traceback.print_exc()

    # Fallback - create error dashboard
    error_text = f"""
    # Dashboard Creation Failed

    **Error:** {e}

    **Arguments:** {getattr(locals(), 'args', 'Not parsed')}

    Please check your input files and try again.
    """

    layout = pn.Column(
        pn.pane.Markdown(error_text),
        pn.widgets.Button(name="Error occurred", button_type="primary", disabled=True)
    )

    layout.servable(title="Maple Dashboard - Error")