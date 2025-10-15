#!/usr/bin/env python3
"""
Test script to create and embed structure viewer in Panel
Combines structure viewer generation and Panel embedding in one script
"""

import sys
import os
import argparse
import base64

# Add parent directory to path to import modules
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'rules', 'utils'))

import pandas as pd
import panel as pn
import py3Dmol
from Bio.PDB import PDBParser
from Bio import Align

from common import cmap_dict

# AA 3-letter to 1-letter mapping
AA_3TO1 = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}


def extract_pdb_sequence(chain):
    """
    Extract amino acid sequence and positions from PDB chain

    Args:
        chain: BioPython chain object

    Returns:
        tuple: (sequence string, list of PDB position numbers)
    """
    residues = []
    positions = []

    for residue in chain:
        if residue.id[0] == ' ' and residue.resname in AA_3TO1:
            residues.append(AA_3TO1[residue.resname])
            positions.append(residue.id[1])

    return ''.join(residues), positions


def align_pdb_to_reference(pdb_file, reference_sequence, min_identity=0.95):
    """
    Align PDB structure sequence to reference sequence and create mapping

    Args:
        pdb_file: path to PDB file
        reference_sequence: reference AA sequence string
        min_identity: minimum sequence identity required (default 0.95)

    Returns:
        dict with keys:
            'chain_id': selected chain ID
            'pdb_sequence': extracted PDB sequence
            'pdb_positions': list of PDB residue numbers
            'mapping': dict mapping reference position -> (chain_id, pdb_position)
            'identity': sequence identity in aligned region

    Raises:
        ValueError: if no chain meets minimum identity threshold
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)

    best_chain = None
    best_identity = 0
    best_coverage = 0
    best_alignment = None
    best_pdb_seq = None
    best_pdb_positions = None
    best_aligned_start = 0
    best_aligned_end = 0


    # Try aligning each chain
    for chain in structure.get_chains():
        # Extract sequence from PDB
        pdb_sequence, pdb_positions = extract_pdb_sequence(chain)

        if not pdb_sequence:
            continue


        # Align sequences
        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'
        # Penalize gaps heavily to favor mismatches over gaps
        aligner.open_gap_score = -10
        aligner.extend_gap_score = -2
        aligner.match_score = 2
        aligner.mismatch_score = -1

        # Get just the best alignment (don't enumerate all optimal alignments)
        alignment = aligner.align(reference_sequence, pdb_sequence)[0]

        # Trim terminal gaps from the alignment (exclude positions where either sequence has only gaps at the ends)
        align_ref = str(alignment[0])
        align_pdb = str(alignment[1])

        # Find first position where both sequences have started (no leading gaps in either)
        start = 0
        for i in range(len(align_ref)):
            if align_ref[i] != '-' and align_pdb[i] != '-':
                start = i
                break

        # Find last position where both sequences are still going (no trailing gaps in either)
        end = len(align_ref)
        for i in range(len(align_ref) - 1, -1, -1):
            if align_ref[i] != '-' and align_pdb[i] != '-':
                end = i + 1
                break

        # Extract aligned region
        aligned_ref = align_ref[start:end]
        aligned_pdb = align_pdb[start:end]

        # Calculate identity in aligned region (gaps count as mismatches)
        matches = sum(1 for a, b in zip(aligned_ref, aligned_pdb) if a == b)
        aligned_length = len(aligned_ref)  # Total length including gaps
        identity = matches / aligned_length if aligned_length > 0 else 0

        # Calculate coverage (fraction of each sequence that has residues in aligned region)
        ref_residues_in_alignment = sum(1 for a in aligned_ref if a != '-')
        pdb_residues_in_alignment = sum(1 for b in aligned_pdb if b != '-')
        ref_coverage = ref_residues_in_alignment / len(reference_sequence) if len(reference_sequence) > 0 else 0
        pdb_coverage = pdb_residues_in_alignment / len(pdb_sequence) if len(pdb_sequence) > 0 else 0
        avg_coverage = (ref_coverage + pdb_coverage) / 2

        # Count positions with residues in both (for reference)
        both_have_residues = sum(1 for a, b in zip(aligned_ref, aligned_pdb) if a != '-' and b != '-')


        # Update if better: prioritize identity, then coverage
        if best_chain is None or identity > best_identity or (identity == best_identity and avg_coverage > best_coverage):
            best_identity = identity
            best_coverage = avg_coverage
            best_chain = chain
            best_alignment = alignment
            best_pdb_seq = pdb_sequence
            best_pdb_positions = pdb_positions
            best_aligned_start = start
            best_aligned_end = end

    if best_chain is None:
        raise ValueError(f"No valid chains found in PDB file")

    if best_identity < min_identity:
        raise ValueError(
            f"Best sequence identity ({best_identity:.1%}) is below minimum threshold ({min_identity:.1%}). "
            f"Chain {best_chain.id} matches reference with {best_identity:.1%} identity. "
            f"This suggests the PDB structure may not match your reference sequence."
        )

    # Create mapping from reference position to PDB position (only for aligned region)
    # Count how many reference positions are skipped before the aligned region starts
    full_align_ref = str(best_alignment[0])
    ref_positions_before_alignment = sum(1 for c in full_align_ref[:best_aligned_start] if c != '-')

    # Slice the alignment to only the aligned region
    align_ref = full_align_ref[best_aligned_start:best_aligned_end]
    align_pdb = str(best_alignment[1])[best_aligned_start:best_aligned_end]

    mapping = {}
    ref_idx = 0
    pdb_idx = 0

    for ref_char, pdb_char in zip(align_ref, align_pdb):
        if ref_char != '-' and pdb_char != '-':
            # Both aligned, create mapping (1-indexed reference positions)
            mapping[ref_positions_before_alignment + ref_idx + 1] = (best_chain.id, best_pdb_positions[pdb_idx])

        if ref_char != '-':
            ref_idx += 1
        if pdb_char != '-':
            pdb_idx += 1


    return {
        'chain_id': best_chain.id,
        'pdb_sequence': best_pdb_seq,
        'pdb_positions': best_pdb_positions,
        'mapping': mapping,
        'identity': best_identity
    }


def frequency_to_color(frequency, colormap='bmy_r', min_freq=0, max_freq=1, frequency_floor=0):
    """
    Map mutation frequency to hex color using colorcet colormaps from common.py

    Args:
        frequency: mutation frequency value
        colormap: colormap name (from common.cmap_dict())
        min_freq: minimum frequency for color scaling
        max_freq: maximum frequency for color scaling
        frequency_floor: minimum frequency threshold for coloring (default 0)

    Returns:
        hex color string
    """
    if frequency <= frequency_floor:
        return '#FFFFFF'  # White for mutations below threshold

    # Get colormap from common.py
    colormap_dict = cmap_dict()
    colors = colormap_dict.get(colormap, colormap_dict['bmy_r'])

    # Normalize frequency to 0-1 range
    if max_freq > min_freq:
        normalized = (frequency - min_freq) / (max_freq - min_freq)
    else:
        normalized = 0

    # Clamp to [0, 1]
    normalized = max(0, min(1, normalized))

    # Map to color index
    color_idx = int(normalized * (len(colors) - 1))

    return colors[color_idx]


def create_hover_callbacks(selected_chain, js_data):
    """
    Create JavaScript callbacks for hover interactions

    Args:
        selected_chain: chain ID to enable hover for
        js_data: dictionary of mutation data for JavaScript

    Returns:
        tuple: (hover_js, unhover_js) callback strings
    """
    hover_js = """
    function(atom, viewer) {
        // Remove previous hover label if it exists
        if(window.hoverLabel) {
            viewer.removeLabel(window.hoverLabel);
            window.hoverLabel = null;
        }

        if(!atom.chain) return;

        // For selected chain
        if(atom.chain === '%s' && atom.resi) {
            var position = String(atom.resi);
            var data = %s;

            if(data[position]) {
                var pdbAA = data[position].pdb_aa;
                var refAA = data[position].ref_aa;
                var refPos = data[position].ref_pos;
                var freq = data[position].freq;

                // Format: PDB_AA or PDB_AA(REF_AA) if mismatch, then position and frequency
                var labelText = pdbAA;
                if(refAA && pdbAA !== refAA) {
                    labelText += "(" + refAA + ")";
                }
                labelText += refPos + ": " + freq.toFixed(3);

                // Position label on CA (alpha carbon) instead of hovered atom
                window.hoverLabel = viewer.addLabel(labelText,
                    {backgroundColor: 'lightgray', fontColor:'black',
                     backgroundOpacity: 0.8, fontSize: 12},
                    {chain: atom.chain, resi: atom.resi, atom: 'CA'});
            }
        }
        // For non-selected chains, just show chain ID
        else if(atom.resi) {
            window.hoverLabel = viewer.addLabel("Chain " + atom.chain,
                {backgroundColor: 'lightgray', fontColor:'black',
                 backgroundOpacity: 0.8, fontSize: 12},
                {chain: atom.chain, resi: atom.resi, atom: 'CA'});
        }
    }
    """ % (selected_chain, str(js_data).replace("'", '"').replace('None', 'null'))

    unhover_js = """
    function(atom, viewer) {
        if(window.hoverLabel) {
            viewer.removeLabel(window.hoverLabel);
            window.hoverLabel = null;
        }
    }
    """

    return hover_js, unhover_js


def create_structure_heatmap(pdb_file, pdb_mutation_data, selected_chain, pdb_sequence, pdb_positions, colormap='bmy_r', width=800, height=600, label_step=None, frequency_floor=0):
    """
    Create py3Dmol structure viewer with mutation frequency coloring

    Args:
        pdb_file: path to PDB file
        pdb_mutation_data: dict mapping PDB position -> {'frequency': float, 'wt': str, 'ref_position': int}
        selected_chain: chain ID to color (others shown but not colored)
        pdb_sequence: full PDB sequence string (for labeling all residues)
        pdb_positions: list of PDB position numbers (for labeling all residues)
        colormap: colormap name
        width: viewer width in pixels
        height: viewer height in pixels
        label_step: if int, label every Nth position (e.g., 10 = label every 10th position)
        frequency_floor: minimum frequency threshold for coloring (default 0)

    Returns:
        py3Dmol view object
    """
    # Load PDB structure
    with open(pdb_file, 'r') as f:
        pdb_data = f.read()

    # Create viewer
    view = py3Dmol.view(width=width, height=height)
    view.addModel(pdb_data, 'pdb')

    # Color by mutation frequency
    max_freq = max(data['frequency'] for data in pdb_mutation_data.values()) if pdb_mutation_data else 1.0
    min_freq = 0  # Always start from 0

    # Set style for all chains to white cartoon with opacity 0.5 (translucent)
    view.setStyle({}, {'cartoon': {'color': '#FFFFFF', 'opacity': 0.5}})

    # Make selected chain opaque white (will be overridden for colored residues below)
    view.setStyle({'chain': selected_chain}, {'cartoon': {'color': '#FFFFFF', 'opacity': 1.0}})

    # Color selected chain residues by mutation frequency (keep them opaque)
    for pdb_pos, data in pdb_mutation_data.items():
        # Get color using colormap from common.py
        color = frequency_to_color(data['frequency'], colormap=colormap, min_freq=min_freq, max_freq=max_freq, frequency_floor=frequency_floor)

        # Style this residue with color (opacity already set to 1.0 above)
        view.setStyle({'chain': selected_chain, 'resi': pdb_pos}, {'cartoon': {'color': color}})

    # Build hover data for ALL residues in the PDB chain
    js_data = {}
    for pdb_pos, aa in zip(pdb_positions, pdb_sequence):
        if pdb_pos in pdb_mutation_data:
            # Has mapping to reference
            data = pdb_mutation_data[pdb_pos]
            js_data[str(pdb_pos)] = {
                'pdb_aa': aa,  # Actual PDB residue
                'ref_aa': data['wt'],  # Reference residue
                'ref_pos': data['ref_position'],
                'freq': data['frequency']
            }
        else:
            # No mapping to reference - use PDB data with 0% frequency
            js_data[str(pdb_pos)] = {
                'pdb_aa': aa,
                'ref_aa': None,
                'ref_pos': pdb_pos,
                'freq': 0.0
            }

    hover_js, unhover_js = create_hover_callbacks(selected_chain, js_data)
    view.setHoverable({}, True, hover_js, unhover_js)

    # Add persistent labels for every Nth position if label_step is specified
    if label_step is not None and isinstance(label_step, int) and label_step > 0:
        # Find the last mapped position
        max_ref_pos = max(data['ref_position'] for data in pdb_mutation_data.values())
        max_ref_pdb_pos = None
        for pdb_pos, data in pdb_mutation_data.items():
            if data['ref_position'] == max_ref_pos:
                max_ref_pdb_pos = pdb_pos
                break

        for pdb_pos, data in pdb_mutation_data.items():
            ref_pos = data['ref_position']
            if ref_pos == 1 or ref_pos % label_step == 0 or pdb_pos == max_ref_pdb_pos:
                wt = data['wt']
                label_text = f"{wt}{ref_pos}"
                view.addLabel(
                    label_text,
                    {'fontSize': 12, 'fontColor': 'black', 'backgroundColor': 'white',
                     'showBackground': True, 'backgroundOpacity': 0.7},
                    {'chain': selected_chain, 'resi': pdb_pos, 'atom': 'CA'}
                )

    # Add some visual enhancements
    view.setBackgroundColor('white')
    # Zoom to selected chain only (center view on it)
    view.zoomTo({'chain': selected_chain})

    return view


def create_structure_pane(pdb_file, agg_df, colormap='bmy_r', width=800, height=600, label_step=None, min_identity=0.95, frequency_floor=0, position_range=''):
    """
    Create structure viewer Panel pane for embedding in dashboard

    Args:
        pdb_file: path to PDB file
        agg_df: aggregated mutations dataframe (from SequenceAnalyzer.aggregate_mutations)
        colormap: colormap name (default: 'bmy_r')
        width: viewer width in pixels (default: 800)
        height: viewer height in pixels (default: 600)
        label_step: if int, label every Nth position with WT_position (e.g., 10 = label every 10th position)
        min_identity: minimum sequence identity required for PDB alignment (default: 0.95)
        frequency_floor: minimum frequency threshold for coloring residues (default: 0)
        position_range: position range(s) to color, e.g., '20-30' or '20-30,50-60' (default: '' = all positions)

    Returns:
        Panel HTML pane with embedded structure viewer

    Raises:
        ValueError: if PDB sequence doesn't match reference with sufficient identity
    """
    # Reconstruct reference sequence from agg_df
    # Note: positions in agg_df are 1-indexed
    wt_map = agg_df.groupby('position')['wt'].first().to_dict()

    # Get min and max positions to build full sequence
    min_pos = min(wt_map.keys())
    max_pos = max(wt_map.keys())

    # Build sequence - fill gaps with 'X' if any positions are missing
    reference_sequence = ''.join(wt_map.get(pos, 'X') for pos in range(min_pos, max_pos + 1))

    # Align PDB to reference and get mapping
    alignment_info = align_pdb_to_reference(
        pdb_file, reference_sequence,
        min_identity=min_identity
    )
    pdb_mapping = alignment_info['mapping']
    selected_chain = alignment_info['chain_id']

    # Calculate mutation frequency per position (fillna(0) converts NaN for WT to 0)
    position_freqs = agg_df.groupby('position')['proportion_of_seqs'].sum().fillna(0).to_dict()

    # Parse position_range if provided
    positions_to_include = None
    if position_range and position_range.strip():
        positions_to_include = set()
        # Parse ranges like "20-30" or "20-30,50-60"
        for range_str in position_range.split(','):
            range_str = range_str.strip()
            if '-' in range_str:
                try:
                    start, end = map(int, range_str.split('-'))
                    positions_to_include.update(range(start, end + 1))
                except ValueError:
                    pass  # Skip invalid ranges

    # Create a unified dict mapping PDB positions to mutation data for ALL mapped residues
    # pdb_position -> {'frequency': float, 'wt': str, 'ref_position': int}
    pdb_mutation_data = {}
    for ref_pos in pdb_mapping.keys():
        chain_id, pdb_pos = pdb_mapping[ref_pos]

        # If position_range is specified, only include positions in range
        # Otherwise set frequency to 0 (will be colored white due to frequency_floor)
        if positions_to_include is not None and ref_pos not in positions_to_include:
            freq = 0.0
        else:
            freq = position_freqs.get(ref_pos, 0.0)

        pdb_mutation_data[pdb_pos] = {
            'frequency': freq,
            'wt': wt_map[ref_pos],
            'ref_position': ref_pos
        }

    print(f"Mapped {len(pdb_mutation_data)} residues to PDB chain {selected_chain}")

    # Create viewer
    view = create_structure_heatmap(
        pdb_file, pdb_mutation_data, selected_chain,
        alignment_info['pdb_sequence'], alignment_info['pdb_positions'],
        colormap=colormap, width=width, height=height, label_step=label_step, frequency_floor=frequency_floor
    )

    # Generate HTML
    structure_html = view._make_html()

    # Create Panel pane with the HTML using an iframe with data URI
    # Encode HTML as base64 for data URI
    encoded_html = base64.b64encode(structure_html.encode('utf-8')).decode('utf-8')
    iframe_html = f'<iframe src="data:text/html;base64,{encoded_html}" width="{width}px" height="{height}px" frameborder="0"></iframe>'

    structure_pane = pn.pane.HTML(
        iframe_html,
        width=width,
        height=height
    )

    return structure_pane


def main():
    """Main entry point for command line and snakemake usage"""
    if 'snakemake' in globals():
        # Running as snakemake script
        args = type('Args', (), {
            'agg_csv': snakemake.input.agg_csv,
            'pdb_file': snakemake.input.pdb,
            'output_html': snakemake.output.html,
            'colormap': snakemake.params.get('colormap', 'bmy_r'),
            'width': snakemake.params.get('width', 800),
            'height': snakemake.params.get('height', 600),
            'label_step': snakemake.params.get('label_step', 50),
            'min_identity': snakemake.params.get('min_identity', 0.95)
        })()
    else:
        # Running as standalone script
        parser = argparse.ArgumentParser(
            description='Generate 3D protein structure viewer with mutation frequency coloring',
            formatter_class=argparse.RawDescriptionHelpFormatter
        )
        parser.add_argument('--agg_csv', required=True,
                          help='Aggregated mutations CSV from SequenceAnalyzer.aggregate_mutations')
        parser.add_argument('--pdb_file', required=True, help='PDB structure file')
        parser.add_argument('--output_html', required=True, help='Output HTML file')
        parser.add_argument('--colormap', default='bmy_r', help='Colormap name (default: bmy_r)')
        parser.add_argument('--width', type=int, default=800, help='Viewer width in pixels')
        parser.add_argument('--height', type=int, default=600, help='Viewer height in pixels')
        parser.add_argument('--label_step', type=int, default=50,
                          help='Label every Nth position (default: 50)')
        parser.add_argument('--min_identity', type=float, default=0.95,
                          help='Minimum sequence identity required (default: 0.95)')
        args = parser.parse_args()

    # Load aggregated mutations CSV
    agg_df = pd.read_csv(args.agg_csv)

    # Generate structure viewer HTML
    structure_pane = create_structure_pane(
        pdb_file=args.pdb_file,
        agg_df=agg_df,
        colormap=args.colormap,
        width=args.width,
        height=args.height,
        label_step=args.label_step,
        min_identity=args.min_identity
    )

    # Save to HTML file
    pn.extension()
    structure_pane.save(args.output_html)
    print(f"Structure viewer saved to {args.output_html}")


if __name__ == '__main__':
    main()
