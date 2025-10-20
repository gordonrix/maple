#!/usr/bin/env python3
"""
Mutation analysis of aligned sequences from BAM files.

This script processes aligned BAM files to identify and quantify mutations at both
nucleotide (NT) and amino acid (AA) levels. It generates comprehensive mutation
statistics including:
- Mutation frequencies at each position
- Mutation distributions across sequences
- Genotype identification and clustering
- Insertion and deletion tracking
- Quality score filtering

Key features:
1. NT and AA mutation analysis: Analyzes substitutions, insertions, and deletions
2. Quality filtering: Filters mutations based on quality scores
3. Genotype clustering: Groups identical mutation patterns
4. Frameshift detection: Identifies and optionally excludes frameshift mutations
5. Representative alignments: Outputs alignments for high-abundance genotypes
6. Barcode integration: Incorporates barcode information when available

The script can run standalone or as part of a Snakemake pipeline.
"""

import os
import sys
import argparse
import tempfile

import numpy as np
import pandas as pd
import pysam
from Bio import SeqIO
from Bio.Seq import Seq

from common import dist_to_DF

# Generate codon table once at module load - faster than Bio.Seq.translate()
CODON_TABLE = {
    ''.join(codon): str(Seq(''.join(codon)).translate())
    for codon in (
        (n1, n2, n3)
        for n1 in 'ATGC'
        for n2 in 'ATGC'
        for n3 in 'ATGC'
    )
}

def load_reference_sequences(reference_fasta, do_aa_analysis):
    """
    Load and process reference sequences from FASTA file.

    Args:
        reference_fasta: Path to reference FASTA file
        do_aa_analysis: Whether to load protein reference sequence

    Returns:
        dict: Dictionary containing reference sequence information
    """
    sequences = list(SeqIO.parse(reference_fasta, 'fasta'))

    # First sequence is alignment reference
    ref_aln = sequences[0]
    ref_aln.seq = ref_aln.seq.upper()
    ref_aln_str = str(ref_aln.seq)

    # Second sequence is trimmed NT reference for analysis
    ref_nt = sequences[1]
    ref_nt.seq = ref_nt.seq.upper()
    ref_nt_str = str(ref_nt.seq)

    # Find trimmed reference in alignment reference
    ref_nt_start = ref_aln_str.find(ref_nt_str)
    use_reverse_complement = False

    if ref_nt_start == -1:
        # Try reverse complement
        use_reverse_complement = True
        ref_nt.seq = ref_nt.seq.reverse_complement()
        ref_nt_str = str(ref_nt.seq)
        ref_nt_start = ref_aln_str.find(ref_nt_str)

    ref_nt_end = len(ref_nt_str) + ref_nt_start

    result = {
        'ref_aln': ref_aln,
        'ref_aln_str': ref_aln_str,
        'ref_nt': ref_nt,
        'ref_nt_str': ref_nt_str,
        'ref_nt_start': ref_nt_start,
        'ref_nt_end': ref_nt_end,
        'use_reverse_complement': use_reverse_complement
    }

    # Third sequence is protein reference (optional)
    if do_aa_analysis and len(sequences) >= 3:
        ref_protein = sequences[2]
        if use_reverse_complement:
            ref_protein.seq = ref_protein.seq.reverse_complement()
        ref_protein_str = str(ref_protein.seq).upper()

        ref_protein_start = ref_nt_str.find(ref_protein_str)
        ref_protein_end = len(ref_protein_str) + ref_protein_start

        if use_reverse_complement:
            ref_protein_str = str(Seq(ref_protein_str).reverse_complement())

        result.update({
            'ref_protein_str': ref_protein_str,
            'ref_protein_start': ref_protein_start,
            'ref_protein_end': ref_protein_end
        })

    return result


def clean_alignment(bam_entry, ref_sequences, quality_minimum, analyze_seqs_with_indels,
                    do_aa_analysis, has_quality_scores):
    """
    Process BAM entry to extract aligned sequences and identify indels.

    Trims alignment to trimmed reference region, identifies insertions and deletions,
    and validates alignment quality. Returns None if alignment fails quality checks.

    Args:
        bam_entry: pysam AlignmentFile entry
        ref_sequences: Dictionary of reference sequence information
        quality_minimum: Minimum quality score threshold
        analyze_seqs_with_indels: Whether to analyze sequences with frameshift indels
        do_aa_analysis: Whether AA analysis is being performed
        has_quality_scores: Whether quality scores are available

    Returns:
        tuple: (ref_aln, align_str, query_aln, query_qualities, insertions, deletions)
               or None if alignment fails validation
    """
    ref_aln_str = ref_sequences['ref_aln_str']
    ref_nt_start = ref_sequences['ref_nt_start']
    ref_nt_end = ref_sequences['ref_nt_end']
    ref_aln_id = ref_sequences['ref_aln'].id

    # Validate alignment reference
    if bam_entry.reference_name != ref_aln_id:
        return None, ('alignment uses wrong reference sequence', 'N/A')

    # Validate alignment spans trimmed reference
    if bam_entry.reference_start > ref_nt_start:
        return None, ('alignment starts past trimmed reference start', 'N/A')

    if bam_entry.reference_end < ref_nt_end:
        return None, ('alignment ends before trimmed reference end', 'N/A')

    # Build alignment strings
    ref_index = bam_entry.reference_start
    query_index = 0
    ref_aln = ' ' * ref_index
    query_aln = ' ' * ref_index
    query_qualities = [-1] * ref_index if has_quality_scores else None
    align_str = ' ' * ref_index
    insertions = []
    deletions = []

    ref_protein_start = ref_sequences.get('ref_protein_start', 0)
    ref_protein_end = ref_sequences.get('ref_protein_end', 0)

    for cigar_op, length in bam_entry.cigartuples:

        if cigar_op == 0:  # Match
            ref_segment = ref_aln_str[ref_index:ref_index + length]
            query_segment = bam_entry.query_alignment_sequence[query_index:query_index + length]
            align_segment = ''.join(['|' if r == q else '.' for r, q in zip(ref_segment, query_segment)])

            ref_aln += ref_segment
            query_aln += query_segment
            if has_quality_scores:
                query_qualities += list(bam_entry.query_alignment_qualities[query_index:query_index + length])
            align_str += align_segment
            ref_index += length
            query_index += length

        elif cigar_op == 1:  # Insertion
            # Check for frameshift if doing AA analysis
            if (do_aa_analysis and not analyze_seqs_with_indels and
                length % 3 != 0 and ref_protein_start <= ref_index < ref_protein_end):
                return None, ('frameshift insertion', query_index)

            # Record insertion if within trimmed reference
            if ref_nt_start <= ref_index < ref_nt_end:
                insertion_seq = bam_entry.query_alignment_sequence[query_index:query_index + length]
                insertions.append((ref_index - ref_nt_start, insertion_seq))

            query_index += length

        elif cigar_op == 2:  # Deletion
            # Check for frameshift if doing AA analysis
            if (do_aa_analysis and not analyze_seqs_with_indels and
                length % 3 != 0 and
                ref_protein_start <= ref_index + length and ref_index < ref_protein_end):
                return None, ('frameshift deletion', query_index)

            ref_aln += ref_aln_str[ref_index:ref_index + length]
            query_aln += '-' * length
            align_str += ' ' * length
            if has_quality_scores:
                query_qualities += [0] * length

            # Record deletion if within trimmed reference
            if ref_nt_start <= ref_index + length and ref_index < ref_nt_end:
                deletions.append((ref_index - ref_nt_start, length))

            ref_index += length

    # Extract trimmed region
    ref_aln_trimmed = ref_aln[ref_nt_start:ref_nt_end]
    align_str_trimmed = align_str[ref_nt_start:ref_nt_end]
    query_aln_trimmed = query_aln[ref_nt_start:ref_nt_end]
    query_qualities_trimmed = query_qualities[ref_nt_start:ref_nt_end] if has_quality_scores else None

    return (ref_aln_trimmed, align_str_trimmed, query_aln_trimmed,
            query_qualities_trimmed, insertions, deletions), None


def clean_alignment_reverse_complement(clean_aln, ref_nt_length):
    """
    Reverse complement aligned sequences and adjust indel positions.

    Args:
        clean_aln: Output from clean_alignment()
        ref_nt_length: Length of trimmed NT reference

    Returns:
        tuple: Reverse complemented alignment data
    """
    ref, align_str, seq, q_scores, insertions, deletions = clean_aln

    # Reverse complement sequences
    ref = str(Seq(ref).reverse_complement())
    align_str = align_str[::-1]
    seq = str(Seq(seq).reverse_complement())
    q_scores = q_scores[::-1] if q_scores else None

    # Adjust insertion positions and sequences
    insertions.reverse()
    insertions_out = []
    for pos, ins_seq in insertions:
        new_pos = ref_nt_length - pos
        new_seq = str(Seq(ins_seq).reverse_complement())
        insertions_out.append((new_pos, new_seq))

    # Adjust deletion positions
    deletions.reverse()
    deletions_out = []
    for pos, length in deletions:
        new_pos = ref_nt_length - pos - length
        deletions_out.append((new_pos, length))

    return ref, align_str, seq, q_scores, insertions_out, deletions_out


def identify_mutations(clean_aln, ref_sequences, quality_minimum,
                       do_aa_analysis, has_quality_scores):
    """
    Identify mutations in an aligned sequence.

    Processes aligned sequence to identify NT and AA substitutions, insertions,
    and deletions. Filters mutations based on quality scores.

    Args:
        clean_aln: Output from clean_alignment()
        ref_sequences: Dictionary of reference sequence information
        quality_minimum: Minimum quality score threshold
        do_aa_analysis: Whether to perform AA mutation analysis
        has_quality_scores: Whether quality scores are available

    Returns:
        tuple: (nt_mut_array, aa_mut_array, genotype_data)
    """
    ref, align_str, seq, q_scores, insertions, deletions = clean_aln

    ref_nt_str = ref_sequences['ref_nt_str']
    nts = 'ATGC'
    aas = 'ACDEFGHIKLMNPQRSTVWY*'

    # Initialize arrays
    nt_mut_array = np.zeros((len(ref_nt_str), len(nts)), dtype=int)
    aa_mut_array = None
    indel_codons = []

    # Process AA indels if doing AA analysis
    if do_aa_analysis:
        ref_protein_start = ref_sequences['ref_protein_start']
        ref_protein_end = ref_sequences['ref_protein_end']
        ref_protein_str = ref_sequences['ref_protein_str']

        # Track codons affected by indels
        for index, insertion_seq in insertions:
            if ref_protein_start <= index < ref_protein_end:
                prot_index = index - ref_protein_start
                if prot_index % 3 != 0:
                    indel_codons.append(int(prot_index / 3))

        for index, length in deletions:
            if (ref_protein_start <= index < ref_protein_end or
                ref_protein_start <= index + length < ref_protein_end):
                prot_index_start = index - ref_protein_start
                prot_index_end = (index + length) - ref_protein_start
                first_codon = int(prot_index_start / 3)
                last_codon = int(prot_index_end / 3)
                indel_codons.extend(range(first_codon, last_codon + 1))

        aa_mut_array = np.zeros((int(len(ref_protein_str) / 3), len(aas)), dtype=int)

    # Find mismatches
    mismatches = [i for i, a in enumerate(align_str) if a == '.']

    nt_substitutions = []
    aa_nonsynonymous = []
    aa_synonymous = []
    codons_checked = []

    # Process each mismatch
    for i in mismatches:
        # Check quality score
        if has_quality_scores and q_scores[i] < quality_minimum:
            continue

        # Record NT mutation
        wt_nt = ref[i]
        mut_nt = seq[i]
        nt_mut_array[i, nts.find(mut_nt)] += 1
        nt_substitutions.append(f'{wt_nt}{i + 1}{mut_nt}')

        # Process AA mutation if applicable
        if do_aa_analysis:
            ref_protein_start = ref_sequences['ref_protein_start']
            ref_protein_end = ref_sequences['ref_protein_end']

            if ref_protein_start <= i < ref_protein_end:
                prot_index = i - ref_protein_start
                codon = int(prot_index / 3)

                # Skip if already checked or affected by indel
                if codon in codons_checked or codon in indel_codons:
                    continue

                codons_checked.append(codon)
                codon_pos = prot_index % 3
                codon_indices = range(i - codon_pos, i + (3 - codon_pos))

                # Check quality scores for entire codon
                if has_quality_scores:
                    codon_q_scores = [q_scores[idx] for idx in codon_indices]
                    if any(qs < quality_minimum for qs in codon_q_scores):
                        continue

                # Translate codons
                wt_codon = ref[codon_indices[0]:codon_indices[-1] + 1]
                mut_codon = seq[codon_indices[0]:codon_indices[-1] + 1]
                wt_aa = CODON_TABLE.get(wt_codon, 'X')
                mut_aa = CODON_TABLE.get(mut_codon, 'X')

                # Build strings
                if wt_aa != mut_aa:
                    aa_mut_array[codon, aas.find(mut_aa)] += 1
                    aa_nonsynonymous.append(f'{wt_aa}{codon + 1}{mut_aa}')
                else:
                    aa_synonymous.append(f'{wt_aa}{codon + 1}')

    # Build genotype data

    ins_output = ', '.join([f'{idx}ins{seq}' for idx, seq in insertions])
    del_output = ', '.join([f'{idx}del{length}' for idx, length in deletions])
    ins_length = sum(len(seq) for _, seq in insertions) if insertions else 0
    del_length = sum(length for _, length in deletions) if deletions else 0

    genotype_data = [
        ', '.join(nt_substitutions),
        len(nt_substitutions),
        ins_output,
        del_output,
        ins_length,
        del_length
    ]

    if do_aa_analysis:
        genotype_data.extend([
            ', '.join(aa_nonsynonymous),
            ', '.join(aa_synonymous),
            len(aa_nonsynonymous)
        ])

    return nt_mut_array, aa_mut_array, genotype_data


def process_bam_file(bam_path, reference_fasta, output_files, quality_minimum,
                     analyze_seqs_with_indels, do_aa_analysis, use_raw_mut_count,
                     highest_abundance_genotypes, desired_genotype_ids,
                     demux_screen_failures, has_barcode_column):
    """
    Process BAM file to identify mutations and generate output files.

    Args:
        bam_path: Path to input BAM file
        reference_fasta: Path to reference FASTA file
        output_files: List of output file paths
        quality_minimum: Minimum quality score threshold
        analyze_seqs_with_indels: Whether to analyze sequences with frameshift indels
        do_aa_analysis: Whether to perform AA mutation analysis
        use_raw_mut_count: Whether to output raw counts vs proportions
        highest_abundance_genotypes: Number of high-abundance genotypes to output
        desired_genotype_ids: Specific genotype IDs to output
        demux_screen_failures: Whether to screen sequences with failed barcodes
        has_barcode_column: Whether barcode information is available

    Returns:
        None (writes output files)
    """
    # Load reference sequences
    ref_sequences = load_reference_sequences(reference_fasta, do_aa_analysis)

    ref_nt_str = ref_sequences['ref_nt_str']
    use_reverse_complement = ref_sequences['use_reverse_complement']

    nts = 'ATGC'
    aas = 'ACDEFGHIKLMNPQRSTVWY*'

    # Initialize data structures
    nt_mut_array = np.zeros((len(ref_nt_str), len(nts)), dtype=int)
    nt_mut_dist = np.zeros(len(ref_nt_str), dtype=int)

    failures_list = []
    genotypes_list = []

    all_seqs_columns = [
        'seq_ID', 'avg_quality_score', 'NT_substitutions', 'NT_substitutions_count',
        'NT_insertions', 'NT_deletions', 'NT_insertion_length', 'NT_deletion_length'
    ]

    if do_aa_analysis:
        ref_protein_str = ref_sequences['ref_protein_str']
        prot_length = int(len(ref_protein_str) / 3)
        aa_mut_array = np.zeros((prot_length, len(aas)), dtype=int)
        aa_mut_dist = np.zeros(prot_length + 1, dtype=int)
        all_seqs_columns.extend([
            'AA_substitutions_nonsynonymous', 'AA_substitutions_synonymous',
            'AA_substitutions_nonsynonymous_count'
        ])

    if has_barcode_column:
        all_seqs_columns.append('barcode(s)')

    genotypes_columns = all_seqs_columns[2:]  # Columns for grouping genotypes

    # Open BAM file
    bam_file = pysam.AlignmentFile(bam_path, 'rb')

    # Detect if quality scores are present
    has_quality_scores = False
    for bam_entry in bam_file:
        if bam_entry.query_alignment_qualities:
            has_quality_scores = True
        bam_file.reset()
        break

    # Process each sequence
    for bam_entry in bam_file:
        # Skip sequences with failed barcodes if screening
        if demux_screen_failures and has_barcode_column:
            if 'fail' in bam_entry.get_tag('BC'):
                continue

        # Clean alignment
        clean_aln, failure = clean_alignment(
            bam_entry, ref_sequences, quality_minimum,
            analyze_seqs_with_indels, do_aa_analysis, has_quality_scores
        )

        if clean_aln is None:
            failures_list.append([bam_entry.query_name, failure[0], failure[1]])
            continue

        # Reverse complement if needed
        if use_reverse_complement:
            clean_aln = clean_alignment_reverse_complement(clean_aln, len(ref_nt_str))

        # Identify mutations
        seq_nt_mut_array, seq_aa_mut_array, genotype_data = identify_mutations(
            clean_aln, ref_sequences, quality_minimum,
            do_aa_analysis, has_quality_scores
        )

        # Update arrays
        nt_mut_array += seq_nt_mut_array
        seq_total_nt_muts = seq_nt_mut_array.sum()
        nt_mut_dist[seq_total_nt_muts] += 1

        if do_aa_analysis:
            aa_mut_array += seq_aa_mut_array
            seq_total_aa_muts = seq_aa_mut_array.sum()
            aa_mut_dist[seq_total_aa_muts] += 1

        # Build record
        avg_q_score = np.average(clean_aln[3]) if has_quality_scores else -1
        record = [bam_entry.query_name, avg_q_score] + genotype_data

        if has_barcode_column:
            record.append(bam_entry.get_tag('BC'))

        genotypes_list.append(record)

    bam_file.close()

    # Build DataFrames
    all_seqs_df = pd.DataFrame(genotypes_list, columns=all_seqs_columns)
    failures_df = pd.DataFrame(failures_list, columns=['seq_ID', 'failure_reason', 'failure_index'])

    # Group identical genotypes
    all_seqs_df['genotype->seq'] = all_seqs_df.groupby(by=genotypes_columns).ngroup()
    all_seqs_df['count'] = all_seqs_df.groupby(by='genotype->seq')['seq_ID'].transform('count')

    # Create condensed genotypes DataFrame
    genotypes_condensed = (all_seqs_df
                          .sort_values('avg_quality_score', ascending=False)
                          .drop_duplicates('genotype->seq', keep='first')
                          [['genotype->seq', 'count'] + all_seqs_columns])

    genotypes_condensed = genotypes_condensed.sort_values(['count', 'NT_substitutions_count'],
                                                          ascending=[False, True])

    # Move wildtype to beginning if present
    wildtype_mask = ((genotypes_condensed['NT_substitutions'] == '') &
                     (genotypes_condensed['NT_insertions'] == '') &
                     (genotypes_condensed['NT_deletions'] == ''))
    wildtype_df = genotypes_condensed[wildtype_mask]

    if len(wildtype_df) > 0:
        genotypes_condensed = genotypes_condensed.drop(index=wildtype_df.index)
        genotypes_condensed = pd.concat([wildtype_df, genotypes_condensed]).reset_index(drop=True)
        genotypes_condensed.rename(index={0: 'wildtype'}, inplace=True)
    else:
        genotypes_condensed.reset_index(drop=True, inplace=True)

    # Make genotype IDs 1-indexed
    if len(genotypes_condensed) != 0 and genotypes_condensed.index[0] == 0:
        genotypes_condensed.index += 1

    # Add genotype_ID column
    genotypes_condensed = genotypes_condensed.reset_index().rename(columns={'index': 'genotype_ID'})
    seq_to_genotype_dict = dict(zip(genotypes_condensed['genotype->seq'],
                                    genotypes_condensed['genotype_ID']))
    all_seqs_df['genotype_ID'] = all_seqs_df['genotype->seq'].map(seq_to_genotype_dict)

    # Write alignments for high-abundance genotypes
    write_genotype_alignments(
        bam_path, genotypes_condensed, ref_sequences,
        highest_abundance_genotypes, desired_genotype_ids,
        quality_minimum, analyze_seqs_with_indels, do_aa_analysis,
        has_quality_scores, output_files[0]
    )

    # Prepare output DataFrames
    total_seqs = int(nt_mut_dist.sum())

    # NT mutation frequencies
    if use_reverse_complement:
        nt_ids = list(str(Seq(ref_nt_str).reverse_complement()))
    else:
        nt_ids = list(ref_nt_str)

    wt_nts = [f'{nt_id}{i + 1}' for i, nt_id in enumerate(nt_ids)]
    nt_mut_df = pd.DataFrame(nt_mut_array, columns=list(nts))
    nt_mut_df['wt_nucleotide'] = pd.Series(wt_nts)
    nt_mut_df.set_index('wt_nucleotide', inplace=True)

    if not use_raw_mut_count and total_seqs > 0:
        nt_mut_df = nt_mut_df.divide(total_seqs)

    # NT mutation distribution
    nt_dist_df = dist_to_DF(np.trim_zeros(nt_mut_dist, 'b'), 'NT mutations', 'sequences')

    # Prepare columns to drop/include
    condensed_drop = ['genotype->seq', 'seq_ID', 'avg_quality_score']
    all_include = ['seq_ID', 'genotype_ID']

    if has_barcode_column:
        condensed_drop.append('barcode(s)')
        all_include.append('barcode(s)')

    # Write output files
    genotypes_condensed.drop(columns=condensed_drop).to_csv(output_files[1], index=False)
    all_seqs_df[all_include].to_csv(output_files[2], index=False)
    failures_df.to_csv(output_files[3], index=False)
    nt_mut_df.to_csv(output_files[4])
    nt_dist_df.to_csv(output_files[5], index=False)

    # AA analysis outputs
    if do_aa_analysis:
        resi_ids = list(str(Seq(ref_protein_str).translate()))
        resi_positions = [str(i) for i in range(1, prot_length + 1)]
        wt_resis = [f'{res_id}{pos}' for res_id, pos in zip(resi_ids, resi_positions)]

        aa_mut_df = pd.DataFrame(aa_mut_array, columns=list(aas))
        aa_mut_df['wt_residues'] = pd.Series(wt_resis)
        aa_mut_df.set_index('wt_residues', inplace=True)

        if not use_raw_mut_count and total_seqs > 0:
            aa_mut_df = aa_mut_df.divide(total_seqs)

        aa_dist_df = dist_to_DF(np.trim_zeros(aa_mut_dist, 'b'), 'AA mutations', 'sequences')

        aa_mut_df.to_csv(output_files[6])
        aa_dist_df.to_csv(output_files[7], index=False)


def write_genotype_alignments(bam_path, genotypes_condensed, ref_sequences,
                               highest_abundance_genotypes, desired_genotype_ids,
                               quality_minimum, analyze_seqs_with_indels,
                               do_aa_analysis, has_quality_scores, output_file):
    """
    Write alignments for selected genotypes to output file.

    Args:
        bam_path: Path to input BAM file
        genotypes_condensed: DataFrame of condensed genotypes
        ref_sequences: Dictionary of reference sequence information
        highest_abundance_genotypes: Number of high-abundance genotypes to output
        desired_genotype_ids: Specific genotype IDs to output
        quality_minimum: Minimum quality score threshold
        analyze_seqs_with_indels: Whether to analyze sequences with frameshift indels
        do_aa_analysis: Whether AA analysis is being performed
        has_quality_scores: Whether quality scores are available
        output_file: Path to output alignments file

    Returns:
        None (writes output file)
    """
    # Select genotypes to output
    genotype_alns_out_df = genotypes_condensed.iloc[0:highest_abundance_genotypes + 1]

    if desired_genotype_ids:
        desired_ids = [int(gid) for gid in str(desired_genotype_ids).split(', ')
                      if int(gid) <= len(genotypes_condensed)]
        genotype_alns_out_df = pd.concat([
            genotype_alns_out_df,
            genotypes_condensed.iloc[desired_ids]
        ])

    with open(output_file, 'w') as txt_out:
        # Create temporary sorted BAM file for efficient lookup
        with tempfile.NamedTemporaryFile(suffix='bam', delete=False) as temp_bam:
            sorted_bam_path = temp_bam.name

        pysam.sort('-n', '-o', sorted_bam_path, bam_path)
        sorted_bam = pysam.AlignmentFile(sorted_bam_path, 'rb')

        name_indexed_bam = pysam.IndexedReads(sorted_bam)
        name_indexed_bam.build()

        use_reverse_complement = ref_sequences['use_reverse_complement']
        ref_nt_length = len(ref_sequences['ref_nt_str'])

        for row in genotype_alns_out_df.itertuples():
            if row.genotype_ID == 'wildtype':
                continue

            seq_id = row.seq_ID
            iterator = name_indexed_bam.find(seq_id)

            for bam_entry in iterator:
                break

            clean_aln, _ = clean_alignment(
                bam_entry, ref_sequences, quality_minimum,
                analyze_seqs_with_indels, do_aa_analysis, has_quality_scores
            )

            if clean_aln:
                if use_reverse_complement:
                    clean_aln = clean_alignment_reverse_complement(clean_aln, ref_nt_length)

                ref, align_string, seq, _, _, _ = clean_aln

                txt_out.write(f'Genotype {row.genotype_ID} representative sequence. Sequence ID: {seq_id}\n')
                for string in [ref, align_string, seq]:
                    txt_out.write(string + '\n')
                txt_out.write('\n')

        sorted_bam.close()
        os.remove(sorted_bam_path)


def main():
    """Main entry point for mutation analysis."""
    # Parse arguments
    if 'snakemake' in globals():
        # Running as snakemake script
        args = type('Args', (), {
            'bam': str(snakemake.input.bam),
            'reference': snakemake.config['runs'][snakemake.wildcards.tag]['reference'],
            'output': list(snakemake.output),
            'quality_minimum': snakemake.config.get('mutation_analysis_quality_score_minimum', 5),
            'analyze_seqs_with_indels': snakemake.params.analyze_seqs_with_indels,
            'do_aa_analysis': snakemake.config['do_AA_mutation_analysis'][snakemake.wildcards.tag],
            'use_raw_mut_count': snakemake.params.mutations_frequencies_raw,
            'highest_abundance_genotypes': snakemake.config.get('highest_abundance_genotypes', 0),
            'desired_genotype_ids': snakemake.config.get('genotype_ID_alignments', 0),
            'demux_screen_failures': snakemake.config.get('demux_screen_failures', False),
            'tag': snakemake.wildcards.tag
        })()

        # Determine if barcode column is needed
        has_barcode_column = False
        if snakemake.config['do_demux'][args.tag]:
            for bc_type in snakemake.config['runs'][args.tag]['barcode_info']:
                if snakemake.config['runs'][args.tag]['barcode_info'][bc_type].get('label_only', False):
                    has_barcode_column = True

        args.has_barcode_column = has_barcode_column
    else:
        # Running as standalone script
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument('--bam', required=True, help='Input BAM file')
        parser.add_argument('--reference', required=True, help='Reference FASTA file')
        parser.add_argument('--output-dir', type=str, default=None,
                          help='Output directory (basename from BAM file will be used for output names)')
        parser.add_argument('--output', nargs='+', default=None,
                          help='Output files (alignments, genotypes, seq-IDs, failures, NT-freqs, NT-dist, [AA-freqs, AA-dist])')
        parser.add_argument('--quality-minimum', type=int, default=5,
                          help='Minimum quality score threshold')
        parser.add_argument('--analyze-seqs-with-indels', action='store_true',
                          help='Analyze sequences with frameshift indels')
        parser.add_argument('--do-aa-analysis', action='store_true',
                          help='Perform AA mutation analysis')
        parser.add_argument('--use-raw-mut-count', action='store_true',
                          help='Output raw counts instead of proportions')
        parser.add_argument('--highest-abundance-genotypes', type=int, default=0,
                          help='Number of high-abundance genotypes to output')
        parser.add_argument('--desired-genotype-ids', type=str, default='',
                          help='Comma-separated list of specific genotype IDs to output')
        parser.add_argument('--demux-screen-failures', action='store_true',
                          help='Screen sequences with failed barcodes')
        parser.add_argument('--has-barcode-column', action='store_true',
                          help='Whether barcode information is available')
        args = parser.parse_args()

        # Handle output paths
        if args.output_dir and not args.output:
            # Generate output paths from BAM filename
            bam_basename = os.path.splitext(os.path.basename(args.bam))[0]
            os.makedirs(args.output_dir, exist_ok=True)
            args.output = [
                os.path.join(args.output_dir, f'{bam_basename}_alignments.txt'),
                os.path.join(args.output_dir, f'{bam_basename}_genotypes.csv'),
                os.path.join(args.output_dir, f'{bam_basename}_seq-IDs.csv'),
                os.path.join(args.output_dir, f'{bam_basename}_failures.csv'),
                os.path.join(args.output_dir, f'{bam_basename}_NT-mutation-frequencies.csv'),
                os.path.join(args.output_dir, f'{bam_basename}_NT-mutation-distribution.csv'),
            ]
            if args.do_aa_analysis:
                args.output.extend([
                    os.path.join(args.output_dir, f'{bam_basename}_AA-mutation-frequencies.csv'),
                    os.path.join(args.output_dir, f'{bam_basename}_AA-mutation-distribution.csv'),
                ])
        elif not args.output:
            parser.error('Either --output-dir or --output must be specified')

    # Process BAM file
    process_bam_file(
        args.bam,
        args.reference,
        args.output,
        args.quality_minimum,
        args.analyze_seqs_with_indels,
        args.do_aa_analysis,
        args.use_raw_mut_count,
        args.highest_abundance_genotypes,
        args.desired_genotype_ids,
        args.demux_screen_failures,
        args.has_barcode_column
    )


if __name__ == '__main__':
    main()
