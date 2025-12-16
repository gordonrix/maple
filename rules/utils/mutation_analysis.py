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

from common import dist_to_DF, load_references_from_csv

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

def load_reference_sequences(references_csv, do_aa_analysis):
    """
    Load and process reference sequences from CSV file.

    Args:
        references_csv: Path to references CSV file
        do_aa_analysis: Whether to load protein reference sequence

    Returns:
        dict: {ref_name: {reference sequence information}} where each reference has:
            - ref_aln_str: alignment reference sequence
            - ref_nt_str: trimmed NT reference sequence
            - ref_nt_start: start position of NT reference in alignment
            - ref_nt_end: end position of NT reference in alignment
            - use_reverse_complement: whether reverse complement is used
            - ref_protein_str: protein reference (if do_aa_analysis)
            - ref_protein_start: start position of protein reference
            - ref_protein_end: end position of protein reference
    """
    # Load from CSV
    references_dict, errors, _, _ = load_references_from_csv(references_csv)

    if errors:
        for error in errors:
            print(error, file=sys.stderr)
        raise ValueError("Failed to load references from CSV")

    # Process each reference
    processed_refs = {}
    for ref_name, ref_data in references_dict.items():
        # Find trimmed reference in alignment reference
        ref_nt_start = ref_data['alignment_seq'].find(ref_data['NT_seq'])
        use_reverse_complement = False
        ref_nt_str = ref_data['NT_seq']

        if ref_nt_start == -1:
            # Try reverse complement
            use_reverse_complement = True
            ref_nt_str = str(Seq(ref_data['NT_seq']).reverse_complement())
            ref_nt_start = ref_data['alignment_seq'].find(ref_nt_str)

        ref_nt_end = len(ref_nt_str) + ref_nt_start

        result = {
            'ref_aln_str': ref_data['alignment_seq'],
            'ref_nt_str': ref_nt_str,
            'ref_nt_start': ref_nt_start,
            'ref_nt_end': ref_nt_end,
            'use_reverse_complement': use_reverse_complement
        }

        # Process protein reference (optional)
        if do_aa_analysis and ref_data.get('coding_seq'):
            ref_protein_str = ref_data['coding_seq']

            if use_reverse_complement:
                ref_protein_str = str(Seq(ref_data['coding_seq']).reverse_complement())

            ref_protein_start = ref_nt_str.find(ref_protein_str)
            ref_protein_end = len(ref_protein_str) + ref_protein_start

            result.update({
                'ref_protein_str': ref_protein_str,
                'ref_protein_start': ref_protein_start,
                'ref_protein_end': ref_protein_end
            })

        processed_refs[ref_name] = result

    return processed_refs


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


def mutation_array_to_tidy_df(mut_array, ref_str, ref_name, NTorAA, use_reverse_complement,
                               total_seqs, use_raw_mut_count):
    """
    Convert mutation array to tidy DataFrame matching SequenceAnalyzer.aggregate_mutations format.

    Args:
        mut_array: Mutation count array (positions x characters)
        ref_str: Reference sequence string
        ref_name: Name of the reference
        NTorAA: 'NT' or 'AA'
        use_reverse_complement: Whether to reverse complement
        total_seqs: Total number of sequences
        use_raw_mut_count: Whether to output raw counts vs proportions

    Returns:
        pd.DataFrame: Tidy format with columns ['wt', 'position', 'mutation', 'total_count',
                                                  'proportion_of_seqs', 'reference_name']
    """
    if NTorAA == 'NT':
        chars = 'ATGC'
    else:
        chars = 'ACDEFGHIKLMNPQRSTVWY*'

    if use_reverse_complement:
        ref_chars = list(str(Seq(ref_str).reverse_complement()))
    else:
        ref_chars = list(ref_str)

    rows = []

    # Loop through all possible mutations
    for posi in range(len(ref_chars)):
        wt = ref_chars[posi]
        position = posi + 1  # 1-indexed

        for char_idx, mut in enumerate(chars):
            count = mut_array[posi, char_idx]

            # Set count to NaN for wt (matching SequenceAnalyzer behavior)
            if wt == mut:
                count = np.nan

            rows.append([wt, position, mut, count])

    df = pd.DataFrame(rows, columns=['wt', 'position', 'mutation', 'total_count'])

    if total_seqs > 0:
        df['proportion_of_seqs'] = df['total_count'] / total_seqs
    else:
        df['proportion_of_seqs'] = 0.0

    df['reference_name'] = ref_name

    return df


def process_single_reference(bam_file, ref_name, ref_sequences, quality_minimum,
                             analyze_seqs_with_indels, do_aa_analysis,
                             demux_screen_failures, has_barcode_column, has_quality_scores):
    """
    Process BAM entries for a single reference.

    Args:
        bam_file: Open pysam AlignmentFile
        ref_name: Name of the reference to process
        ref_sequences: Dictionary with reference sequence information
        quality_minimum: Minimum quality score threshold
        analyze_seqs_with_indels: Whether to analyze sequences with frameshift indels
        do_aa_analysis: Whether to perform AA mutation analysis
        demux_screen_failures: Whether to screen sequences with failed barcodes
        has_barcode_column: Whether barcode information is available
        has_quality_scores: Whether BAM file has quality scores

    Returns:
        tuple: (nt_mut_array, nt_mut_dist, aa_mut_array, aa_mut_dist, genotypes_list, failures_list)
    """
    ref_nt_str = ref_sequences['ref_nt_str']
    use_reverse_complement = ref_sequences['use_reverse_complement']

    nts = 'ATGC'
    aas = 'ACDEFGHIKLMNPQRSTVWY*'

    # Initialize data structures
    nt_mut_array = np.zeros((len(ref_nt_str), len(nts)), dtype=int)
    nt_mut_dist = np.zeros(len(ref_nt_str), dtype=int)
    aa_mut_array = None
    aa_mut_dist = None

    if do_aa_analysis and 'ref_protein_str' in ref_sequences:
        ref_protein_str = ref_sequences['ref_protein_str']
        prot_length = int(len(ref_protein_str) / 3)
        aa_mut_array = np.zeros((prot_length, len(aas)), dtype=int)
        aa_mut_dist = np.zeros(prot_length + 1, dtype=int)

    failures_list = []
    genotypes_list = []

    # Fetch BAM entries for this reference
    for bam_entry in bam_file.fetch(reference=ref_name):
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

        if do_aa_analysis and aa_mut_array is not None:
            aa_mut_array += seq_aa_mut_array
            seq_total_aa_muts = seq_aa_mut_array.sum()
            aa_mut_dist[seq_total_aa_muts] += 1

        # Build record
        avg_q_score = np.average(clean_aln[3]) if has_quality_scores else -1
        record = [bam_entry.query_name, avg_q_score, ref_name] + genotype_data

        if has_barcode_column:
            record.append(bam_entry.get_tag('BC'))

        genotypes_list.append(record)

    # Convert mutation arrays to tidy DataFrames
    total_seqs = int(nt_mut_dist.sum())

    nt_freq_df = mutation_array_to_tidy_df(
        nt_mut_array, ref_nt_str, ref_name, 'NT',
        use_reverse_complement, total_seqs, use_raw_mut_count=False
    )

    # NT mutation distribution
    nt_dist_df = dist_to_DF(np.trim_zeros(nt_mut_dist, 'b'), 'NT mutations', 'sequences')
    nt_dist_df['reference_name'] = ref_name

    # AA mutation frequencies if applicable
    aa_freq_df = None
    aa_dist_df = None

    if do_aa_analysis and aa_mut_array is not None:
        ref_protein_str = ref_sequences['ref_protein_str']
        total_aa_seqs = int(aa_mut_dist.sum())

        # Translate nucleotide coding sequence to amino acids
        ref_protein_translated = str(Seq(ref_protein_str).translate())

        aa_freq_df = mutation_array_to_tidy_df(
            aa_mut_array, ref_protein_translated, ref_name, 'AA',
            use_reverse_complement, total_aa_seqs, use_raw_mut_count=False
        )

        # AA mutation distribution
        aa_dist_df = dist_to_DF(np.trim_zeros(aa_mut_dist, 'b'), 'AA mutations', 'sequences')
        aa_dist_df['reference_name'] = ref_name

    return nt_freq_df, nt_dist_df, aa_freq_df, aa_dist_df, genotypes_list, failures_list


def process_bam_file(bam_path, references_csv, output_files, quality_minimum,
                     analyze_seqs_with_indels, do_aa_analysis, use_raw_mut_count,
                     highest_abundance_genotypes, desired_genotype_ids,
                     demux_screen_failures, has_barcode_column):
    """
    Process BAM file to identify mutations and generate output files.

    Args:
        bam_path: Path to input BAM file
        references_csv: Path to references CSV file
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
    references_dict = load_reference_sequences(references_csv, do_aa_analysis)

    # Open BAM file
    bam_file = pysam.AlignmentFile(bam_path, 'rb')

    # Detect if quality scores are present
    has_quality_scores = False
    for bam_entry in bam_file:
        if bam_entry.query_alignment_qualities:
            has_quality_scores = True
        bam_file.reset()
        break

    # Count sequences per reference to determine order
    ref_counts = {}
    for ref_name in references_dict.keys():
        ref_counts[ref_name] = bam_file.count(reference=ref_name)

    # Sort references by abundance (most abundant first)
    sorted_refs = sorted(ref_counts.items(), key=lambda x: x[1], reverse=True)
    most_abundant_ref = sorted_refs[0][0] if sorted_refs else None

    # Process each reference in order of abundance
    all_nt_freq_dfs = []
    all_nt_dist_dfs = []
    all_aa_freq_dfs = []
    all_aa_dist_dfs = []
    all_genotypes_lists = []
    all_failures_lists = []

    for ref_name, _ in sorted_refs:
        ref_sequences = references_dict[ref_name]
        nt_freq_df, nt_dist_df, aa_freq_df, aa_dist_df, genotypes_list, failures_list = process_single_reference(
            bam_file, ref_name, ref_sequences, quality_minimum,
            analyze_seqs_with_indels, do_aa_analysis,
            demux_screen_failures, has_barcode_column, has_quality_scores
        )

        # Collect DataFrames and lists
        all_nt_freq_dfs.append(nt_freq_df)
        all_nt_dist_dfs.append(nt_dist_df)

        if aa_freq_df is not None:
            all_aa_freq_dfs.append(aa_freq_df)
        if aa_dist_df is not None:
            all_aa_dist_dfs.append(aa_dist_df)

        all_genotypes_lists.extend(genotypes_list)
        all_failures_lists.extend(failures_list)

    bam_file.close()

    # Concatenate all DataFrames
    nt_freq_combined = pd.concat(all_nt_freq_dfs, ignore_index=True)
    nt_dist_combined = pd.concat(all_nt_dist_dfs, ignore_index=True)

    aa_freq_combined = None
    aa_dist_combined = None
    if all_aa_freq_dfs:
        aa_freq_combined = pd.concat(all_aa_freq_dfs, ignore_index=True)
    if all_aa_dist_dfs:
        aa_dist_combined = pd.concat(all_aa_dist_dfs, ignore_index=True)

    # Define column names
    all_seqs_columns = [
        'seq_ID', 'avg_quality_score', 'reference_name', 'NT_substitutions', 'NT_substitutions_count',
        'NT_insertions', 'NT_deletions', 'NT_insertion_length', 'NT_deletion_length'
    ]

    if do_aa_analysis:
        all_seqs_columns.extend([
            'AA_substitutions_nonsynonymous', 'AA_substitutions_synonymous',
            'AA_substitutions_nonsynonymous_count'
        ])

    if has_barcode_column:
        all_seqs_columns.append('label_barcodes')

    # Columns for grouping genotypes (skip seq_ID, avg_quality_score, reference_name, and label_barcodes)
    genotypes_columns = [col for col in all_seqs_columns[3:] if col != 'label_barcodes']

    # Build DataFrames
    all_seqs_df = pd.DataFrame(all_genotypes_lists, columns=all_seqs_columns)
    failures_df = pd.DataFrame(all_failures_lists, columns=['seq_ID', 'failure_reason', 'failure_index'])

    # Process genotypes per reference
    genotypes_condensed_list = []

    for ref_name in all_seqs_df['reference_name'].unique():
        ref_mask = all_seqs_df['reference_name'] == ref_name
        ref_df = all_seqs_df[ref_mask].copy()

        # Group identical genotypes within this reference
        ref_df['genotype->seq'] = ref_df.groupby(by=genotypes_columns).ngroup()
        ref_df['count'] = ref_df.groupby(by='genotype->seq')['seq_ID'].transform('count')

        # Create condensed genotypes DataFrame
        ref_genotypes_condensed = (ref_df
                              .sort_values('avg_quality_score', ascending=False)
                              .drop_duplicates('genotype->seq', keep='first')
                              [['genotype->seq', 'count'] + all_seqs_columns])

        ref_genotypes_condensed = ref_genotypes_condensed.sort_values(['count', 'NT_substitutions_count'],
                                                              ascending=[False, True])

        # Move wildtype to beginning if present
        wildtype_mask = ((ref_genotypes_condensed['NT_substitutions'] == '') &
                         (ref_genotypes_condensed['NT_insertions'] == '') &
                         (ref_genotypes_condensed['NT_deletions'] == ''))
        wildtype_df = ref_genotypes_condensed[wildtype_mask]

        if len(wildtype_df) > 0:
            ref_genotypes_condensed = ref_genotypes_condensed.drop(index=wildtype_df.index)
            ref_genotypes_condensed = pd.concat([wildtype_df, ref_genotypes_condensed]).reset_index(drop=True)
            ref_genotypes_condensed.rename(index={0: 'wildtype'}, inplace=True)
        else:
            ref_genotypes_condensed.reset_index(drop=True, inplace=True)

        # Make genotype IDs 1-indexed
        if len(ref_genotypes_condensed) != 0 and ref_genotypes_condensed.index[0] == 0:
            ref_genotypes_condensed.index += 1

        # Add genotype_ID column
        ref_genotypes_condensed = ref_genotypes_condensed.reset_index().rename(columns={'index': 'genotype_ID'})
        seq_to_genotype_dict = dict(zip(ref_genotypes_condensed['genotype->seq'],
                                        ref_genotypes_condensed['genotype_ID']))
        ref_df['genotype_ID'] = ref_df['genotype->seq'].map(seq_to_genotype_dict)

        # Update all_seqs_df
        all_seqs_df.loc[ref_mask, 'genotype_ID'] = ref_df['genotype_ID'].values
        all_seqs_df.loc[ref_mask, 'count'] = ref_df['count'].values
        all_seqs_df.loc[ref_mask, 'genotype->seq'] = ref_df['genotype->seq'].values

        genotypes_condensed_list.append(ref_genotypes_condensed)

    # Concatenate all reference genotypes
    genotypes_condensed = pd.concat(genotypes_condensed_list, ignore_index=True)

    # Write alignments for high-abundance genotypes (from most abundant reference only)
    if most_abundant_ref:
        most_abundant_ref_sequences = references_dict[most_abundant_ref]
        most_abundant_genotypes = genotypes_condensed[genotypes_condensed['reference_name'] == most_abundant_ref]

        write_genotype_alignments(
            bam_path, most_abundant_genotypes, most_abundant_ref_sequences,
            highest_abundance_genotypes, desired_genotype_ids,
            quality_minimum, analyze_seqs_with_indels, do_aa_analysis,
            has_quality_scores, output_files[0]
        )

    # Prepare columns to drop/include
    condensed_drop = ['genotype->seq', 'seq_ID', 'avg_quality_score']
    all_include = ['seq_ID', 'reference_name', 'genotype_ID']

    if has_barcode_column:
        condensed_drop.append('label_barcodes')
        all_include.append('label_barcodes')

    # Write output files
    genotypes_condensed.drop(columns=condensed_drop).to_csv(output_files[1], index=False)
    all_seqs_df[all_include].to_csv(output_files[2], index=False)
    failures_df.to_csv(output_files[3], index=False)
    nt_freq_combined.to_csv(output_files[4], index=False)
    nt_dist_combined.to_csv(output_files[5], index=False)

    # AA analysis outputs
    if do_aa_analysis and aa_freq_combined is not None:
        aa_freq_combined.to_csv(output_files[6], index=False)
        aa_dist_combined.to_csv(output_files[7], index=False)


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
            'references_csv': snakemake.config['runs'][snakemake.wildcards.tag]['reference_csv'],
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
        parser.add_argument('--references-csv', required=True, help='References CSV file')
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
                os.path.join(args.output_dir, f'{bam_basename}_NT-mutations-aggregated.csv'),
                os.path.join(args.output_dir, f'{bam_basename}_NT-mutation-distribution.csv'),
            ]
            if args.do_aa_analysis:
                args.output.extend([
                    os.path.join(args.output_dir, f'{bam_basename}_AA-mutations-aggregated.csv'),
                    os.path.join(args.output_dir, f'{bam_basename}_AA-mutation-distribution.csv'),
                ])
        elif not args.output:
            parser.error('Either --output-dir or --output must be specified')

    # Process BAM file
    process_bam_file(
        args.bam,
        args.references_csv,
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
