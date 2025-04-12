#!/usr/bin/env python3
"""
UMI_extract.py

Extracts UMI barcodes from each read in an aligned BAM file using a reference and provided UMI contexts.
For each read where all UMIs are identified:
  - In BAM mode: Appends the concatenated UMI to the read's query name, sets a GN tag, and writes the record to an output BAM file (which is then indexed).
  - In FASTQ mode: Writes two FASTQ files: one with full-length sequences (with the UMI appended to the read name)
    and one with only the extracted UMI sequence.
Reads for which no UMI can be extracted are omitted (though failures are logged).
"""

import argparse
import pandas as pd
import numpy as np
import pysam
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

try:
    snakemake
    snakemake_mode = True
except NameError:
    snakemake_mode = False

# Set up arguments for snakemake or command line
if snakemake_mode:
    class Args: pass
    args = Args()
    args.input_bam = snakemake.input.bam
    args.mode = snakemake.params.mode
    if args.mode == "bam":
        args.output_bam = snakemake.output.bam
    else:
        args.output_fastq_full = snakemake.output.sequences
        args.output_fastq_umi = snakemake.output.UMIs
    args.log = snakemake.output.log
    args.reference = snakemake.params.reference
    args.UMI_contexts = snakemake.params.UMI_contexts

else:
    parser = argparse.ArgumentParser(description="Extract UMIs from an aligned BAM file")
    parser.add_argument("--input_bam", required=True)
    parser.add_argument("--mode", choices=["bam", "fastq"], default="bam")
    parser.add_argument("--output_bam", help="Output BAM file (if mode is bam)")
    parser.add_argument("--output_fastq_full", help="Output full-length FASTQ file (if mode is fastq)")
    parser.add_argument("--output_fastq_umi", help="Output UMI FASTQ file (if mode is fastq)")
    parser.add_argument("--log", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--UMI_contexts", required=True, help="Comma-separated UMI contexts")
    args = parser.parse_args()

class UMI_Extractor:
    def __init__(self, BAMin, reference, UMI_contexts):
        self.BAMin = BAMin
        self.reference = list(SeqIO.parse(reference, 'fasta'))[0]
        self.referenceSequence = str(self.reference.seq).upper()
        self.UMI_contexts = [c.strip().upper() for c in UMI_contexts.split(',')]
        for context in self.UMI_contexts:
            contextIndex = self.referenceSequence.find(context)
            assert contextIndex != -1, f'UMI context not found: {context}'
            assert self.referenceSequence[contextIndex+1:].find(context) == -1, f'UMI context appears multiple times: {context}'

    def find_N_start_end(self, sequence, context, i):
        location = sequence.find(context)
        if location == -1:
            self.logFailure[i] += 1
            return 'fail', 'fail'
        N_start = location + context.find('N')
        N_end = location + len(context) - context[::-1].find('N')
        if sequence[N_end:].find(context) == -1:
            return N_start, N_end
        else:
            self.logFailure[i] += 1
            return 'fail', 'fail'

    def align_reference(self, BAMentry):
        index = BAMentry.reference_start
        refAln = ''
        for op, length in BAMentry.cigartuples:
            if op == 0:
                refAln += self.reference.seq.upper()[index:index+length]
                index += length
            elif op == 1:
                refAln += '-' * length
            elif op == 2:
                index += length
        return refAln

    def id_UMIs(self, refAln, BAMentry):
        UMItag = ''
        for i, context in enumerate(self.UMI_contexts):
            start, stop = self.find_N_start_end(refAln, context, i)
            if start == 'fail' or stop == 'fail':
                return None
            try:
                UMItag += BAMentry.query_alignment_sequence[start:stop]
            except Exception:
                return None
        return UMItag

    def extract_UMIs(self, mode, log_out, **kwargs):
        bam_in = pysam.AlignmentFile(self.BAMin, 'rb')
        self.logList = []
        cols = ['read_id','umi','success','failure'] + [f'umi_{i+1}_failure' for i in range(len(self.UMI_contexts))]
        if mode == "bam":
            bam_out_path = kwargs.get("bam_out")
            bam_out = pysam.AlignmentFile(bam_out_path, 'wb', template=bam_in)
        else:
            fastq_full = []
            fastq_umi = []
        for BAMentry in bam_in.fetch(self.reference.id):
            refAln = self.align_reference(BAMentry)
            self.logFailure = np.zeros(len(self.UMI_contexts))
            UMIs = self.id_UMIs(refAln, BAMentry)
            if UMIs:
                self.logList.append([BAMentry.qname, UMIs, 1, 0] + list(map(int, self.logFailure)))
                new_qname = BAMentry.qname + '_' + UMIs
                if mode == "bam":
                    BAMentry.qname = new_qname
                    BAMentry.set_tag('GN','0','H')
                    bam_out.write(BAMentry)
                else:
                    strand = '+' if not BAMentry.is_reverse else '-'
                    full = SeqRecord(Seq(BAMentry.query_sequence), id=new_qname, description="strand={strand}")
                    full.letter_annotations["phred_quality"] = (BAMentry.query_qualities 
                                                                 if BAMentry.query_qualities 
                                                                 else [40]*len(BAMentry.query_sequence))
                    fastq_full.append(full)
                    umi_rec = SeqRecord(Seq(UMIs), id=new_qname, description=f"strand={strand}")
                    umi_rec.letter_annotations["phred_quality"] = [35]*len(UMIs)
                    fastq_umi.append(umi_rec)
            else:
                self.logList.append([BAMentry.qname, '', 0, 1] + list(map(int, self.logFailure)))
        bam_in.close()
        pd.DataFrame(self.logList, columns=cols).to_csv(log_out, index=False)
        if mode == "bam":
            bam_out.close()
            pysam.index(bam_out_path)
        else:
            SeqIO.write(fastq_full, kwargs.get("fastq_full_out"), "fastq")
            SeqIO.write(fastq_umi, kwargs.get("fastq_umi_out"), "fastq")

def main(args):
    extractor = UMI_Extractor(args.input_bam, args.reference, args.UMI_contexts)
    if args.mode == "bam":
        extractor.extract_UMIs("bam", log_out=args.log, bam_out=args.output_bam)
    else:
        extractor.extract_UMIs("fastq", log_out=args.log,
                                 fastq_full_out=args.output_fastq_full,
                                 fastq_umi_out=args.output_fastq_umi)

if __name__ == '__main__':
    main(args)
