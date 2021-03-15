"""part of nanopype-MACE pipeline, written by Gordon Rix
UMI_extract.py
identifies UMI barcodes in each sequence of a .bam file based on sequence context provided
in config file, appends this sequence to the sequence name, and writes these to a new .bam file,
including a GN tag that is used by the UMI_group rule to group reads"""

import gzip
import re
import pandas as pd
import numpy as np
import pysam
from Bio.Seq import reverse_complement
from Bio import SeqIO
from demux import BarcodeParser

def main():

    ### Asign variables from config file and input
    config = snakemake.config
    tag = snakemake.wildcards.tag
    BAMin = str(snakemake.input.bam)

    BAMout = snakemake.output.extracted
    logOut = snakemake.output.log

    xUMIs = UMI_Extractor(config['runs'], tag, BAMin, BAMout, logOut)
    xUMIs.extract_UMIs()

class UMI_Extractor:

    def __init__(self, runsConfig, tag, BAMin, BAMout, logOut):
        """
        arguments:

        runsConfig      - snakemake config dictionary for all runs
        tag             - tag for sequences to be UMI extracted, defined in config file
        BAMin           - BAM input file
        BAMout          - BAM file containing each sequence for which all UMIs could be identified
                            with each UMI concatenated together and added to the end of the read ID
        """
        self.tag = tag
        self.BAMin = BAMin
        self.BAMout = BAMout
        self.logOut = logOut
        self.refSeqfasta = runsConfig[tag]['reference']
        self.reference = list(SeqIO.parse(self.refSeqfasta, 'fasta'))[0]
        self.referenceSequence = str(self.reference.seq).upper()
        self.UMI_contexts = [context.upper() for context in runsConfig[tag]['UMI_contexts']]
        for context in self.UMI_contexts:
            contextIndex = self.referenceSequence.find(context)
            assert contextIndex != -1, f'UMI context not found in reference sequence. Modify context or reference sequence to ensure an exact match is present.\n\nRun tag: `{self.tag}`\nsequence context: `{context}`\nreference sequence: `{self.reference.id}`\nreference sequence fasta file: `{self.refSeqfasta}`'
            assert self.referenceSequence[contextIndex+1:].find(context) == -1, f'UMI context found in reference sequence more than once. Modify context or reference sequence to ensure only one exact match is present.\n\nRun tag: `{self.tag}`\nsequence context: `{context}`\nreference sequence: `{self.reference.id}`\nreference sequence fasta file: `{self.refSeqfasta}`'
            

    def find_N_start_end(self, sequence, context, i):
        """given a sequence with barcode locations marked as Ns, and a barcode sequence context (e.g. ATCGNNNNCCGA),
        this function will return the beginning and end of the Ns within the appropriate context if it exists only once.
        If the context does not exist or if the context appears more than once, then will assign the failure mode
        to `self.failureReason`, including the index of the context `i`, and will return 'fail', 'fail'"""

        location = sequence.find(context)
        if location == -1:
            self.logFailure[i] += 1
            return 'fail', 'fail'
        N_start = location + context.find('N')
        N_end = location + len(context) - context[::-1].find('N')

        if sequence[N_end:].find(context) == -1:
            return N_start, N_end
        else:
            self.logDict[f'context_failure_{i+1}_appears_in_alignment_more_than_once'][0] += 1
            return 'fail', 'fail'


    def align_reference(self, BAMentry):
        """given a pysam.AlignmentFile BAM entry,
        builds the reference alignment string with indels accounted for"""
        index = BAMentry.reference_start
        refAln = ''
        for cTuple in BAMentry.cigartuples:
            if cTuple[0] == 0: #match
                refAln += self.reference.seq.upper()[index:index+cTuple[1]]
                index += cTuple[1]
            elif cTuple[0] == 1: #insertion
                refAln += '-'*cTuple[1]
            elif cTuple[0] == 2: #deletion
                index += cTuple[1]
        return refAln


    def id_UMIs(self, refAln, BAMentry):
        """Inputs:
            refAln:         aligned reference string from align_reference()
            BAMentry:       pysam.AlignmentFile entry
        
        Returns the combined UMIs if all can be identified, or the reason for failure
        if they cannot"""

        UMItag = ''

        for i, context in enumerate(self.UMI_contexts):

            start,stop = self.find_N_start_end(refAln, context, i)

            try:
                UMItag += BAMentry.query_alignment_sequence[ start:stop ]
            except (TypeError or SyntaxError) as Err:
                UMItag = None

        return UMItag

    def extract_UMIs(self):
        """ loops through BAM file and uses alignments to identify sequences aligned to UMI contexts,
        and appends these sequences to the query name for use by UMI_tools group,
        and if a UMI is identified then the alignment will be written to self.BAMout
        """

        BAMin = pysam.AlignmentFile(self.BAMin, 'rb')
        BAMout = pysam.AlignmentFile(self.BAMout, 'wb', template=BAMin)

        self.logList = []
        columns = ['read_id', 'umi', 'success', 'failure'] + [f'umi__{i+1}_failure' for i,a in enumerate(self.UMI_contexts)]

        count = 0
        for BAMentry in BAMin.fetch(self.reference.id):
            refAln = self.align_reference(BAMentry)
            self.logFailure = np.zeros(len(self.UMI_contexts))
            UMIs = self.id_UMIs(refAln, BAMentry)

            if UMIs:
                self.logList.append([BAMentry.qname, UMIs, 1, 0] + list(map(int, self.logFailure)))
                BAMentry.qname = BAMentry.qname + '_' + UMIs
                BAMentry.set_tag('GN', '0', 'H')
                BAMout.write(BAMentry)
            else:
                self.logList.append([BAMentry.qname, '', 0, 1] + list(map(int, self.logFailure)))
        
        BAMin.close()
        BAMout.close()
        pysam.index(self.BAMout)

        logDF = pd.DataFrame(self.logList, columns=columns)
        logDF.to_csv(self.logOut, index=False)


if __name__ == '__main__':
    main()