from Bio.Seq import Seq
from Bio import SeqIO
import Bio
import numpy as np
import pandas as pd
import re
import pysam

from common import dist_to_DF

### Asign variables from config file
config = snakemake.config
BAMin = str(snakemake.input.bam)
tag = BAMin.strip('.bam').split('/')[1].split('_')[0]
###

outputDir = 'mutation_data'

def main():
    x = MutationAnalysis(config, tag, BAMin, snakemake.output)
    x.process_seqs()

class MutationAnalysis:

    def __init__(self, config, tag, BAMin, output):
        """
        arguments:

        config          - snakemake config dictionary for all runs
        tag             - tag for which all BAM files will be demultiplexed, defined in config file
        BAMin           - BAM file input
        output          - list of output file names
        """
        refSeqfasta = config['runs'][tag]['reference']
        self.ref = list(SeqIO.parse(refSeqfasta, 'fasta'))[0]
        self.ref.seq = self.ref.seq.upper()
        self.refTrimmed = list(SeqIO.parse(refSeqfasta, 'fasta'))[1]
        self.refTrimmed.seq = self.refTrimmed.seq.upper()
        self.refStr = str(self.ref.seq)
        self.refTrimmedStr = str(self.refTrimmed.seq)
        self.config = config
        self.highestAbundanceGenotypes = config.get('highest_abundance_genotypes', 0)
        self.desiredGenotypeIDs = config.get('genotype_ID_alignments', 0)
        self.BAMin = BAMin
        self.outputList = output
        self.refTrimmedStart = self.refStr.find(self.refTrimmedStr)
        self.useReverseComplement = False
        self.use_raw_mut_count = snakemake.params.mutations_frequencies_raw
        if self.refTrimmedStart == -1:
            self.useReverseComplement = True
            self.refTrimmed.seq = self.refTrimmed.seq.reverse_complement()
            self.refTrimmedStr = str(self.refTrimmed.seq)
            self.refTrimmedStart = self.refStr.find(self.refTrimmedStr)
        self.refTrimmedEnd = len(self.refTrimmedStr) + self.refTrimmedStart
        self.analyze_seqs_with_indels = snakemake.params.analyze_seqs_with_indels

        if 'mutation_analysis_quality_score_minimum' in config:
            self.QSminimum = config['mutation_analysis_quality_score_minimum']
        self.doAAanalysis = config['do_AA_mutation_analysis'][tag]

        if self.doAAanalysis: # define protein sequence
            refProtein = list(SeqIO.parse(refSeqfasta, 'fasta'))[2]
            if self.useReverseComplement: # sequence only needs to be reversed temporarily to find start/end positions
                refProtein.seq = refProtein.seq.reverse_complement()
            refProteinStr = str(refProtein.seq).upper()
            self.refProtein = refProteinStr
            self.refProteinStart = self.refTrimmedStr.find(self.refProtein)
            self.refProteinEnd = len(self.refProtein) + self.refProteinStart
            if self.useReverseComplement:
                self.refProtein = str(Seq(refProteinStr).reverse_complement())

        self.AAs = "ACDEFGHIKLMNPQRSTVWY*"
        self.NTs = "ATGC"


    def clean_alignment(self, BAMentry):
        """given a pysam.AlignmentFile BAM entry,
        trims the ends off the query and reference sequences according to the trimmed reference,
        aligns these two strings as well as the quality scores, and creates an alignment string between the
        two trimmed sequences  where '|'=match, and '.'=mismatch, and returns these three strings, the list of quality
        scores, a list of all insertions, and a list of all deletions
        """

        if BAMentry.reference_name != self.ref.id:
            self.alignmentFailureReason = ('alignment uses wrong reference sequence', 'N/A')
            return None
        
        if BAMentry.reference_start > self.refTrimmedStart:
            self.alignmentFailureReason = ('alignment starts past trimmed reference start', 'N/A')
            return None

        if BAMentry.reference_end < self.refTrimmedEnd:
            self.alignmentFailureReason = ('alignment ends before trimmed reference end', 'N/A')
            return None

        refIndex = BAMentry.reference_start
        queryIndex = 0
        refAln = ' ' * refIndex  #pad beginning of alignment when alignment begins after beginning of reference
        queryAln = ' ' * refIndex
        queryQualities = [-1] * refIndex if self.fastq else None
        alignStr = ' ' * refIndex
        insertions = [] # list of tuples where first element is index within trimmed reference, second element is sequence inserted
        deletions = []  # list of tuples where first element is index within trimmed reference, second element is number of bases deleted

        for cTuple in BAMentry.cigartuples: # build pretty alignment through combining consecutive segments, which are determined by the cigar tuples, in which the first value determines the operation and the second value determines the length

            if cTuple[0] == 0: #no indel
                refSegment = self.refStr[refIndex:refIndex+cTuple[1]]
                querySegment = BAMentry.query_alignment_sequence[queryIndex:queryIndex+cTuple[1]]
                alignSegment = ''.join(['|' if r==q else '.' for r,q in zip(refSegment,querySegment)])
                refAln += refSegment
                queryAln += querySegment
                if self.fastq:
                    queryQualities += BAMentry.query_alignment_qualities[queryIndex:queryIndex+cTuple[1]]
                alignStr += alignSegment
                refIndex += cTuple[1]
                queryIndex += cTuple[1]

            elif cTuple[0] == 1: #insertion, not added to sequence to maintain alignment to reference
                if self.doAAanalysis and not self.analyze_seqs_with_indels and cTuple[1]%3 != 0 and self.refProteinStart <= refIndex < self.refProteinEnd: # frameshift, discard sequence if protein sequence analysis is being done and indel sequences are being ignored
                    self.alignmentFailureReason = ('frameshift insertion', queryIndex)
                    return None
                if self.refTrimmedStart <= refIndex < self.refTrimmedEnd: # record insertions as tuples of position and sequence
                    insertions.append((refIndex-self.refTrimmedStart, BAMentry.query_alignment_sequence[queryIndex:queryIndex+cTuple[1]]))
                queryIndex += cTuple[1]

            elif cTuple[0] == 2: #deletion, '-' added to sequence to maintain alignment to reference
                if self.doAAanalysis and not self.analyze_seqs_with_indels and cTuple[1]%3 != 0 and ( self.refProteinStart <= refIndex + cTuple[1] ) and ( refIndex < self.refProteinEnd ): # frameshift, discard sequence if protein sequence analysis is being done and indel sequences are being ignored
                    self.alignmentFailureReason = ('frameshift deletion', queryIndex)
                    return None
                refAln += self.refStr[refIndex:refIndex+cTuple[1]]
                queryAln += '-'*cTuple[1]
                alignStr += ' '*cTuple[1]
                if self.fastq:
                    queryQualities += [0]*cTuple[1]
                # record deletions that are present within the nucleotide analysis window
                if ( self.refTrimmedStart <= refIndex + cTuple[1] ) and ( refIndex < self.refTrimmedEnd ): # record deletions as tuples of position and length
                    deletions.append((refIndex-self.refTrimmedStart, cTuple[1]))
                refIndex += cTuple[1]

        return  [refAln[self.refTrimmedStart:self.refTrimmedEnd],
                alignStr[self.refTrimmedStart:self.refTrimmedEnd],
                queryAln[self.refTrimmedStart:self.refTrimmedEnd],
                queryQualities[self.refTrimmedStart:self.refTrimmedEnd] if self.fastq else None,
                insertions,
                deletions]

    def clean_alignment_reverse_complement(self, cleanAlignment):
        """
        given the output from `clean_alignment`, will remake each of the components using the reverse
        (for alignment string and deletions), and reverse complement (for reference, query sequence,
        and insertions) of the sequence
        """
        
        ref, alignStr, seq, qScores, insertions, deletions = cleanAlignment

        ref = str(Seq(ref).reverse_complement())
        alignStr = alignStr[::-1]
        seq = str(Seq(seq).reverse_complement())
        qScores = qScores[::-1] if qScores else None

        insertions.reverse()
        deletions.reverse()

        insertionsOut = []
        for ins in insertions:
            insertionsOut.append( (len(self.refTrimmedStr) - ins[0], str(Seq(ins[1]).reverse_complement())) )

        deletionsOut = []
        for dels in deletions:
            deletionsOut.append( (len(self.refTrimmedStr) - dels[0] - dels[1], dels[1]) )

        return [ref, alignStr, seq, qScores, insertionsOut, deletionsOut]

    def ID_muts(self, cleanAlignment):
        """ Identify mutations in an aligned sequence

        Input: output of clean_alignment(). Mutations only counted if quality score
        is above globally set threshold.

        Outputs:
            ntMutArray  - x by y array where x is the length of self.refTrimmed, representing each
                            nucleotide position within this sequence, and y is four, representing the
                            possible nucleotides that a base could be mutated to. values in the array
                            will be set to 1 if a mutation at position x is mutated to nucleotide that
                            corresponds to position y. If no mutations are found, an array of 0s is returned 
            aaMutArray  - x by y array where x is the length of self.refORF, representing each
                            amino acid position within this sequence, and y is the length of self.AAs,
                            representing the number of
                            possible nucleotides that a base could be mutated to. values in the array
                            will be set to 1 if a mutation at position x is mutated to amino acid that
                            corresponds to position y. If no mutations are found, an array of 0s is returned
            genotype    - list of strings where each string is a list of different types of mutations separated by ', '.
                            Includes AA mutations in list only if self.doAAanalysis is True
        """
        ref, alignStr, seq, qScores, insertions, deletions = cleanAlignment
        AA_insertions = []
        AA_deletions = []

        mismatches = [i for i,a in enumerate(alignStr) if a=='.']

        if self.doAAanalysis:
            indelCodons = []    # list of amino acid positions that are affected by indel (for insertion, insertion is within a codon; for deletion, at least one base of codon deleted)
            for index, insertion_seq in insertions:
                if self.refProteinStart <= index < self.refProteinEnd:
                    protIndex = index-self.refProteinStart
                    if protIndex%3 == 0:
                        if len(insertion_seq)%3 == 0:
                            AA_insertions.append((protIndex, Seq.translate(insertion_seq)))
                        continue # don't add to indelCodons if insertion occurs between codons
                    else: indelCodons.append( int(protIndex/3) )

            for index, length in deletions:
                if (self.refProteinStart <= index < self.refProteinEnd) or (self.refProteinStart <= index+length < self.refProteinEnd):
                    protIndexStart = index-self.refProteinStart
                    protIndexEnd = (index+length)-self.refProteinStart
                    firstCodon = int(protIndexStart/3)
                    lastCodon = int(protIndexEnd/3)
                    indelCodons.extend([i for i in range(firstCodon,lastCodon+1)])
                    if length%3 == 0:
                        AA_deletions.append( (firstCodon, length/3) )
            self.indelCodons = indelCodons
            AAmutArray = np.zeros((int(len(self.refProtein)/3), len(self.AAs)), dtype=int)
        else:
            AAmutArray = None

        insOutput = ', '.join([str(index)+'ins'+NTs for index,NTs in insertions])               # string of all insertions for genotype output
        delOutput = ', '.join([str(index)+'del'+str(length) for index,length in deletions])     # string of all deletions for genotype output
        AAins_output = ', '.join([str(index)+'ins'+AAs for index,AAs in AA_insertions])
        AAdel_output = ', '.join([str(index)+'del'+str(length) for index,length in deletions])

        NTmutArray = np.zeros((int(len(self.refTrimmedStr)), len(self.NTs)), dtype=int)
        codonsChecked = []
        NTsubstitutions = []
        AAnonsynonymous = []
        AAsynonymous = []

        for i in mismatches:

            if self.fastq:
                if qScores[i] < self.QSminimum: continue

            wtNT = ref[i]
            mutNT = seq[i]

            NTmutArray[i,self.NTs.find(mutNT)] += 1
            NTsubstitutions.append(wtNT+str(i+1)+mutNT) # genotype output 1-index

            if self.doAAanalysis and self.refProteinStart <= i < self.refProteinEnd:

                protIndex = i-self.refProteinStart
                codon = int(protIndex/3) # 0-index amino acid position
                if codon in codonsChecked: continue
                codonsChecked.append(codon)
                codonPosi = protIndex%3
                codonIndices = list(range(i-codonPosi, i+(3-codonPosi)))

                # check that codon doesn't contain any bases influenced by an indel
                if codon in indelCodons: continue

                #check that all three quality scores in codon are above threshold
                if qScores:
                    QStooLow = False
                    codonQS = qScores[codonIndices[0]:codonIndices[2]]
                    for qs in codonQS:
                        if qs < self.QSminimum:
                            QStooLow = True
                    if QStooLow: continue

                wtAA = str(Seq(ref[codonIndices[0]:codonIndices[2]+1]).translate())
                mutAA = str(Seq(seq[codonIndices[0]:codonIndices[2]+1]).translate())

                if wtAA!=mutAA:
                    AAmutArray[codon, self.AAs.find(mutAA)] += 1
                    AAnonsynonymous.append(wtAA+str(codon+1)+mutAA) # genotype output 1-index
                else:
                    AAsynonymous.append(wtAA+str(codon+1))

        genotype = [', '.join(NTsubstitutions)]
        genotype.append(len(NTsubstitutions))
        genotype.append(insOutput)
        genotype.append(delOutput)
        # ins_total = ''
        # del_total = ''
        # if insOutput:
        ins_total = sum([len(i_ins.split('ins')[1]) for i_ins in insOutput.split(',')]) if insOutput else 0
        del_total = sum([int(i_L.split('del')[1]) for i_L in delOutput.split(',')]) if delOutput else 0
        genotype.append(ins_total)
        genotype.append(del_total)

        if self.doAAanalysis:
            for subType in [AAnonsynonymous, AAsynonymous]:
                genotype.append(', '.join(subType))
            genotype.append(len(AAnonsynonymous))
            genotype.append(AAins_output)
            genotype.append(AAdel_output)
            
        return NTmutArray, AAmutArray, genotype

    def process_seqs(self):
        """loops through a BAM file and produces appropriate .csv files to describe mutation data.
        If config['do_AA_analysis']==False, will produce only files for NT mutation data, otherwise
        will also produce AA mutation data"""

        trimmedSeqLength = int(len(self.refTrimmedStr))
        NTmutArray = np.zeros((trimmedSeqLength, len(self.NTs)), dtype=int)         # see ID_muts() docstring
        NTmutDist = np.zeros(trimmedSeqLength, dtype=int)                           # 1D array distribution where position (x) is number of mutations and value is number of sequences with x mutations 
        failuresList = []                                                           # list of data for failed sequences that will be used to generate DataFrame
        failuresColumns = ['seq_ID', 'failure_reason', 'failure_index']             # columns for failures DataFrame
        genotypesList = []                                                          # list of genotypes to be converted into a DataFrame
        genotypesColumns = ['seq_ID', 'avg_quality_score', 'NT_substitutions', 'NT_substitutions_count', 'NT_insertions', 'NT_deletions', 'NT_insertion_length', 'NT_deletion_length'] # columns for genotypes DataFrame
        wildTypeCount = 0
        wildTypeRow = [wildTypeCount, 0, '', 0, '', '']

        if self.doAAanalysis:
            protLength = int( len(self.refProtein) / 3 )
            AAmutArray = np.zeros((protLength, len(self.AAs)), dtype=int)
            AAmutDist = np.zeros(protLength, dtype=int)
            genotypesColumns.extend(['AA_substitutions_nonsynonymous', 'AA_substitutions_synonymous', 'AA_substitutions_nonsynonymous_count', 'AA_insertions', 'AA_deletions'])
            wildTypeRow.extend(['', '', 0])

        # if any barcodes are not used to demultiplex, add a column that shows what these barcodes are
        self.barcodeColumn = False
        if self.config['do_demux'][tag]:
            for bcType in self.config['runs'][tag]['barcodeInfo']:
                if self.config['runs'][tag]['barcodeInfo'][bcType].get('noSplit', False):
                    self.barcodeColumn = True
        if self.barcodeColumn:
            genotypesColumns.append('barcode(s)')
            wildTypeRow.append('')

        ## should use SequenceAnalyzer for this
        # # if there are any mutations of interest for this tag, add genotype columns for these
        # if self.config['runs'][tag].get('NT_muts_of_interest', False):
        #     genotypesColumns.append('NT_muts_of_interest')
        #     wildTypeRow.append('')
        #     self.NT_muts_of_interest = self.config['runs'][tag]['NT_muts_of_interest'].split(', ')
        #     for mut in self.NT_muts_of_interest:
        #         genotypesColumns.append(mut)
        #         wildTypeRow.append(0)
        # if self.doAAanalysis and self.config['runs'][tag].get('AA_muts_of_interest', False):
        #     genotypesColumns.append('AA_muts_of_interest')
        #     wildTypeRow.append('')
        #     self.AA_muts_of_interest = self.config['runs'][tag]['AA_muts_of_interest'].split(', ')
        #     for mut in self.AA_muts_of_interest:
        #         genotypesColumns.append(mut)
        #         wildTypeRow.append(0)

        bamFile = pysam.AlignmentFile(self.BAMin, 'rb')

        # set whether to use quality score features based on whether or not quality scores are present
        for bamEntry in bamFile:
            self.fastq = False
            if bamEntry.query_alignment_qualities:
                self.fastq = True
            bamFile.reset()
            break
        
        for bamEntry in bamFile:
            if self.config.get('demux_screen_failures', False):
                if self.barcodeColumn and ('fail' in bamEntry.get_tag('BC')):
                    continue
            cleanAln = self.clean_alignment(bamEntry)
            if cleanAln:
                if self.useReverseComplement:
                    cleanAln = self.clean_alignment_reverse_complement(cleanAln)
                seqNTmutArray, seqAAmutArray, seqGenotype = self.ID_muts(cleanAln)                    
            else:
                failuresList.append([bamEntry.query_name, self.alignmentFailureReason[0], self.alignmentFailureReason[1]])
                continue

            seqTotalNTmuts = sum(sum(seqNTmutArray))
            NTmutArray += seqNTmutArray
            NTmutDist[seqTotalNTmuts] += 1
            if self.doAAanalysis:
                AAmutArray += seqAAmutArray
                seqTotalAAmuts = sum(sum(seqAAmutArray))
                AAmutDist[seqTotalAAmuts] += 1

            if not self.fastq:
                avgQscore = -1
            else:
                avgQscore = np.average(np.array(cleanAln[3]))
            seqGenotype = [bamEntry.query_name, avgQscore] + seqGenotype

            if self.barcodeColumn:
                seqGenotype.append(bamEntry.get_tag('BC'))
            
            genotypesList.append(seqGenotype)

        genotypesDF = pd.DataFrame(genotypesList, columns=genotypesColumns)
        failuresDF = pd.DataFrame(failuresList, columns=failuresColumns)
        
        genotypesDF['genotype->seq'] = genotypesDF.groupby(by=genotypesColumns[2:]).ngroup()        # identify duplicate genotypes, assign a unique value. This value will not be the genotype_ID, which will be assigned after further manipulation, but will be used to correlate seq_IDs to genotype_IDs
        genotypesDF['count'] = genotypesDF.groupby(by='genotype->seq')['seq_ID'].transform('count') # count duplicate genotypes, add to new column
        genotypesDFcondensed = genotypesDF.sort_values('avg_quality_score', ascending=False).drop_duplicates('genotype->seq', keep='first')[ ['genotype->seq','count']+genotypesColumns ] # remove duplicate genotype rows, keep the seq_ID with the highest average quality score
        sort = ['count','NT_substitutions_count']
        ascendBool = [False,True]
        if self.barcodeColumn:
            sort.append('barcode(s)')
            ascendBool.append(True)
        genotypesDFcondensed = genotypesDFcondensed.sort_values(sort, ascending=ascendBool)

        # move wildtype row(s) to the beginning, if they exist. rename as wild type only if there aren't any barcodes in the genotype, as this would result in many different 'wildtype' rows
        wildtypeDF = genotypesDFcondensed.loc[(genotypesDFcondensed['NT_substitutions']=='')&(genotypesDFcondensed['NT_insertions']=='')&(genotypesDFcondensed['NT_deletions']=='')]
        if len(wildtypeDF) > 0:
            wildtype_in_df = True
        else: wildtype_in_df = False
        if wildtype_in_df:
            genotypesDFcondensed = genotypesDFcondensed.drop(index=wildtypeDF.index)
            genotypesDFcondensed = pd.concat([wildtypeDF, genotypesDFcondensed]).reset_index(drop=True)
            if not self.barcodeColumn:
                genotypesDFcondensed.rename(index={0:'wildtype'}, inplace=True)
        else:
            genotypesDFcondensed.reset_index(drop=True, inplace=True)
        if (len(genotypesDFcondensed)!=0) and (genotypesDFcondensed.index[0] == 0) : # make barcode IDs 1-indexed if necessary
                genotypesDFcondensed.index += 1

        # now that genotype IDs are established, add column that correlates every sequence ID with a genotype ID from the condensed genotypes DF
        genotypesDFcondensed = genotypesDFcondensed.reset_index().rename(columns={'index':'genotype_ID'})
        seq_to_genotype_dict = dict(zip(genotypesDFcondensed['genotype->seq'], genotypesDFcondensed['genotype_ID']))
        genotypesDF['genotype_ID'] = genotypesDF['genotype->seq'].map(seq_to_genotype_dict)

        # iterate through x genotypes with highest counts and genotypes of specific ID # (both defined in config file) , get a representative sequence for each (that w highest avg_quality_score, or essentially random if there are no quality scores), and write alignments to file
        genotypeAlignmentsOutDF = genotypesDFcondensed.iloc[0:self.highestAbundanceGenotypes+1,]
        if self.desiredGenotypeIDs:
            desiredGenotypeIDs = [int(ID) for ID in str(self.desiredGenotypeIDs).split(', ') if int(ID) <= len(genotypesDFcondensed)]
            genotypeAlignmentsOutDF = pd.concat( [genotypeAlignmentsOutDF, genotypesDFcondensed.iloc[desiredGenotypeIDs,]] )
        with open(self.outputList[0], 'w') as txtOut:
            nameIndexedBAM = pysam.IndexedReads(bamFile)
            nameIndexedBAM.build()
            for row in genotypeAlignmentsOutDF.itertuples():
                if row.genotype_ID=='wildtype':
                    continue
                seqID = row.seq_ID
                iterator = nameIndexedBAM.find(seqID)
                for BAMentry in iterator:
                    break
                x = self.clean_alignment(BAMentry)
                if self.useReverseComplement:
                    x = self.clean_alignment_reverse_complement(x)
                ref, alignString, seq, _, _, _ = x
                txtOut.write(f'Genotype {row.genotype_ID} representative sequence. Sequence ID: {seqID}\n')
                for string in [ref, alignString, seq]:
                    txtOut.write(string+'\n')
                txtOut.write('\n')
            txtOut.write('')

        # output to files
        if self.useReverseComplement:
            ntIDs = list(str(Seq(self.refTrimmedStr).reverse_complement()))
        else:
            ntIDs = list(self.refTrimmedStr)
        WTnts = [f'{ID}{i+1}' for i, ID in enumerate(ntIDs)]
        NTmutDF = pd.DataFrame(NTmutArray, columns=list(self.NTs))
        NTmutDF['wt_nucleotide'] = pd.Series(WTnts)
        NTmutDF.set_index('wt_nucleotide', inplace=True)

        totalSeqs = int(NTmutDist.sum())

        NTdistDF = dist_to_DF(np.trim_zeros(NTmutDist,'b'), 'NT mutations', 'sequences')

        genotypesDFcondensed.drop(columns=['genotype->seq', 'seq_ID', 'avg_quality_score']).to_csv(self.outputList[1], index=False)
        genotypesDF.drop(columns=genotypesDF.columns.difference(['seq_ID', 'genotype_ID'])).to_csv(self.outputList[2], index=False)

        failuresDF.to_csv(self.outputList[3], index=False)
        if not self.use_raw_mut_count and totalSeqs>0:
            NTmutDF = NTmutDF.divide(totalSeqs)
        NTmutDF.to_csv(self.outputList[4])
        NTdistDF.to_csv(self.outputList[5], index=False)
        
        if self.doAAanalysis:
            resiIDs = list(str(Seq(self.refProtein).translate()))
            protLength = int(len(self.refProtein)/3)
            resiPositions = [str(i) for i in range(1, int((len(self.refProtein)/3)+1) )]
            WTresis = [ID+posi for ID,posi in zip(resiIDs,resiPositions)]
            AAmutDF = pd.DataFrame(AAmutArray, columns=list(self.AAs))
            AAmutDF['wt_residues'] = pd.Series(WTresis)
            AAmutDF.set_index('wt_residues', inplace=True)
            if not self.use_raw_mut_count and totalSeqs > 0:
                AAmutDF = AAmutDF.divide(totalSeqs)
            AAmutDF.to_csv(self.outputList[6])

            AAdistDF = dist_to_DF(np.trim_zeros(AAmutDist,'b'), 'AA mutations', 'sequences')
            AAdistDF.to_csv(self.outputList[7], index=False)

if __name__ == '__main__':
    main()