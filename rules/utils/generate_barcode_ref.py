"""
Automated generation of barcode fasta files
    For each barcode type, if generate:x option is set, will generate a fasta file containing
    x barcode sequences appearing in the specified sequence context in the provided BAM file.
    Sequences are chosen based upon frequency of appearance, with only the most frequently appearing
    sequences written to the file. If hammingDistance:h option is also set, all sequences written
    will be of hamming distance > h from each other

Input: BAM file

Output: barcodes fasta file, name defined in config file for each barcodeType
"""

import datetime
import multiprocessing as mp
import os
import re
import shutil
import time
import gzip
import itertools
import pprint
import pandas as pd

import pysam
from Bio import Align, pairwise2, Seq
from Bio.pairwise2 import format_alignment
from Bio import SeqIO

from demux import BarcodeParser

def main():

    ### Asign variables from config file and input
    config = snakemake.config
    tag = snakemake.wildcards.tag
    BAMin = snakemake.input[0]
    reference = list(SeqIO.parse(config['runs'][tag]['reference'], 'fasta'))[0]

    ### Output variables
    output = snakemake.output[0]

    bcp = BarcodeParser(config, tag)

    barcodeDict = {} # nested dictionary where keys are types of barcodes for which barcode file should be generated, subkeys are barcode sequences that are more than the set hamming distance away from all other subkey barcode sequences,
                        # subsubkeys are barcodes encountered that are all at most the set hamming distance away from each other, and values are the number of times those specific sequences were found
    
    for barcodeType in config['runs'][tag]['barcodeInfo']:
        bcTypeDict = config['runs'][tag]['barcodeInfo'][barcodeType]
        if 'generate' in bcTypeDict:
            if os.path.exists(bcTypeDict['fasta']):
                print((f"[WARNING] Barcode type `{barcodeType}` in run tag `{tag}` `generate` option is set, but barcode fasta file `{bcTypeDict['fasta']}` already exists. Not generating fasta file, delete the pre-existing fasta file if a new one should be generated"), file=sys.stderr)
            else:
                barcodeDict[barcodeType] = {}

    if len(barcodeDict.keys()) == 0:
        with open(output, 'w') as f:
            f.write('flag file for fasta file generation')
        sys.exit()

    # hammingDistanceListDict # nested dictionary where keys are barcodeType, subkeys are and values are a list of all barcode sequences that are within the set hamming distance from any barcode sequence previously encountered
        
    bamfile = pysam.AlignmentFile(BAMin, 'rb')
    for BAMentry in bamfile.fetch(reference.id):
        refAln = bcp.align_reference(BAMentry)
        for barcodeType in barcodeDict:
            barcodeName = None
            start, stop = bcp.find_N_start_end(refAln, config['runs'][tag]['barcodeInfo'][barcodeType]['context'])

            try:
                barcode = BAMentry.query_alignment_sequence[ start:stop ]
            except TypeError:
                barcodeName = 'fail'
                pass

            if barcodeName != 'fail':
                if barcode in barcodeDict[barcodeType]:
                    barcodeDict[barcodeType][barcode] += 1
                else:
                    barcodeDict[barcodeType][barcode] = 1

    barcodeDFs = {}
    # generate DataFrame of counts for all barcodes, sort by count
    for barcodeType in barcodeDict:
        barcodeRowList = []
        for barcode in barcodeDict[barcodeType]:
            barcodeRowList.append([barcode, barcodeDict[barcodeType][barcode]])
        df = pd.DataFrame(barcodeRowList, columns=['barcode', 'count'])
        df = df[df['count']>1] # ignore barcodes that only appear once
        df = df.sort_values('count', ascending=False).reset_index(drop=True)
        barcodeDFs[barcodeType] = df

    # use DataFrame of barcodes, to write fasta file
    for barcodeType, df in barcodeDFs.items():
        HDbarcodesList = [] # list to be populated with all sequences that are within the set hamming distance from any barcode that was already written to the file, to prevent hamming distance overlap
        maxBCs = config['runs'][tag]['barcodeInfo'][barcodeType]['generate']
        fileName = config['runs'][tag]['barcodeInfo'][barcodeType]['fasta']
        with open(fileName, 'w') as f:
            count = 0
            for row in df.itertuples():
                if str(row.barcode) not in HDbarcodesList: # write barcode to fasta file if it's not within the set hamming distance of any barcodes previously written to the fasta file
                    count += 1
                    f.write(f'>bc{str(count)}\n')
                    f.write(row.barcode + '\n')

                    if count >= maxBCs:
                        break

                    # if hamming distance is set, find all sequences that are within the set hamming distance of the observed barcode and add them to a list
                    if 'hammingDistance' in config['runs'][tag]['barcodeInfo'][barcodeType]:
                        hammingDistanceList = list(bcp.hamming_distance_dict(str(row.barcode), config['runs'][tag]['barcodeInfo'][barcodeType]['hammingDistance']).keys())
                    else:
                        hammingDistanceList = [str(row.barcode)]

                    HDbarcodesList.extend(hammingDistanceList)

    with open(output, 'w') as f:
        f.write('flag file for fasta file generation')            

if __name__ == '__main__':
    main()