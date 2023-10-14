""" script for maple pipeline

Script to convert amino acid data from *_AA-mutation-frequencies.csv files into table that is usable by
dms-view, a tool for visualization of mutation data onto pdb files: https://dms-view.github.io/
"""

import pandas as pd
import os

config = snakemake.config
inputList = snakemake.input

countorfreq = 'count' if config['mutations_frequencies_raw'] else 'frequency'
cols = ['site', 'label_site', 'wildtype', 'mutation', 'condition', 'protein_chain', 'protein_site', f'mut_{countorfreq}', f'site_{countorfreq}']
proteinChain = config['dms_view_chain']
numberingDifference = config['dms_view_chain_numbering_difference']

def dmsviewDF_from_mut_data(filename):
    rows = []
    tag, barcodes = filename.split('/')[-1].split('_')[-3:-1]
    
    condition = tag+'_'+barcodes

    mutDataDF = pd.read_csv(filename, index_col=0)
    mutDataDF.loc[:, f'site_{countorfreq}'] = mutDataDF.sum(axis=1)
    for index, row in mutDataDF.iterrows():
        wtAA = index[0]
        posi = int(index[1:])
        site = row[f'site_{countorfreq}']
        mutDataCols = list(mutDataDF.columns)
        mutDataCols.remove(f'site_{countorfreq}')
        for mutAA in mutDataCols:
            mut = row[mutAA]
            if mut == 0:
                continue
            rows.append([posi, wtAA+str(posi), wtAA, mutAA, condition, proteinChain, posi, mut, site])
    df = pd.DataFrame(rows, columns=cols)
    return df

dfList = []
for f in inputList:
    df = dmsviewDF_from_mut_data(f)
    if len(df) == 0:
        continue
    dfList.append(dmsviewDF_from_mut_data(f))

dfAll = pd.concat(dfList)
dfAll.to_csv(snakemake.output[0], index=False)