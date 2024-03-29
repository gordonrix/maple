"""part of nanopype-MACE pipeline, written by Gordon Rix
Plots the number of sequences output in each step
as well as the time spent for each step of the pipeline
"""

import pandas as pd
import itertools
import os
import pysam
import gzip
import math
from bokeh.plotting import figure, output_file, show, save
from bokeh.layouts import column
from bokeh.models import (BasicTicker, ColorBar, ColumnDataSource, FactorRange,
                          HoverTool, Legend, LinearColorMapper,
                          PrintfTickFormatter)
from bokeh.models.ranges import DataRange1d
from Bio import SeqIO

# generate list to be populated with rows for pipeline step, # of sequences, and time
# between pipeline steps in minutes, which will be converted to dataframe later
outList = []
tag = snakemake.input.initial.split('/')[1].split('.fastq.gz')[0]

# Get reference id
refID = list(SeqIO.parse(snakemake.config['runs'][tag]['reference'], 'fasta'))[0].id
if snakemake.config['do_UMI_analysis'][tag]:
    UMItag = snakemake.config['consensusCopyDict'][tag]
    UMIrefID = list(SeqIO.parse(snakemake.config['runs'][UMItag]['reference'], 'fasta'))[0].id


# retrieve data from each pipeline step that was performed

def linecount(filename):
    f = open(filename, 'rb')
    out = 0
    for _ in f:
        out += 1
    f.close()
    return out/4

# get time from minimap2 log file
def get_runtime(logfile):
    with open(logfile, 'r') as log:
        loglines = ''.join(log.readlines())
    before, after = 'Real time: ',  'sec; CPU'
    start = loglines.find(before) + len(before) -1
    end = loglines.find(after)
    return float(loglines[start:end])

# initial fastq file
outList.append(['initial', linecount(snakemake.input.initial), 0])

if snakemake.config['do_UMI_analysis'][tag]:

    # UMI_preconsensus_alignment
    BAMin = pysam.AlignmentFile(snakemake.input.UMI_preconsensus_alignment, 'rb')
    count = 0
    for BAMentry in BAMin.fetch(UMIrefID):
        count += 1
    time = get_runtime(snakemake.input.UMI_preconsensus_log)
    outList.append(['preconsensus alignment', count, time/60])

    # UMI extraction
    extractCSV = pd.read_csv(snakemake.input.UMI_extract, dtype={'umi':str})
    count = extractCSV['success'].sum()
    previousTimestamp = os.path.getmtime(snakemake.input.UMI_preconsensus_alignment)
    timestamp = os.path.getmtime(snakemake.input.UMI_extract)
    outList.append(['UMI extraction', count, (timestamp-previousTimestamp)/60])

    # UMI group, reports total number of raw reads that are in UMI groups above the read per UMI threshold
    groupCSV = pd.read_csv(snakemake.input.UMI_group)
    UMIreadMin = snakemake.config['UMI_consensus_minimum']
    count = groupCSV['reads_in_UMI_groups_with_n_reads'].iloc[UMIreadMin-1:].sum()
    previous_timestamp = timestamp
    timestamp = os.path.getmtime(snakemake.input.UMI_group)
    outList.append([f'UMI group (UMI read minimum {UMIreadMin})', count, (timestamp-previousTimestamp)])

    # UMI consensus, total number of consensus sequences
    previous_timestamp = timestamp
    timestamp = os.path.getmtime(snakemake.input.UMI_consensus)
    count = 0
    for line in gzip.open(snakemake.input.UMI_consensus):
        if line.decode('utf-8')[0] == '>':
            count += 1
    outList.append(['consensus', count, timestamp-previousTimestamp])

# alignment
BAMin = pysam.AlignmentFile(snakemake.input.alignment, 'rb')
count = 0
for BAMentry in BAMin.fetch(refID):
    count += 1
time = get_runtime(snakemake.input.alignment_log)
outList.append(['alignment', count, time/60])

def fail_check(output_barcodes):
    return 'fail' not in output_barcodes

# demultiplexing
if snakemake.config['do_demux'][tag]:
    demuxCSV = pd.read_csv(snakemake.input.demux)
    demuxCSVfiltered = demuxCSV[demuxCSV.apply(lambda row:
        fail_check(row['output_file_barcodes']), axis=1)]
    count = demuxCSVfiltered['demuxed_count'].sum()
    time = os.path.getmtime(snakemake.input.demux)-os.path.getmtime(snakemake.input.alignment)
    outList.append(['demux', count, time/60])    
                

outDF = pd.DataFrame(outList, columns=['pipeline_step', 'sequences', 'time'])

outDF.to_csv(snakemake.output.csv, index=False)


# plot sequence counts as bar plot

source = ColumnDataSource(outDF)
plotTitle = f'{tag}: sequence count after rule'
yLabel = 'sequences'
xCategories = list(outDF['pipeline_step'])
TOOLTIPS = [('count', '@sequences')]
seqsPlot = figure(title=plotTitle, x_range=FactorRange(*xCategories), tooltips=TOOLTIPS,
    plot_width=600, plot_height=400)
seqsPlot.vbar(x='pipeline_step', top='sequences', source=source,
    width=0.5, bottom=0, color='black')
seqsPlot.yaxis.axis_label = yLabel
seqsPlot.xaxis.major_label_orientation = math.pi/6

# plot rule time as bar plot

plotTitle = f'{tag}: time spent to complete each rule'
xLabel, yLabel = 'rule', 'time (min)'
xCategories = list(outDF['pipeline_step'])
TOOLTIPS = [('time (min)', '@time')]
timePlot = figure(title=plotTitle, x_range=FactorRange(*xCategories), tooltips=TOOLTIPS,
    plot_width=600, plot_height=400)
timePlot.vbar(x='pipeline_step', top='time', source=source,
    width=0.5, bottom=0, color='black')
timePlot.xaxis.axis_label = xLabel
timePlot.yaxis.axis_label = yLabel
timePlot.xaxis.major_label_orientation = math.pi/6


output_file(snakemake.output.plot)
save(column([seqsPlot, timePlot]))