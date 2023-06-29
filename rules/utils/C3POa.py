#!/usr/bin/env python3
# Roger Volden

import os
import sys
import numpy as np
import argparse
import multiprocessing as mp
import mappy as mm
from conk import conk
import gc
import gzip
import time

PATH = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/bin/'
sys.path.append(os.path.abspath(PATH))

C3POaPath = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/'
abpoa=C3POaPath+'/abPOA-v1.4.1/bin/abpoa'
racon='racon'
blat='blat'

from preprocess import preprocess
from call_peaks import call_peaks
from determine_consensus import determine_consensus

VERSION = "v3 - It's over Anakin. I have the consensus!"

def parse_args():
    '''Parses arguments.'''
    parser = argparse.ArgumentParser(description='Makes consensus sequences from R2C2 reads.',
                                     add_help=True,
                                     prefix_chars='-')
    parser.add_argument('--reads', '-r', type=str, action='store',
                          help='FASTQ file that contains the long R2C2 reads or a folder containing multiple of these FASTQ files.')
    parser.add_argument('--splint_file', '-s', type=str, action='store',
                          help='Path to the splint FASTA file.')
    parser.add_argument('--out_path', '-o', type=str, action='store', default=os.getcwd(),
                        help='''Directory where all the files will end up.
                                Defaults to your current directory.''')
    parser.add_argument('--lencutoff', '-l', type=int, action='store', default=1000,
                        help='''Sets the length cutoff for your raw sequences. Anything
                                shorter than the cutoff will be excluded. Defaults to 1000.''')
    parser.add_argument('--mdistcutoff', '-d', type=int, action='store', default=500,
                        help='''Sets the median distance cutoff for consensus sequences.
                                Anything shorter will be excluded. Defaults to 500.''')
    parser.add_argument('--zero', '-z', action='store_false', default=True,
                        help='Use to exclude zero repeat reads. Defaults to True (includes zero repeats).')
    parser.add_argument('--numThreads', '-n', type=int, default=1,
                        help='Number of threads to use during multiprocessing. Defaults to 1.')
    parser.add_argument('--groupSize', '-g', type=int, default=100000,
                        help='Number of reads processed by each thread in each iteration. Defaults to 100000.')
    parser.add_argument('--blatThreads', '-b', action='store_true', default=False,
                        help='''Use to chunk blat across the number of threads instead of by groupSize (faster).''')
    parser.add_argument('--compress_output', '-co', action='store_true', default=False,
                        help='Use to compress (gzip) both the consensus fasta and subread fastq output files.')

    parser.add_argument('--resume', '-u', action='store_true', default=False,
                        help='''If set, C3POa will look for c3poa.log file in output directory. 
                                If c3poa.log exists, files marked as processed in the file will be skipped. 
                                Output will be appended to existing output files.''')
    parser.add_argument('--peakFinderSettings', '-p', action='store', default='20,3,41,2',
                        help='Only set this if you have a really short splint (<50nt) and all your reads are discarded. Defaults to "20,3,41,2". Try "30,3,15,2" for a short splint. No promises though. We only tested C3POa for splints >100nt')

    parser.add_argument('--version', '-v', action='version', version=VERSION, help='Prints the C3POa version.')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    return parser.parse_args()

def getFileList(query_path,done):
    file_list=[]
    '''Takes a path and returns a list of fast5 files excluding the most recent fast5.'''
    for file in os.listdir(query_path):
        if 'fastq' in file and 'temp' not in file:
            if os.path.abspath(query_path+file) not in done:
                file_list.append(os.path.abspath(query_path+file))
    if file_list:
        exclude = max(file_list, key=os.path.getctime)
        if time.time()-os.path.getctime(exclude)<600:
            file_list.remove(exclude)
        file_list.sort(key=lambda x:os.path.getctime(x))
    return file_list



def rounding(x, base):
    '''Rounds to the nearest base, we use 50'''
    return int(base * round(float(x) / base))
def analyze_reads(args, read, splint, read_adapter, racon, tmp_dir,abpoa):

    peakFinderSettings=args.peakFinderSettings.split(',')
    penalty, iters, window, order = int(peakFinderSettings[0]),int(peakFinderSettings[1]),int(peakFinderSettings[2]),int(peakFinderSettings[3])
#        penalty, iters, window, order = 30, 3, 15, 2 #shorter splint
    name, seq, qual = read[0], read[1], read[2]
    seq_len = len(seq)
    use=False
    read_consensus = ''
    subs=[]
    scores = conk.conk(splint, seq, penalty)
    start=time.time()
    peaks = call_peaks(scores, args.mdistcutoff, iters, window, order)
    if  list(peaks):
        peaks = list(peaks + len(splint) // 2)
        for i in range(len(peaks) - 1, -1, -1):
            if peaks[i] >= seq_len:
                del peaks[i]
        if peaks:
            use=True

        # check for outliers in subread length
    if use:
        subreads, qual_subreads, dangling_subreads, qual_dangling_subreads = [], [], [], []
        if len(peaks) > 1:
            subread_lens = np.diff(peaks)
            subread_lens = [rounding(x, 50) for x in subread_lens]
            median_subread_len = np.median(subread_lens)
            for i in range(len(subread_lens)):
                bounds = [peaks[i], peaks[i+1]]
                if median_subread_len*0.8 <= subread_lens[i] <= median_subread_len*1.2:
                    subreads.append(seq[bounds[0]:bounds[1]])
                    qual_subreads.append(qual[bounds[0]:bounds[1]])
            if peaks[0] > 100:
                dangling_subreads.append(seq[:peaks[0]])
                qual_dangling_subreads.append(qual[:peaks[0]])
            if seq_len - peaks[-1] > 100:
                dangling_subreads.append(seq[peaks[-1]:])
                qual_dangling_subreads.append(qual[peaks[-1]:])
        else:
            dangling_subreads.append(seq[:peaks[0]])
            qual_dangling_subreads.append(qual[:peaks[0]])
            dangling_subreads.append(seq[peaks[0]:])
            qual_dangling_subreads.append(qual[peaks[0]:])


        consensus, repeats, subs = determine_consensus(
                args, read, subreads, qual_subreads, dangling_subreads, qual_dangling_subreads,
                racon, tmp_dir,abpoa
            )
        if consensus:
            numbers=[]
            for Q in qual:
                numbers.append(ord(Q)-33)
            avg_qual=str(round(np.average(numbers),1))
            cons_len = len(consensus)
            read_consensus ='>'+name+'_'+str(avg_qual)+'_'+str(seq_len)+'_'+str(repeats)+'_'+str(cons_len)+'\n'+consensus+'\n'
    return read_consensus,read_adapter,subs,peaks


def create_files(adapter,args,outDict,outSubDict,outCountDict,previous,resume):
    outCountDict[adapter]=0
    if not os.path.exists(args.out_path + adapter):
        os.mkdir(args.out_path + adapter)
    outCountDict[adapter]=0
    if args.compress_output:
        if resume:
            writeMode='ab+'
        elif adapter not in previous:
            writeMode='wb+'
        else:
            writeMode='ab+'
        outDict[adapter]=gzip.open(args.out_path + adapter +'/R2C2_Consensus.fasta.gz',writeMode)
        outSubDict[adapter]=gzip.open(args.out_path + adapter +'/R2C2_Subreads.fastq.gz',writeMode)
    else:
        if resume:
            writeMode='a'
        elif adapter not in previous:
            writeMode='w'
        else:
            writeMode='a'
        outDict[adapter]=open(args.out_path + adapter +'/R2C2_Consensus.fasta',writeMode)
        outSubDict[adapter]=open(args.out_path + adapter +'/R2C2_Subreads.fastq',writeMode)
    previous.add(adapter)
    return outDict,outSubDict,outCountDict,previous


def main(args):
    print(f'C3POa {VERSION} \nGenerating consensus sequences from R2C2 read data')
    splint_dict = {}
    for splint in mm.fastx_read(args.splint_file, read_comment=False):
        splint_dict[splint[0]] = [splint[1]]
        splint_dict[splint[0]].append(mm.revcomp(splint[1]))


    if not args.out_path.endswith('/'):
        args.out_path += '/'
    if not os.path.exists(args.out_path):
        os.mkdir(args.out_path)


    resume=args.resume
    done=set()
    if resume:
        print(f'--resume option is True: Looking for existing log file in {args.out_path}')
        if os.path.isfile(args.out_path + 'c3poa.log'):
            print('log file found')
            for line in open(args.out_path + 'c3poa.log'):
                if line.startswith('processed '):
                    processed = line.strip()[10:]
                    done.add(os.path.abspath(processed))
        print(f'{len(done)} processed files found in log file. They will be skipped')

    log_file = open(args.out_path + 'c3poa.log', 'a+')
    log_file.write(f'C3POa version: {VERSION}\n')
    iterate=True
    timeAtLastFile=time.time()
    timeSinceLastFile=0
    first=True
    previous=set()
    pool = mp.Pool(args.numThreads)
    consNumberTotal=0
    file_list = []
    while iterate:
        fileTimes=[]
        fileStart=time.time()
        log_file.write('new iteration\n')
        tmp_dir = args.out_path + 'tmp/'
        if not os.path.isdir(tmp_dir):
            os.mkdir(tmp_dir)
        else:
            os.system(f'rm -r {tmp_dir}')
            os.mkdir(tmp_dir)

        outDict={}
        outSubDict={}
        outCountDict={}

        print('Starting consensus calling iteration - if input is directory it will check for files that are new since last iteration')


        input_path=args.reads
        if os.path.isdir(input_path):
             print('\tRead input is a folder')
             file_list=getFileList(input_path,done)

        elif os.path.isfile(input_path):
            print('\tRead input is a file')
            file_list.append(input_path)
            iterate=False
        else:
            print('\tno file provided')
            iterate=False
        print(f'\t{len(file_list)} file(s) provided')

        if len(file_list)==0:
            timeSinceLastFile=time.time()-timeAtLastFile
            print(f'\t{round(timeSinceLastFile/60,2)} minutes since last file was provided. Will terminate if more than 30 minutes')
            time.sleep(30)
            if timeSinceLastFile>1800:
                iterate=False
            continue
        else:
            timeAtLastFile=time.time()
            total_reads=0
            short_reads = 0
            no_splint_reads=0

            outDict={}
            outSubDict={}
            outCountDict={}
            fileCounter=0
            totalFileCount=len(file_list)
            log_file.write(f'Total files to process: {totalFileCount}\n')
            for reads in file_list:
                fileCounter+=1
                adapter_dict, adapter_set, no_splint = preprocess(blat, args, tmp_dir, reads)
                for adapter in adapter_set:
                        outDict,outSubDict,outCountDict,previous=create_files(adapter,args,outDict,outSubDict,outCountDict,previous,resume)
                pool = mp.Pool(args.numThreads)
                results={}
                for name,seq,q in mm.fastx_read(reads, read_comment=False):
                    total_reads+=1
                    if len(seq) < args.lencutoff:
                        short_reads+=1
                    elif name not in adapter_dict:
                        no_splint_reads+=1
                    else:
                        read_adapter_info=adapter_dict[name]
                        strand = read_adapter_info[1]
                        if strand == '-':
                        # use reverse complement of the splint
                            splint = splint_dict[read_adapter_info[0]][1]
                        else:
                            splint = splint_dict[read_adapter_info[0]][0]
                        results[name]=pool.apply_async(analyze_reads,[args, [name,seq,q], splint, read_adapter_info[0], racon, tmp_dir,abpoa])

#                print(f'\tfinished file {fileCounter} of {totalFileCount}',' '*60,end='\r')
                start=time.time()
                if fileCounter%10==0:
                    pool.close()
                    pool.join()
                    gc.collect()
                    print('\t','-'*20,'restarting multithreading pool','-'*20,' '*70,end='\r')
                    pool = mp.Pool(args.numThreads)
#                print('garbage collection took ',time.time()-start)


                consNumber=0
                adapters_with_reads=set()
                for index,result in results.items():
                    consensus,adapter,subs,peaks = result.get()
                    if consensus:
                        if args.compress_output:
                             consensus=consensus.encode()
                        outDict[adapter].write(consensus)
                        outCountDict[adapter]+=1
                        consNumber+=1
                        consNumberTotal+=1
                        for subname,subseq,subq in subs:
                            entry=f'@{subname}\n{subseq}\n+\n{subq}\n'
                            if args.compress_output:
                                entry=entry.encode()
                            outSubDict[adapter].write(entry)

                log_file.write(f'Too short reads: {short_reads}'+' ({:.2f}%)'.format((short_reads/total_reads)*100)+'\n')
                log_file.write(f'No splint reads: {no_splint_reads}'+' ({:.2f}%)'.format((no_splint/total_reads)*100)+'\n')
                log_file.write(f'Successful consensus reads: {consNumber}'+' ({:.2f}%)'.format((no_splint/total_reads)*100)+'\n')

                fileEnd=time.time()
                fileTimes.append(fileEnd-fileStart)
                fileStart=fileEnd
                averageTime=round(np.average(fileTimes)/60,1)
                projectedTime = round(averageTime*(totalFileCount-fileCounter),1)
                unit='m'
                if projectedTime > 60:
                    projectedTime= round(projectTime/60,1)
                    unit ='h'
                print(f'\tFinished generating consensus sequences for file {fileCounter} of {totalFileCount} ({consNumberTotal} consensus reads total). Avg. time per file: {round(np.average(fileTimes)/60,1)}m. Estimated {projectedTime}{unit} remaining', ' '*10,end='\r')
                for adapter in adapter_set:
                    outDict[adapter].close()
                    outSubDict[adapter].close()
                    log_file.write(f'\t{outCountDict[adapter]} consensus reads generated for {adapter}\n')
#                    print(f'\t{outCountDict[adapter]} consensus reads generated for {adapter}')
                log_file.write(f'processed {os.path.abspath(reads)}\n')
                done.add(reads)
            print('\n')
    log_file.close()
    print('\n')

if __name__ == '__main__':
    args = parse_args()
    if not args.reads or not args.splint_file:
        print('Reads (--reads/-r) and splint (--splint_file/-s) are required', file=sys.stderr)
        sys.exit(1)
    mp.set_start_method("spawn")
    main(args)

