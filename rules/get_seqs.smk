#
#  DESCRIPTION   : Supplementary snakefile in the Maple snakemake pipeline.
#                   Rules for retrieval of sequence files from locations outside of the working directory
#
#  RESTRICTIONS  : none
#
#  AUTHOR(S)     : Gordon Rix
#

import os

def retrieve_fastqs(rootFolder, folderList, subfolderString, select=''):
    """
    Search for a uniquely named folder within the root folder, and retrieve the file names ending in '.fastq.gz' for all files within a list of subfolders.
    
    Parameters:
        rootFolder (str):       The path to the root folder to search in.
        folderName (list(str)): A list of folders to search for files within the root directory.
                                    Each one must appear only once within the root folder
        subfolderNames (str):   A string of comma separated subfolder names to search for within the uniquely named folder.
                                    At least one must be present, but all don't need to be
        select (str):           An optional string to select a specific file name within the subfolders. If not specified, all files ending in '.fastq.gz' will be retrieved.
    
    Returns:
        List[str]: A list of the full file paths, including the root folder, for all files ending in '.fastq.gz' or 'fq.gz' within all of the subfolders found within the folders.
        If the folder or any of the subfolders were not found, or if the folder was found more than once, or if none of the subfolders were found within the folder, returns an empty list.
    """
    subfolders = subfolderString.replace(' ','').split(',')
    filePaths = []
    select_matches = 0

    for folder in folderList:

        folders_with_seqs = []
        folderFilePaths = []

        if not os.path.isdir(rootFolder):
            print(f'[WARNING] Provided root folder "{rootFolder}" does not exist.')

        # Search for the uniquely named folder within the root folder
        for root, dirs, _ in os.walk(rootFolder):
            if folder in dirs:
                # Search for each subfolder
                folderPath = os.path.join(root, folder)
                subfolderCount = 0
                for subfolder in subfolders:
                    for _, subdirs, _ in os.walk(folderPath):
                        if subfolder in subdirs:
                            # Found the subfolder, now retrieve the file names
                            subfolderPath = os.path.join(folderPath, subfolder)
                            for _, _, files in os.walk(subfolderPath):
                                for file in files:
                                    if file.endswith('.fastq.gz') or file.endswith('.fq.gz'):
                                        if os.path.join(root,folder) not in folders_with_seqs:
                                            folders_with_seqs.append(os.path.join(root,folder))
                                        subfolderCount += 1 # only care about folders / subfolders that contain the sequences in the right place
                                        if select: # only grab the file if it matches provided select argument
                                            if file == select:
                                                select_matches += 1
                                            else:
                                                continue
                                        folderFilePaths.append(os.path.join(root, folder, subfolder, file))

                if subfolderCount == 0:
                    # None of the subfolders were found, print warning
                    print(f'[NOTICE] None of the subfolders "{subfolder}" were found in directory "{os.path.join(root, folder)}".')

        if len(folders_with_seqs) > 1:
            print(f'[NOTICE] Found folder "{folder}" more than once in root folder "{rootFolder}". Merging all sequences within the designated subfolders.')

        if len(folders_with_seqs) == 0:
            print(f'[WARNING] Folder "{folder}" not found in root folder "{rootFolder}".')
        elif len(folderFilePaths) == 0:
            print(f'[WARNING] No .fastq.gz files found within subfolders {subfolderString} within folder "{folder}" within root folder "{rootFolder}".')

        filePaths.extend(folderFilePaths)
    if select:
        if select_matches == 0:
            print('[ERROR] More than one fastq file found that matches user input:\n'+''.join([f'{file}\n' for file in input]))
        elif select_matches > 1:
            print(f'[WARNING] More than one file found that matches user input "{select}":\n'+''.join([f'{file}\n' for file in folderFilePaths]))

    if len(filePaths) == 0:
        print('[NOTICE] Did not find any sequences to import. This is likely because no runname for a tag was specified. If sequence import is necessary then you may have change that.')

    return filePaths

rule merge_paired_end:
    input:
        fwd = lambda wildcards: retrieve_fastqs(config['sequences_dir'], [config['runs'][wildcards.tag].get('runname','')[0]], config['fastq_dir'], select=config['runs'][wildcards.tag]['fwdReads'])[0],
        rvs = lambda wildcards: retrieve_fastqs(config['sequences_dir'], [config['runs'][wildcards.tag].get('runname','')[0]], config['fastq_dir'], select=config['runs'][wildcards.tag]['rvsReads'])[0]
    output:
        merged = temp("sequences/paired/{tag, [^\/_]*}.fastq.gz"),
        log = "sequences/paired/{tag, [^\/_]*}_NGmerge.log",
        failedfwd = "sequences/paired/{tag, [^\/_]*}_failed-merge_1.fastq.gz",
        failedrvs = "sequences/paired/{tag, [^\/_]*}_failed-merge_2.fastq.gz"
    params:
        flags = config['NGmerge_flags'],
        failed = "sequences/paired/{tag, [^\/_]*}_failed-merge"
    shell:
        """
        NGmerge -1 {input.fwd} -2 {input.rvs} -o {output.merged} -l {output.log} -f {params.failed} -z {params.flags}
        """

rule basecaller_combine_tag: # combine batches of basecalled reads into a single file
    input:
        lambda wildcards: retrieve_fastqs(config['sequences_dir'], config['runs'][wildcards.tag].get('runname',''), config['fastq_dir'])
    output:
        temp("sequences/{tag, [^\/_]*}_combined.fastq.gz")
    run:
        with open(output[0], 'wb') as fp_out:
            if len(input)==0:
                raise RuntimeError(f"Basecalled sequence batches not found for tag `{wildcards.tag}`.")
            for f in input:
                with open(f, 'rb') as fp_in:
                    fp_out.write(fp_in.read())

rule move_seqs: # allows for merging batches of sequences or merging paired end reads depending on the tag definition using the above rules
    input:
        lambda wildcards: f'sequences/{wildcards.tag}_combined.fastq.gz' if config['merge_paired_end'][wildcards.tag]==False else f'sequences/paired/{wildcards.tag}.fastq.gz'
    output:
        'sequences/{tag, [^\/_]*}.fastq.gz'
    shell:
        """
        mv {input} {output}
        """
