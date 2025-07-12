#
#  DESCRIPTION   : Supplementary snakefile in the Maple snakemake pipeline.
#                   Rules for retrieval of sequence files from locations outside of the working directory
#
#  RESTRICTIONS  : none
#
#  AUTHOR(S)     : Gordon Rix
#

import os
import glob
import pandas as pd
import subprocess
import re

# Helper functions for retrieve_fastqs
def _is_fastq_file(filename):
    """Check if a filename is a FASTQ file."""
    return filename.endswith(('.fastq', '.fastq.gz', '.fq.gz'))

def _find_fastq_files_in_path(base_path, extensions=('*.fastq.gz', '*.fq.gz', '*.fastq')):
    """Find all FASTQ files in a given path."""
    files = []
    for ext in extensions:
        pattern = os.path.join(base_path, ext)
        files.extend(glob.glob(pattern, recursive=True))
    return files

def _search_for_specific_file(root_folder, filename):
    """Search for a specific FASTQ file in the root folder structure."""
    pattern = os.path.join(root_folder, '**', filename)
    matches = glob.glob(pattern, recursive=True)
    
    if len(matches) == 0:
        raise FileNotFoundError(f"File '{filename}' not found in '{root_folder}'")
    elif len(matches) > 1:
        raise ValueError(f"Multiple files found matching '{filename}': {matches}")
    
    return matches

def _validate_folder_structure(root_folder, folder, subfolders):
    """Validate and report on folder structure."""
    warnings = []
    
    # Check if folder exists
    folder_pattern = os.path.join(root_folder, '**', folder)
    if not glob.glob(folder_pattern, recursive=True):
        warnings.append(f"Folder '{folder}' not found in root folder '{root_folder}'")
        return warnings, False
    
    # Check for each subfolder
    for sub in subfolders:
        sub_pattern = os.path.join(root_folder, '**', folder, sub)
        if not glob.glob(sub_pattern, recursive=True):
            warnings.append(f"Subfolder '{sub}' not found under folder '{folder}' in '{root_folder}'")
    
    return warnings, True

def retrieve_fastqs(root_folder, list_of_folders_or_fqs, subfolder_string, select=''):
    """
    Search for a uniquely named folder within the root folder, and retrieve the file names 
    ending in '.fastq.gz' for all files within a list of subfolders.
    
    Parameters:
        root_folder (str):          The path to the root folder to search in.
        list_of_folders_or_fqs:     Either a list of folders to search within the root directory,
                                   or a list of fastq filenames to find directly in the root folder.
                                   Can also be a single string (folder or fastq filename).
        subfolder_string (str):     A string of comma separated subfolder names to search for within 
                                   the uniquely named folder. At least one must be present, but all 
                                   don't need to be.
        select (str):              An optional string to select a specific file name within the 
                                   subfolders. If not specified, all files ending in '.fastq.gz' 
                                   will be retrieved.
    
    Returns:
        List[str]: A list of the full file paths for all matching FASTQ files.
    """
    # Early validation
    if not os.path.isdir(root_folder):
        print(f"[WARNING] Provided root folder '{root_folder}' does not exist.")
        return []
    
    # Convert single string to list for uniform handling
    if isinstance(list_of_folders_or_fqs, str):
        list_of_folders_or_fqs = [list_of_folders_or_fqs]
    
    # Check if all items are fastq files
    if all(_is_fastq_file(item) for item in list_of_folders_or_fqs):
        # Handle list of fastq files
        found_files = []
        for fq_file in list_of_folders_or_fqs:
            pattern = os.path.join(root_folder, '**', fq_file)
            matches = glob.glob(pattern, recursive=True)
            
            if len(matches) == 0:
                raise FileNotFoundError(f"File '{fq_file}' not found in '{root_folder}'")
            elif len(matches) > 1:
                raise ValueError(f"Multiple files found matching '{fq_file}': {matches}")
            
            found_files.extend(matches)
        
        # Ensure we found exactly the same number of files as requested
        if len(found_files) != len(list_of_folders_or_fqs):
            raise ValueError(f"Expected {len(list_of_folders_or_fqs)} files but found {len(found_files)}")
        
        return found_files
    
    # If not fastq files, treat as folder names
    folder_list = list_of_folders_or_fqs
    
    # Parse subfolders
    subfolders = [s.strip() for s in subfolder_string.split(',')]
    file_paths = []
    select_matches = 0
    all_warnings = []

    # Process each folder
    for folder in folder_list:
        # Validate folder structure
        warnings, folder_exists = _validate_folder_structure(root_folder, folder, subfolders)
        for warning in warnings:
            if "Subfolder" in warning:
                print(f"[NOTICE] {warning}")
            else:
                print(f"[WARNING] {warning}")
        
        if not folder_exists:
            continue
        
        # Collect FASTQ files from all subfolders
        folder_hits = []
        for sub in subfolders:
            # Build search patterns
            search_paths = glob.glob(os.path.join(root_folder, '**', folder, sub), recursive=True)
            
            for search_path in search_paths:
                files = _find_fastq_files_in_path(search_path)
                
                # Apply select filter if specified
                if select:
                    files = [f for f in files if os.path.basename(f) == select]
                    select_matches += len(files)
                
                folder_hits.extend(files)
        
        # Report findings for this folder
        if folder_hits:
            if len(folder_hits) > 1 and select:
                print(f"[WARNING] More than one file found that matches select '{select}':")
                for f in folder_hits:
                    print(f"  {f}")
                print(f"Using all {len(folder_hits)} matching files.")
            file_paths.extend(folder_hits)
        else:
            print(f"[WARNING] No .fastq.gz files found within subfolders {subfolders} "
                  f"within folder '{folder}' in root folder '{root_folder}'.")
    
    # Final reporting
    if select and select_matches == 0:
        print(f"[ERROR] No fastq file matching '{select}' was found.")
    if not file_paths:
        print('[NOTICE] Did not find any sequences to import.')
    
    return file_paths
    

rule combine_tag: # combine batches of reads into a single file
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

checkpoint basespace_retrieve_project:
    output:
        flag = "sequences/bs{bs_project_ID}/.done"
    shell:
        """
        bs download project -i {wildcards.bs_project_ID} -o sequences/bs{wildcards.bs_project_ID} --extension fastq.gz
        touch {output.flag}
        """

def get_bs_paired_end(tag):
    """
    Retrieve paired-end fastq files from BaseSpace for a given tag.
    Assumes paired end files are located in the following directory structure after a bs project download:
        sequences/bs{bs_project_ID}/{sample_ID}*/{tag}*R1*.fastq.gz
        sequences/bs{bs_project_ID}/{sample_ID}*/{tag}*R2*.fastq.gz
    """
    project_ID = config['runs'][tag].get('bs_project_ID', None)
    if not project_ID:
        raise ValueError(f"No BaseSpace project ID found for tag {tag}")
    flag = checkpoints.basespace_retrieve_project.get(bs_project_ID=project_ID).output[0]
    sample_ID = config['runs'][tag].get('sample_ID', None)
    if not sample_ID:
        raise ValueError(f"[ERROR] Sample ID not provided for tag {tag}, but bs_project_ID is provided.")
    project_dir = os.path.dirname(flag)
    fwd = glob.glob(os.path.join(project_dir, f"{config['runs'][tag]['sample_ID']}*", f"*R1*.fastq.gz"))
    rvs = glob.glob(os.path.join(project_dir, f"{config['runs'][tag]['sample_ID']}*", f"*R2*.fastq.gz"))
    if fwd and rvs:
        return fwd[0], rvs[0]
    else:
        raise FileNotFoundError(f"Paired-end files for sample {tag} not found in {project_dir}")

rule merge_paired_end:
    input:
        fwd = lambda wildcards: retrieve_fastqs(config['sequences_dir'], [config['runs'][wildcards.tag].get('runname','')[0]], config['fastq_dir'], select=config['runs'][wildcards.tag]['fwdReads'])[0],
        rvs = lambda wildcards: retrieve_fastqs(config['sequences_dir'], [config['runs'][wildcards.tag].get('runname','')[0]], config['fastq_dir'], select=config['runs'][wildcards.tag]['rvsReads'])[0]
    output:
        merged = "sequences/paired/{tag, [^\/_]*}.fastq.gz",
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

rule merge_paired_end_from_bs:
    input:
        fwd = lambda wildcards: get_bs_paired_end(wildcards.tag)[0],
        rvs = lambda wildcards: get_bs_paired_end(wildcards.tag)[1]
    output:
        merged = "sequences/bs{bs_project_ID, [^\/_]*}_merged/{tag, [^\/_]*}.fastq.gz",
        log = "sequences/bs{bs_project_ID, [^\/_]*}_merged/{tag, [^\/_]*}_NGmerge.log",
        failedfwd = "sequences/bs{bs_project_ID, [^\/_]*}_merged/{tag, [^\/_]*}_failed-merge_1.fastq.gz",
        failedrvs = "sequences/bs{bs_project_ID, [^\/_]*}_merged/{tag, [^\/_]*}_failed-merge_2.fastq.gz"
    params:
        flags = config['NGmerge_flags'],
        failed = "sequences/bs{bs_project_ID, [^\/_]*}_merged/{tag, [^\/_]*}_failed-merge"
    shell:
        """
        NGmerge -1 {input.fwd} -2 {input.rvs} -o {output.merged} -l {output.log} -f {params.failed} -z {params.flags}
        """

def get_merged_seqs(tag):

    if any([x in config['runs'][tag] for x in ['fwdReads', 'rvsReads']]):
        return f'sequences/paired/{tag}.fastq.gz'

    elif 'runname' in config['runs'][tag]:
        return f'sequences/{tag}_combined.fastq.gz'

    elif 'bs_project_ID' in config['runs'][tag]:
        project_ID = config['runs'][tag]['bs_project_ID']
        return f'sequences/bs{project_ID}_merged/{tag}.fastq.gz'
    
    else:
        return ''

rule move_seqs: # allows for merging batches of sequences or merging paired end reads depending on the tag definition
    input:
        lambda wildcards: get_merged_seqs(wildcards.tag)
    output:
        'sequences/{tag, [^\/_]*}.fastq.gz'
    shell:
        """
        mv {input} {output}
        """