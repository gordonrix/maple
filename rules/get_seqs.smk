#
#  DESCRIPTION   : Supplementary snakefile in the Maple snakemake pipeline.
#                   Rules for retrieval of sequence file from locations outside of the working directory
#
#  RESTRICTIONS  : none
#
#  AUTHOR(S)     : Gordon Rix
#

def get_batch_sequences(runpath, fastq_dir):
    batches_fastqgz, = glob_wildcards("{runpath}/{fastq_dir}/{{id}}.fastq.gz".format(runpath=runpath, fastq_dir=fastq_dir))
    return batches_fastqgz

# get batches of sequences for each provided runname to be merged together
def get_batches(wildcards):
    output = []
    for runname in config['runs'][wildcards.tag]['runname']:

        # search through all directories two directories deep within the minknowDir for the provided runname, and stop once a matching directory is found
        experimentDirs = [d for d in os.listdir(config['minknowDir']) if os.path.isdir(os.path.join(config['minknowDir'],d))]
        expt_sample = False
        for expt in experimentDirs:
            if expt_sample: break # stop once a matching directory is found
            sampleDirs = [d for d in os.listdir(os.path.join(config['minknowDir'], expt)) if os.path.isdir(os.path.join(config['minknowDir'], expt,d))]
            for sample in sampleDirs:
                if runname in os.listdir(os.path.join(config['minknowDir'], expt, sample)):
                    expt_sample = (expt, sample)
                    break

        if expt_sample: # add all batches of sequences to a list to be merged together
            outputs = []
            for fastq_dir in config['fastq_dir'].replace(' ','').split(','):
                batch = get_batch_sequences(os.path.join(config['minknowDir'].rstrip('/'), expt_sample[0], expt_sample[1], runname), fastq_dir)
                outputs.extend(expand("{minknowDir}/{expt}/{sample}/{runname}/{fastq_dir}/{batch}.fastq.gz",
                                        minknowDir = config['minknowDir'].rstrip('/'),
                                        expt = expt_sample[0],
                                        sample = expt_sample[1],
                                        runname = runname,
                                        fastq_dir = fastq_dir,
                                        batch = batch))
        else:
            print('[WARNING] No folders matching the provided runname was found. This is fine if you have already combined the sequences you want to combine but if not, then it is not fine and you should refer to the documentation.')
        output.extend(outputs)
    return output

rule merge_paired_end:
    input:
        fwd = lambda wildcards: os.path.join('sequences', 'paired', config['runs'][wildcards.tag]['fwdReads']),
        rvs = lambda wildcards: os.path.join('sequences', 'paired', config['runs'][wildcards.tag]['rvsReads'])
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
        lambda wildcards: get_batches(wildcards)
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
        lambda wildcards: f'sequences/{wildcards.tag}_combined.fastq.gz' if config['merge_paired_end'][tag]==False else f'sequences/paired/{wildcards.tag}.fastq.gz'
    output:
        'sequences/{tag, [^\/_]*}.fastq.gz'
    shell:
        """
        mv {input} {output}
        """
