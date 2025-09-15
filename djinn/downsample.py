import random
import os
import re
import rich_click as click
import subprocess
from djinn.extract import extract_barcodes_sam, extract_barcodes_fq
from djinn.utils import which_linkedread

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/djinn/downsample")
@click.option('-b', '--barcode-tag', type = str, default = "BX", show_default = True, help = 'Tag that contains the barcode')
@click.option('-d', '--downsample', type = click.FloatRange(min = 0.0000001), help = 'Number/fraction of barcodes to retain')
@click.option('-i', '--invalid', default = 1, show_default = True, type=click.FloatRange(min=0,max=1), help = "Proportion of invalid barcodes to sample")
@click.option('-p', '--prefix', type = click.Path(exists = False), default = "downsampled", show_default = True, help = 'Prefix for output file(s)')
@click.option('--random-seed', type = click.IntRange(min = 1), help = "Random seed for sampling")
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(1,999, clamp = True), help = 'Number of threads to use')
@click.argument('inputs', required=True, type=click.Path(exists=True, readable=True, dir_okay=False, resolve_path=True), nargs=-1)
def downsample(inputs, invalid, output_dir, prefix, downsample, random_seed, threads):
    """
    Downsample data by barcode
    
    Downsamples FASTQ or BAM file(s) by barcode in the `BX` (default) tag to keep all reads
    containing `-d` randomly sampled barcodes. Downsamples by keeping `-d` random barcodes if `d >= 1`,
    otherwise, samples a fraction of the total barcodes if `0 < d < 1` (e.g. `-d 0.5` retains 50% of all barcodes).
    The SAM tag TYPE is expected to be `Z`, i.e. `--barcode-type BX` will search through the `BX:Z:` tags. 
    Use `--invalid/-i` to specify the proportion of invalid barcodes to consider for sampling. Input can be:
    - one BAM file
    - two FASTQ files (R1 and R2 of a paired-end read set)
    
    | `--invalid` | effect                                           |
    |:---:|:---------------------------------------------------|
    | `0` | removes all invalid barcodes from the sampling pool |
    | `1` | adds all invalid barcodes to the sampling pool |
    | 0<`i`<1| keeps `i` proprotion of invalids in the sampling pool |
    """
    ## checks and validations ##
    if len(inputs) > 2:
        raise click.BadParameter('inputs must be 1 BAM file or 2 FASTQ files.', param_hint="INPUT")
    if len(inputs) == 1:
        if not inputs[0].lower().endswith(".bam"):
            raise click.BadParameter('inputs must be 1 BAM (.bam) file or 2 FASTQ (.fastq|.fq) files. The FASTQ files can be gzipped.', param_hint="INPUT")
        bamfile = inputs[0]
        from_ = which_linkedread(inputs[0])
        barcodes = extract_barcodes_fq(from_, bamfile)

    if len(inputs) == 2:
        if inputs[0] == inputs[1]:
            raise click.BadParameter('the two input files cannot be identical', param_hint="INPUT")
        re_ext = re.compile(r"\.(fq|fastq)(?:\.gz)?$", re.IGNORECASE)
        for i in inputs:
            if not re_ext.search(i):
                raise click.BadParameter('inputs must be 1 BAM (.bam) file or 2 FASTQ (.fastq|.fq) files. The FASTQ files can be gzipped.', param_hint="INPUT")  
        subprocess.run(f'samtools import -@ {threads} -T "*" {inputs} -O BAM > {prefix}.bam'.split())
        bamfile = f"{prefix}.bam"    
        #TODO unclear what the SAM logic should be
        barcodes = extract_barcodes_sam(from_, bamfile)

    n_bc = len(barcodes)
    if downsample < 1:
        downsample = int(n_bc * downsample)
    else:
        downsample = int(downsample)
        if n_bc < downsample:
            raise ValueError(f"The input has fewer barcodes ({n_bc}) than the requested downsampling amount ({downsample})")

    # downsample the barcodes
    #TODO ADD THE INVALID LOGIC HERE BEFORE THE MAIN BARCODE SAMPLING
    for i in rng.sample(sorted(barcodes), downsample):
        sys.stdout.write(f"{i}\n")

    if not os.path.exists(f"{bamfile}.bai"):
        subprocess.run("samtools index {input}".split())

    #TODO THESE NEED TO BE FIXED
    subprocess.run("samtools view -O BAM -h -D {params}:{input.bc_list} {input.bam} > {output.bam}".split())
    subprocess.run("samtools fastq -@ {threads} -n -c 6 -T \"*\" -1 {output.R1} -2 {output.R2} {input}".split())