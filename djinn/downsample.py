from itertools import zip_longest
import random
import pysam
import rich_click as click
from djinn.extract import extract_barcodes_sam, extract_barcodes_fq
from djinn.utils.file_ops import compress_fq, print_error, validate_fq_sam, which_linkedread
from djinn.utils.fq_tools import FQRecord

def downsample_fastq(fq1: str, fq2: str, prefix: str, downsample: int|float, keep_invalid: bool, randseed: None|int|float) -> None:
    from_ = which_linkedread(fq1)

    if randseed:
        random.seed(randseed)

    if keep_invalid:
        barcodes = extract_barcodes_fq(from_, fq1, fq2, separate_invalid = False)
    else:
        barcodes, invalid = extract_barcodes_fq(from_, fq1, fq2, separate_invalid = True)
        # remove invalid b/c it's not being used
        del invalid
    barcodes = list(barcodes)

    n_bc = len(barcodes)
    random.shuffle(barcodes)

    if downsample < 1:
        downsample = int(n_bc * downsample)
    else:
        downsample = int(downsample)
        if n_bc < downsample:
            print_error("not enough barcodes", f"The input has fewer barcodes ({n_bc}) than the requested downsampling amount ({downsample})")
    barcodes = barcodes[:downsample]
    with open(f"{prefix}.bc", "w") as bc_out:
        bc_out.write("\n".join(barcodes))

    with (
        pysam.FastxFile(fq1, persist=False) as R1,
        pysam.FastxFile(fq2, persist=False) as R2,
        open(f"{prefix}.R1.fq", "w") as R1_out,
        open(f"{prefix}.R2.fq", "w") as R2_out,
    ):
        for r1,r2 in zip_longest(R1,R2):
            if r1:
                _r1 = FQRecord(r1, True, from_, 0)
                if _r1.barcode in barcodes:
                    R1_out.write(str(_r1.convert(from_, _r1.barcode)))
            if r2:
                _r2 = FQRecord(r2, False, from_, 0)
                if _r2.barcode in barcodes:
                    R2_out.write(str(_r2.convert(from_, _r2.barcode)))

    compress_fq(f"{prefix}.R1.fq", f"{prefix}.R2.fq")

def downsample_sam(bam: str, prefix: str, downsample: int|float, keep_invalid: bool, randseed: None|int|float, threads: float) -> None:
    if randseed:
        random.seed(randseed)

    if keep_invalid:
        barcodes = extract_barcodes_sam(bam, separate_invalid = False)
    else:
        barcodes, invalid = extract_barcodes_sam(bam, separate_invalid = True)
        # rm invalid bc it's not being used
        del invalid
    barcodes = list(barcodes)

    n_bc = len(barcodes)
    random.shuffle(barcodes)

    if downsample < 1:
        downsample = int(n_bc * downsample)
    else:
        downsample = int(downsample)
        if n_bc < downsample:
            print_error("not enough barcodes", f"The input has fewer barcodes ({n_bc}) than the requested downsampling amount ({downsample})")
    barcodes = barcodes[:downsample]
    with open(f"{prefix}.bc", "w") as bc_out:
        bc_out.write("\n".join(barcodes))

    try:
        pysam.view("-O", "BAM", "-@", f"{threads-1}", "-o", f"{prefix}.bam", "-h", "-D", f"BX:{prefix}.bc", bam, catch_stdout=False)
    except pysam.SamtoolsError as e:
        print_error("samtools experienced an error", f"Filtering the input alignment file using samtools view resulted in an error. See the samtools error information below:\n{e}")

@click.command(panel = "Other Tools", no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/djinn/downsample")
@click.option('-d', '--downsample', type = click.FloatRange(min = 0.0000001), help = 'Number/fraction of barcodes to retain')
@click.option("-i", "--invalid", is_flag=True, default=True, help = "Include invalid barcodes in downsampling")
#@click.option('-i', '--invalid', default = 1, show_default = True, type=click.FloatRange(min=0,max=1), help = "Proportion of invalid barcodes to sample")
@click.option("--threads", "-t", type = click.IntRange(min = 4, max_open=True), default=10, show_default=True, help = "Number of threads to use (BAM only)")
@click.option('--random-seed', type = click.IntRange(min = 1), help = "Random seed for sampling")
@click.argument('prefix', type = click.Path(exists = False))
@click.argument('inputs', required=True, type=click.Path(exists=True, readable=True, dir_okay=False, resolve_path=True), callback = validate_fq_sam, nargs=-1)
def downsample(prefix, inputs, invalid, downsample, random_seed, threads):
    """
    Downsample data by barcode
    
    Downsamples a FASTQ read pair or BAM file by barcode to keep all records containing `-d` randomly sampled barcodes.
    If `d >= 1`, the downsampling is a fixed number of barcodes, whereas `d < 1` would indicate a fraction of the total
    number of barcodes `(e.g. `-d 0.5` retains 50% of all barcodes). Use `--invalid/-i` to specify if invalid barcodes
    should be included in downsampling.

    **FASTQ Input**
    - expects two files, R1 and R2 (can be gzipped)
    - fastq can be in haplotagging, stlfr, or tellseq formats

    ** BAM Input**
    - expects one sam/bam file
    - barcode must be in `BX:Z` SAM tag

    | `--invalid` | effect                                           |
    |:---:|:---------------------------------------------------|
    | `0` | removes all invalid barcodes from the sampling pool |
    | `1` | adds all invalid barcodes to the sampling pool |
    | 0<`i`<1| keeps `i` proprotion of invalids in the sampling pool |
    """
    if len(inputs) == 1:
        downsample_sam(inputs[0], prefix, downsample, invalid, random_seed, threads)
    else:
        downsample_fastq(inputs[0], inputs[1], prefix, downsample, invalid, random_seed)
