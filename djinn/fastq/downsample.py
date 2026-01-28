import random
import pysam
from rich_click import click
from djinn.fastq.extract import extract_barcodes_fq
from djinn.utils.file_ops import print_error, which_linkedread, validate_fq_sam
from djinn.utils.fq_tools import FQRecord, CachedFQWriter

@click.command(panel = "FASTQ", no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/djinn/downsample")
@click.option("-c", "--cache-size", hidden=True, type=click.IntRange(min=1000, max_open=True), default=5000, help = "Number of cached reads for write operations")
@click.option('-d', '--downsample', type = click.FloatRange(min = 0.0000001), help = 'Number/fraction of barcodes to retain')
@click.option("-i", "--invalid", is_flag=True, default=True, help = "Include invalid barcodes in downsampling")
@click.option('-i', '--invalid', default = 1, show_default = True, type=click.FloatRange(min=0,max=1), help = "Proportion of invalid barcodes to sample")
@click.option("--threads", "-t", type = click.IntRange(min = 4, max_open=True), default=10, show_default=True, help = "Number of threads to use (BAM only)")
@click.option('--random-seed', type = click.IntRange(min = 1), help = "Random seed for sampling")
@click.argument('prefix', type = str, callback=make_dir)
@click.argument('input', required=True, type=click.Path(exists=True, readable=True, dir_okay=False, resolve_path=True), callback = validate_fq_sam, nargs=-1)
@click.help_option('--help', hidden = True)
def downsample(prefix, input, invalid, downsample, random_seed, threads, cache_size):
    """
    Downsample data by barcode
    
    Downsamples a FASTQ read pair or BAM file by barcode to keep all records containing `-d` randomly sampled barcodes.
    If `d >= 1`, the downsampling is a fixed number of barcodes, whereas `d < 1` would indicate a fraction of the total
    number of barcodes `(e.g. `-d 0.5` retains 50% of all barcodes). Use `--invalid/-i` to specify if invalid barcodes
    should be included in downsampling. Inputs can be single-end or paired-end reads and must be in haplotagging, stlfr,
    or tellseq formats.

    | `--invalid` | effect                                           |
    |:---:|:---------------------------------------------------|
    | `0` | removes all invalid barcodes from the sampling pool |
    | `1` | adds all invalid barcodes to the sampling pool |
    | 0<`i`<1| keeps `i` proprotion of invalids in the sampling pool |
    """
    if len(input) > 2:
        print_error('invalid input files', 'Inputs can be one single-ended or 2 paired-end FASTQ files.')

    if random_seed:
        random.seed(random_seed)

    from_ = which_linkedread(input[0])

    if invalid:
        barcodes = extract_barcodes_fq(from_, input, separate_invalid = False)
    else:
        barcodes, invalid = extract_barcodes_fq(from_, input, separate_invalid = True)
        # remove invalid b/c it's not being used
        del invalid
    barcodes = list(barcodes)

    n_bc = len(barcodes)
    downsample = int(n_bc * downsample) if downsample < 1 else int(downsample)
    if n_bc < downsample:
        print_error("not enough barcodes", f"The input has fewer barcodes ({n_bc}) than the requested downsampling amount ({downsample})")

    with open(f"{prefix}.bc", "w") as bc_out:
        random.shuffle(barcodes)
        barcodes = barcodes[:downsample]
        bc_out.write("\n".join(barcodes))

    with CachedFQWriter(prefix, cache_size) as writer:
        for i,j in enumerate(input):
            _fw = i == 0
            with pysam.FastxFile(j, persist=False) as FQ:
                for _read in FQ:
                    _read = FQRecord(_read, _fw, from_, 0)
                    if _read.barcode in barcodes:
                        writer.queue(_read.convert(from_, _read.barcode), None)
