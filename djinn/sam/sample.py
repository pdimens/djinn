import os
import random
import pysam
from djinn.sam.extract import extract_barcodes
from djinn.utils.file_ops import print_error, validate_sam
import rich_click as click

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/djinn/downsample")
@click.option('-d', '--downsample', type = click.FloatRange(min = 0.00000001), help = 'Number/fraction of barcodes to retain')
@click.option('-i', '--invalid', default = 0, show_default = True, type=click.FloatRange(min=0,max=1), help = "Proportion of invalid barcodes to sample")
@click.option("--threads", "-t", type = click.IntRange(min = 4, max_open=True), default=10, show_default=True, help = "Number of threads to use")
@click.option('--random-seed', type = click.IntRange(min = 1), help = "Random seed for sampling")
@click.argument('input', required=True, type=click.Path(exists=True, readable=True, dir_okay=False, resolve_path=True), callback = validate_sam)
@click.help_option('--help', hidden = True)
def sample(input, invalid, downsample, random_seed, threads):
    """
    Downsample data by barcode
    
    Downsamples a SAM/BAM file by barcode to keep all records containing `-d` randomly sampled barcodes.
    If `d >= 1`, the downsampling is a fixed number of barcodes, whereas `d < 1` would indicate a fraction of the total
    number of barcodes `(e.g. `-d 0.5` retains 50% of all barcodes). Use `--invalid/-i` to specify if invalid barcodes
    should be included in downsampling. Barcode must be in `BX:Z` SAM tag. Writes to stdout.

    | `--invalid` | effect                                           |
    |:---:|:---------------------------------------------------|
    | `0` | removes all invalid barcodes from the sampling pool |
    | `1` | adds all invalid barcodes to the sampling pool |
    | 0<`i`<1| keeps `i` proprotion of invalids in the sampling pool |
    """
    if random_seed:
        random.seed(random_seed)

    barcodes, _invalid = extract_barcodes(input, separate_invalid = True)
    # rm invalid bc it's not being used
    if invalid == 1:
        barcodes = barcodes.union(_invalid)
    elif invalid > 0:
        _invalid = random.sample(list(_invalid), round(len(_invalid) * invalid))
        barcodes = barcodes.union(_invalid)
    else:
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
    _input = os.path.basename(input)
    with open(f'{_input}.bc', 'w') as bc_file:
        bc_file.write("\n".join(barcodes))

    try:
        pysam.view("-O", "BAM", "-@", f"{threads-1}", "-h", "-D", f"BX:{_input}.bc", input, catch_stdout=False)
    except pysam.SamtoolsError as e:
        print_error("Samtools experienced an error", f"Filtering the input alignment file using samtools view resulted in an error. See the samtools error information below:\n{e}")
