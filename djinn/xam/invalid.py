#! /usr/bin/env python

import rich_click as click
import pysam
from djinn.utils.file_ops import print_error, validate_sam, which_linkedread_sam

@click.command(no_args_is_help = True, epilog = "Documentation: https://pdimens.github.io/djinn/filter/")
@click.option("-i", "--invalid", type = str, default=False, help = "Output records with invalid barcodes to this file")
@click.option("-t", "--threads", type = click.IntRange(min = 2, max_open=True), default=4, show_default=True, help = "Number of threads to use")
@click.argument('input', nargs=1, required=True, type=click.Path(exists = True,dir_okay=False,readable=True,resolve_path=True), callback = validate_sam)
def filter_invalid(input, invalid, threads):
    '''
    Retain only valid-barcoded reads

    Use `--invalid` to separately output reads with invalid barcodes. Barcode must be in `BX:Z` SAM tag.
    Writes to stdout.
    '''
    lr_type = which_linkedread_sam(input[0])
    if lr_type == "none":
        print_error("unrecognized barcode format", f"The values associated with BX tags in {input[0]} do not conform to haplotagging, stlfr, or tellseq/10X formats.")
    if lr_type == "haplotagging":
        invalid_rx = '[BX]!~"[ABCD]0{2,4}"'
    elif lr_type == "stlfr":
        invalid_rx = '[BX]!~"^0_|_0_|_0$"'
    else:
        invalid_rx = '[BX]!~"[N]"'

    _inv = f"--unoutput {invalid}" if invalid else ""
    try:
        pysam.view("-O", "BAM", "-@", f"{threads-1}", *_inv.split(), "-h", "-e", invalid_rx, input[0], catch_stdout=False)
    except pysam.SamtoolsError as e:
        print_error("Samtools error", str(e))
