#! /usr/bin/env python

import rich_click as click
import pysam
from djinn.utils.file_ops import make_dir, print_error, validate_fq_sam, which_linkedread, which_linkedread_sam

@click.command(panel = "Other Tools", no_args_is_help = True, epilog = "Documentation: https://pdimens.github.io/djinn/filter/")
@click.option("-c", "--cache-size", hidden=True, type=click.IntRange(min=1000, max_open=True), default=5000, help = "Number of cached reads for write operations")
@click.option("-i", "--invalid", is_flag=True, default=False, help = "Separately output records with invalid barcodes")
@click.option("-t", "--threads", type = click.IntRange(min = 2, max_open=True), default=4, show_default=True, help = "Number of threads to use (BAM only)")
@click.argument('prefix', required=True, type = str, callback=make_dir)
@click.argument('inputs', required=True, type=click.Path(exists = True,dir_okay=False,readable=True,resolve_path=True), callback = validate_fq_sam, nargs=-1)
def filter_invalid(prefix, inputs, invalid, threads):
    '''
    Retain only valid-barcoded reads

    **FASTQ Input**
    - expects two files, R1 and R2 (can be gzipped)
    - fastq can be in haplotagging, stlfr, or tellseq formats

    **SAM/BAM Input**
    - expects one sam/bam file
    - barcode must be in `BX:Z` SAM tag

    Use `--invalid` to separately output reads with invalid barcodes.
    '''
    lr_type = which_linkedread_sam(inputs[0]) if len(inputs) == 1 else which_linkedread(inputs[0])

    if len(inputs) == 1:
        if lr_type == "none":
            print_error("unrecognized barcode format", f"The values associated with BX tags in {inputs[0]} do not conform to haplotagging, stlfr, or tellseq/10X formats.")
        if lr_type == "haplotagging":
            invalid_rx = '[BX]!~"[ABCD]0{2,4}"'
        elif lr_type == "stlfr":
            invalid_rx = '[BX]!~"^0_|_0_|_0$"'
        else:
            invalid_rx = '[BX]!~"[N]"'

        _inv = f"--unoutput {prefix}.invalid.bam" if invalid else ""
        try:
            pysam.view("-O", "BAM", "-@", f"{threads-1}", *_inv.split(), "-o", f"{prefix}.bam", "-h", "-e", invalid_rx, inputs[0], catch_stdout=False)
        except pysam.utils.SamtoolsError as e:
            print_error("Samtools error", str(e))
