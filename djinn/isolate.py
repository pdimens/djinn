#! /usr/bin/env python

import os
import subprocess
import sys
import rich_click as click
import subprocess
from djinn.utils.file_ops import print_error, validate_fq_sam, which_linkedread_sam

@click.command(panel = "Other Tools", no_args_is_help = True, epilog = "Documentation: https://pdimens.github.io/djinn/isolate/")
@click.option("--threads", "-t", type = click.IntRange(min = 2, max_open=True), default=10, show_default=True, help = "Number of threads to use")
@click.argument('prefix', required=True, type = str)
@click.argument('inputs', required=True, type=click.Path(exists = True,dir_okay=False,readable=True,resolve_path=True), callback = validate_fq_sam, nargs=-1)
def ncbi(prefix, inputs, threads):
    '''
    Split a BAM file or FASTQ file pair to isolate valid barcodes from invalid barcodes.
    '''
    if len(inputs) == 1:
        lr_type = which_linkedread_sam(inputs[0])
        if lr_type == "none":
            print_error("unrecognized barcode format", f"The values associated with BX tags in {inputs[0]} do not conform to haplotagging, stlfr, or tellseq/10X formats.")
        if lr_type == "haplotagging":
            invalid = "[BX]!~\"[ABCD]0{2,4}\""
        elif lr_type == "stlfr":
            invalid = "[BX]!~\"^0_|_0_|_0$\""
        elif lr_type == "tellseq":
            invalid = "[BX]!~\"[N]\""

        subprocess.run(
            (f"samtools view -@ {threads-1} -O BAM -e '{invalid}' --unoutput " + f"{prefix}.invalid.bam {inputs[0]}").split()            
        )
