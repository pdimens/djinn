#! /usr/bin/env python

from collections import Counter
import pysam
import rich_click as click
import sys
from djinn.utils.barcodes import is_invalid
from djinn.utils.file_ops import validate_sam
from djinn.utils.xam_tools import bam_barcode

def count_barcodes(bamfile: str, invalid: bool = False) -> Counter:
    '''
    Return a Counter of the unique barcodes in bamfile. If invalid is True,
    the Counter includes invalid barcodes as well.
    '''
    barcodes = Counter()
    with pysam.AlignmentFile(bamfile, check_sq=False) as infile:
        for record in infile.fetch(until_eof=True):
            _bc = bam_barcode(record)
            if _bc:
                if invalid:
                    if is_invalid(_bc):
                        barcodes[_bc] += 1
                else:
                    barcodes[_bc] += 1
    return barcodes

@click.command(no_args_is_help = True, epilog = "Documentation: https://pdimens.github.io/djinn/filter/")
@click.option("-s", "--singletons", type = str, help = "Print valid singleton records to this file")
@click.argument('input', nargs=1, required=True, type=click.Path(exists = True,dir_okay=False,readable=True,resolve_path=True), callback = validate_sam)
@click.help_option('--help', hidden = True)
def filter_singletons(input, singletons):
    '''
    Retain only non-singleton reads

    This method also filters out invalid barcodes, since they are not considered linked. Barcode must be in `BX:Z`
    SAM tag. Use `--singletons` to optionally output reads with singleton barcodes into a separate BAM file whose name
    is provided for this option. Writes to stdout.
    '''
    bc_counts = count_barcodes(input[0])
    linked = list(filter(lambda x: bc_counts[x] > 2, bc_counts.keys()))
    if not linked:
        print("There are no barcodes with >2 reads. Nothing to do.")
        return

    with (
        pysam.AlignmentFile(input[0], check_sq=False) as xam,
        pysam.AlignmentFile(sys.stdout.buffer, "wb", template = xam) as xam_out,
    ):
        if singletons:
            singleton_out = pysam.AlignmentFile(singletons, "wb", template = xam)
        for record in xam.fetch(until_eof=True):
            #SEARCH FOR BARCODE
            bc = bam_barcode(record)
            if not bc:
                continue
            # IF BARCODE IN GOOD LIST, OUTPUT RECORD
            if not is_invalid(bc):
                if bc in linked:
                    xam_out.write(record)
            # IF SINGLETONS AND BC IN SINGLETON LIST, OUTPUT TO OTHER FILE
                elif singletons:
                    singleton_out.write(record)
    if singletons:
        singleton_out.close()
