from collections import Counter
import pysam
import rich_click as click
import sys
from djinn.utils.fq_tools import FQRecord
from djinn.utils.file_ops import which_linkedread, validate_fq

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/djinn/extract")
@click.option("-i", "--invalid", is_flag=True, default=False, help = "Include invalid barcodes")
@click.argument('input', nargs = -1, required=True, type=click.Path(exists=True, readable=True, dir_okay=False, resolve_path=True), callback = validate_fq)
@click.help_option('--help', hidden = True)
def count(input, invalid):
    '''
    Get counts for unique barcodes.
    
    R2 is optional and only the barcodes not present in R1 will be counted. Writes to stdout.
    '''
    lr = which_linkedread(input[0])
    barcodes = Counter()
    with pysam.FastxFile(input[0], persist=False) as FQ:
        for _read in FQ:
            _read = FQRecord(_read, lr, 0)
            if _read.valid or (not _read.valid and invalid):
                barcodes.update([_read.barcode])
    if len(input) > 1:
        with pysam.FastxFile(input[1], persist=False) as FQ:
            for _read in FQ:
                _read = FQRecord(_read, lr, 0)
                if _read.barcode not in barcodes:
                    if _read.valid or (not _read.valid and invalid):
                        barcodes.update([_read.barcode])

    for k,v in barcodes.items():
        sys.stdout.buffer.write(f"{k}\t{v}\n".encode())
