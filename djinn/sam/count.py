import pysam
from collections import Counter
import sys
from djinn.utils.file_ops import validate_sam
import rich_click as click

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/djinn/extract")
@click.argument('input', required=True, type=click.Path(exists=True, readable=True, dir_okay=False, resolve_path=True), callback = validate_sam)
@click.help_option('--help', hidden = True)
def count(input):
    '''
    Get counts for unique barcodes
    
    Writes to stdout.
    '''
    barcodes = Counter()
    with pysam.AlignmentFile(input, check_sq=False) as infile:
        for record in infile.fetch(until_eof=True):
            if not record.has_tag("BX"):
                continue
            _bc = record.get_tag("BX")
            barcodes.update([_bc])

    for k,v in barcodes.items():
        sys.stdout.write(f"{k}\t{v}\n")
