from collections import Counter
import pysam
import rich_click as click
import sys
from djinn.utils.file_ops import validate_sam
from djinn.utils.barcodes import ANY_INVALID
import signal
if hasattr(signal, "SIGPIPE"):
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/djinn/extract")
@click.option("-i", "--invalid", is_flag=True, default=False, help = "Include invalid barcodes")
@click.argument('input', required=True, type=click.Path(exists=True, readable=True, dir_okay=False, resolve_path=True), callback = validate_sam)
@click.help_option('--help', hidden = True)
def count(input, invalid):
    '''
    Get counts for unique barcodes
    
    To avoid double-counting, all Read-1 alignments are counted, whereas only Read-2 alignments that are unpaired or contain
    a novel barcode are counted. Writes to stdout.
    '''
    barcodes = Counter()
    with pysam.AlignmentFile(input, check_sq=False) as infile:
        for record in infile.fetch(until_eof=True):
            if not record.has_tag("BX"):
                continue
            _bc = record.get_tag("BX")
            is_invalid = ANY_INVALID.search(_bc)
            if record.is_read1:
                if not is_invalid or (is_invalid and invalid):
                    barcodes.update([_bc])
            elif record.is_read2 and (not record.is_paired or _bc not in barcodes):
                if not is_invalid or (is_invalid and invalid):
                    barcodes.update([_bc])

    for k,v in barcodes.items():
        sys.stdout.write(f"{k}\t{v}\n")
