import pysam
import sys
from djinn.utils.barcodes import ANY_INVALID, TELLSEQ_STLFR_RX
from djinn.utils.file_ops import validate_sam
import rich_click as click

def extract_barcodes(bamfile: str, separate_invalid: bool = False):
    '''
    Return a set of the unique barcodes in bamfile. If return_invalid is True,
    returns a tuple of set(valid),set(invalid) barcodes.
    '''
    barcodes = set()
    invalid = set()
    with pysam.AlignmentFile(bamfile, check_sq=False) as infile:
        for record in infile.fetch(until_eof=True):
            _bc = None
            if record.has_tag("BX"):
                _bc = record.get_tag("BX")
            else:
                _check =  TELLSEQ_STLFR_RX.search(record.name)
                if _check:
                    _bc = _check.group(0)[1:]

            if _bc:
                if separate_invalid and ANY_INVALID.search(_bc):
                    invalid.add(_bc)
                barcodes.add(_bc)

    if separate_invalid:
        return barcodes, invalid
    return barcodes


@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/djinn/extract")
@click.argument('input', required=True, type=click.Path(exists=True, readable=True, dir_okay=False, resolve_path=True), callback = validate_sam)
@click.help_option('--help', hidden = True)
def extract(input):
    '''
    Extract all barcodes

    Inputs must be one SAM/BAM file or two FASTQ files (R1 and R2, can be gzipped). Both FASTQ and SAM/BAM
    inputs expect barcodes to follow the standard haplotagging (BX tag), stlfr (@seq_id#barcode), or tellseq
    (@seq_id:barcode) formats.  Writes to stdout.
    '''
    bc = extract_barcodes(input)
    sys.stdout.write("\n".join(bc))
