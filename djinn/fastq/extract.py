import pysam
import sys
from djinn.utils.fq_tools import FQRecord
from djinn.utils.file_ops import print_error, which_linkedread, validate_fq
import rich_click as click

def extract_barcodes_fq(barcode_type: str, fq: list[str], separate_invalid: bool = False):
    '''
    Return a set of the unique barcodes in the fastq files. If separate_invalid is True,
    returns a list of [set(valid), set(invalid)] barcodes.
    '''
    barcodes = set()
    invalid = set()
    for i,j in enumerate(fq):
        _fw = i == 0
        with pysam.FastxFile(j, persist=False) as FQ:
            for _read in FQ:
                _read = FQRecord(_read, _fw, barcode_type, 0)
                if separate_invalid and not _read.valid:
                    invalid.add(_read.barcode)
                else:
                    barcodes.add(_read.barcode)
   
    if separate_invalid:
        return barcodes, invalid
    return barcodes


@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/djinn/extract")
@click.argument('input', required=True, type=click.Path(exists=True, readable=True, dir_okay=False, resolve_path=True), callback = validate_fq, nargs=-1)
@click.help_option('--help', hidden = True)
def extract(input):
    '''
    Extract all barcodes

    Inputs must be one or two FASTQ files (R1 and optional R2, can be gzipped). Input
    expect barcodes to follow the standard haplotagging (BX tag), stlfr (@seq_id#barcode), or tellseq
    (@seq_id:barcode) formats. Writes to stdout.
    '''
    from_ = which_linkedread(input[0])
    bc = extract_barcodes_fq(from_, input)
    sys.stdout.write("\n".join(bc))