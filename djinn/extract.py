import pysam
import sys
import rich_click as click
from djinn.utils.file_ops import validate_fq_sam, which_linkedread
from djinn.utils.barcodes import ANY_INVALID, TELLSEQ_STLFR_RX
from djinn.utils.fq_tools import FQRecord

def extract_barcodes_sam(bamfile: str, separate_invalid: bool = False):
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


def extract_barcodes_fq(barcode_type: str, fq1: str, fq2: str, separate_invalid: bool = False):
    '''
    Return a set of the unique barcodes in the fastq files. If separate_invalid is True,
    returns a list of [set(valid), set(invalid)] barcodes.
    '''
    barcodes = set()
    invalid = set()
    with (
        pysam.FastxFile(fq1, persist=False) as R1,
        pysam.FastxFile(fq2, persist=False) as R2,
    ):
        for _r1 in R1:
            _r1 = FQRecord(_r1, True, barcode_type, 0)
            if separate_invalid and not _r1.valid:
                invalid.add(_r1.barcode)
            else:
                barcodes.add(_r1.barcode)
        for _r2 in R2:
            _r2 = FQRecord(_r2, False, barcode_type, 0)
            if separate_invalid and not _r2.valid:
                invalid.add(_r2.barcode)
            barcodes.add(_r2.barcode)
    
    if separate_invalid:
        return barcodes, invalid
    return barcodes

@click.command(panel = "Other Tools", no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/djinn/extract")
@click.argument('inputs', required=True, type=click.Path(exists=True, readable=True, dir_okay=False, resolve_path=True), callback = validate_fq_sam, nargs=-1)
def extract(inputs):
    '''
    Extract all barcodes in BAM/FASTQ file(s)

    Inputs must be one SAM/BAM file or two FASTQ files (R1 and R2, can be gzipped). Both FASTQ and SAM/BAM
    inputs expect barcodes to follow the standard haplotagging (BX tag), stlfr (@seq_id#barcode), or tellseq
    (@seq_id:barcode) formats.  Writes to stdout.
    '''
    if len(inputs) == 1:
        for i in extract_barcodes_sam(inputs[0]):
            sys.stdout.write(f"{i}\n")
    else:
        from_ = which_linkedread(inputs[0])
        for i in extract_barcodes_fq(from_, inputs[0], inputs[1]):
            sys.stdout.write(f"{i}\n")
