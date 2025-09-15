import pysam
import sys
import rich_click as click
from djinn.utils import which_linkedread
from djinn.common import FQRecord, TELLSEQ_STLFR_RX


def extract_barcodes_sam(barcode_type: str, bamfile: str) -> set[str]:
    barcodes = set()
    with pysam.AlignmentFile(bamfile, check_sq=False) as infile:
        for record in infile:
            if record.has_tag("BX"):
                barcodes.add(record.get_tag("BX"))
            elif record.has_tag("BC"):
                barcodes.add(record.get_tag("BC"))
            elif TELLSEQ_STLFR_RX.search(record.name):
                barcodes.add(TELLSEQ_STLFR_RX.search(record.name).group(0)[1:])
    return barcodes

def extract_barcodes_fq(barcode_type: str, fq1: str, fq2: str) -> set[str]:
    with (
        pysam.FastxFile(fq1, persist=False) as R1,
        pysam.FastxFile(fq2, persist=False) as R2,
    ):
        barcodes = set()
        for _r1 in R1:
            _r1 = FQRecord(R1, True, barcode_type, 0)
            barcodes.add(_r1.barcode)
        for _r2 in R2:
            _r2 = FQRecord(R2, False, barcode_type, 0)
            barcodes.add(_r2.barcode)
    return barcodes

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/harpy/workflows/extract")
@click.argument('barcode-tag', type = str, required = True)
@click.argument('inputs', required=True, type=click.Path(exists=True, readable=True, dir_okay=False, resolve_path=True), nargs=-1)
def extract(barcode_tag, inputs):
    '''
    Extracts all the barcodes present in a SAM/BAM file. Can optionally subsample the barcodes. The invalid parameter specifies a proportion of invalid barcodes to output.
    '''
    if len(inputs) == 2:
        from_ = which_linkedread(inputs[0])
        for i in extract_barcodes_fq(from_, inputs[0], inputs[1]):
            sys.stdout.write(f"{i}\n")
    else:
        for i in extract_barcodes_sam(barcode_tag, inputs):
            sys.stdout.write(f"{i}\n")
