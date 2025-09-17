import pysam
import re
import sys
import rich_click as click
from djinn.utils.file_ops import which_linkedread
from djinn.utils.barcodes import TELLSEQ_STLFR_RX
from djinn.utils.fq_tools import FQRecord

def extract_barcodes_sam(bamfile: str) -> set[str]:
    barcodes = set()
    with pysam.AlignmentFile(bamfile, check_sq=False) as infile:
        for record in infile:
            if record.has_tag("BX"):
                barcodes.add(record.get_tag("BX"))
            else:
                _check =  TELLSEQ_STLFR_RX.search(record.name)
                if _check:
                    barcodes.add(_check.group(0)[1:])
    return barcodes

def extract_barcodes_fq(barcode_type: str, fq1: str, fq2: str) -> set[str]:
    with (
        pysam.FastxFile(fq1, persist=False) as R1,
        pysam.FastxFile(fq2, persist=False) as R2,
    ):
        barcodes = set()
        for _r1 in R1:
            _r1 = FQRecord(_r1, True, barcode_type, 0)
            barcodes.add(_r1.barcode)
        for _r2 in R2:
            _r2 = FQRecord(_r2, False, barcode_type, 0)
            barcodes.add(_r2.barcode)
    return barcodes


@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/djinn/extract")
@click.argument('inputs', required=True, type=click.Path(exists=True, readable=True, dir_okay=False, resolve_path=True), nargs=-1)
def extract(inputs):
    '''
    Extract all barcodes in BAM/FASTQ file(s)

    Inputs must be one SAM/BAM file or two FASTQ files (R1 and R2, can be gzipped). Both FASTQ and SAM/BAM
    inputs expect barcodes to follow the standard haplotagging (BX tag), stlfr (@seq_id#barcode), or tellseq
    (@seq_id:barcode) formats.  Writes to stdout.
    '''
    if len(inputs) > 2:
        raise click.BadParameter('inputs must be 1 BAM file or 2 FASTQ files.', param_hint="INPUT")
    if len(inputs) == 1:
        if not inputs[0].lower().endswith(".bam"):
            raise click.BadParameter('inputs must be 1 BAM (.bam) file or 2 FASTQ (.fastq|.fq) files. The FASTQ files can be gzipped.', param_hint="INPUT")
        
        for i in extract_barcodes_sam(inputs[0]):
            sys.stdout.write(f"{i}\n")

    else:
        if inputs[0] == inputs[1]:
            raise click.BadParameter('the two input files cannot be identical', param_hint="INPUT")
        re_ext = re.compile(r"\.(fq|fastq)(?:\.gz)?$", re.IGNORECASE)
        for i in inputs:
            if not re_ext.search(i):
                raise click.BadParameter('inputs must be 1 BAM (.bam) file or 2 FASTQ (.fastq|.fq) files. The FASTQ files can be gzipped.', param_hint="INPUT")  

        from_ = which_linkedread(inputs[0])
        for i in extract_barcodes_fq(from_, inputs[0], inputs[1]):
            sys.stdout.write(f"{i}\n")
