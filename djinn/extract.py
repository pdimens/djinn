import sys
import rich_click as click
from djinn.utils.file_ops import validate_fq_sam, which_linkedread

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
