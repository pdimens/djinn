import rich_click as click
import os
import pysam
from djinn.utils.file_ops import compress_fq, print_error, which_linkedread
from djinn.utils.fq_tools import FQRecord, FQPool

@click.command(no_args_is_help = True, epilog = "Documentation: https://pdimens.github.io/djinn/ncbi/")
@click.argument('prefix', required=True, type = str)
@click.argument('inputs', required=True, type=click.Path(dir_okay=False,readable=True,resolve_path=True), nargs=2)
def hic_spoof(prefix, inputs):
    """
    Convert linked-read fastq into fake HI-C data

    Reads with the same barcode will have their forward/reverse reads combinatorally
    rearranged to mimic the long-range data captured with HI-C. The resulting
    fastq files will be in TELLseq\-_ish_ format (original barcode appended
    to sequence ID). See the documentation for more details.
    
    **Input FASTQ files must be sorted by barcode**
    
    For example, if a barcode has two read pairs [1,1'] and [2,2'], the resulting "HI-C" data will be:
    [1,1'], [1,2'], [2,1'], [2,2'].
    """
    from_ = which_linkedread(inputs[0])
    if from_ == "none":
        print_error("Error: unknown input format\nThe input FASTQ files were not recognized as either haplotagging, stlfr, or tellseq formats. Please check that they conform to format standards for those chemistries.")

    # create the output directory in case it doesn't exist
    if os.path.dirname(prefix):
        os.makedirs(os.path.dirname(prefix), exist_ok=True)

    with (
        pysam.FastxFile(inputs[0], persist=True) as R1,
        pysam.FastxFile(inputs[1], persist=True) as R2,
        open(f"{prefix}.R1.fq", "w") as R1_out,
        open(f"{prefix}.R2.fq", "w") as R2_out
    ):
        _fqpool = FQPool()
        for r1,r2 in zip(R1,R2):
            _r1 = FQRecord(r1, True, from_, 0)
            _r2 = FQRecord(r2, False, from_, 0)

            if not _fqpool.barcode or _r1.barcode == _fqpool.barcode:
                _fqpool.add(_r1, _r2)
            elif _fqpool.barcode and _r1.barcode != _fqpool.barcode:
                # it's a new barcode, do the spoofing and writing, reset the barcode pool
                _fqpool.spoof_hic(R1_out, R2_out)
                _fqpool = FQPool()
        # convert last record pool manually, since it would be missed by the loop
        _fqpool.spoof_hic(R1_out, R2_out)

    compress_fq(f"{prefix}.R1.fq", f"{prefix}.R2.fq")