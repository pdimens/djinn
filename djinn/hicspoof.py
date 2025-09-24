import subprocess
import os
import pysam
import rich_click as click
from djinn.utils.file_ops import compress_fq, print_error, which_linkedread
from djinn.utils.fq_tools import FQRecord, FQPool

@click.command(panel = "Conversion Commands", no_args_is_help = True, epilog = "Documentation: https://pdimens.github.io/djinn/ncbi/")
@click.option("-i", "--invalid", is_flag=True, default=False, help = "Include invalid barcodes in the output")
@click.option("-s", "--singletons", is_flag=True, default=False, help = "Include singleton barcodes in the output")
@click.option("-m", "--max-pairs", type=click.IntRange(min=1, max_open=True), default=1, show_default=True, help = "Maximum number of R2 reads per R1 per barcode")
@click.argument('prefix', required=True, type = str)
@click.argument('inputs', required=True, type=click.Path(exists=True,dir_okay=False,readable=True,resolve_path=True), nargs=2)
def hic_spoof(prefix, inputs, invalid, singletons, max_pairs):
    """
    Convert linked-read fastq into fake HI-C data

    Reads with the same barcode will have their forward reads paired with up to `-m` random 
    reverse reads to mimic the long-range data captured with HI-C. A large `-m` value may
    significantly increase output file size. The resulting fastq files will be in TELLseq-ish
    format (original barcode appended to sequence ID). See the documentation for more details.
    
    **Input FASTQ pair must be sorted by barcode and properly paired**
    """
    from_ = which_linkedread(inputs[0])
    if from_ == "none":
        print_error("unknown input format", "The input FASTQ files were not recognized as either haplotagging, stlfr, or tellseq formats. Please check that they conform to format standards for those chemistries.")

    # create the output directory in case it doesn't exist
    if os.path.dirname(prefix):
        os.makedirs(os.path.dirname(prefix), exist_ok=True)

    with (
        pysam.FastxFile(inputs[0], persist=True) as R1,
        pysam.FastxFile(inputs[1], persist=True) as R2,
        open(f"{prefix}.R1.fq.gz", "wb") as R1_out,
        open(f"{prefix}.R2.fq.gz", "wb") as R2_out,
        subprocess.Popen("gzip -c".split(), stdout= R1_out, stdin=subprocess.PIPE) as gz_r1,
        subprocess.Popen("gzip -c".split(), stdout= R2_out, stdin=subprocess.PIPE) as gz_r2,
    ):
        _fqpool = FQPool()
        for r1,r2 in zip(R1,R2):
            _r1 = FQRecord(r1, True, from_, 0)
            _r2 = FQRecord(r2, False, from_, 0)
            # if invalid barcode, do not add to pool, just convert and write
            if (not _r1.valid or not _r2.valid):
                if invalid:
                    gz_r1.stdin.write(str(_r1.convert("tellseq", _r1.barcode)).encode("utf-8"))
                    gz_r2.stdin.write(str(_r2.convert("tellseq", _r1.barcode)).encode("utf-8"))
            elif not _fqpool.barcode or _r1.barcode == _fqpool.barcode:
                # barcode pool is empty/new, so add new read pair and barcode
                _fqpool.barcode = _r1.barcode
                _fqpool.add(_r1, _r2)
            elif _fqpool.barcode and _r1.barcode != _fqpool.barcode:
                # it's a new barcode, do the spoofing, writing to the out files and resetting the pool
                _fqpool.spoof_hic(gz_r1, gz_r2, singletons, max_pairs)
                _fqpool.add(_r1, _r2)
        # convert last record pool manually, since it would be missed by the loop
        _fqpool.spoof_hic(gz_r1, gz_r2, singletons, max_pairs)

#    compress_fq(f"{prefix}.R1.fq", f"{prefix}.R2.fq")