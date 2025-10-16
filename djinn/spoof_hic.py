import pysam
import rich_click as click
from djinn.utils.file_ops import make_dir, print_error, which_linkedread
from djinn.utils.fq_tools import FQRecord, FQBarcodePool

@click.command(panel = "File Conversions", no_args_is_help = True, epilog = "Documentation: https://pdimens.github.io/djinn/ncbi/")
@click.option("-c", "--cache-size", hidden=True, type=click.IntRange(min=1000, max_open=True), default=5000, help = "Number of cached reads for write operations")
@click.option("-i", "--invalid", is_flag=True, default=False, help = "Include invalid barcodes in the output")
@click.option("-s", "--singletons", is_flag=True, default=False, help = "Include singleton barcodes in the output")
@click.option("-m", "--max-pairs", type=click.IntRange(min=1, max_open=True), default=1, show_default=True, help = "Maximum number of R2 reads per R1 per barcode")
@click.argument('prefix', required=True, type = str, callback=make_dir)
@click.argument('inputs', required=True, type=click.Path(exists=True,dir_okay=False,readable=True,resolve_path=True), nargs=2)
def spoof_hic(prefix, inputs, invalid, singletons, max_pairs, cache_size):
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

    with (
        pysam.FastxFile(inputs[0], persist=False) as R1,
        pysam.FastxFile(inputs[1], persist=False) as R2,
        FQBarcodePool(prefix, singletons, max_pairs, cachemax=cache_size) as _fqpool
    ):
        for r1,r2 in zip(R1,R2):
            _r1 = FQRecord(r1, True, from_, 0)
            _r2 = FQRecord(r2, False, from_, 0)
            # if invalid barcode, do not add to pool, just convert and add directly to writer cache
            if (not _r1.valid or not _r2.valid):
                if invalid:
                    _fqpool.writer.queue(
                        _r1.convert2("tellseq", _r1.barcode),
                        _r2.convert2("tellseq", _r2.barcode)
                    )
                continue
            _fqpool.add(_r1,_r2)
