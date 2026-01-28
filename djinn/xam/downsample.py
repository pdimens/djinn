import random
import pysam
from djinn.extract import extract_barcodes_sam
from djinn.utils.file_ops import print_error

def downsample_sam(bam: str, prefix: str, downsample: int|float, keep_invalid: bool, randseed: None|int|float, threads: float) -> None:
    if randseed:
        random.seed(randseed)

    if keep_invalid:
        barcodes = extract_barcodes_sam(bam, separate_invalid = False)
    else:
        barcodes, invalid = extract_barcodes_sam(bam, separate_invalid = True)
        # rm invalid bc it's not being used
        del invalid
    barcodes = list(barcodes)

    n_bc = len(barcodes)
    random.shuffle(barcodes)

    if downsample < 1:
        downsample = int(n_bc * downsample)
    else:
        downsample = int(downsample)
        if n_bc < downsample:
            print_error("not enough barcodes", f"The input has fewer barcodes ({n_bc}) than the requested downsampling amount ({downsample})")
    barcodes = barcodes[:downsample]
    with open(f"{prefix}.bc", "w") as bc_out:
        bc_out.write("\n".join(barcodes))

    try:
        pysam.view("-O", "BAM", "-@", f"{threads-1}", "-o", f"{prefix}.bam", "-h", "-D", f"BX:{prefix}.bc", bam, catch_stdout=False)
    except pysam.SamtoolsError as e:
        print_error("samtools experienced an error", f"Filtering the input alignment file using samtools view resulted in an error. See the samtools error information below:\n{e}")
