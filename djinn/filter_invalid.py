#! /usr/bin/env python

from itertools import zip_longest
import os
import rich_click as click
import pysam
from djinn.utils.file_ops import make_dir, print_error, validate_fq_sam, which_linkedread, which_linkedread_sam
from djinn.utils.fq_tools import FQRecord, CachedFQWriter

@click.command(panel = "Other Tools", no_args_is_help = True, epilog = "Documentation: https://pdimens.github.io/djinn/filter/")
@click.option("-c", "--cache-size", hidden=True, type=click.IntRange(min=1000, max_open=True), default=5000, help = "Number of cached reads for write operations")
@click.option("-i", "--invalid", is_flag=True, default=False, help = "Separately output records with invalid barcodes")
@click.option("-t", "--threads", type = click.IntRange(min = 2, max_open=True), default=4, show_default=True, help = "Number of threads to use (BAM only)")
@click.argument('prefix', required=True, type = str, callback=make_dir)
@click.argument('inputs', required=True, type=click.Path(exists = True,dir_okay=False,readable=True,resolve_path=True), callback = validate_fq_sam, nargs=-1)
def filter_invalid(prefix, inputs, cache_size, invalid, threads):
    '''
    Retain only valid-barcoded reads

    **FASTQ Input**
    - expects two files, R1 and R2 (can be gzipped)
    - fastq can be in haplotagging, stlfr, or tellseq formats

    **SAM/BAM Input**
    - expects one sam/bam file
    - barcode must be in `BX:Z` SAM tag

    Use `--invalid` to separately output reads with invalid barcodes.
    '''
    lr_type = which_linkedread_sam(inputs[0]) if len(inputs) == 1 else which_linkedread(inputs[0])

    if len(inputs) == 1:
        if lr_type == "none":
            print_error("unrecognized barcode format", f"The values associated with BX tags in {inputs[0]} do not conform to haplotagging, stlfr, or tellseq/10X formats.")
        if lr_type == "haplotagging":
            invalid_rx = '[BX]!~"[ABCD]0{2,4}"'
        elif lr_type == "stlfr":
            invalid_rx = '[BX]!~"^0_|_0_|_0$"'
        else:
            invalid_rx = '[BX]!~"[N]"'

        _inv = f"--unoutput {prefix}.invalid.bam" if invalid else ""
        try:
            pysam.view("-O", "BAM", "-@", f"{threads-1}", *_inv.split(), "-o", f"{prefix}.bam", "-h", "-e", invalid_rx, inputs[0], catch_stdout=False)
        except pysam.utils.SamtoolsError as e:
            print_error("Samtools error", str(e))

    else:
        with (
            pysam.FastxFile(inputs[0], persist=False) as R1,
            pysam.FastxFile(inputs[1], persist=False) as R2,
            CachedFQWriter(prefix, cache_size) as writer,
            CachedFQWriter(f"{prefix}.invalid", cache_size) as writer_inv,
        ):
            for r1,r2 in zip_longest(R1,R2):
                if r1:
                    _r1 = FQRecord(r1, True, lr_type, 0)
                    _r1.convert(lr_type, _r1.barcode)
                    if not _r1.valid:
                        if invalid:
                            writer_inv.queue(_r1, None)
                        # spoof the record as None so it doesn't get written to valid
                        _r1 = None
                else:
                    _r1 = None
                
                if r2:
                    _r2 = FQRecord(r2, False, lr_type, 0)
                    _r2.convert(lr_type, _r2.barcode)
                    if not _r2.valid:
                        if invalid:
                            writer_inv.queue(None, _r2)
                        _r2 = None
                else:
                    _r2 = None

                writer.queue(_r1, _r2)

            # dump remaining
            writer.write()
            if invalid:
                writer_inv.write()
        
        if not invalid:
            # remove the empty invalid files if those weren't requested
            os.unlink(f"{prefix}.invalid.R1.fq.gz")
            os.unlink(f"{prefix}.invalid.R2.fq.gz")

