#! /usr/bin/env python

from itertools import zip_longest
import os
import rich_click as click
import pysam
from djinn.utils.file_ops import make_dir, print_error, validate_fq_sam, which_linkedread
from djinn.utils.fq_tools import FQRecord, CachedFQWriter

@click.command(panel = "FASTQ", no_args_is_help = True, epilog = "Documentation: https://pdimens.github.io/djinn/filter/")
@click.option("-c", "--cache-size", hidden=True, type=click.IntRange(min=1000, max_open=True), default=10000, help = "Number of cached reads for write operations")
@click.option("-i", "--invalid", is_flag=True, default=False, help = "Separately output records with invalid barcodes")
#@click.option("-t", "--threads", type = click.IntRange(min = 2, max_open=True), default=4, show_default=True, help = "Number of threads to use (BAM only)")
@click.argument('prefix', required=True, type = str, callback=make_dir)
@click.argument('input', required=True, type=click.Path(exists = True,dir_okay=False,readable=True,resolve_path=True), callback = validate_fq_sam, nargs=-1)
@click.help_option('--help', hidden = True)
def filter_invalid(prefix, input, cache_size, invalid):
    '''
    Retain only valid-barcoded reads

    Use `--invalid` to separately output reads with invalid barcodes.
    Barcodes can be in haplotagging, stlfr, or tellseq formats.
    '''
    if len(input) > 2:
        print_error('invalid input files', 'Inputs can be one single-ended or 2 paired-end FASTQ files.')

    lr_type = which_linkedread(input[0])

    if len(input) == 1:
        with (
            pysam.FastxFile(input[0], persist=False) as R1,
            CachedFQWriter(prefix, cache_size) as writer,
            CachedFQWriter(f"{prefix}.invalid", cache_size) as writer_inv,
        ):
            for r1 in R1:
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
                
                writer.queue(_r1, None)

            # dump remaining
            writer.write()
            if invalid:
                writer_inv.write()

        if not invalid:
            os.unlink(f"{prefix}.invalid.R1.fq.gz")
        return

    with (
        pysam.FastxFile(input[0], persist=False) as R1,
        pysam.FastxFile(input[1], persist=False) as R2,
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

