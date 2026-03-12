#! /usr/bin/env python

from itertools import zip_longest
import os
import rich_click as click
import pysam
from djinn.utils.file_ops import make_dir, validate_fq, which_linkedread
from djinn.utils.fq_tools import FQRecord, CachedFQWriter

@click.command(no_args_is_help = True, epilog = "Documentation: https://pdimens.github.io/djinn/filter/")
@click.option("-c", "--cache-size", hidden=True, type=click.IntRange(min=1000, max_open=True), default=10000, help = "Number of cached reads for write operations")
@click.option("-i", "--invalid", is_flag=True, default=False, help = "Separately output records with invalid barcodes")
@click.option("-t", "--threads", type = click.IntRange(min = 1, max_open=True), default=4, show_default=True, help = "Number of compression threads to use for output files")
@click.argument('prefix', required=True, type = str, callback=make_dir)
@click.argument('input', nargs=-1, required=True, type=click.Path(exists = True,dir_okay=False,readable=True,resolve_path=True), callback = validate_fq)
@click.help_option('--help', hidden = True)
def filter_invalid(prefix, input, cache_size, invalid, threads):
    '''
    Retain only valid-barcoded reads

    Use `--invalid` to separately output reads with invalid barcodes.
    Barcodes can be in haplotagging, stlfr, or tellseq formats.
    '''
    lr_type = which_linkedread(input[0])

    with CachedFQWriter(prefix, cache_size,len(input), threads = threads) as writer:
        if invalid:
            writer_inv = CachedFQWriter(f"{prefix}.invalid", cache_size, len(input))
        for i in input:
            with pysam.FastxFile(i, persist=False) as FQ:
                for record in FQ:
                    _read = FQRecord(record, lr_type, 0)
                    _read.convert(lr_type, _read.barcode)
                    if not _read.valid:
                        if invalid:
                            writer_inv.queue(_read)
                    else:
                        writer.queue(_read)
        if invalid:
            writer_inv.close()

