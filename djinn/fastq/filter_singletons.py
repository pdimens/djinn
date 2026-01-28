#! /usr/bin/env python

from collections import Counter
import os
from itertools import zip_longest
import pysam
import rich_click as click
from djinn.utils.file_ops import make_dir, print_error, validate_fq, which_linkedread
from djinn.utils.fq_tools import FQRecord, CachedFQWriter

def count_barcodes_fq(barcode_type: str, fq: list[str], invalid: bool = False) -> Counter:
    '''
    Return a Counter of the unique barcodes in the fastq files. If invalid is True,
    the Counter includes invalid barcodes as well.
    '''
    barcodes = Counter()
    for i,j in enumerate(fq):
        _fw = i == 0
        with pysam.FastxFile(j, persist=False) as FQ:
            for _read in FQ:
                _read = FQRecord(_read, _fw, barcode_type, 0)
                if invalid:
                    if not _read.valid:
                        barcodes[_read.barcode] += 1
                else:
                    barcodes[_read.barcode] += 1
    return barcodes

@click.command(no_args_is_help = True, epilog = "Documentation: https://pdimens.github.io/djinn/filter/")
@click.option("-c", "--cache-size", hidden=True, type=click.IntRange(min=1000, max_open=True), default=5000, help = "Number of cached reads for write operations")
@click.option("-s", "--singletons", is_flag=True, default=False, help = "Separately output records with valid singleton barcodes")
#@click.option("-t", "--threads", type = click.IntRange(min = 2, max_open=True), default=4, show_default=True, help = "Number of threads to use")
@click.argument('prefix', required=True, type = str, callback=make_dir)
@click.argument('input', required=True, type=click.Path(exists = True,dir_okay=False,readable=True,resolve_path=True), callback = validate_fq, nargs=-1)
@click.help_option('--help', hidden = True)
def filter_singletons(prefix, input, cache_size, singletons):
    '''
    Retain reads with non-singleton barcodes

    Use `--singletons` to output reads with singleton barcodes into a separate file(s). This method also
    filters out invalid barcodes, since they are not considered linked. Expects barcode to be in haplotagging,
    stlfr, or tellseq formats.
    '''
    from_ = which_linkedread(input[0])
    bc_counts = count_barcodes_fq(from_, input)
    linked = list(filter(lambda x: bc_counts[x] > 2, bc_counts.keys()))
    if not linked:
        print("There are no barcodes with >2 reads. Nothing to do.")
        return

    if len(input) == 1:
        with (
            pysam.FastxFile(input[0], persist=False) as R1,
            CachedFQWriter(prefix, cache_size) as writer,
            CachedFQWriter(f"{prefix}.singletons", cache_size) as writer_single,
        ):
            for r1 in R1:
                _r1 = FQRecord(r1, True, from_, 0)
                if _r1.valid:
                    if _r1.barcode in linked:
                        writer.queue(_r1.convert2(from_, _r1.barcode), None)
                    elif singletons:
                        writer_single.queue(_r1.convert2(from_, _r1.barcode), None)

        if not singletons:
            os.unlink(f"{prefix}.singletons.R1.fq.gz")
        return

    with (
        pysam.FastxFile(input[0], persist=False) as R1,
        pysam.FastxFile(input[1], persist=False) as R2,
        CachedFQWriter(prefix, cache_size) as writer,
        CachedFQWriter(f"{prefix}.singletons", cache_size) as writer_single,
    ):
        for r1,r2 in zip_longest(R1,R2):
            if r1:
                _r1 = FQRecord(r1, True, from_, 0)
                if _r1.valid:
                    if _r1.barcode in linked:
                        writer.queue(_r1.convert2(from_, _r1.barcode), None)
                    elif singletons:
                        writer_single.queue(_r1.convert2(from_, _r1.barcode), None)

            if r2:
                _r2 = FQRecord(r2, False, from_, 0)
                if _r2.valid:
                    if _r2.barcode in linked:
                        writer.queue(None, _r2.convert2(from_, _r2.barcode))
                    elif singletons:
                        writer_single.queue(None, _r2.convert2(from_, _r2.barcode))
    if not singletons:
        os.unlink(f"{prefix}.singletons.R1.fq.gz")
        os.unlink(f"{prefix}.singletons.R2.fq.gz")