#! /usr/bin/env python

from collections import Counter
from itertools import zip_longest
import pysam
import rich_click as click
from djinn.utils.file_ops import make_dir, validate_fq, which_linkedread
from djinn.utils.fq_tools import FQRecord, CachedFQWriter

def count_barcodes_fq(barcode_type: str, fq: list[str], invalid: bool = False) -> Counter:
    '''
    Return a `Counter` of the unique barcodes in the fastq files. If invalid is `True`,
    the `Counter` includes invalid barcodes as well.
    '''
    barcodes = Counter()
    if len(fq) == 1:
        with pysam.FastxFile(fq[0], persist=False) as FQ:
            for _read in FQ:
                _read = FQRecord(_read, barcode_type, 0)
                if invalid or _read.valid:
                    barcodes.update(_read.barcode)
    else:
        with(
            pysam.FastxFile(fq[0], persist=False) as R1,
            pysam.FastxFile(fq[1], persist=False) as R2
        ):
            for r1,r2 in zip_longest(R1,R2):
                bc_r1 = FQRecord(r1, barcode_type, 0).barcode if r1 else None
                bc_r2 = FQRecord(r2, barcode_type, 0).barcode if r2 else None
                # if paired reads have same barcode, update barcode once
                if bc_r1 == bc_r2:
                    barcodes.update(bc_r1)
                else:
                    # update barcode counts separately
                    [barcodes.update(i) for i in [bc_r1, bc_r2] if i]
    return barcodes


@click.command(no_args_is_help = True, epilog = "Documentation: https://pdimens.github.io/djinn/filter/")
@click.option("-c", "--cache-size", hidden=True, type=click.IntRange(min=1000, max_open=True), default=10000, help = "Number of cached reads for write operations")
@click.option("-s", "--singletons", is_flag=True, default=False, help = "Separately output records with valid singleton barcodes")
@click.option("-t", "--threads", type = click.IntRange(min = 1, max_open=True), default=4, show_default=True, help = "Number of compression threads to use for output files")
@click.argument('prefix', required=True, type = str, callback=make_dir)
@click.argument('input', nargs = -1, required=True, type=click.Path(exists = True,dir_okay=False,readable=True,resolve_path=True), callback = validate_fq)
@click.help_option('--help', hidden = True)
def filter_singletons(prefix, input, cache_size, singletons, threads):
    '''
    Retain reads with non-singleton barcodes

    Use `--singletons` to output reads with singleton barcodes into a separate file(s). This method also
    filters out invalid barcodes, since they are not considered linked. Expects barcode to be in haplotagging,
    stlfr, or tellseq formats. Specify `--threads` if `pigz` is available in your PATH (the value will be divided
    between the number of input files).
    '''
    from_ = which_linkedread(input[0])
    bc_counts = count_barcodes_fq(from_, input)
    linked = list(filter(lambda x: bc_counts[x] > 2, bc_counts.keys()))
    if not linked:
        print("There are no barcodes with >2 reads. Nothing to do.")
        return

    with CachedFQWriter(prefix, cache_size, len(input), threads = threads) as writer:
        if singletons:
            writer_single = CachedFQWriter(f"{prefix}.singletons", cache_size, len(input))
        for i in input:
            with pysam.FastxFile(i, persist=False) as FQ:
                for record in FQ:
                    _read = FQRecord(record, from_, 0)
                    if _read.valid:
                        if _read.barcode in linked:
                            writer.queue(_read.convert2(from_, _read.barcode))
                        elif singletons:
                            writer_single.queue(_read.convert2(from_, _read.barcode))
        if singletons:
            writer_single.close()
