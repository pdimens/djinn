#! /usr/bin/env python

from collections import Counter
import os
import re
import subprocess
import sys
import pysam
import rich_click as click
from djinn.utils.barcodes import ANY_INVALID, TELLSEQ_STLFR_RX
from djinn.utils.file_ops import print_error, validate_fq_sam, which_linkedread, which_linkedread_sam
from djinn.utils.fq_tools import FQRecord, CachedWriter

def count_barcodes_sam(bamfile: str, invalid: bool = False) -> Counter:
    '''
    Return a Counter of the unique barcodes in bamfile. If invalid is True,
    the Counter includes invalid barcodes as well.
    '''
    barcodes = Counter()
    with pysam.AlignmentFile(bamfile, check_sq=False) as infile:
        for record in infile:
            _bc = None
            if record.has_tag("BX"):
                _bc = record.get_tag("BX")
            else:
                _check =  TELLSEQ_STLFR_RX.search(record.name)
                if _check:
                    _bc = _check.group(0)[1:]

            if _bc:
                if invalid:
                    if ANY_INVALID.search(_bc):
                        barcodes[_bc] += 1
                else:
                    barcodes[_bc] += 1
    return barcodes

def count_barcodes_fq(barcode_type: str, fq1: str, fq2: str, invalid: bool = False) -> Counter:
    '''
    Return a Counter of the unique barcodes in the fastq files. If invalid is True,
    the Counter includes invalid barcodes as well.
    '''
    barcodes = Counter()
    with (
        pysam.FastxFile(fq1, persist=False) as R1,
        pysam.FastxFile(fq2, persist=False) as R2,
    ):
        for _r1 in R1:
            _r1 = FQRecord(_r1, True, barcode_type, 0)
            if invalid:
                if not _r1.valid:
                    barcodes[_r1.barcode] += 1
            else:
                barcodes[_r1.barcode] += 1
        for _r2 in R2:
            _r2 = FQRecord(_r2, False, barcode_type, 0)
            if invalid:
                if not _r2.valid:
                    barcodes[_r2.barcode] += 1
            barcodes[_r2.barcode] += 1
    
    return barcodes

@click.command(panel = "Other Tools", no_args_is_help = True, epilog = "Documentation: https://pdimens.github.io/djinn/filter/")
@click.option("-c", "--cache-size", hidden=True, type=click.IntRange(min=1000, max_open=True), default=5000, help = "Number of cached reads for write operations")
@click.option("-s", "--singletons", is_flag=True, default=False, help = "Separately output records with valid singleton barcodes")
@click.option("-t", "--threads", type = click.IntRange(min = 2, max_open=True), default=4, show_default=True, help = "Number of threads to use")
@click.argument('prefix', required=True, type = str)
@click.argument('inputs', required=True, type=click.Path(exists = True,dir_okay=False,readable=True,resolve_path=True), callback = validate_fq_sam, nargs=-1)
def filter_singletons(prefix, inputs, cache_size, singletons, threads):
    '''
    Isolate linked reads from singletons

    Use `--singleton` to separately output reads with singleton barcodes.
    '''
    if len(inputs) == 1:
        bc_counts = count_barcodes_sam(inputs[0])
        linked = filter(lambda x: bc_counts[x] > 2, bc_counts.keys())
        if singletons:
            _singles = list(filter(lambda x: bc_counts[x] <= 2, bc_counts.keys()))
        
        with pysam.AlignmentFile(inputs[0]) as xam, pysam.AlignmentFile(f"{prefix}.bam", "wb", template = xam) as xam_out:
            for record in xam:
                #SEARCH FOR BARCODE

                # IF BARCODE IN GOOD LIST, OUTPUT RECORD

                # IF SINGLETONS AND BC IN SINGLETON LIST, OUTPUT TO OTHER FILE


    else:
        from_ = which_linkedread(inputs[0])
        bc_counts = count_barcodes_fq(from_, inputs[0], inputs[1])
        linked = list(filter(lambda x: bc_counts[x] > 2, bc_counts.keys()))
        if singletons:
            _singles = list(filter(lambda x: bc_counts[x] <= 2, bc_counts.keys()))

        with (
            pysam.FastxFile(inputs[0], persist=False) as R1,
            pysam.FastxFile(inputs[1], persist=False) as R2,
            open(f"{prefix}.R1.fq.gz", "wb") as R1_out,
            open(f"{prefix}.R2.fq.gz", "wb") as R2_out,
            subprocess.Popen("gzip", stdout= R1_out, stdin=subprocess.PIPE) as gz_r1,
            subprocess.Popen("gzip", stdout= R2_out, stdin=subprocess.PIPE) as gz_r2,
            open(f"{prefix}.singleton.R1.fq.gz", "wb") as R1_single,
            open(f"{prefix}.singleton.R2.fq.gz", "wb") as R2_single,
            subprocess.Popen("gzip", stdout= R1_single, stdin=subprocess.PIPE) as gz_r1_single,
            subprocess.Popen("gzip", stdout= R2_single, stdin=subprocess.PIPE) as gz_r2_single,
        ):
            writer = CachedWriter(gz_r1, gz_r2, cache_size)
            writer_inv = CachedWriter(gz_r1_single, gz_r2_single, cache_size)

