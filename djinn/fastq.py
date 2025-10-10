from itertools import zip_longest
import os
import pysam
import subprocess
import rich_click as click
from djinn.utils.file_ops import print_error, validate_barcodefile, which_linkedread
from djinn.utils.fq_tools import FQRecord, CachedFQWriter
from djinn.utils.barcodes import haplotagging, tellseq, stlfr, tenx

@click.command(panel = "File Conversions", no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/djinn/fastq/")
@click.option('-b','--barcodes', type = click.Path(exists=True, readable=True, dir_okay=False), help='barcodes file [10x input only]', required=False)
@click.option("-c", "--cache-size", hidden=True, type=click.IntRange(min=1000, max_open=True), default=5000, help = "Number of cached reads for write operations")
@click.argument('prefix', metavar = "PREFIX", type = str,  required=True, nargs = 1)
@click.argument('target', metavar = "TARGET", type = click.Choice(["10x", "haplotagging", "stlfr", "tellseq"], case_sensitive=False), nargs = 1)
@click.argument('fq1', metavar="R1_FASTQ", type = click.Path(exists=True,dir_okay=False,readable=True,resolve_path=True), required = True, nargs=1)
@click.argument('fq2', metavar="R2_FASTQ", type = click.Path(exists=True,dir_okay=False,readable=True,resolve_path=True), required=True, nargs= 1)
def fastq(target,fq1,fq2,prefix, barcodes, cache_size):
    """
    Convert between linked-read FASTQ formats
    
    Auto-detects the input data format and takes the positional argument `TARGET` specifying the target data format.
    10X data as input requires a `--barcodes` file (often called a barcode whitelist) so djinn can identify the inline
    barcodes. In all cases, a file will be created with the barcode conversion map. Requires 2 threads.
    
    | from/to      | barcode format                                     | example                     |
    |:-------------|:---------------------------------------------------|:----------------------------|
    | 10x          | the first N base pairs of R1, given `--barcodes`   |                             |
    | haplotagging | a `BX:Z:ACBD` SAM tag in the sequence header       | `@SEQID BX:Z:A01C93B56D11`  |
    | stlfr        | `#1_2_3` format appended to the sequence ID        | `@SEQID#1_2_3`              |
    | tellseq      | `:ATCG` format appended to the sequence ID         | `@SEQID:GGCAAATATCGAGAAGTC` |
    """
    from_ = which_linkedread(fq1)
    if from_ == "none" and not barcodes:
        print_error("no barcodes provided", "The input file was inferred to be 10X format, which requires a list of --barcodes so Djinn knows how to identify legitimate barcodes from the sequences.")

    if from_ == target:
        print_error("identical conversion target", f"The input file was inferred to be {from_}, which is identical to the conversion target {target}. The formats must be different from each other. If the input data is not {from_}, then it is formatted incorrectly for whatever technology it was generated with.")
    to_ = target.lower()

    if to_ == "tellseq":
        BX = tellseq()
    elif to_ == "stlfr":
        BX = stlfr()
    elif to_ == "10x":
        BX = tenx()
    else:
        BX = haplotagging()

    if from_ == "none":
        from_ = "10x"
        barcodelist = validate_barcodefile(barcodes)
        BX.length = 0
        for e in barcodelist:
            BX.length += len(e)
            break

    # create the output directory in case it doesn't exist
    if os.path.dirname(prefix):
        os.makedirs(os.path.dirname(prefix), exist_ok=True)

    with (
        pysam.FastxFile(fq1, persist=False) as R1,
        pysam.FastxFile(fq2, persist=False) as R2,
        CachedFQWriter(prefix, cache_size) as writer,
        open(f"{prefix}.bc", "w") as bc_out
    ):
        for r1,r2 in zip_longest(R1,R2):
            if r1:
                _r1 = FQRecord(r1, True, from_, BX.length)
                if _r1.barcode not in BX.inventory:
                    if from_ == "10x":
                        _r1.valid = _r1.barcode in barcodelist
                    if _r1.valid:
                        try:
                            BX.inventory[_r1.barcode] = BX.next()
                        except StopIteration:
                            print_error("too many barcodes", f"There are more {from_} barcodes in the input data than it is possible to generate {to_} barcodes from.")
                    else:
                        BX.inventory[_r1.barcode] = BX.invalid
                    bc_out.write(f"{_r1.barcode}\t{BX.inventory[_r1.barcode]}\n")
                converted_bc = BX.inventory[_r1.barcode]
                writer.queue(_r1.convert2(to_, converted_bc), None)
            if r2:
                if r1 and from_ == "10x":
                    _bc = _r1.barcode
                elif not r1 and from_ == "10x":
                    _bc = "N" * BX.length
                else:
                    _bc = from_
                # if input format is 10x, copy the barcode to R2
                _r2 = FQRecord(r2, False, _bc, 0)
                # check the inventory for existing barcode match
                if _r2.barcode not in BX.inventory:
                    # if it's just tellseq<->10x, keep the existing nucleotide barcode
                    if _r2.valid:
                        try:
                            BX.inventory[_r2.barcode] = BX.next()
                        except StopIteration:
                            print_error("too many barcodes", f"There are more {from_} barcodes in the input data than it is possible to generate {to_} barcodes from.")
                    else:
                        BX.inventory[_r2.barcode] = BX.invalid
                    bc_out.write(f"{_r2.barcode}\t{BX.inventory[_r2.barcode]}\n")
                converted_bc = BX.inventory[_r2.barcode]
                writer.queue(None, _r2.convert2(to_, converted_bc))
