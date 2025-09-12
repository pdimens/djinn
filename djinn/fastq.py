from concurrent.futures import ThreadPoolExecutor
from itertools import zip_longest
import os
import pysam
import rich_click as click
from djinn.utils import compress_fq, FQRecord, print_error, which_linkedread
from djinn.common import haplotagging, tellseq, stlfr, tenx

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/harpy/convert")
@click.option('-o','--output', type = str, metavar= "PREFIX", help='file prefix for output fastq files', required=True)
@click.argument('target', metavar = "TARGET", type = click.Choice(["10x", "haplotagging", "stlfr", "tellseq"], case_sensitive=False), nargs = 1)
@click.argument('fq1', metavar="R1_FASTQ", type = click.Path(dir_okay=False,readable=True,resolve_path=True), required = True, nargs=1)
@click.argument('fq2', metavar="R2_FASTQ", type = click.Path(dir_okay=False,readable=True,resolve_path=True), required=True, nargs= 1)
def fastq(target,fq1,fq2,output):
    """
    Convert between linked-read FASTQ formats
    
    Autodetects the input data format and takes the positional argument `TARGET` specifying the target data format.
    10X input data requires a `--barcodes` file containing one nucleotide barcode per line to
    determine which barcodes are valid/invalid. In all cases, a file will be created with
    the barcode conversion map. Requires 2 threads.
    
    | from/to      | barcode format                                     | example                     |
    |:-------------|:---------------------------------------------------|:----------------------------|
    | 10x          | the first N base pairs of R1, given `--barcodes`   |                             |
    | haplotagging | a `BX:Z:ACBD` SAM tag in the sequence header       | `@SEQID BX:Z:A01C93B56D11`  |
    | stlfr        | `#1_2_3` format appended to the sequence ID        | `@SEQID#1_2_3`              |
    | tellseq      | `:ATCG` format appended to the sequence ID         | `@SEQID:GGCAAATATCGAGAAGTC` |
    """
    from_ = which_linkedread(fq1)
    if from_ == target:
        print_error(f"Error: identical conversion target\nThe input file was inferred to be {from_}, which is identical to the conversion target {target}. The formats must be different from each other. If the input data is not {from_}, then it is formatted incorrectly for whatever technology it was generated with.")
    to_ = target.lower()

    # for barcodes, use sample() so the barcodes don't all start with AAAAAAAAAAAAA (or 1)
    # it's not functionally important, but it does make the barcodes *look* more distinct
    if to_ == "tellseq":
        BX = tellseq()
    elif to_ == "stlfr":
        BX = stlfr()
    elif to_ == "10x":
        BX = tenx()
    else:
        BX = haplotagging()

    # create the output directory in case it doesn't exist
    if os.path.dirname(output):
        os.makedirs(os.path.dirname(output), exist_ok=True)

    with (
        pysam.FastxFile(fq1, persist=False) as R1,
        pysam.FastxFile(fq2, persist=False) as R2,
        open(f"{output}.R1.fq", "w") as R1_out,
        open(f"{output}.R2.fq", "w") as R2_out,
        open(f"{output}.bc", "w") as bc_out
    ):
        for r1,r2 in zip_longest(R1,R2):
            if r1:
                _r1 = FQRecord(r1, True, from_, BX.length)
                if _r1.barcode not in BX.inventory:
                    if _r1.valid:
                        try:
                            BX.inventory[_r1.barcode] = BX.next()
                        except StopIteration:
                            print_error(f"Error: too many barcodes\nThere are more {from_} barcodes in the input data than it is possible to generate {to_} barcodes from.")
                    else:
                        BX.inventory[_r1.barcode] = BX.invalid
                    bc_out.write(f"{_r1.barcode}\t{BX.inventory[_r1.barcode]}\n")
                converted_bc = BX.inventory[_r1.barcode]
                R1_out.write(str(_r1.convert(to_, converted_bc)))
            if r2:
                _bc = from_
                # if input format is 10x, copy the barcode to R2
                _r2 = FQRecord(r2, False, _bc, BX.length)
                # check the inventory for existing barcode match
                if _r2.barcode not in BX.inventory:
                    # if it's just tellseq<->10x, keep the existing nucleotide barcode
                    if _r2.valid:
                        try:
                            BX.inventory[_r2.barcode] = BX.next()
                        except StopIteration:
                            print_error(f"Error: too many barcodes\nThere are more {from_} barcodes in the input data than it is possible to generate {to_} barcodes from.")
                    else:
                        BX.inventory[_r2.barcode] = BX.invalid
                    bc_out.write(f"{_r2.barcode}\t{BX.inventory[_r2.barcode]}\n")
                converted_bc = BX.inventory[_r2.barcode]
                R2_out.write(str(_r2.convert(to_, converted_bc)))
    
    # bgzip compress the output, one file per thread
    with ThreadPoolExecutor(max_workers=2) as executor:
        executor.submit(compress_fq, f"{output}.R1.fq")
        executor.submit(compress_fq, f"{output}.R2.fq")
