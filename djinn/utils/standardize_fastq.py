
from itertools import zip_longest
import os
import pysam
import subprocess
from djinn.utils.file_ops import print_error, which_linkedread
from djinn.utils.barcodes import haplotagging, tellseq, stlfr, tenx
from djinn.utils.fq_tools import FQRecord, CachedFQWriter

def std_fastq(prefix: str, r1_fastq: str, r2_fastq: str, style: str, cache_size: int):
    """
    Move barcodes to `BX`/`VX` sequence header tags

    This conversion moves the barcode to the `BX:Z` tag in fastq records, maintaining the same barcode type by default (auto-detected).
    See `harpy convert fastq` for the location and format expectations for different linked-read technologies.
    Also writes a `VX:i` tag to describe barcode validation `0` (invalid) or `1` (valid). Use `harpy convert fastq`
    if your data is in 10X format, as this command will not work on 10X format (i.e. barcode is the first 16 bases of read 1).
    Use `--style` to also convert the barcode to a different style (`haplotagging`, `stlfr`, `tellseq`, `10X`).

    | Option         | Style                                        |
    |:---------------|:---------------------------------------------|
    | `haplotagging` | AxxCxxBxxDxx                                 |
    | `stlfr`        | 1_2_3                                        |
    | `tellseq`      | 18-base nucleotide (e.g. AGCCATGTACGTATGGTA) |
    | `10X`          | 16-base nucleotide (e.g. GGCTGAACACGTGCAG)   |
    """
    BC_TYPE = which_linkedread(r1_fastq)
    if not BC_TYPE:
        print_error("undertermined file type", f"Unable to determine the linked-read barcode type after scanning the first 100 records of {os.path.basename(r1_fastq)}. Please make sure the format is one of haplotagging, stlfr, or tellseq. 10X-style with the barcode as the first 16 nucleotides of read 1 is not supported here.")
    
    if style:
        bc_out = open(f"{prefix}.bc", "w")
        style = style.lower()
        if style == "tellseq":
            BX = tellseq()
        elif style == "10x":
            BX = tenx()
        elif style == "stlfr":
            BX = stlfr()
        else:
            BX = haplotagging()

    with (
        pysam.FastxFile(r1_fastq, persist=False) as R1,
        pysam.FastxFile(r2_fastq, persist=False) as R2,
        CachedFQWriter(prefix, cache_size) as writer,
    ):
        for r1,r2 in zip_longest(R1,R2):
            if r1:
                _r1 = FQRecord(r1, True, BC_TYPE, 0)
                if style:
                    if _r1.barcode not in BX.inventory:
                        if _r1.valid:
                            try:
                                BX.inventory[_r1.barcode] = BX.next()
                            except StopIteration:
                                print_error("too many barcodes", f"There are more {BC_TYPE} barcodes in the input data than it is possible to generate {style} barcodes from.")
                        else:
                            BX.inventory[_r1.barcode] = BX.invalid
                        # write the barcode to file
                        bc_out.write(f"{_r1.barcode}\t{BX.inventory[_r1.barcode]}\n")
                    # overwrite the record's barcode
                    _r1.barcode = BX.inventory[_r1.barcode]
                writer.queue(_r1.convert2("standard", _r1.barcode), None)
            if r2:
                _r2 = FQRecord(r2, False, BC_TYPE, 0)
                if style:
                    if _r2.barcode not in BX.inventory:
                        if _r2.valid:
                            try:
                                BX.inventory[_r2.barcode] = BX.next()
                            except StopIteration:
                                print_error("too many barcodes", f"There are more {BC_TYPE} barcodes in the input data than it is possible to generate {style} barcodes from.")
                        else:
                            BX.inventory[_r2.barcode] = BX.invalid
                        # write the barcode to file
                        bc_out.write(f"{_r2.barcode}\t{BX.inventory[_r2.barcode]}\n")
                    # overwrite the record's barcode
                    _r2.barcode = BX.inventory[_r2.barcode]
                writer.queue(None, _r2.convert2("standard", _r2.barcode))
    if style:
        bc_out.close()
