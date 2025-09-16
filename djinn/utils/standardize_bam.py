import os
import pysam
from djinn.utils.barcodes import haplotagging, tellseq, stlfr, tenx, TELLSEQ_STLFR_RX, TELLSEQ_HAPLOTAGGING_INVALID_RX
from djinn.utils.file_ops import print_error

def std_bam(prefix, sam, style):
    """
    Move barcodes to `BX`/`VX` alignment tags

    If barcodes are present in the sequence name (stlfr, tellseq), this method moves the barcode to the `BX:Z`
    tag of the alignment, maintaining the same barcode style by default. If moved to or already in a `BX:Z` tag,
    will then write a complementary `VX:i` tag to describe barcode validation `0` (invalid) or `1` (valid).
    Use `--style` to also convert the barcode to a different style (`haplotagging`, `stlfr`, `tellseq`, `10X`),
    which also writes a `conversion.bc` file to the working directory mapping the barcode conversions. Writes to `stdout`.

    | Option         | Style                                        |
    |:---------------|:---------------------------------------------|
    | `haplotagging` | AxxCxxBxxDxx                                 |
    | `stlfr`        | 1_2_3                                        |
    | `tellseq`      | 18-base nucleotide (e.g. AGCCATGTACGTATGGTA) |
    | `10X`          | 16-base nucleotide (e.g. GGCTGAACACGTGCAG)   |
    """
    convert = None
    if style:
        bc_out = open(f"{prefix}.bc", "w")
        convert = style.lower()
        if convert == "tellseq":
            BX = tellseq()
        if convert == "10x":
            BX = tenx()
        elif convert == "stlfr":
            BX = stlfr()
        else:
            BX = haplotagging()

    with (
        pysam.AlignmentFile(sam, require_index=False, check_sq=False) as samfile, 
        pysam.AlignmentFile(f"{prefix}.bam", "wb", template=samfile) as outfile,
    ):
        for record in samfile.fetch(until_eof=True):
            if record.has_tag("BX"):
                bx_sanitized = str(record.get_tag("BX"))
                if record.has_tag("VX"):
                    print_error(f"Error: BX/VX tags present\nThe BX:Z and VX:i tags are already present in {os.path.basename(sam)} and does not need to be standardized.")
            if record.has_tag("BX"):
                bx_sanitized = str(record.get_tag("BX"))
                # try to split by "_" (stlfr) and if any of the ints are zero, it's invalid
                # otherwise look for the tellseq N or haplotag 00
                if "0" in bx_sanitized.split("_") or TELLSEQ_HAPLOTAGGING_INVALID_RX.search(bx_sanitized):
                    record.set_tag("VX", 0, "i")
                    vx = 0
                else:
                    record.set_tag("VX", 1, "i")
                    vx = 1
            else:
                # matches either tellseq or stlfr   
                bx = TELLSEQ_STLFR_RX.search(record.query_name)
                if bx:
                    # the 1:0 ignores the first character, which will either be : or #
                    bx_sanitized = bx[0][1:]
                    record.query_name = record.query_name.remove_suffix(bx[0])
                    if "0" in bx_sanitized.split("_") or "N" in bx_sanitized:
                        record.set_tag("VX", 0, "i")
                        vx = 0
                    else:
                        record.set_tag("VX", 1, "i")
                        vx = 1
                    record.set_tag("BX", bx_sanitized, "Z")
                else:
                    outfile.write(record)
            if convert:
                if bx_sanitized not in BX.inventory:
                    if bool(vx):
                        try:
                            BX.inventory[bx_sanitized] = BX.next()
                        except StopIteration:
                            print_error(f"Error: too many barcodes\nThere are more barcodes in the input data than it is possible to generate {convert} barcodes from.")
                    else:
                        BX.inventory[bx_sanitized] = BX.invalid
                    bc_out.write(f"{bx_sanitized}\t{BX.inventory[bx_sanitized]}\n")
                    record.set_tag("BX", BX.inventory[bx_sanitized], "Z")
            outfile.write(record)
    if convert:
        bc_out.close()