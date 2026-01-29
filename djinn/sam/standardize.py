import os
import pysam
import sys
import rich_click as click
from djinn.utils.file_ops import print_error, validate_sam
from djinn.utils.barcodes import haplotagging, tellseq, stlfr, tenx, TELLSEQ_STLFR_RX, TELLSEQ_HAPLOTAGGING_INVALID_RX

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/djinn/standardize/#fastq")
@click.option('-s', '--style', type = click.Choice(["haplotagging", "stlfr", "tellseq", "10x"], case_sensitive=False), help = 'Change the barcode style')
@click.argument('input', type = click.Path(exists=True,dir_okay=False,readable=True,resolve_path=True), required = True, callback = validate_sam)
def standardize(input, style):
    """
    Move barcodes to `BX`+`VX` sequence header tags

    This conversion moves the barcode to the `BX:Z` tag in sam/bam records, maintaining the same barcode type by default (auto-detected).
    See the documentation for a deeper look into the location and format expectations for different linked-read technologies.
    Also writes a `VX:i` tag to describe barcode validation `0` (invalid) or `1` (valid).
    Use `--style` to also convert the barcode to a different style (`haplotagging`, `stlfr`, `tellseq`, `10X`).

    | Option         | Style                                        |
    |:---------------|:---------------------------------------------|
    | `haplotagging` | AxxCxxBxxDxx                                 |
    | `stlfr`        | 1_2_3                                        |
    | `tellseq`      | 18-base nucleotide (e.g. AGCCATGTACGTATGGTA) |
    | `10X`          | 16-base nucleotide (e.g. GGCTGAACACGTGCAG)   |
    """
    convert = None
    if style:
        bc_out = open(f"{input}.bc", "w")
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
        pysam.AlignmentFile(input, require_index=False, check_sq=False) as samfile,
        pysam.AlignmentFile(sys.stdout.buffer, "wb", template = samfile) as outfile
    ):
        for record in samfile.fetch(until_eof=True):
            if record.has_tag("BX"):
                if record.has_tag("VX"):
                    outfile.write(record)
                    continue
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
                bx = TELLSEQ_STLFR_RX.search(str(record.query_name))
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
                    # write record that is missing BX tag to output
                    outfile.write(record)

            if convert:
                if bx_sanitized not in BX.inventory:
                    if bool(vx):
                        try:
                            BX.inventory[bx_sanitized] = BX.next()
                        except StopIteration:
                            print_error(
                                "too many barcodes",
                                f"There are more barcodes in the input file than it is possible to generate {convert} barcodes from."
                            )
                    else:
                        BX.inventory[bx_sanitized] = BX.invalid
                    bc_out.write(f"{bx_sanitized}\t{BX.inventory[bx_sanitized]}\n")
                    record.set_tag("BX", BX.inventory[bx_sanitized], "Z")
            outfile.write(record)
    if convert:
        bc_out.close()