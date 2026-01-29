import os
import pysam
from djinn.utils.file_ops import print_error, which_linkedread, make_dir, validate_fq
from djinn.utils.barcodes import haplotagging, tellseq, stlfr, tenx
from djinn.utils.fq_tools import FQRecord, CachedFQWriter
import rich_click as click

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/djinn/standardize/#fastq")
@click.option("-c", "--cache-size", hidden=True, type=click.IntRange(min=1000, max_open=True), default=10000, help = "Number of cached reads for write operations")
@click.option('-s', '--style', type = click.Choice(["haplotagging", "stlfr", "tellseq", "10x"], case_sensitive=False), help = 'Change the barcode style')
@click.option("-t", "--threads", type = click.IntRange(min = 1, max_open=True), default=4, show_default=True, help = "Number of compression threads to use per output file")
@click.argument('prefix', metavar="output_prefix", type = str, required = True, nargs=1, callback=make_dir)
@click.argument('input', nargs = -1, type = click.Path(exists=True,dir_okay=False,readable=True,resolve_path=True), required = True, callback = validate_fq)
@click.help_option('--help', hidden = True)
def standardize(prefix, input, style, cache_size, threads):
    """
    Move barcodes to `BX`/`VX` sequence header tags

    This conversion moves the barcode to the `BX:Z` tag in fastq records, maintaining the same barcode type by default (auto-detected).
    See the documentation for a deeper look into the location and format expectations for different linked-read technologies.
    Also writes a `VX:i` tag to describe barcode validation `0` (invalid) or `1` (valid). Use `djinn fastq`
    if your fastq data is in 10X format, as this command will not work on 10X format (i.e. barcode is the first 16 bases of read 1).
    Use `--style` to also convert the barcode to a different style (`haplotagging`, `stlfr`, `tellseq`, `10X`).
    Specify `--threads` if `pigz` is available in your PATH (the value will be divided
    between the number of input files).

    | Option         | Style                                        |
    |:---------------|:---------------------------------------------|
    | `haplotagging` | AxxCxxBxxDxx                                 |
    | `stlfr`        | 1_2_3                                        |
    | `tellseq`      | 18-base nucleotide (e.g. AGCCATGTACGTATGGTA) |
    | `10X`          | 16-base nucleotide (e.g. GGCTGAACACGTGCAG)   |
    """
    BC_TYPE = which_linkedread(input[0])
    if not BC_TYPE:
        print_error(
            "undertermined file type",
            f"Unable to determine the linked-read barcode type after scanning the first 100 records of {os.path.basename(input[0])}. "
            "Please make sure the format is one of haplotagging, stlfr, or tellseq. 10X-style with the barcode as the first 16 nucleotides"
            " of read 1 is not supported here."
        )
    
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

    with CachedFQWriter(prefix, cache_size, len(input), threads=threads) as writer:
        for i in input:
            with pysam.FastxFile(i, persist=False) as FQ:
                for _read in FQ:
                    _record = FQRecord(_read, BC_TYPE, 0)
                    if style:
                        if _record.barcode not in BX.inventory:
                            if _record.valid:
                                try:
                                    BX.inventory[_record.barcode] = BX.next()
                                except StopIteration:
                                    print_error(
                                        "insufficient {style} barcodes",
                                        f"There are more {BC_TYPE} barcodes in the input data than it is possible to generate {style} barcodes from."
                                    )
                            else:
                                BX.inventory[_record.barcode] = BX.invalid
                            # write the barcode to file
                            bc_out.write(f"{_record.barcode}\t{BX.inventory[_record.barcode]}\n")
                        # overwrite the record's barcode
                        _record.convert("standard", BX.inventory[_record.barcode])
                    writer.queue(_record)
    if style:
        bc_out.close()

