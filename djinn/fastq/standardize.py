import os
import pysam
from djinn.utils.file_ops import print_error, which_linkedread, make_dir, validate_fq
from djinn.utils.barcodes import haplotagging, tellseq, stlfr, tenx
from djinn.utils.fq_tools import FQRecord, CachedFQWriter
import rich_click as click

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/djinn/standardize/#fastq")
@click.option("-c", "--cache-size", hidden=True, type=click.IntRange(min=1000, max_open=True), default=5000, help = "Number of cached reads for write operations")
@click.option('-s', '--style', type = click.Choice(["haplotagging", "stlfr", "tellseq", "10x"], case_sensitive=False), help = 'Change the barcode style')
@click.argument('prefix', metavar="output_prefix", type = str, required = True, nargs=1, callback=make_dir)
@click.argument('input', nargs = -1, type = click.Path(exists=True,dir_okay=False,readable=True,resolve_path=True), required = True, callback = validate_fq)
def standardize(prefix, input, style, cache_size):
    """
    Move barcodes to `BX`/`VX` sequence header tags

    This conversion moves the barcode to the `BX:Z` tag in fastq records, maintaining the same barcode type by default (auto-detected).
    See the documentation for a deeper look into the location and format expectations for different linked-read technologies.
    Also writes a `VX:i` tag to describe barcode validation `0` (invalid) or `1` (valid). Use `djinn fastq`
    if your fastq data is in 10X format, as this command will not work on 10X format (i.e. barcode is the first 16 bases of read 1).
    Use `--style` to also convert the barcode to a different style (`haplotagging`, `stlfr`, `tellseq`, `10X`).

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

    with CachedFQWriter(prefix, cache_size, len(input)) as writer:
        for i,j in enumerate(input):
            _fw = i == 0
            with pysam.FastxFile(j, persist=False) as FQ:
                for _read in FQ:
                    _record = FQRecord(_read, _fw, BC_TYPE, 0)
                    if style:
                        if _record.barcode not in BX.inventory:
                            if _record.valid:
                                try:
                                    BX.inventory[_record.barcode] = BX.next()
                                except StopIteration:
                                    print_error("too many barcodes", f"There are more {BC_TYPE} barcodes in the input data than it is possible to generate {style} barcodes from.")
                            else:
                                BX.inventory[_record.barcode] = BX.invalid
                            # write the barcode to file
                            bc_out.write(f"{_record.barcode}\t{BX.inventory[_record.barcode]}\n")
                        # overwrite the record's barcode
                        _record.barcode = BX.inventory[_record.barcode]
                    if _fw:
                        writer.queue(_record.convert2("standard", _record.barcode), None)
                    else:
                        writer.queue(None, _record.convert2("standard", _record.barcode))

    if style:
        bc_out.close()

