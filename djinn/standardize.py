import re
import rich_click as click
from djinn.utils.standardize_bam import std_bam
from djinn.utils.standardize_fastq import std_fastq
from djinn.utils.file_ops import validate_fq_sam
@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/djinn/standardize/#fastq")
@click.option('-s', '--style', type = click.Choice(["haplotagging", "stlfr", "tellseq", "10x"], case_sensitive=False), help = 'Change the barcode style')
@click.argument('prefix', metavar="output_prefix", type = str, required = True, nargs=1)
@click.argument('inputs', type = click.Path(exists=True,dir_okay=False,readable=True,resolve_path=True), required = True, callback = validate_fq_sam, nargs=-1)
def standardize(prefix, inputs, style):
    """
    Move barcodes to `BX`/`VX` sequence header tags

    This conversion moves the barcode to the `BX:Z` tag in fastq/bam records, maintaining the same barcode type by default (auto-detected).
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
    if len(inputs) == 1:
        std_bam(prefix, inputs[0], style)
    else:
        std_fastq(prefix, inputs[0], inputs[0], style)

