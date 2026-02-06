import rich_click as click
import subprocess
import sys
from djinn.utils.file_ops import print_error, validate_sam

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/djinn/sort/")
@click.option("--threads", "-t", type = click.IntRange(min = 6, max_open=True), default=10, show_default=True, help = "Number of threads to use")
@click.option('-S', '--sam', is_flag = True, default = False, help = 'Output as SAM instead of BAM')
@click.argument('samtag', metavar="SAM_tag", type = str, required = True, nargs=1)
@click.argument('input', type = click.Path(dir_okay=False,readable=True,resolve_path=True, exists = True), required = True, callback = validate_sam)
@click.help_option('--help', hidden = True)
def sort(samtag,input,threads, sam):
    """
    Sort by barcode

    The barcode **must** be in a SAM tag (e.g. `BX`, `BC`).
    """
    if len(samtag) != 2:
        print_error('incorrect SAM tag', 'The SAM TAG is expected to be exactly 2 letters (e.g. BX).')

    outfmt = "SAM" if sam else "BAM"
    sam_sort = subprocess.Popen(
        f"samtools sort -@ {threads-1} -O {outfmt} -t {samtag} {input}".split(),
        stdout = sys.stdout,
        stderr = subprocess.PIPE
    )
    if sam_sort.returncode == 1:
        print_error("samtools failure", f"Samtools was unable to process your input file. See the error log from samtools sort:\n\033[31m{sam_sort.stderr.decode()}\033[0m")
