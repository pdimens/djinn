import rich_click as click
import subprocess
from djinn.utils.file_ops import make_dir, print_error, validate_fq_sam

@click.command(panel = "FASTQ", no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/djinn/sort/")
@click.option("--threads", "-t", type = click.IntRange(min = 6, max_open=True), default=10, show_default=True, help = "Number of threads to use")
@click.argument('samtag', metavar="SAM_tag", type = str, required = True, nargs=1)
@click.argument('prefix', metavar="output_prefix", type = str, required = True, nargs=1, callback=make_dir)
@click.argument('input', type = click.Path(dir_okay=False,readable=True,resolve_path=True, exists = True), required = True, nargs=-1, callback = validate_fq_sam)
@click.help_option('--help', hidden = True)
def sort(samtag,prefix,input,threads):
    """
    Sort FASTQ/BAM by barcode

    The barcode **must** be in a SAM tag (e.g. `BX`, `BC`) whether in
    FASTQ or SAM/BAM format.
    """
    if len(input) > 2:
        print_error('invalid input files', 'Inputs can be one single-ended or 2 paired-end FASTQ files.')
    if len(samtag) != 2:
        print_error('incorrect SAM tag', 'The SAM TAG is expected to be exactly 2 letters (e.g. BX).')

    quotient, remainder = divmod(threads - 2, 2)
    threads_sort = quotient + remainder
    threads_fastq = quotient
    _outfiles = f"-1 {prefix}.R1.fq.gz"
    if len(input) == 2:
        _outfiles += " -2 {prefix}.R2.fq.gz"

    sam_import = subprocess.Popen(
        f'samtools import -@ 1 -T *'.split() + input,
        stdout = subprocess.PIPE
    )

    sam_sort = subprocess.Popen(
        f"samtools sort -@ {threads_sort -1} -O SAM -t {samtag}".split(),
        stdout = subprocess.PIPE,
        stdin = sam_import.stdout
    )


    sam_fastq = subprocess.run(
        f'samtools fastq -@ {threads_fastq -1} -N -c 6 -T * {_outfiles}'.split(),
        stdin = sam_sort.stdout,
        stderr = subprocess.PIPE
    )

    if sam_fastq.returncode == 1:
        print_error("samtools failure", f"Samtools was unable to process your input file(s). See the error log from samtools fastq:\n\033[31m{sam_fastq.stderr.decode()}\033[0m")
