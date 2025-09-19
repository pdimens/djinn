import rich_click as click
import subprocess
from djinn.utils.file_ops import print_error, validate_fq_sam

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/djinn/sort/")
@click.argument('samtag', metavar="SAM_tag", type = str, required = True, nargs=1)
@click.argument('prefix', metavar="output_prefix", type = str, required = True, nargs=1)
@click.argument('inputs', type = click.Path(dir_okay=False,readable=True,resolve_path=True, exists = True), required = True, nargs=-1, callback = validate_fq_sam)
def sort(samtag,prefix,inputs):
    """
    Sort FASTQ/BAM by barcode

    The barcode **must** be in a SAM tag (e.g. `BX`, `BC`) whether in
    FASTQ or SAM/BAM format.
    """
    ## checks and validations ##
    if len(samtag) != 2:
        print_error('incorrect SAM tag', 'The SAM TAG is expected to be exactly 2 letters (e.g. BX).')

    if len(inputs) == 1:
        sam_sort = subprocess.Popen(
            f"samtools sort -@ 2 -O BAM -o {prefix}.bam -t {samtag} {inputs[0]}".split(),
            stdout = subprocess.PIPE,
            stderr = subprocess.PIPE
        )
        if sam_sort.returncode == 1:
            print_error("samtools failure", f"Samtools was unable to process your input file. See the error log from samtools sort:\n\033[31m{sam_sort.stderr}\033[0m")

    else:
        sam_import = subprocess.Popen(
            f'samtools import -@ 1 -T * {inputs[0]} {inputs[1]}'.split(),
            stdout = subprocess.PIPE
        )

        sam_sort = subprocess.Popen(
            f"samtools sort -@ 1 -O SAM -t {samtag}".split(),
            stdout = subprocess.PIPE,
            stdin = sam_import.stdout
        )

        sam_fastq = subprocess.run(
            f'samtools fastq -@ 1 -N -c 6 -T * -1 {prefix}.R1.fq.gz -2 {prefix}.R2.fq.gz'.split(),
            stdin = sam_sort.stdout,
            stderr = subprocess.PIPE
        )

        if sam_fastq.returncode == 1:
            print_error("samtools failure", f"Samtools was unable to process your input file. See the error log from samtools fastq:\n\033[31m{sam_fastq.stderr}\033[0m")
