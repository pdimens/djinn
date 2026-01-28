import rich_click as click
import subprocess
from djinn.utils.file_ops import make_dir, print_error, validate_sam

@click.command(no_args_is_help = True, epilog = "Documentation: https://pdimens.github.io/djinn/ncbi/")
@click.option("--threads", "-t", type = click.IntRange(min = 2, max_open=True), default=10, show_default=True, help = "Number of threads to use")
@click.argument('prefix', required=True, type = str, callback=make_dir)
@click.argument('inputs',nargs=1, required=True, type=click.Path(exists = True,dir_okay=False,readable=True,resolve_path=True), callback = validate_sam)
@click.help_option('--help', hidden = True)
def ncbi(prefix, inputs, threads):
    """
    BAM → FASTQ conversion from NCBI

    Converts an unmapped SAM/BAM file to FASTQ sequences, losslessly, preserving barcode tags.
    """
    fq = subprocess.run(
        f'samtools fastq -@ {threads-1} -N -c 6 -T * -1 {prefix}.R1.fq.gz -2 {prefix}.R2.fq.gz {inputs[0]}'.split(),
        stderr = subprocess.PIPE
    )
    if fq.returncode == 1:
        print_error("samtools failure", f"Samtools was unable to process your input file. See the error log from samtools fastq:\n\033[31m{fq.stderr.decode()}\033[0m")
