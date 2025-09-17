import re
import rich_click as click
import subprocess

#BAM
@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/djinn/sort/")
@click.argument('samtag', metavar="SAM_tag", type = str, required = True, nargs=1)
@click.argument('prefix', metavar="output_prefix", type = str, required = True, nargs=1)
@click.argument('inputs', metavar="inputs(s)", type = click.Path(dir_okay=False,readable=True,resolve_path=True), required = True, nargs=-1)
def sort(prefix, inputs, samtag):
    """
    Sort FASTQ/BAM by barcode

    The barcode **must** be in a SAM tag (e.g. `BX`, `BC`) whether in
    FASTQ or SAM/BAM format.
    """
    ## checks and validations ##
    if len(samtag) != 2:
        raise click.BadParameter('The SAM TAG is expected to be exactly 2 letters (e.g. BX).', param_hint="SAM_TAG")

    if len(inputs) > 2:
        raise click.BadParameter('inputs must be 1 BAM file or 2 FASTQ files.', param_hint="INPUT")

    if len(inputs) == 1:
        if not inputs[0].lower().endswith(".bam"):
            raise click.BadParameter('inputs must be 1 BAM (.bam) file or 2 FASTQ (.fastq|.fq) files. The FASTQ files can be gzipped.', param_hint="INPUT")
        sam_sort = subprocess.Popen(
            f"samtools sort -@ 2 -O BAM -o {prefix}.bam -t {samtag} {inputs[0]}".split(),
            stdout = subprocess.PIPE,
            stderr = subprocess.PIPE
        )
        if sam_sort.returncode == 1:
            print(f"Error: samtools failure\nSamtools was unable to process your input file. See the error log from samtools sort:\n\033[33m{sam_sort.stderr}\033[0m")

    else:
        if inputs[0] == inputs[1]:
            raise click.BadParameter('the two input files cannot be identical', param_hint="INPUT")
        re_ext = re.compile(r"\.(fq|fastq)(?:\.gz)?$", re.IGNORECASE)
        for i in inputs:
            if not re_ext.search(i):
                raise click.BadParameter('inputs must be 1 BAM (.bam) file or 2 FASTQ (.fastq|.fq) files. The FASTQ files can be gzipped.', param_hint="INPUT")  

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
            print(f"Error: samtools failure\nSamtools was unable to process your input file. See the error log from samtools fastq:\n\033[33m{sam_fastq.stderr}\033[0m")
