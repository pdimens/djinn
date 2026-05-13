from itertools import zip_longest
import rich_click as click
import subprocess
from pysam import FastxFile
from djinn.utils.file_ops import print_error, which_linkedread, make_dir, validate_fq
from djinn.utils.fq_tools import FQRecord, CachedFQWriter
from djinn.fastq.singletons import count_barcodes_fq
from djinn.__main__ import config

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/djinn/arachne/")
@click.option("--threads", "-t", type = click.IntRange(min = 5, max_open=True), default=4, show_default=True, help = "Number of threads to use (min: 4)")
@click.argument('prefix', type = str, required = True, nargs=1, callback=make_dir)
@click.argument('input', nargs=2, type = click.Path(dir_okay=False,readable=True,resolve_path=True, exists = True), required = True, callback = validate_fq)
@click.rich_config(config)
@click.help_option('--help', hidden = True)
def arachne(input, prefix, threads):
    """
    Prepare FASTQ pair for arachne input

    This command shortcuts performing the necessary file processing to make a pair of FASTQ files compliant with
    the Arachne linked-read sequence aligner, which includes standardizing linked-read barcode format, sorting by barcode,
    and filtering out invalid and singleton barcodes. You **must** provide a prefix for the output FASTQ files.
    """
    BC_TYPE = which_linkedread(input[0])
    bc_counts = count_barcodes_fq(BC_TYPE, input)
    linked = list(filter(lambda x: bc_counts[x] > 1, bc_counts.keys()))
    samheader = (
        b"@HD\tVN:1.6\tSO:unsorted\tGO:query\n"
        b"@CO\tReverse with: samtools fastq -1 R1.fastq -2 R2.fastq\n"
        b"@PG\tID:samtools\tPN:samtools\tVN:1.23.1\tCL:samtools import\n"
    )
    sam_sort = subprocess.Popen(
        f"samtools sort --no-PG -O SAM -t BX -".split(),
        stdout=subprocess.PIPE,
        stdin=subprocess.PIPE,
        stderr=subprocess.PIPE
    )

    sam_fastq = subprocess.Popen(
        f'samtools fastq -@ {threads - 3} -N -c 6 -T * -1 {prefix}.R1.fq.gz -2 {prefix}.R2.fq.gz'.split(),
        stdin=sam_sort.stdout,
        stderr=subprocess.PIPE
    )

    sam_sort.stdin.write(samheader)

    with (
        FastxFile(input[0], persist=False) as FQF,
        FastxFile(input[1], persist=False) as FQR,
        CachedFQWriter(f"{prefix}.filtered_out", 10000, 2) as invalids
    ):
        for _read1, _read2 in zip_longest(FQF, FQR):
            if _read1 is not None:
                _record1 = FQRecord(_read1, BC_TYPE, 0)
                _record1.convert("standard", _record1.barcode)
                _record1.id = _record1.id[:-2] 
                if _record1.valid and _record1.barcode in linked:
                    sam_sort.stdin.write(_record1.asBam())
                else:
                    invalids.queue(_record1)

            if _read2 is not None:
                _record2 = FQRecord(_read2, BC_TYPE, 0)
                _record2.convert("standard", _record2.barcode)
                _record2.id = _record2.id[:-2] 
                if _record2.valid and _record2.barcode in linked:
                    sam_sort.stdin.write(_record2.asBam())
                else:
                    invalids.queue(_record2)

    # Close stdin to signal EOF to samtools sort
    sam_sort.stdin.close()

    # Wait for the pipeline to finish in order
    sam_sort.wait()
    sam_fastq.wait()

    errs = []
    for proc, name in [
        (sam_sort,   "samtools sort"),
        (sam_fastq,  "samtools fastq"),
    ]:
        if proc.returncode != 0:
            errs.append(f"[red]{name}[/]\n" + proc.stderr.read().decode() + "\n")

    if errs:
        print_error(
            "samtools failure",
            f"Samtools was unable to process your input. See the error logs:\n" + "\n".join(errs)
    )