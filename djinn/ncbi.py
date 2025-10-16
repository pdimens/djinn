import os
import rich_click as click
import subprocess
from djinn.utils.file_ops import make_dir, print_error, validate_fq_sam

@click.command(panel = "File Conversions", no_args_is_help = True, epilog = "Documentation: https://pdimens.github.io/djinn/ncbi/")
@click.option("--threads", "-t", type = click.IntRange(min = 2, max_open=True), default=10, show_default=True, help = "Number of threads to use")
@click.argument('prefix', required=True, type = str, callback=make_dir)
@click.argument('inputs', required=True, type=click.Path(exists = True,dir_okay=False,readable=True,resolve_path=True), callback = validate_fq_sam, nargs=-1)
def ncbi(prefix, inputs, threads):
    """
    FASTQ ⇆ BAM conversion for/from NCBI

    The input FASTQ files must have their barcode in an auxilary tag (e.g. `BX:Z:`), otherwise
    you run the risk of NCBI removing any barcode information stored in the sequence header.
    If the barcodes are in default tellseq/stlfr format, use `djinn standardize` to move the barcode into
    the BX tag. If given a single input SAM/BAM file, will instead convert it back to
    two FASTQ files.

    """
    ## checks and validations ##
    if len(inputs) == 1:
        fq = subprocess.run(
            f'samtools fastq -@ {threads-1} -N -c 6 -T * -1 {prefix}.R1.fq.gz -2 {prefix}.R2.fq.gz {inputs[0]}'.split(),
            stderr = subprocess.PIPE
        )
        if fq.returncode == 1:
            print_error("samtools failure", f"Samtools was unable to process your input file. See the error log from samtools fastq:\n\033[31m{fq.stderr.decode()}\033[0m")

    else:
        fq = subprocess.run(
            f'samtools import -@ {threads-1} -O BAM -o {prefix}.bam -T * -1 {inputs[0]} -2 {inputs[1]}'.split(),
            stderr = subprocess.PIPE      
        )
        if fq.returncode == 1:
            print_error("samtools failure", f"Samtools was unable to process your input files. See the error log from samtools import:\n\033[31m{fq.stderr.decode()}\033[0m")


#BC_QUAL = "I"*20
#SPACER_NUC = "N"*10
#SPACER_QUAL = "!"*10

#@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "https://pdimens.github.io/djinn/ncbi/")
#@click.option('-m', '--barcode-map',  is_flag = True, default = False, help = 'Write a map of the barcode-to-nucleotide conversion')
#@click.option('-p', '--prefix', required=True, type = str, help = "Output file name prefix")
#@click.option('-i', '--preserve-invalid',  is_flag = True, default = False, help = 'Retain the uniqueness of invalid barcodes')
#@click.option('-s', '--scan', show_default = True, default = 100, type = click.IntRange(min=1), help = 'Number of reads to scan to identify barcode location and format')
#@click.argument('r1_fq', required=True, type=click.Path(dir_okay=False,readable=True,resolve_path=True), nargs=1)
#@click.argument('r2_fq', required=True, type=click.Path(dir_okay=False,readable=True,resolve_path=True), nargs=1)
#def ncbi_fastq(prefix, r1_fq, r2_fq, scan, preserve_invalid, barcode_map):
#    """
#    Convert FASTQ files for NCBI submission
#
#    Converts the linked-read barcodes into nucleotide format (if necessary) and adds it to the beginning
#    of the sequence, retaining the linked-read barcodes after NCBI/SRA reformats the FASTQ after submission.
#    The barcode will be stored as the first 18 bases in both R1 and R2 reads, followed by a 5bp spacer of "NNNNN",
#    then the actual sequence. Requires a `--prefix` to name the output files. Use `--barcode-map`/`-m` to write a file with
#    the barcode conversion map if you intend to keep the same barcodes after downloading sequences from NCBI and
#    demultiplexing with `harpy demultiplex ncbi`. Invalid barcodes will be generalized to 18bp of `N`, but you can
#    use `--preserve-invalid`/`-i` to keep invalid barcodes unique (likely not useful for most applications).
#    """
#    from_ = which_linkedread(r1_fq, scan)
#    if from_ == "none":
#        print_error("Error: format not recognized\nThe provided FASTQ files were not identified as haplotagging, stlfr, or tellseq data. Please check that your input files conform to the standards of those technologies.")
#    
#    parser = generic_parser(from_)
#    BX = _ncbi(preserve_invalid)
#
#    for i,fq in enumerate([r1_fq, r2_fq],1):
#        with pysam.FastxFile(fq) as in_fq, open(f"{prefix}.R{i}.fq.gz", "wb") as out_fq:
#            gzip = subprocess.Popen(["gzip"], stdin = subprocess.PIPE, stdout = out_fq)
#            try:
#                for record in in_fq:
#                    parser.process_barcode(record)
#                    # if the barcode is missing or flagged as invalid
#                    if not parser.valid:
#                        inline_bc = BX.next_invalid()
#                        if preserve_invalid and parser.barcode not in BX.inventory:
#                            BX.inventory[parser.barcode] = inline_bc
#                    else:
#                        if parser.barcode not in BX.inventory:
#                            inline_bc = BX.next()
#                            BX.inventory[parser.barcode] = inline_bc
#                        else:
#                            inline_bc = BX.inventory[parser.barcode]
#
#                    record.sequence = inline_bc + SPACER_NUC + record.sequence
#                    record.quality  = BC_QUAL + SPACER_QUAL + record.quality
#                    gzip.stdin.write(str(record).encode("utf-8") + b"\n")
#            finally:
#                gzip.stdin.close()
#                retcode = gzip.wait()
#                if retcode != 0:
#                    click.echo(f"Error: gzip exited with status {retcode}", err=True)
#                    sys.exit(retcode)
#    if barcode_map:
#        with open(f"{prefix}.barcode.map", "w") as bc_out:
#            for bx,nuc in BX.inventory.items():
#                bc_out.write(f"{nuc}\t{bx}\n")
