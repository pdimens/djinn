from concurrent.futures import ThreadPoolExecutor
import gzip
import os
import pysam
import re
import sys
from djinn.utils.barcodes import HAPLOTAGGING_RX, STLFR_RX, TELLSEQ_RX
from djinn.utils.barcodes import HAPLOTAGGING_SIMPLE, STLFR_SIMPLE, TELLSEQ_SIMPLE

def _compress_fq(fq: str):
    """use pysam bgzip to compress fastq and delete the original"""
    try:
        pysam.tabix_compress(fq, f"{fq}.gz", force=True)
        os.remove(fq)
    except Exception as e:
        print(f"Failed to compress {fq}: {str(e)}")
        sys.exit(1)

def compress_fq(fq1: str, fq2: str):
    '''bgzip compress the output, one thread per file'''
    with ThreadPoolExecutor(max_workers=2) as executor:
        executor.submit(_compress_fq, fq1)
        executor.submit(_compress_fq, fq2)

class generic_parser():
    def __init__(self, bx_type: str):
        self.bx_type = bx_type
        self.regex = re.compile(r"(?:\:[ATCGN]+$|#\d+_\d+_\d+$)")

    def process_barcode(self, record):
        if self.bx_type == "haplotagging":
            bx = [i for i in record.comment.split() if i.startswith("BX:Z:")]
            if bx:
                self.barcode = bx[0].removeprefix("BX:Z:")
            else:
                self.barcode = None
        else:
            bx = self.regex.search(record.name)
            if bx:
                self.barcode = bx[0][1:]
            else:
                self.barcode = None
        if self.barcode:
            self.valid = "0" in self.barcode.split("_") or bool(re.search(r"(?:N|[ABCD]00)", self.barcode))
        else:
            self.valid = 0

def print_error(title:str, text: str):
    """Print the error text to stderr and exit with code 1"""
    print(f"\033[33mError: {title}\033[0m", file = sys.stderr)
    try:
        print(text.encode('ascii', 'ignore').decode('unicode_escape'), file = sys.stderr)
    except:
        print(text, file = sys.stderr)
    sys.exit(1)

def safe_read(file_path: str):
    """returns the proper file opener for reading if a file_path is gzipped"""
    try:
        with gzip.open(file_path, 'rt') as f:
            f.read(10)
        return gzip.open(file_path, 'rt')
    except gzip.BadGzipFile:
        return open(file_path, 'r')

def validate_barcodefile(infile: str, limit: int = 60) -> set[str]:
    """
    Does validations to make sure it's one length, within a length limit, one per line, and nucleotides. Returns
    a list of the barcodes.
    """
    barcodes = set()
    lengths = set()
    nucleotides = {'A','C','G','T'}
    
    def validate(line_num, bc_line):
        barcode = bc_line.rstrip()
        if len(barcode.split()) > 1:
            print_error("incorrect barcode format", f"There must be one barcode per line, but multiple entries were detected on line {line_num} in {infile}")
        if not set(barcode).issubset(nucleotides) or barcode != barcode.upper():
            print_error("incorrect barcode format", f"Invalid barcode format on line {line_num }: {barcode}.\nBarcodes in {infile} must be captial letters and only contain standard nucleotide characters ATCG.")
        return len(barcode)
    with safe_read(infile) as bc_file:
        for line,bc in enumerate(bc_file, 1):
            length = validate(line, bc)
            if length > limit:
                print_error("barcodes too long", f"Barcodes in {infile} are {length}bp and cannot exceed a length of {limit}. Please use shorter barcodes.")
            lengths.add(length)
            if len(lengths) > 1:
                str_len = ", ".join(str(_length) for _length in lengths)
                print_error("inconsistent length", f"Barcodes in {infile} must all be a single length, but multiple lengths were detected: {str_len}")
            barcodes.add(bc)
    if not lengths:
        print_error("no barcodes detected", f"No barcodes were found in {infile}. Please check the input file.")
    return barcodes

def which_linkedread(fastq: str, n: int = 100) -> str:
    """
    Scans the first 100 records of a FASTQ file and tries to determine the barcode technology
    Returns one of: "haplotagging", "stlfr", "tellseq", or "none"
    """
    with pysam.FastxFile(fastq, persist=False) as fq:
        for i,record in enumerate(fq, 1):
            if i > n:
                break
            if record.comment and HAPLOTAGGING_RX.search(record.comment):
                return "haplotagging"
            if STLFR_RX.search(str(record.name)):
                return "stlfr"
            if TELLSEQ_RX.search(str(record.name)):
                return "tellseq"
    return "none"

def which_linkedread_sam(sam: str, n: int = 100) -> str:
    """
    Scans the first 100 records of a SAM/BAM file and tries to determine the barcode technology
    Returns one of: "haplotagging", "stlfr", "tellseq", or "none"
    """
    with pysam.AlignmentFile(sam, check_sq=False) as _sam:
        i = 1
        for record in _sam.fetch(until_eof=True):
            if i > n:
                break
            if record.has_tag("BX"):
                bc = str(record.get_tag("BX"))
                if HAPLOTAGGING_SIMPLE.search(bc):
                    return "haplotagging"
                if STLFR_SIMPLE.search(bc):
                    return "stlfr"
                if TELLSEQ_SIMPLE.search(bc):
                    return "tellseq"
            i += 1
    return "none"

def validate_fq(ctx, param, value):
    """
    Take input fastq files or sam/bam file and do quick checks. Either errors or returns None.
    """
    if len(value) > 2:
        print_error("incorrect number of input files","Inputs must be 1 or 2 FASTQ (.fastq|.fq) files, which can be gzipped.")
        sys.exit(1)

    if len(value) == 2 and value[0] == value[1]:
        print_error('identical fastq files', 'The two input fastq files cannot be the same file')
    
    re_ext = re.compile(r"\.(fq|fastq)(?:\.gz)?$", re.IGNORECASE)
    for i in value:
        if not re_ext.search(i):
            print_error('unrecognized format', "Inputs must be 1 or 2 FASTQ (.fastq|.fq) files, which can be gzipped.")
    return value

def validate_sam(ctx, param, value):
    """
    Take input fastq files or sam/bam file and do quick checks. Either errors or returns None.
    """
    if not value[0].lower().endswith(".bam") or value[0].lower().endswith(".sam"):
        print_error('unrecognized format','Input must be 1 SAM (.sam|.bam) file.')

    return value

def make_dir(ctx, param, value):
    """
    CLI callback method create the output directory preceding the prefix in case it doesn't exist.
    e.g. `prefix = this/that` will create `this` folder, whereas `prefix = that` won't create anything.
    """
    if os.path.dirname(value):
        os.makedirs(os.path.dirname(value), exist_ok=True)
    return value