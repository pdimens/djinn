from itertools import product
from random import sample
import re

STLFR_INVALID_RX = re.compile("^0_|_0_|_0$")
TELLSEQ_HAPLOTAGGING_INVALID_RX = re.compile(r"(?:N|[ABCD]00)")
HAPLOTAGGING_RX = re.compile(r'\s?BX:Z:(A[0-9]{2}C[0-9]{2}B[0-9]{2}D[0-9]{2})')
STLFR_RX = re.compile(r'#([0-9]+_[0-9]+_[0-9]+)(\s|$)')
TELLSEQ_RX = re.compile(r':([ATCGN]+)(\s|$)')
TELLSEQ_STLFR_RX = re.compile(r"(?:\:[ATCGN]+$|#\d+_\d+_\d+$)")

class haplotagging():
    def __init__(self):
        """Initialize a haplotagging object that generates formatted barcodes on the fly"""
        self.length = 0
        self.inventory = {}
        self.invalid = "A00C00B00D00"
        self.barcodes = product(
            ["A" + str(i).zfill(2) for i in sample(range(1,97), 96)],
            ["C" + str(i).zfill(2) for i in sample(range(1,97), 96)],
            ["B" + str(i).zfill(2) for i in sample(range(1,97), 96)],
            ["D" + str(i).zfill(2) for i in sample(range(1,97), 96)]
        )

    def next(self):
        return "".join(next(self.barcodes))

class stlfr():
    def __init__(self):
        """Initialize a stlfr object that generates formatted barcodes on the fly"""
        self.length = 0
        self.inventory = {}
        self.invalid = "0_0_0"
        self.barcodes = product(*[sample(range(1,1538), 1537) for i in range(3)])

    def next(self):
        return "_".join(str(i) for i in next(self.barcodes))

class tellseq():
    def __init__(self):
        """Initialize a tellseq object that generates formatted barcodes on the fly"""
        self.length = 0
        self.inventory = {}
        self.invalid = "N" * 18
        self.barcodes = product(*[sample("ATCG", 4) for i in range(18)])

    def next(self):
        return "".join(next(self.barcodes))

class tenx():
    def __init__(self):
        """Initialize a tenx object that generates formatted barcodes on the fly"""
        self.length = 16
        self.inventory = {}
        self.invalid =  "N" * 16
        self.barcodes = product(*[sample("ATCG", 4) for i in range(16)])

    def next(self):
        return "".join(next(self.barcodes))

class _ncbi():
    def __init__(self, preserve_invalid: bool):
        """Initialize an ncbi object that generates formatted barcodes on the fly"""
        self.length = 20
        self.inventory = {}
        if preserve_invalid:
            self.invalid = product(*(["N"] + ["ATCGN" for i in range(19)]))
        else:
            self.invalid =  "N" * 20
        self.barcodes = product(*[sample("ATCG", 4) for i in range(20)])

    def next(self):
        return "".join(next(self.barcodes))
    
    def next_invalid(self):
        try:
            return "".join(next(self.invalid))
        except TypeError:
            return self.invalid



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
            self.valid = "0" in bx.split("_") or bool(re.search(r"(?:N|[ABCD]00)", bx))
        else:
            self.valid = 0