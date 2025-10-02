import copy
import random
import subprocess
from djinn.utils.barcodes import STLFR_INVALID_RX

class FQRecord():
    def __init__(self, pysamfq, FORWARD: bool, bc: str, length: int):
        """Initialize a FASTQ record. FORWARD denotes if it's a forward read, bc is the barcode type, length is the length"""
        self.forward = FORWARD
        fr = "1" if self.forward else "2"
        comments = pysamfq.comment.strip().split()
        designation = [i for i in comments if i.startswith(fr)]
        self.illumina_new = designation[0] if designation else f"{fr}:N:0:CAGATC"
        self.illumina_old = f"/{fr}"
        self.id = pysamfq.name.removesuffix(self.illumina_old)
        self.comment = "\t".join(i for i in comments if not i.startswith(fr))
        self.seq = pysamfq.sequence
        self.qual = pysamfq.quality
        self.valid = True
        # account for /1 and /2 formatting at the end of the read (if present)
        if bc == "stlfr":
            # identify the trailing #1_2_3 barcode and remove it from the ID
            if self.id.endswith("#"):
                # is invalid
                self.id = self.id.rstrip("#")
                self.barcode = "0_0_0"
                self.valid = False
            else:
                _id = pysamfq.name.split("#")
                self.barcode = _id.pop(-1)
                self.id = "#".join(_id)
                self.valid = not bool(STLFR_INVALID_RX.search(self.barcode))
        elif bc == "10x":
            # identify the first N bases and remove it from the seq and qual of R1
            if self.forward:
                self.seq = pysamfq.sequence[length:]
                self.qual = pysamfq.quality[length:]
                self.barcode = pysamfq.sequence[:length]
        elif bc == "tellseq":
            # identify the trailing :ATCG barcode and remove it from the ID
            if self.id.endswith(":"):
                # is invalid
                self.id = self.id.rstrip(":")
                self.barcode = "N" * 18
                self.valid = False
            else:
                _id = pysamfq.name.split(":")
                self.barcode = _id.pop(-1)
                self.id = ":".join(_id)
                self.valid = "N" not in self.barcode
        elif bc in ["haplotagging", "standard"]:
            # identify the BX:Z SAM tag and remove it from the comment
            bc = [i for i in self.comment.split() if i.startswith("BX:Z")]
            if len(bc) < 1:
                # is invalid
                self.barcode = "A00C00B00D00"
            else:
                self.barcode = bc.pop().removeprefix("BX:Z:")
            self.comment = "\t".join(i for i in self.comment.split() if not i.startswith("BX:Z"))
            self.valid = "00" not in self.barcode
        else:
            # the barcode was provided outright, so just use it
            # this is used in cases where the R2/R1 doesn't have retrievable barcode info (10X R2)
            self.barcode = bc
            # clear out existing BX or # tags
            self.comment = "\t".join(i for i in self.comment.split() if not i.startswith("BX:Z"))
            if "#" in self.id:
                _id = pysamfq.name.split("#")
                bc = _id.pop(-1)
                self.id = "#".join(_id)
            self.id = self.id.rstrip(":")
    def copy(self):
        return copy.deepcopy(self)

    def __str__(self):
        """Default string method returns a formatted FASTQ record."""
        return f"@{self.id}\t{self.comment}\n{self.seq}\n+\n{self.qual}\n"

    def convert(self, _type: str, BC: str):
        """In-place mutating conversion"""
        if _type == "10x":
            if self.forward:
                self.seq = BC + self.seq
                self.qual = "F"*len(BC) + self.qual
            self.id += f" {self.illumina_new}"
        elif _type == "tellseq":
            self.id += f":{BC} {self.illumina_new}"
        elif _type == "stlfr":
            self.id += f"#{BC} {self.illumina_new}"
        elif _type in ["haplotagging", "standard"]:
            self.id += self.illumina_old
            if not self.comment:
                self.comment = f"VX:i:{int(self.valid)}\tBX:Z:{BC}"
            else:
                _comments = [i for i in self.comment.split() if not i.startswith("BX") and not i.startswith("VX")] + [f"VX:i:{int(self.valid)}", f"BX:Z:{BC}"]
                self.comment = "\t".join(_comments)
        return self

    def convert2(self, _type: str, BC: str):
        """Non-mutating conversion that just returns a new, converted record"""
        new_rec = self.copy()
        if _type == "10x":
            if new_rec.forward:
                new_rec.seq = BC + new_rec.seq
                new_rec.qual = "F"*len(BC) + new_rec.qual
            new_rec.id += f" {new_rec.illumina_new}"
        elif _type == "tellseq":
            new_rec.id += f":{BC} {new_rec.illumina_new}"
        elif _type == "stlfr":
            new_rec.id += f"#{BC} {new_rec.illumina_new}"
        elif _type in ["haplotagging", "standard"]:
            new_rec.id += new_rec.illumina_old
            if not new_rec.comment:
                new_rec.comment = f"VX:i:{int(new_rec.valid)}\tBX:Z:{BC}"
            else:
                _comments = [i for i in new_rec.comment.split() if not i.startswith("BX") and not i.startswith("VX")] + [f"VX:i:{int(new_rec.valid)}", f"BX:Z:{BC}"]
                new_rec.comment = "\t".join(_comments)
        return new_rec


class CachedWriter():
    def __init__(self, gz_R1: subprocess.Popen, gz_R2: subprocess.Popen, cachemax: int = 5000):
        """
        A cache for R1 and R2 reads that writes to the stdin of an open gzip subprocesses when the cache exceeds cachemax entries
        """
        self.r1_cache: list[bytes] = []
        self.r2_cache: list[bytes] = []
        self.writer_r1 = gz_R1.stdin
        self.writer_r2 = gz_R2.stdin
        self.cachemax: int = cachemax

    def add(self, r1: FQRecord|None, r2: FQRecord|None):
        """add r1 and r2 reads as bytestrings to the cache, write and empty cache when cache exceeds seld.max reads"""
        if r1:
            self.r1_cache.append(str(r1).encode("utf-8"))
        if r2:
            self.r2_cache.append(str(r2).encode("utf-8"))

        if len(self.r1_cache) >= self.cachemax or len(self.r2_cache) >= self.cachemax:
            self.write()

    def write(self):
        """writes r1 and r2 cache to self.writer_r* and empties caches"""
        self.writer_r1.write(b"".join(self.r1_cache))
        self.writer_r2.write(b"".join(self.r2_cache))
        self.r1_cache = []
        self.r2_cache = []


class FQPool():
    def __init__(self, gz1, gz2, singletons: bool, max_pairs: int, cachemax: int = 5000):
        """
        Initialize a FASTQ record pool for a given barcode. The pools track barcode and the reads therein.
        Passes things off to a CachedWriter() when spoofing all the reads of a barcode
        """
        self.barcode: str = ""
        self.forward: list[FQRecord] = []
        self.reverse: list[FQRecord] = []
        self.writer = CachedWriter(gz1, gz2, cachemax)
        self.singletons = singletons
        self.max_pairs = max_pairs

    def add(self, fq1: FQRecord, fq2: FQRecord) -> None:
        '''add a read pair to the pool'''
        self.forward.append(fq1)
        self.reverse.append(fq2)

    def randomize_id(self, rec: FQRecord) -> FQRecord:
        '''
        Replace sequence ID from rec with 3 random integers between 1000 and 99999 replacing the ones in the last 3 positions
        Returns a new FQRecord
        '''
        seq_id = rec.id.split(":")[:4]
        for _ in range(3):
            seq_id.append(str(random.randint(1000, 99999)))
        seq_id = ":".join(seq_id)
        rec.id = seq_id
        return rec

    def spoof_hic(self) -> None:
        '''
        Create all self.max_pairs unique combinations of forward and reverse reads and write to open file
        connections r1_filecon and r2_filecon. Randomizes the last three numbers in the sequence ID
        in the process to make sure reads don't have identical read headers. Resets the self.forward,
        self.reverse and self.barcode when done.
        '''
        n_reads = len(self.forward)
        n_choice = min(self.max_pairs, n_reads)
        if n_reads == 1:
            if self.singletons:
                self.writer.add(
                    self.forward[0].convert2("tellseq", self.barcode),
                    self.reverse[0].convert2("tellseq", self.barcode)
                )
        elif n_choice == 1:
            # only 1 pair requested, make sure it doesn't pair with itself
            all_idx = set(range(n_reads))
            for idx in range(len(self.forward)):
                options = list(all_idx.difference([idx]))
                r1 = self.randomize_id(self.forward[idx])
                r2 = self.reverse[random.sample(options, k = 1)[0]]
                r2.id = r1.id
                self.writer.add(
                    r1.convert2("tellseq", self.barcode),
                    r2.convert2("tellseq", self.barcode)
                )
        else:
            for _r1 in self.forward:
                for r2 in random.sample(self.reverse, k = n_choice):
                    r1 = self.randomize_id(_r1)
                    r2.id = _r1.id
                    self.writer.add(
                        r1.convert2("tellseq", self.barcode),
                        r2.convert2("tellseq", self.barcode)
                    )
        # reset the pool, keeping the writer open
        self.barcode = ""
        self.forward = []
        self.reverse = []

