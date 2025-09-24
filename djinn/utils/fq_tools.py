import random
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
    
    def __str__(self):
        """Default string method returns a formatted FASTQ record."""
        return f"@{self.id}\t{self.comment}\n{self.seq}\n+\n{self.qual}\n"

    def convert(self, _type: str, BC: str):
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


class FQPool():
    def __init__(self):
        """Initialize a FASTQ record pool for a specific barcode bc"""
        self.barcode: str = ""
        self.forward = []
        self.reverse = []

    def add(self, fq1, fq2) -> None:
        '''add a read pair to the pool'''
        self.forward.append(fq1)
        self.reverse.append(fq2)

    def randomize_id(self, idx: int) -> None:
        '''replace sequence ID from self.forward[idx] with 3 random integers between 1000 and 99999 replacing the ones in the last 3 positions'''
        seq_id = self.forward[idx].id.split(":")[:4]
        for _ in range(3):
            seq_id.append(str(random.randint(1000, 99999)))
        seq_id = ":".join(seq_id)
        self.forward[idx].id = seq_id
        self.reverse[idx].id = seq_id

    def spoof_hic(self, r1_filecon, r2_filecon, singletons:bool, n: int) -> None:
        '''
        Create all possible unique combinations of forward and reverse reads and write to open file
        connections r1_filecon and r2_filecon. Randomizes the last three numbers in the sequence ID
        in the process to make sure reads don't have identical read headers. Resets the self.forward,
        self.reverse and self.barcode when done.
        '''
        n_reads = len(self.forward)
        n_choice = min(n, n_reads)
        if n_reads == 1 and singletons:
            r1_filecon.write(str(self.forward[0].convert("tellseq", self.forward[0].barcode)))
            r2_filecon.write(str(self.reverse[0].convert("tellseq", self.forward[0].barcode)))
        elif n_choice == 1:
            # only 1 pair requested, make sure it doesn't pair with itself
            all_idx = set(range(n_reads))
            for idx,r1 in enumerate(self.forward):
                options = list(all_idx.difference([idx]))
                r2 = self.reverse[random.sample(options, k = 1)[0]]
                self.randomize_id(idx)
                r1_filecon.write(str(r1.convert("tellseq", r1.barcode)))
                r2_filecon.write(str(r2.convert("tellseq", r1.barcode)))
        else:
            for idx,r1 in enumerate(self.forward):
                for r2 in random.sample(self.reverse, k = n_choice):
                    self.randomize_id(idx)
                    r1_filecon.write(str(r1.convert("tellseq", r1.barcode)))
                    r2_filecon.write(str(r2.convert("tellseq", r1.barcode)))
        self.__init__()
