import pysam
from djinn.utils.barcodes import ANY_INVALID, TELLSEQ_STLFR_RX

def extract_barcodes_sam(bamfile: str, separate_invalid: bool = False):
    '''
    Return a set of the unique barcodes in bamfile. If return_invalid is True,
    returns a tuple of set(valid),set(invalid) barcodes.
    '''
    barcodes = set()
    invalid = set()
    with pysam.AlignmentFile(bamfile, check_sq=False) as infile:
        for record in infile.fetch(until_eof=True):
            _bc = None
            if record.has_tag("BX"):
                _bc = record.get_tag("BX")
            else:
                _check =  TELLSEQ_STLFR_RX.search(record.name)
                if _check:
                    _bc = _check.group(0)[1:]

            if _bc:
                if separate_invalid and ANY_INVALID.search(_bc):
                    invalid.add(_bc)
                barcodes.add(_bc)

    if separate_invalid:
        return barcodes, invalid
    return barcodes