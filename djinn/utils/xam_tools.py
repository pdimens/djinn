import pysam
from djinn.utils.barcodes import TELLSEQ_STLFR_RX

def bam_barcode(record: pysam.AlignedSegment) -> str:
    if record.has_tag("BX"):
        return str(record.get_tag("BX"))
    else:
        _check =  TELLSEQ_STLFR_RX.search(str(record.query_name))
        if _check:
            return _check.group(0)[1:]
    return ""