---
label: Count
description: Count linked-read barcodes
category: [fastq,bam]
tags: [fastq,bam]
icon: fold-down
order: 10
---

# :icon-fold-down: Count linked-read barcodes

You might want a table of the counts (frequency) of all the unique barcodes in your data and that's what [!badge corners="pill" text="count"] is for. Records
with missing barcodes are ignored.

### FASTQ
You can use one file of a pair (e.g. R1 or R2), but including both files will prioritize counting the first file and only count occurences of _new_ barcodes in the second file that aren't present in the first file. This is done to avoid double-counting barcodes.
```bash usage
djinn fastq count INPUT...
```

```bash example 
djinn fastq count sample1.F.fq.gz sample1.R.fq.gz > sample1.bxcount
```

### SAM
Similarly, counting barcodes in SAM/BAM files avoids double-counting by counting forward reads if they are unpaired or unpaired, and reverse reads if they are unpaired or the barcode has not been observed yet.
```bash usage
djinn sam extract INPUT
```

```bash example
djinn sam count sample1.bam > sample1.bxcount
```
