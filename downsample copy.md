---
label: Extract
description: Extract barcodes from data
category: [fastq,bam]
tags: [fastq,bam]
icon: codescan
order: 10
---

# :icon-codescan: Extract barcodes from data
Sometimes you just want to get a list of the barcodes in your data.
Inputs must be one SAM/BAM file or two FASTQ files (R1 and R2, can be gzipped). Both FASTQ and SAM/BAM
inputs expect barcodes to follow the standard haplotagging (BX tag), stlfr (@seq_id#barcode), or tellseq
(@seq_id:barcode) formats. Writes to `stdout`.

```bash usage
djinn extract INPUTS > output.bc
```

```bash example | pull out all the barcodes from a bam file 
djinn extract sample1.bam > sample1.bc
```

```bash example | pull out all the barcodes from a fastq file pair
djinn extract sample2.F.fq.gz sample2.R.fq.gz > sample2.barcodes
```

## :icon-terminal: Running Options
{.compact}
| argument                | description                                                                                          |
|:---------------------|:-----------------------------------------------------------------------------------------------------|
| `INPUTS`             | [!badge variant="info" text="required"] One BAM file or both read files from a paired-end FASTQ pair |
