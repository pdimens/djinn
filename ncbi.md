---
label: NCBI Storage
description: Convert linked-read data to unaligned BAM for NCBI subission
category: [fastq]
tags: [fastq]
icon: package-dependencies
order: 10
---

# :icon-package-dependencies: Convert to/from NCBI format

When submitting sequences to NCBI, [they reformat the read headers](https://www.ncbi.nlm.nih.gov/sra/docs/sra-data-formats/),
which means any useful information in the read headers disappears. This applies to haplotagging, TELLseq, and stLFR FASTQ
formats, where the barcode is encoded in the sequence header, and thus vanishes from the public archive when uploaded to
NCBI. Obviously this isn't great, so we propose a simple approach to solving this problem: uploading sequence data as
unaligned BAM files (recommended).

```bash usage
djinn ncbi PREFIX INPUTS...
```

## :icon-terminal: Running Options
{.compact}
| argument          | description                                                                   |
|:------------------|:------------------------------------------------------------------------------|
| `PREFIX`          | [!badge variant="info" text="required"] output filename prefix                |
| `INPUTS`          | [!badge variant="info" text="required"] FASTQ pair or SAM/BAM file            |
| `-t` `--threads`  | Number of threads to use (default: 10)                                        |


## Convert to unaligned BAM
If you didn't already know, the BAM format is quite flexible and contains all the fields one would already use in FASTQ format.
BAM (or SAM) files can also have unaligned records in them, meaning you can quite easily convert a paired-end set of FASTQ
files into a single unaligned BAM file without any data loss (and also free up disk space). The conversion is a simple
`samtools` command in each direction (`samtools import` and `samtools fastq`), but as a convenience, Djinn provides a
[!badge corners="pill" text="ncbi"] wrapper to accomplish this.

!!! Barcode Placement
**For this to work as intended** the barcodes should be stored in the `BX:Z` tag (or some other SAM-compliant tag e.g. `BC:Z`).
It's possible NCBI will still strip the read name from alignment records.
!!!

```bash losslessly convert to unaligned BAM
djinn ncbi file.R1.fq file.R2.fq > out.bam

# is equivalent to #
samtools import -O BAM -T "*" -1 file.R1.fq file.R2.fq > out.bam
```

## Convert to fastq from BAM
The reverse of this process is conspicuously named `from-ncbi`.

```bash losslessly convert to fastq from unaligned BAM
djinn ncbi PREFIX infile.bam

# is equivalent to #
samtools fastq -N -c 6 -T "*" -1 PREFIX.R1.fq.gz -2 PREFIX.R2.fq.gz infile.bam
```