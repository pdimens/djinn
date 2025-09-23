---
label: Sort
description: Sort linked-read data by barcode
category: [fastq, bam]
tags: [fastq, bam]
icon: sort-desc
order: 10
---

# :icon-sort-desc: Sort linked-read data by barcode

Sometimes you want/need linked-read data to be sorted by barcode.
In order to accomplish this with Djinn, the barcode **must** be in a SAM tag (e.g. `BX`, `BC`) FASTQ or SAM/BAM format. If your 
barcodes aren't in a SAM tag, consider using [!badge corners="pill" text="fastq"](/convert_fastq.md) or
[!badge corners="pill" text="standardize"](/standardize.md) to convert it into a compliant format.

```bash usage
djinn sort SAMTAG PREFIX INPUTS...
```

```bash example | sort fastq with barcodes in BC tag
djinn sort BC pebbles_sorted data/pebbles.R1.fq.gz data/pebbles.R2.fq.gz
```

## :icon-terminal: Running Options
{.compact}
| argument          | description                                                                   |
|:------------------|:------------------------------------------------------------------------------|
| `PREFIX`          | [!badge variant="info" text="required"] output filename prefix                |
| `SAMTAG`          | [!badge variant="info" text="required"] 2-letter name of the SAM tag storing the barcode (e.g. `BC` or `BX`)  |
| `INPUTS`          | [!badge variant="info" text="required"] FASTQ pair or SAM/BAM file           |
| `--threads` `-t`     | Number of threads to use (default: 10)                                                  |

