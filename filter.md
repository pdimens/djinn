---
label: Filter
description: Remove singletons or invalid barcodes
category: [fastq, bam]
tags: [fastq, bam]
icon: filter
order: 10
---

# :icon-filter: Remove singletons or invalid barcodes

You may be interested in retaining only "useful" linked reads, that is, reads with a valid barcode or reads
that are not singletons (i.e. there are more than 2 reads with the same barcode). Like their names suggests, [!badge corners="pill" text="filter-invalid"]() will remove invalid barcodes and [!badge corners="pill" text="filter-singletons"]() will remove singletons.
Each of these commands can optionally retain the removed reads in separate files.


## FASTQ

### filter-invalid

```bash usage
djinn fastq filter-invalid [--invalid] PREFIX INPUTS...
```

```bash example | remove invalid barcodes and retain them in separate files
djinn fastq filter-invalid --invalid echidna_valid echidna.R1.fq.gz echidna.R2.fq.gz
```

#### :icon-terminal: Running Options
{.compact}
| argument          | description                                                                   |
|:------------------|:------------------------------------------------------------------------------|
| `PREFIX`          | [!badge variant="info" text="required"] output filename prefix                |
| `INPUTS`          | [!badge variant="info" text="required"] FASTQ or file pair                    |
| `-c` `--cache-size` [!badge variant="ghost" text="hidden"] number of reads to store before writing (bigger is faster, default: `10000`) |
| `-i` `--invalid` | write reads with invalid barcodes to `prefix.invalid.bam` or `prefix.invalid.R[12].fq.gz` |

### filter-singletons

```bash usage
djinn fastq filter-singletons [--singletons] PREFIX INPUT...
```

```bash example | remove singleton barcodes and retain them in separate files
djinn fastq filter-singletons --singletons hyena_filt hyena.R1.fq.gz hyena.R2.fq.gz
```

#### :icon-terminal: Running Options
{.compact}
| argument          | description                                                                   |
|:------------------|:------------------------------------------------------------------------------|
| `PREFIX`          | [!badge variant="info" text="required"] output filename prefix                |
| `INPUT`          | [!badge variant="info" text="required"] FASTQ file or file pair                |
| `-c` `--cache-size` | [!badge variant="ghost" text="hidden"] number of reads to store before writing (bigger is faster, default: `10000`) |
| `-s` `--singletons` | write reads with singleton barcodes to `prefix.singletons.R[12].fq.gz` |


## SAM

### filter-invalid

```bash usage
djinn sam filter-invalid [--invalid] INPUT > output.bam
```

```bash example | remove invalid barcodes and retain them in a separate file
djinn sam filter-invalid --invalid echidnda.invalid.bam echidna.bam > echidna.filtered.bam
```

#### :icon-terminal: Running Options
{.compact}
| argument          | description                                                                   |
|:------------------|:------------------------------------------------------------------------------|
| `INPUT`          | [!badge variant="info" text="required"] SAM/BAM file                           |
| `-i` `--invalid`  | output records with invalid barcodes to this file                    |

### filter-singletons

```bash usage
djinn sam filter-singletons [--singletons] PREFIX INPUT > output.bam
```

```bash example | remove singleton barcodes and retain them in a separate file
djinn sam filter-singletons --singletons hyena.singles.bam hyena.bam > hyena.filter.bam
```

#### :icon-terminal: Running Options
{.compact}
| argument          | description                                                                   |
|:------------------|:------------------------------------------------------------------------------|
| `INPUTS`          | [!badge variant="info" text="required"] SAM/BAM file           |
| `-s` `--singletons` | write valid singleton records to this file |
