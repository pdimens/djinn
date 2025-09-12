---
label: Standard Format
description: Convert linked-read data to the Standard format
category: [fastq, bam]
tags: [fastq,bam]
icon: arrow-switch
order: 10
---

In the effort of making it painless to have your data in the preferred standard format, use `djinn std-*`
to quickly standardize FASTQ and BAM files. By default, standardization just moves the barcode (wherever it may be)
into a `BX:Z` SAM tag as-is and does a technology-appropriate validation of the barcode value, which it writes to the
`VX:i` tag. However, you can use `--style` to also convert the barcode style between formats. Keep in mind that each
barcode style has a different upper limit as to how many unique barcodes it can support, which may prevent successful conversions.
The styles are given as:

{.compact}
| Style          | Maximum Unique | What they look like | Example            |
|:---------------|:--------------:|:--------------------|:-------------------|
| `haplotagging` |     $96^4$     | AxxCxxBxxDxx        | A41C22B70D93       |
| `stlfr`        |    $1537^3$    | 1_2_3               | 901_3_1121         |
| `tellseq`      |    $4^{18}$    | 18-base nucleotide  | AGCCATGTACGTATGGTA |
| `10X`          |    $4^{16}$    | 16-base nucleotide  | GGCTGAACACGTGCAG   |

### BAM
If barcodes are present in the sequence name (stlfr, tellseq), this method moves the barcode to the `BX:Z`
tag of the alignment, maintaining the same barcode style by default (auto-detected). If moved to or already in a `BX:Z` tag,
will then write a complementary `VX:i` tag to describe barcode validation `0` (invalid) or `1` (valid).
Use `--style` to also convert the barcode to a different style (`haplotagging`, `stlfr`, `tellseq`, `10X`),
which also writes a `conversion.bc` file to the working directory mapping the barcode conversions. Writes to `stdout`.

```bash usage
djinn standardize-bam [--style] SAM > output.bam
```

```bash example | standardize a bam and change the barcodes to stLFR style
djinn standardize-bam --style stflr yucca.bam > yucca.std.stlfr.bam
```

#### :icon-terminal: Running Options
{.compact}
| argument       | default | description                                                                        |
|:---------------|:-------:|:-----------------------------------------------------------------------------------|
| `SAM`          |         | [!badge variant="info" text="required"] input SAM/BAM alignment file               |
| `--quiet`      |   `0`   | `0` and `1` (all) or `2` (no) output                                               |
| `-s`/`--style` |         | change barcode style in the output BAM: [`10x`,`haplotagging`, `stlfr`, `tellseq`] |

### FASTQ
This conversion moves the barcode to the `BX:Z` tag in fastq records, maintaining the same barcode type by default (auto-detected).
See [this section](https://pdimens.github.io/harpy/getting_started/linked_read_data/#linked-read-data-types) of the Harpy documentation for the location and format expectations for different linked-read technologies.
Also writes a `VX:i` tag to describe barcode validation `0` (invalid) or `1` (valid).
Use `--style` to also convert the barcode to a different style (`haplotagging`, `stlfr`, `tellseq`, `10X`),
which will also write a `conversion.bc` file to the working directory mapping the barcode conversions.

!!!warning Incompatible with 10X data
Standardization will **not** work with the 10X FASTQ format, where the barcodes are the first 16 bases of read 1.
You will first need to demultiplex the barcodes from the sequences into the read headers.
!!!

```bash usage
djinn std-fastq [--style] PREFIX R1.fq R2.fq
```

```bash example | standardize a fastq pair and change the barcodes to stLFR style
djinn std-fastq --style stflr myotis.stlfr myotis.R1.fq.gz myotis.R2.fq.gz
```

#### :icon-terminal: Running Options
{.compact}
| argument       | description                                                                          |
|:---------------|:-------------------------------------------------------------------------------------|
| `PREFIX`       | [!badge variant="info" text="required"] prefix for output filenames                  |
| `R1.fq`        | [!badge variant="info" text="required"] input FASTQ forward-read file                |
| `R1.fq`        | [!badge variant="info" text="required"] input FASTQ reverse-read file                |
| `-s`/`--style` | change barcode style in the output FASTQ: [`10x`,`haplotagging`, `stlfr`, `tellseq`] |


----

:::info Useless trivia
The original version of this command was written while I was waiting at a mechanic shop for car repairs. During development, it was called `lr-switcheroo`.
:::