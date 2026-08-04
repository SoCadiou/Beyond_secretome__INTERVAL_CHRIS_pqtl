# pQTL Epitope Annotation

## Overview

This R script identifies pQTL associations that may be influenced by an epitope effect.

It uses VEP annotations for independent COJO variants and nearby proxy variants to detect moderate- or high-impact consequences in protein-coding genes. Results are produced at both the individual COJO association and regional association levels.

## Requirements

Required R packages:

```r
library(tidyverse)
library(data.table)
library(future)
library(furrr)
```

## Inputs

The input paths and filenames must be updated at the beginning of the script.

### Regional pQTL associations

Default filename:

```text
mapped_LB_gp_ann_va_ann_bl_ann_collapsed_hf_ann.csv
```

Update:

```r
path_lb_cistrans <- "mapped_LB_gp_ann_va_ann_bl_ann_collapsed_hf_ann.csv"
```

### Independent COJO SNP associations

Default filename:

```text
16-Dec-24_collected_independent_snps.csv
```

Update:

```r
path_cojo <- "16-Dec-24_collected_independent_snps.csv"
```

The script currently uses `path_freez` when reading this file. Define the directory:

```r
path_freez <- "/path/to/cojo/files/"
```

Alternatively, replace:

```r
fread(paste0(path_freez, path_cojo))
```

with:

```r
fread(path_cojo)
```

### VEP annotation files

Default directory:

```text
/exchange/healthds/pQTL/pQTL_workplace/annotations/VEP/data/unzipped/
```

Update:

```r
path_vep_extract <- "/path/to/extracted/VEP/files/"
```

The extracted annotation files must be located in:

```text
<path_vep_extract>/snps_ld_in_meta_annot/
```

Expected VEP columns include:

* `Consequence`
* `BIOTYPE`
* `SYMBOL`
* `SNPID`

## Annotation logic

For **cis associations**, an epitope effect is flagged when a moderate- or high-impact variant affects the gene encoding the measured protein.

For **trans associations**, the script flags high-impact variants affecting any protein-coding gene.

The script also records affected variants and gene symbols.

## Outputs

### COJO-level output

```text
cojo_epitope_symbol_matching.tsv
```

Contains epitope annotations for each independent COJO association.

### Regional output

```text
mapped_LB_gp_ann_va_ann_bl_ann_collapsed_hf_ann_epitope_symbol_matching.tsv
```

Contains regional summaries, including the number and proportion of epitope-positive COJO signals and the variants and genes implicated by VEP.

Output filenames can be updated using:

```r
path_cojo_epitop <- "cojo_epitope_symbol_matching.tsv"

out_lb_epitop_cojo <- paste0(
  "mapped_LB_gp_ann_va_ann_bl_ann_collapsed_",
  "hf_ann_epitope_symbol_matching.tsv"
)
```

## Running the script

```bash
Rscript annotate_pqtl_epitope.R
```

Replace `annotate_pqtl_epitope.R` with the actual script filename.

## Parallel processing

The script uses 32 parallel workers:

```r
future::plan(multicore, workers = 32)
```

Adjust this value according to the available resources. On Windows, use `multisession` instead of `multicore`.
