# pQTL Epitope Annotation

## Overview

This R script identifies pQTL associations that may be influenced by an epitope effect.

It uses VEP annotations for independent COJO variants and nearby proxy variants to detect moderate- or high-impact consequences in protein-coding genes. 
Results are produced at both the level of COJO associations, i.e flagging each independant COJO association as possibly driven by epitope, and at the level of 
the regional association, i.e if the locus association is possibly driven by epitope (see annotation logic and output below).

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
Demonstration datasets are provided under /data/

### Regional pQTL associations

Default filename:

```text
supplementary_table_2.xlsx
```
Default file is the supplementary table 2 of the corresponding article.

### Independent COJO SNP associations

Default filename:

```text
supplementary_table_3.xlsx
```
Default file is the supplementary table 3 of the corresponding article.


### VEP annotation files

Default directory:

```text
/VEP/snps_ld_in_meta_annot/
```
By default, a subset restricted to chromosome 20 of VEP annotations from the current article is provided.

The extraction of lead SNP, locus and aptamer ID (SeqID) from the VEP name must be adapted depending on the naming convention of VEP files.

Expected VEP columns must include:

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

Contains the input regional associations file, with additional epitope annotation, including the number and proportion of epitope-positive COJO signals per regional association and the variants and genes implicated by VEP.


## Running the script

```bash
Rscript annotate_pqtl_epitope.R
```


## Parallel processing

The script uses 32 parallel workers:

```r
future::plan(multicore, workers = 32)
```

Adjust this value according to the available resources.

Computing time for demonstration subset (VEP annotation restricted to chromosome 20) is monitored in the script, thus allowing estimation of computing time on user-machine

