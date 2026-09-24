# Coloc_networks 

## Overview

This repository contains all scripts necessary to analyze colocalization results through network analysis. 

The 2 main steps are:

1. Preparing input files starting from annotated regional associations and colocalization results, throught the script Summary_coloc_results.R. 
Especially, this script extracts:
  - Colocalization_results_filtered: Filtered results containing only significant colocalizations.
  - Some files to study cis-colocalizing signals and identify those linked to new somamers and new Uniprot entries (according to literature review)
  - Hotspots_df: A table listing hotspots found in the pQTL results along with the number of colocalizations at each hotspot.
  - LB_in_hotsposts: A subset of Lonespot Binding (LB) signals that fall within defined hotspots.
2. Perform network analysis, once the summary script has been run and inputs are prepared. 
The Coloc_networks_top_cond folder contains the scripts to run the network analysis, with a dedicated README.

## Input Files for Script 1

The input paths and filenames must be updated at the beginning of the script.
Default inputs are stored within /data/.

### Regional pQTL annotation file

Default file:
Supplementary table 2 of corresponding manuscript: regional associations annotated.

Default filename:

```text
supplementary_table_2.xlsx
```

Update this line:

```r
mapped_LB_gp_ann_va_ann_bl_ann_collapsed_hf_ann <- read_excel("data/supplementary_table_2.xlsx", sheet = "ST2", skip = 1)
```

This file is expected to contain regional pQTL and hotspot annotations, including:

* `CHR`
* `locus_START_END_37`
* `SeqID`
* `cis_or_trans`
* `UniProt_ID`
* `hotspot`
* `full_hotspot_gene_window`
* `uniprot_new_in_somascan7k_vs5k`
* `uniprot_match`

### Colocalization results file

Default file:
Mock dataset subset (first 297231 rows) of the colocalization pairs for hotspot from the corresponding manuscript. 

Default filename:

```text
combined_colocalization_results.csv
```

Update this line:

```r
colocalization_results <- fread(
  "data/combined_colocalization_results_mock_dataset.csv"
)
```

This file is expected to contain paired trait and locus information, including:

* `trait_a`
* `trait_b`
* `locus_a`
* `locus_b`
* `target_a`
* `target_b`
* `PP.H4.abf`
* `top_cond_a`
* `top_cond_b`
* `top_freq_a`
* `top_freq_b`
* `top_freq_geno_a`
* `top_freq_geno_b`
* `top_mlog10pC_a`
* `top_mlog10pC_b`

