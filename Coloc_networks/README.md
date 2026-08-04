# Coloc_networks overview

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

# Input Files for Script 1

The input paths and filenames must be updated at the beginning of the script.

## Regional pQTL annotation file

Default filename:

```text
mapped_LB_gp_ann_va_ann_bl_ann_collapsed_hf_ann.csv
```

Update this line:

```r
mapped_LB_gp_ann_va_ann_bl_ann_collapsed_hf_ann <- read_delim(
  "mapped_LB_gp_ann_va_ann_bl_ann_collapsed_hf_ann.csv",
  delim = ";",
  escape_double = FALSE,
  trim_ws = TRUE
)
```

This file is expected to contain regional pQTL and hotspot annotations, including:

* `chr`
* `start`
* `end`
* `phenotype_id`
* `cis_or_trans`
* `UniProt_ID`
* `hotspot`
* `full_hotspot_gene_window`
* `new_somamer`
* `uniprot_match`

## Colocalization results file

Default filename:

```text
14-Apr-25_combined_colocalization_results.csv
```

Update this line:

```r
colocalization_results <- fread(
  "14-Apr-25_combined_colocalization_results.csv"
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
