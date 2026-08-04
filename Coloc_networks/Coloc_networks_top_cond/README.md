This folder contains scripts to analyze colocalization results within individual hotspots defined in Hotspots_df.csv.

## colocalization_in_hotspots_parallel_target.R
This script takes a single row from the Hotspots_df.csv (i.e., one hotspot) and performs the following steps:
  1. Constructs a colocalization graph for the signals within the hotspot.
  2. Applies a clustering algorithm (defined by the type_of_clustering parameter).
  3. Generates a network plot where nodes are colored by community and node shapes indicate cis vs trans. This plot is mainly for investigation reasons since each node is labeled and a corresponding mapping with the node name is available.
  4. Generates two network plots where nodes are colored by community and node shapes indicate cis vs trans. These plots are mainly for display. 
  5. Saves an annotated table that includes all signals in the hotspot, lists the community assignment, and indicates the type of community (cis, trans, multi_cis, multi_trans, cis_multi_trans, multi_cis_multi_trans).
  6. Creates a positional plot of the genomic region where signals are colored by community and shapes distinguish cis and trans.
  7. Computes and saves the median genomic distance between signals in each community.
  8. Computes and saves the median protein correlation within and between communities.

#### Key Parameters to Edit (inside colocalization_in_hotspots_parallel_target.R):
  - Line 1: Set your working directory – where output files will be saved.
  - Line 16: Path to the filtered colocalization results.
  - Line 17: Path to the hotspot summary dataframe (Hotspots_df.csv).
  - Line 18: Path to the LB results filtered for hotspots.
  - Line 19: Specify the clustering type you want to use.

## Sbatch_parallel.sbatch and Sbatch_parallel_luncher.sh
To process all hotspots efficiently, we use these two scripts to run the analysis in parallel.
Sbatch_parallel.sbatch: Template job script to run a single instance of the R script.
Sbatch_parallel_luncher.sh: Launches 22 parallel jobs, each processing one row (hotspot) from Hotspots_df.csv.

#### Key Parameters to Edit (inside Sbatch_parallel.sbatch):
- Line 2: Set the log directory where SLURM output/error logs will be saved.

## enrichment_proteins_community_hotspots.py
The script enrichment_proteins_community_hotspots.py perform the enrichment analysis for each of the communities in a hotspot with at least 4 seq.ids and save the results
