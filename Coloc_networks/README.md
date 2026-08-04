# Coloc_networks
This repository contains all scripts necessary to analyze colocalization results through network analysis. The scripts and generated results are also available in the following directory on the shared system:
"/group/diangelantonio/users/alessia_mapelli/pQTL/INTERVAL/Hotspots and lonespots/New_coloc_results/"

To begin the analysis, start with the script: Summary_coloc_results.R. 
Starting from the colocalization results extarct:
  - Colocalization_results_filtered: Filtered results containing only significant colocalizations.
  - Some files to study cis-colocalizing signals and identify those linked to new somamers and new Uniprot entries (according to literature review)
  - Hotspots_df: A table listing hotspots found in the pQTL results along with the number of colocalizations at each hotspot.
  - LB_in_hotsposts: A subset of Lonespot Binding (LB) signals that fall within defined hotspots.
    
Once the summary script has been run and inputs are prepared, you can proceed to the network analysis. The Coloc_networks_top_cond folder contains the scripts to run the network analysis. 
