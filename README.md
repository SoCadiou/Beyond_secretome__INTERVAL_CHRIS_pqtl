# Repository Overview

This repository contains paper-specific code accompanying: **Beyond the classical plasma secretome: genetic architecture and disease associations of the expanded human plasma proteome in 13,445 Europeans**. 

Detailed descriptions of all methods are available in the Supplementary Methods provided with the preprint.

## Code provided in this repository

Paper-specific scripts are provided for:

- **[Epitope annotation of regional and conditional associations](./Epitope_definition)**  
  Scripts used to derive VEP-based functional and epitope-effect annotations of conditionally independent pQTL variants and their LD proxies.

- **[Colocalization network analysis](./Coloc_networks)**  
  Scripts used to construct and analyse networks of robust pQTL-pQTL colocalization signals within and outside pQTL hotspots.

Each folder contains a description of the analysis, software requirements, and a small reproducible example including example input data, instructions to run the analysis, expected output and approximate runtime.

## Software and computational workflows

The analyses reported in the manuscript also used the following key publicly available software and computational resources. Additional software used for cohort-specific genotype processing and data preparation is described in the Supplementary Methods.

- **REGENIE v3.3** — genome-wide association analyses  
  https://github.com/rgcgithub/regenie

- **Nextflow v22.10.1** — workflow execution for the REGENIE analyses  
  The workflow used in this study is available at:  
  https://github.com/HTGenomeAnalysisUnit/nf-pipeline-regenie

- **PLINK 1.9 / PLINK 2.0** — genotype processing and quality control  
  https://www.cog-genomics.org/plink/

- **METAL** — inverse-variance weighted GWAS meta-analysis  
  https://genome.sph.umich.edu/wiki/METAL_Documentation

- **GCTA v1.94.0** — COJO conditional analyses (`cojo-slct` and `cojo-cond`)  
  https://yanglab.westlake.edu.cn/software/gcta/

- **coloc v5.2.3** — fine-mapping and colocalization analyses  
  https://github.com/chr1swallace/coloc

- **TwoSampleMR** — Mendelian randomization analyses  
  https://github.com/MRCIEU/TwoSampleMR

- **Ensembl Variant Effect Predictor (VEP) v113** — functional variant annotation  
  https://www.ensembl.org/info/docs/tools/vep/

- **bcftools / bcftools-liftover** — variant normalization and genome-build conversion  
  https://github.com/samtools/bcftools

Additional R packages used for data processing and annotation included `bigsnpr`, `biomaRt`, `AnnotationDbi` v1.68.0, `org.Hs.eg.db` v3.20.0 and `SomaScan.db`. Additional methodological details and software usage are provided in the Supplementary Methods.

### Conditional analysis and pQTL-pQTL colocalization workflow

Locus definition, GCTA-COJO conditional analyses, fine-mapping and pQTL-pQTL colocalization were performed using established methods described in detail in Supplementary Methods SM6 and SM14–15.

These analyses were orchestrated using a Snakemake workflow developed for the pQTL analysis project. The exact workflow version used for the analyses in this manuscript is provided for computational provenance at:

https://github.com/ht-diva/pqtl_conditional/commit/889a015802db374f0b8baa8e64618a9915494746

The workflow combines established software and methods including GCTA-COJO and coloc and is provided here as a record of the implementation used in the study rather than as a general-purpose supported software package.

## Reproducible examples

Minimal reproducible examples are provided for the two paper-specific analyses contained in this repository.

## Data availability and full-scale reproduction

The example datasets provided here are intended to demonstrate the functionality of the paper-specific code.

Full reproduction of all analyses additionally requires controlled-access individual-level INTERVAL and CHRIS data and external GWAS resources. Access to these datasets is described in the Data Availability statement of the manuscript.

GWAS summary statistics generated in this study are available for scientific and non-commercial use at:

https://pgwas-chris-interval.gm.eurac.edu

## License

Code contained directly in this repository is licensed under the [MIT License](./LICENSE).

External software and workflows linked above remain subject to their respective licences.

## Repository versioning

This repository contains the code and examples made available during peer review of Nature Genetics manuscript NG-A73928.

If the manuscript is accepted for publication, the paper-specific code and reproducible examples will be transferred or archived under the institutional/Centre GitHub repository and a frozen release corresponding to the published version will be created.


## Citation

Preprint: [Beyond the classical plasma secretome: genetic architecture and disease associations of the expanded human plasma proteome in 13,445 Europeans.](https://www.medrxiv.org/content/10.64898/2026.07.23.26358667v1.full)
