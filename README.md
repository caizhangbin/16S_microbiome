# 16S_microbiome

**Public repository** for downstream analysis of 16S rRNA gene sequencing data processed through the mothur pipeline.

## 📖 Overview

This collection of R scripts takes your mothur‐generated OTU table, taxonomy assignments, and sample metadata and produces:

- **OTU rarefaction** (`rarefy_OTU.R`)  
- **Alpha diversity** metrics (Observed OTUs, Shannon, Simpson…) (`alpha_diversity.R`)  
- **Beta diversity** distance matrices & ordination plots (PCoA, NMDS…) (`beta_diversity.R`)  
- **Good’s coverage** calculations & rarefaction curves (`goods_coverage.R` & `goods_coverage_and_rarefaction.R`)  
- **Relative abundance** barplots at any taxonomic level (`relative_abundance.R`)  


---


---

## ⚙️ Requirements

- **R** (≥ 4.0)  
- **Packages** (install with `install.packages()` or Bioconductor):
  ```r
  install.packages(c(
    "vegan", "ggplot2", "dplyr", "tidyr", 
    "readr", "QsRutils"
  ))
📚 References
mothur MiSeq SOP
https://mothur.org/wiki/miseq_sop/

Phyloseq documentation & tutorials
https://joey711.github.io/phyloseq/

## Examples images
![D3_chao1](https://github.com/user-attachments/assets/d0e983db-e67a-4ca8-ad12-99a5fb31c670)
![D3_beta](https://github.com/user-attachments/assets/14528f9f-0415-425b-a9f7-3b23a09892a2)

