# Chatterjee et al. — Nature Communications 2026

Code and analysis scripts for:

> **Chatterjee R, Simoni-Nieves A, et al.** *Distinct gain-of-function p53
> mutations drive immunosuppression in ovarian cancer.*
> *Nature Communications* (2026). DOI: [to be added]

## Overview

This study investigates how different gain-of-function p53 mutations remodel
the immune landscape of ovarian cancer. We combine bulk RNA-seq, pathway
enrichment, PDO–PBMC co-culture assays, high-dimensional flow cytometry, and
cytokine profiling to show that mutant p53 rewires NRF2/ROS and EMT programs,
and to test combination therapies that convert immune-cold tumors into
responsive malignancies.

## Repository structure

- `analyses/` — one R script per figure panel
- `data/raw-data/` — small reference files (full dataset on Zenodo, DOI below)
- `figures/` — generated plots
- `outputs/` — generated tables

## Reproducing the figures

```r
install.packages(c("here", "DESeq2", "ggplot2", "dplyr", "ggrepel",
                   "pheatmap", "RColorBrewer", "scales", "tidyverse",
                   "uwot", "multcomp", "rstatix", "ggpubr", "patchwork",
                   "writexl", "viridis", "pals", "khroma", "wesanderson",
                   "randomcoloR"))
