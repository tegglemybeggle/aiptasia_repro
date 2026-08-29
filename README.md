## Overview

This repository contains code from an independent reproduction of computational analyses from **[Insights into coral bleaching under heat stress from analysis of gene expression in a sea anemone model system](https://doi.org/10.1073/pnas.2015737117)**. 

The original paper compared gene expression in symbiotic and aposymbiotic anemones in the 96 hours following initiation of heat stress. The dataset includes samples collected at 0, 3, 12, 24, 48, and 96 hours from both populations.

I reproduced major portions of the analysis starting from the published sequencing data, including read alignment, clustering of genes of interest, various queries relating to differential expression, and identification and analysis of enriched promoter motifs.

This project was undertaken independently as a way to develop practical experience with RNA-seq analysis and computational biology workflows.

**Project presentation:** [View slides](presentation/aiptasia_presentation.pdf)

## Analysis Workflow

Major steps in the analysis included:


1. **Read processing and quantification** — RNA-seq data were downloaded, aligned to the Aiptasia reference genome with STAR, and quantified with HTSeq.
2. **Differential expression analysis** — DESeq2 was used to identify early stress-responsive genes and compare symbiotic and aposymbiotic animals across the time course.
3. **Clustering and visualization** — early-response genes were clustered by temporal expression pattern and visualized with R and ComplexHeatmap.
4. **Functional and regulatory analysis** — GO enrichment was performed with topGO, and promoter sequences were analyzed with MEME Suite to find recurring motifs

Snakemake was used to automate substantial portions of the workflow.

## Key Results

The analysis successfully reproduced several major findings of the paper, including:

* Two distinct temporal expression patterns among stress-responsive genes
* Similar lists of enriched GO terms for the two clusters
* NFkB and HSF among the most enriched motifs in the promoters of stress-responsive genes
* Moderate correlation between cluster 1/2 and NFkB/HSF enrichment
* Convergence of initially differentially-expressed genes between the two populations over the time course


## Reference Genome and Annotation

The reproduction uses Aiptasia genome version 1.1:

* https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_001417965.1/

The published study used an earlier genome/annotation version, so gene identifiers and some results differ from those reported in the original paper.


## Repository Structure

The repository contains:

* `Snakefile` — workflow definition
* `config/` — workflow configuration files
* `scripts/` — R and other analysis scripts
* `data/` — analysis inputs and intermediate data
* `figures/` — generated figures

Some directories may contain intermediate files generated during development.

## Reproducibility Notes

This project was developed as an independent research reproduction rather than as a packaged software pipeline.

Snakemake automates the majority of the processing and analysis steps, but the repository does not currently provide a fully containerized, one-command reproduction of the entire project.

Several reference files and intermediate inputs were obtained or prepared manually during development. Most of the figures underwent minor cosmetic editing (removing legends, creating some labels, etc.) Exact software versions were also not frozen at the beginning of the project.

