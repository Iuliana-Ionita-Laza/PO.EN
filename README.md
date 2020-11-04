# PO.EN

<!-- An Elastic-Net Regularized Presence-Only Model -->
This repository contains implementations and R package to accompany our PO-EN paper. The subdirectory [manuscript](https://github.com/Iuliana-Ionita-Laza/PO.EN/tree/master/manuscript) contains the code and data to reproduce the analysis in our paper.

If you use the code here, please cite our publication:

> Zikun Yang et al. A robust semi-supervised model to predict regulatory effects of genetic variants at single nucleotide resolution using massively parallel reporter assays, with applications to fine-mapping of GWAS loci.

## Introduction

Presence-only model with Elastic Net penalty is a regularized generalized linear model training on the presence-absence response. This package provides functions for tuning and fitting the presence-only model. The presence-only model can be used to predict regulatory effects of genetic variants at sequence-level resolution by integrating a large number of epigenetic features and massively parallel reporter assays (MPRAs).

## Manual

Reference manual can be found at https://github.com/Iuliana-Ionita-Laza/PO.EN/blob/master/PO.EN.pdf.

## Vignettes

Vignettes can be found at https://github.com/Iuliana-Ionita-Laza/PO.EN/blob/master/PO.EN_demonstration.pdf.

## Installation

Install PO.EN in R: 

```
  devtools::install_github("Iuliana-Ionita-Laza/PO.EN")
```

# Citation

Zikun Yang et al. A robust semi-supervised model to predict regulatory effects of genetic variants at single nucleotide resolution using massively parallel reporter assays, with applications to fine-mapping of GWAS loci.
