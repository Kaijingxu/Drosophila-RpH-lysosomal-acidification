# Drosophila-RpH-lysosomal-acidification

Reproducible R code, data and representative images accompanying the manuscript:

> **Quantitative imaging of lysosomal acidification in Drosophila using a genetically encoded pH sensor**

## Overview

This repository contains the data and R scripts used to quantify Lamp1-positive endolysosomal acidification in *Drosophila* using the Lamp1-pHluorin-mCherry ratiometric reporter, termed RpH.

The repository reproduces the principal analyses presented in the manuscript, including:

- pHlys Red and Lamp1-GFP analysis in the larval fat body
- Bafilomycin A1 validation of RpH in the larval fat body
- Reproducibility across independent larval fat-body experiments
- RP107 pharmacological application in the larval fat body
- Bafilomycin A1 titration in the adult eye
- RP107 pharmacological application in the adult eye
- Image-quality and processing-sensitivity analyses

The scripts generate the statistical models, estimated marginal means, contrasts, quality-control summaries and **ggplot2** objects underlying the manuscript figures.

---

## Repository structure

```text
Drosophila-RpH-lysosomal-acidification/
│
├── data/
│   ├── GMR_RpH/
│   ├── cg_RpH/
│   └── pHlys_Red/
│
├── scripts/
│   ├── GMR_RpH_BafA1.R
│   ├── GMR_RpH_RP107.R
│   ├── cg_RpH_BafA1.R
│   ├── cg_RpH_RP107.R
│   ├── cg_RpH_vehicle_control.R
│   └── pHlys_Red.R
│
├── example_images/
│
├── sessionInfo.txt
├── README.md
├── LICENSE
└── .gitignore