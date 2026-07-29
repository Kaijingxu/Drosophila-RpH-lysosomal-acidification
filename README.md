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
```

---

## Requirements

The analyses were performed in **R 4.3.1** using packages including:

- tidyverse
- readxl
- glmmTMB
- emmeans
- sandwich
- ggpubr

The complete R version, operating-system information and package versions are provided in `sessionInfo.txt`.

---

## Running the analyses

Each script is independent and can be run separately.

1. `pHlys_Red.R`
2. `cg_RpH_BafA1.R`
3. `cg_RpH_vehicle_control.R`
4. `cg_RpH_RP107.R`
5. `GMR_RpH_BafA1.R`
6. `GMR_RpH_RP107.R`

The scripts assume that the working directory is the repository root and read the corresponding data files from the `data/` directory.

---

## Outputs

The scripts generate:

- Regression models
- Estimated marginal means
- Genotype and treatment contrasts
- Multiplicity-adjusted statistical comparisons
- Image-quality summaries
- Processing-sensitivity summaries
- **ggplot2** objects corresponding to the manuscript figures

The plots and statistical summaries are generated within the R session. Final figure assembly and export were performed during manuscript preparation.

---

## Image analysis

Microscopy images were quantified in FIJI.

Detailed FIJI image-analysis procedures, including tissue-specific segmentation, background measurement, fluorescence correction and calculation of RpH ratios, are described in the Supplementary Methods of the associated manuscript.

Representative microscopy images are provided in the `example_images/` directory.

---

## Interpretation of RpH measurements

RpH provides a relative fluorescence readout of acidification within Lamp1-positive endolysosomal compartments.

Because in situ pH calibration was not performed, RpH ratios should not be interpreted as absolute lysosomal pH values. Comparisons are intended to be made within the corresponding tissue-specific acquisition and analysis workflow.

---

## Data availability

The minimally processed datasets used in the manuscript are included in the `data/` directory.

Because of file size, the complete raw microscopy dataset is not hosted in this repository but is available from the corresponding author upon reasonable request.

---

## Code availability

All statistical analyses were performed in R and are reproducible using the scripts provided in this repository.

The FIJI image-analysis workflows are described in the Supplementary Methods of the associated manuscript.

---

## Citation

If you use this repository, please cite the associated manuscript once published.