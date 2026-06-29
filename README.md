# Illumina 450K DNA Methylation Analysis

![R](https://img.shields.io/badge/R-Analysis-276DC3?logo=r&logoColor=white)
![Bioconductor](https://img.shields.io/badge/Bioconductor-minfi-87B13F)
![Platform](https://img.shields.io/badge/Platform-Illumina%20450K-F28C28)
![Status](https://img.shields.io/badge/Status-In%20Progress-yellow)

An R-based workflow for quality control, normalization, exploratory analysis, and differential methylation testing of Illumina HumanMethylation450 array data.

**Workflow:** Raw fluorescence data → Quality control → SWAN normalization → PCA → Differential methylation analysis

---

## Overview

This repository contains the code developed for a university group project focused on the analysis of DNA methylation data generated with the Illumina HumanMethylation450 BeadChip.

The analysis compares control (`CTRL`) and disease (`DIS`) samples and includes the inspection of raw fluorescence intensities, sample quality assessment, normalization, exploratory data analysis, and probe-level differential methylation testing.

---

## Analysis Workflow

The project includes the following steps:

1. Import raw methylation data using the `minfi` package.
2. Extract Red and Green fluorescence intensity matrices.
3. Inspect the assigned probe address and identify its probe chemistry.
4. Create a raw methylation dataset.
5. Perform quality control using:
   - median methylated and unmethylated signal intensities;
   - negative control probes;
   - detection p-values.
6. Calculate raw Beta values and M-values.
7. Compare methylation distributions between `CTRL` and `DIS` samples.
8. Normalize the data using SWAN.
9. Compare raw and normalized methylation distributions.
10. Perform principal component analysis to investigate variation associated with:
    - experimental group;
    - biological sex;
    - Sentrix ID.
11. Identify candidate differentially methylated probes using the Mann–Whitney test.
12. Apply multiple-testing correction and visualize the results.

---

## Assigned Analysis Parameters

| Parameter | Assigned value |
| --- | --- |
| Probe address | `10804411` |
| Detection p-value threshold | `0.01` |
| Normalization method | `preprocessSWAN` |
| Differential methylation test | Mann–Whitney test |

---

## Main Tools

The analysis is performed in R using packages including:

- `minfi`
- `ggplot2`
- `gplots`
- `qqman`
- `ggrepel`

---

## Repository Contents

```text
illumina-450k-methylation-analysis/
├── README.md
├── methylation_analysis_workflow.R
├── methylation_analysis_report.html
└── .gitignore
```

- `methylation_analysis_workflow.R` contains the complete analysis code.
- `methylation_analysis_report.html` contains the knitted report, figures, results, and interpretation.

---

## Data and Environment Setup

The raw methylation data and the preconfigured R library files are not stored directly in this repository.

The analysis script downloads the required files from Google Drive:

- [Raw data files](https://drive.google.com/drive/u/1/folders/1KumBGpRTrO_EAoWVIwdYvo33fE5Cu0Bd)
- [Preconfigured R libraries](https://drive.google.com/file/d/1jxRXgic6slClamJrgKFRrJtALZBRzFfT)

The Google Drive links must be accessible to anyone with the link for the automated download steps in the R script to work.

The raw data files should be placed in the directory expected by the script, or the corresponding file paths should be updated before running the analysis.

---

## Reproducibility

To run the analysis:

1. Install R and the required CRAN and Bioconductor packages.
2. Download the required raw IDAT files and sample metadata.
3. Update the data paths in `methylation_analysis_workflow.R`.
4. Run the script in sequential order.

The analysis uses a fixed random seed where required to support reproducible results.

---

## Academic Context

This repository was developed as a collaborative university project for educational purposes.

The analysis uses a small sample subset and should be interpreted as an exploratory bioinformatics workflow rather than a clinical or diagnostic study.

---

## Authors

Developed collaboratively by:

- [Maria Kotsikori](https://github.com/cloudy00000)
- [Paolo F. Zarcone](https://github.com/Aikon003)
- [Henrik Celaj](https://github.com/jhenrik14) 
- [Simona Ramella](https://github.com/simoram04-wq) 
- [Name](https://github.com/username) #add
- [Name](https://github.com/username) #add
- [Name](https://github.com/username) #add

---

## License

This project is shared for educational purposes. Please contact the authors before reusing substantial parts of the analysis.
