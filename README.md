# scRNA-seq Analysis Pipeline RUNX1 Human Data (Phs003508.v1.p1) - Platform: 10x Genomics

This repository contains a complete pipeline for analyzing single-cell RNA-seq data and single-cell transposable elements.

Data downloaded from Phs003508.v1.p1 - DOI: 10.1126/scitranslmed.adn9832

This dataset has private access the user has to reach out to the authors to receive the dbGaP repository key approval.

Human genome reference 10x Genomics: refdata-gex-GRCh38-2024-A

The analysis uses multiple programming languages:
- Bash (for preprocessing and pipeline steps)
- Python (for handling .h5ad files)
- R (for downstream analysis and visualization)

All scripts used in this project are already included in the folders of this repository.

------------------------------------------------------------

## Main Report

The primary deliverable is a rendered R Markdown report:

RUNX1_code.html – rendered analysis report (recommended starting point)
RUNX1_code_code.Rmd – R Markdown source containing all code and narrative
Please open the RUNX1_code_code.html file in a web browser to view the full analysis.

All data and repositories are stored in "Z:\home\common\Data\NM\RUNX1_Phs003508"

------------------------------------------------------------

## Repository Structure

- `counts/` – input data objects
  *(Not included in this repository due to dataset size and restrictions)*
- `csv/` – CSV files generated during analysis
  *(Not included in this repository due to dataset size and restrictions)*
  *(Not included in this repository due to file size)*
- `gencode/` – files related to GENCODE annotation
  *(Not included in this repository due to dataset size and restrictions)*
- `hscscore/` – files and outputs from HSC score analysis
  *(Not included in this repository due to dataset size and restrictions)*
- `plots/` – figures and plots generated during the analysis  
- `rda/` – saved `.rda` objects  
  *(Not included in this repository due to file size)*

All code used in the analysis is contained directly in the R Markdown document.

Large data files were not added to their respective repository.

------------------------------------------------------------

## Notes on data files

A large serialized R workspace object (workspace.rda) generated during analysis
is not included in the repository to keep it lightweight. All figures and results
are fully reproducible from the provided R Markdown document.

- Paths in scripts may need to be adjusted to your system
- Some steps may require an HPC environment (cluster)
- Make sure all required software is installed before running

------------------------------------------------------------

## scTE Important Note

To convert .h5ad files into .h5seurat format:

- You MUST use Seurat version lower than 5.0
- After the conversion is completed, upgrade Seurat to version higher than 5.0

This step is required for compatibility between file formats.

------------------------------------------------------------

## Required Files

Links to download all required data are provided inside the RMarkdown files.

------------------------------------------------------------

## How to Run the Pipeline

The main steps of the analysis include:

1. FASTQ preprocessing and merging (Bash scripts)
	- 0_cellranger.sh
2. BAM file processing using samtools
	- included in "0_cellranger.sh"
3. Transposable element quantification using scTE
	- 1_scte.sh
4. Downstream analysis using R, Bash and Python
	- 2_hscscore.py
	
------------------------------------------------------------

## Package Versions

All R package versions used in this project are recorded using:

sessionInfo()

Please check the RUNX1_code.html file in this repository to see the exact versions.

------------------------------------------------------------

## Author

Prepared by: **Silva Morosini, Natalia**  
Date: **June 03, 2026**