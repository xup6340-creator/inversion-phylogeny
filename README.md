# Inversion-Based Phylogeny Analysis

This repository contains the analysis pipeline used in the thesis:

"Inversion-based phylogeny under a Dollo-style evolutionary framework"

The project includes:

- Empirical analyses (cotton and radish datasets)
- Simulation studies (Simulation 1 and Simulation 2)

------------------------------------------------------------

## R Environment

All analyses were conducted using:

R version 4.5.1 (2025-06-13 ucrt)  
Platform: x86_64-w64-mingw32/x64  
Operating system: Windows 11 x64  
Time zone: America/Toronto  

Matrix products: default  
LAPACK version: 3.12.1  

Locale:
Chinese (Simplified)_China.utf8  

------------------------------------------------------------

## Core Package Versions

The following packages were attached during analysis:

ape        5.8-1  
phangorn   2.12.1  
dplyr      1.1.4  
ggplot2    3.5.2  
ggtree     3.16.3  
patchwork  1.3.1  

Additional dependencies were loaded automatically.

To install required packages:

    install.packages(c(
      "ape",
      "phangorn",
      "dplyr",
      "ggplot2",
      "patchwork",
      "reshape2"
    ))

For tree visualization:

    install.packages("ggtree")

------------------------------------------------------------

## Directory Structure

data/       Inversion matrices (cotton and radish)  
scripts/    Analysis scripts  
figures/    Generated figures  
outputs/    Simulation outputs  

------------------------------------------------------------

## Reproducing Analyses

### Cotton Dataset

    source("scripts/cotton_analysis.R")

Generates:

- NJ tree
- Reference comparison
- Bipartition overlap
- Randomization test
- Bootstrap MAST analysis
- Noise perturbation analysis

------------------------------------------------------------

### Radish Dataset

    source("scripts/radish_analysis.R")

Generates:

- NJ tree
- Reference comparison
- Bipartition overlap
- Randomization test
- Bootstrap MAST analysis

------------------------------------------------------------

### Simulation 1

    source("scripts/simulation1.R")

Generates:

- 20-bin MAST sensitivity curves
- Cophenetic correlation curves
- Sister-pair preservation analysis

------------------------------------------------------------

### Simulation 2

    source("scripts/simulation2.R")

Generates:

- Example trees with controlled trait ratios
- MAST sensitivity curves under controlled state distributions

------------------------------------------------------------

## Notes on Reproducibility

Results may vary slightly depending on:

- R version
- Package versions
- Platform differences

All major simulations use fixed random seeds to ensure reproducibility.

