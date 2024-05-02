# How marine are Marine Stramenopiles (MAST)? A cross-system evaluation 

This repository contains code and data included in:

--------
Obiol, A., del Campo, J., de Vargas, C., Mahé, F., and Massana R. How marine are Marine Stramenopiles (MAST)? A cross-system evaluation 

--------

## How to use this code

### 1. Clone this repository

Download all scripts and data by using:

```
git clone https://github.com/aleixop/mast_eukbank.git
```

### 2. Install R packges

The following packages are needed:

```
Biostrings
ComplexHeatmap
ggtext
ggtree
patchwork
phytools
speedyseq
tidyverse
treeio
viridis
```

### 3. Run the code

Open the R scripts (`.Rmd`) and run the code.

## Data and code explained

### Code

Code is divided into three parts, which create the following objects:

- [`1-overall_analysis.Rmd`](1-overall_analysis.Rmd): Figure 2, Table S1, Figure S7
- [`2-mast_analysis.Rmd`](2-mast_analysis.Rmd): Rest of figures and tables (except Figure 1)
- [`3-mast_trees.Rmd`](3-mast_trees.Rmd): Figure S1

Additionaly, [`aux_functions.R`](aux_functions.R) loads packages and functions needed for the scripts above to work.

### Data

Data is divided into the following directories:

- [`data/aux_files/`](data/aux_files/): output of analyses that take long to run (code is also available inside the R scripts, commented) and auxiliary data.
- [`data/phyloseq/`](data/phyloseq/): phyloseq objects for all Stramenopiles ASVs and non-ochrophyta ASVs in EukBank.
- [`data/trees/`](data/trees/): trees, fasta files and alignments for references trees alone and references trees with ASVs.
    * [`data/trees/mast/`](data/trees/mast/): files for separated trees per clade (with and without ASVs), used for Figure S1.
    * [`data/trees/stramenopiles/`](data/trees/stramenopiles/): reference trees of Stramenopiles. Tree from Figure 1 is located in [`main`](data/trees/stramenopiles/main/) and tree used for phylogenetic placement is located in [`epa`](data/trees/stramenopiles/epa/).
- [`data/vegan/`](data/vegan/): NMDS analysis data used for Figure 2.

Additionaly, all currently available MAST references (used in this study) can be found [here](data/trees/tableS3.tsv).

## Original data

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7804946.svg)](https://doi.org/10.5281/zenodo.7804946)

Original data used in this study can be downloaded by clicking the above DOI.
