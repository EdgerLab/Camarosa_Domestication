# Purpose:
This repository contains data, analyses, and scripts related to the investigation of TEs, domestication, and genome diversity of wild and cultivated strawberries.
The project is at its core, a comparative genomics project that examines the presence of TEs relative to genes, and how TEs influence gene expression, genome architecture, and evolutionary traits in strawberry.

# Abstract:
Transposable elements (TEs) are a major source of mutation and genetic diversity.
In plants, TEs have been linked to a variety of traits, including fruit color, disease resistance, stress response, crop domestication, and more.
However, the role of TEs in the domestication of octoploid strawberry *Fragaria x ananassa* (Royal Royce) remains poorly understood.
In this study we leverage a new high-quality reference genome for domesticated strawberry's wild octoploid progenitor *Fragaria chiloensis* (Del Norte), and a newly constructed TE-pangenome to examine the influences of TEs on domestication-related phenotypes.

Using comparative genomics, we show significant differences in TE content near positionally conserved orthologs between wild and domesticated strawberry.
Functional enrichment analyses point to distinct patterns of TE-associated genes between the two genomes, particularly in pathways relevant to strawberry breeding.
We identified a number of fruit, defense, transcription factor, and development genes enriched for novel TE presence in domesticated strawberry, demonstrating the potential for TEs to influence a variety of important domestication related traits.
Together, these findings point to a dynamic and important role for TEs in shaping phenotypes in domesticated strawberry, laying the groundwork for future functional validation studies.

# Table of Contents:
- Overview
- Repository Structure
- Installation
- Usage
- Acknowledgements
- Data Availability
- Future Directions

# Overview:
The strawberry domestication projects explores the relationship between TEs and genes in cultivated and wild strawberry, focusing on:
- Construction of a TE-pangenome for select *Fragaria* species
- Identification of orthologs between strawberry genomes and *Arabidopsis*
- Comparative analyses of TE distributions in wild and domesticated strawberry
- Identification of genes with unusual or high levels of TE Density
	- Identification of syntelogs with TE differences between cultivated and wild strawberry
- Functional enrichments of these TE-dense genes
- Overlap the TE-dense and/or interesting gene ontology genes with selective sweep regions from previous research.

# Repository Structure:
```
Strawberry_Domestication/
├── data/             # Raw and processed genomic data
├── doc/              # Miscellaneous documentation
├── EDTA/             # EDTA Git Submodule
├── JCVI_Strawberry/  # JCVI Git Submodule
├── logs/             # Miscellaneous logs from HPC jobs
├── requirements/     # Contains pip freeze output for Python package
├── results/          # Output files and figures
├── src/              # Analysis and visualization scripts
├── TE_Density/       # TE Density Git Submodule
├── README.md         # Project overview (this file)
└── Makefile          # Organized commands to recreate analyses
```

# Installation:
This project primarily uses Python 3.10.4, with other software such as EDTA described further below.
To create and activate a new, blank slate Python environment to install packages, assuming you have the correct version of Python active:
```bash
python -m venv /path/to/new/virtual/environment
source <venv>/bin/activate
```

To clone this repository and set up the required environment:
```bash
git clone https://github.com/sjteresi/Strawberry_Domestication.git
cd Strawberry_Domestication
git submodule init
git submodule update
pip install -r requirements/python_requirements.txt
```

## Other Installation Requirements:
I used [TE Density](https://github.com/sjteresi/TE_Density), [EDTA](https://github.com/oushujun/EDTA), and [JCVI](https://github.com/tanghaibao/jcvi).
I listed all 3 as Git Submodules to keep track of the commit I used for my project.
Please take the time to read more about [submodules](https://git-scm.com/book/en/v2/Git-Tools-Submodules) if you are unfamiliar.
Install EDTA according to its instructions, do not use the same environment to hold this software **or** TE Density.

## A Note on Conda:
Please do not use Conda.
Just use pip.
If you must conda, please use miniconda to minimize unforeseen issues.

# Usage:
The `Makefile` is the best authority for recreating the analyses.
If all you do is create a Python environment, install the correct packages, read this README, and read the Makefile, you should be in a good spot.
I tried to organize the Makefile chronologically in the order that I performed analyses, but that wasn't always feasible due to the circular nature of some analyses.

I am not the best at creating Makefiles, and it is likely that a large amount of the rules do not resemble traditional C programming language rules.
Most of the Makefile rules involve transforming some data input via a Python script, and outputting a large number of output text files or figures.
Most of the time this could be confusing as I am running the script on two parallel datasets, cultivated and wild strawberry.


# Acknowledgements:
- Pat Edger
- Adrian Platts
- Nick Panchy
- Ning Jiang
- Jordan Brock

# Data Availability:
- Royal Royce genome: [bioRxiv](https://www.biorxiv.org/content/10.1101/2021.11.03.467115v1.full.pdf)
- Del Norte genome: Sourced from Pat, unpublished
- H4 genome: [Giga Science](https://academic.oup.com/gigascience/article/7/2/gix124/4739363?login=true#supplementary-data)

# Future Directions:
Currently, generating a list of potential TE-impacted genes can be quite manual towards the end of the pipeline.
As described in the [Overview section](#overview), there is a lot of merging and subsetting.

- Basically, genes with high levels of TE presence are identified (this cutoff should be explored more, and I purposefully used a high TE presence cutoff).
- Then, these genes with high TE presence are subsetted against a GO table, so that only genes with Arabidopsis orthologs with functional annotations remain (this cuts down the data by a lot!)
- Next, some analyses proceed, the UpSet plot and tables of the GO terms that are unique to each genome or are interesting from a biological perspective are generated.
- Finally, this table is further subsetted by a list of selective sweep regions (this dataset was obtained from [The Plant Cell](https://academic.oup.com/plcell/article/36/5/1622/7479895)
- Then, I manually inspected this list of genes for a final visual check of interesting GO terms and manually BLAST'ed the genes on NCBI to see if my strawberry gene was at least found in other Rosaceae and to "spot-check" any functional information I was getting from my GO enrichment. Oftentimes I would find that the GO enrichment did not encompass the whole body of knowledge around a gene, and that NCBI BLAST was useful here for giving me some of the informal names for genes such as COP1.

Future work could consider using 
