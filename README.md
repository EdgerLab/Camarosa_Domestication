# Purpose:
This repository contains data, analyses, and scripts related to the investigation of TEs, domestication, and genome diversity of wild and cultivated strawberries.
The project is at its core, a comparative genomics project that examines how TEs influence gene expression, genome architecture, and evolutionary traits in strawberry.

# Abstract:
Transposable elements (TEs) are a major source of mutation and genetic diversity.
In plants, TEs have been linked to a variety of traits, including fruit color, disease resistance, stress response, crop domestication, and more.
However, the role of TEs in the domestication of octoploid strawberry (*Fragaria x ananassa*) remains poorly understood.
In this study we leverage a new high-quality reference genome for domesticated strawberry's wild octoploid progenitor *Fragaria chiloensis* (`Del Norte'), and a newly constructed TE-pangenome to examine the influences of TEs on domestication-related phenotypes.
Using comparative genomics, we show significant differences in TE content near positionally conserved orthologs between wild and domesticated strawberry.
Functional enrichment analyses point to distinct patterns of TE-associated genes between the two genomes, particularly in pathways relevant to strawberry breeding.
We identified a number of fruit, defense, transcription factor, and development genes enriched for novel TE presence in domesticated strawberry, demonstrating the potential for TEs to influence a variety of important domestication related traits.
Together, these findings point to a dynamic and important role for TEs in shaping phenotypes in domesticated strawberry, laying the groundwork for future functional validation studies.

# Table of Contents:
- Overview
- Installation
- Usage
- Repository Structure
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

# Installation:
This project primarily uses Python 3.10.4, with other softwares such as EDTA described further below.
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

# Other Installation Requirements:
I used [TE Density](https://github.com/sjteresi/TE_Density), [EDTA](https://github.com/oushujun/EDTA), and [JCVI](https://github.com/tanghaibao/jcvi).
I listed all 3 as Git Submodules to keep track of the commit I used for my project.
Please take the time to read more about [submodules](https://git-scm.com/book/en/v2/Git-Tools-Submodules) if you are unfamiliar.


# TODO
