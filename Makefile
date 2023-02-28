# scripts for development
# __file__ Makefile
# __author__ Scott Teresi

ROOT_DIR:=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))
DATA_DIR := $(ROOT_DIR)/data
GENOMES_DIR := $(DATA_DIR)/Genomes
RESULTS_DIR := $(ROOT_DIR)/results

# TEs
DEV_H4_UNCLEAN_TEs := $(RESULTS_DIR)/Pangenome_Annotation_Output/Fvesca_H4.fasta.mod.EDTA.TEanno.gff3
DEV_DN_UNCLEAN_TEs := $(RESULTS_DIR)/Pangenome_Annotation_Output/Del_Norte_NewNames.fasta.mod.EDTA.TEanno.gff3
DEV_RR_UNCLEAN_TEs := $(RESULTS_DIR)/Pangenome_Annotation_Output/Royal_Royce.fasta.mod.EDTA.TEanno.gff3

DEV_FVI_UNCLEAN_TEs := $(RESULTS_DIR)/Pangenome_Annotation_Output/FVI_NewNames.fasta.mod.EDTA.TEanno.gff3
DEV_FNI_UNCLEAN_TEs := $(RESULTS_DIR)/Pangenome_Annotation_Output/FNI_NewNames.fasta.mod.EDTA.TEanno.gff3
DEV_FII_UNCLEAN_TEs := $(RESULTS_DIR)/Pangenome_Annotation_Output/FII.fasta.mod.EDTA.TEanno.gff3

DEV_H4_UNCLEAN_GENES := $(DATA_DIR)/Genomes/Fvesca_H4/Fvesca_H4_Gene_Annotation.gff
DEV_DN_UNCLEAN_GENES := $(DATA_DIR)/Genomes/Del_Norte/Del_Norte_Gene_Annotation.gff
DEV_RR_UNCLEAN_GENES := $(DATA_DIR)/Genomes/Royal_Royce/Royal_Royce_Gene_Annotation.gff

DEV_FVI_UNCLEAN_GENES := $(DATA_DIR)/Genomes/F_viridis/FVI_maker_annotation.gff
DEV_FNI_UNCLEAN_GENES := $(DATA_DIR)/Genomes/F_nipponica/FNI_maker_annotation.gff
DEV_FII_UNCLEAN_GENES := $(DATA_DIR)/Genomes/F_iinumae/FII_gene_annotation.gff3

setup:
	mkdir -p requirements doc results src data

sync_hpcc_to_onedrive:
	# MUST be standing in root folder for project
	ml Rclone
	rclone sync . remote:HPCC_Mirror/Strawberry_Domestication/ --exclude=.git/** -P -L

#sync_onedrive_to_hpcc:
	# MUST be standing in root folder for project
	#ml Rclone
	#rclone sync remote:HPCC_Mirror/Strawberry_Domestication/ . --exclude=.git/** -P

filter_TEs:
	@echo Filtering strawberry TEs into appropriate format for TE Density
	python $(ROOT_DIR)/src/import_strawberry_EDTA.py $(DEV_H4_UNCLEAN_TEs) $(RESULTS_DIR) H4
	python $(ROOT_DIR)/src/import_strawberry_EDTA.py $(DEV_DN_UNCLEAN_TEs) $(RESULTS_DIR) DN
	python $(ROOT_DIR)/src/import_strawberry_EDTA.py $(DEV_RR_UNCLEAN_TEs) $(RESULTS_DIR) RR
	python $(ROOT_DIR)/src/import_strawberry_EDTA.py $(DEV_FVI_UNCLEAN_TEs) $(RESULTS_DIR) FVI
	python $(ROOT_DIR)/src/import_strawberry_EDTA.py $(DEV_FNI_UNCLEAN_TEs) $(RESULTS_DIR) FNI
	python $(ROOT_DIR)/src/import_strawberry_EDTA.py $(DEV_FII_UNCLEAN_TEs) $(RESULTS_DIR) FII

filter_genes:
	@echo Filtering strawberry genes into appropriate format for TE Density
	python $(ROOT_DIR)/src/import_strawberry_gene_anno.py $(DEV_H4_UNCLEAN_GENES) $(RESULTS_DIR) H4
	python $(ROOT_DIR)/src/import_strawberry_gene_anno.py $(DEV_DN_UNCLEAN_GENES) $(RESULTS_DIR) DN
	python $(ROOT_DIR)/src/import_strawberry_gene_anno.py $(DEV_RR_UNCLEAN_GENES) $(RESULTS_DIR) RR
	python $(ROOT_DIR)/src/import_strawberry_gene_anno.py $(DEV_FVI_UNCLEAN_GENES) $(RESULTS_DIR) FVI
	python $(ROOT_DIR)/src/import_strawberry_gene_anno.py $(DEV_FNI_UNCLEAN_GENES) $(RESULTS_DIR) FNI
	python $(ROOT_DIR)/src/import_strawberry_gene_anno.py $(DEV_FII_UNCLEAN_GENES) $(RESULTS_DIR) FII
