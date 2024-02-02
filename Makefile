# scripts for development
# __file__ Makefile
# __author__ Scott Teresi

ROOT_DIR:=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))
DATA_DIR := $(ROOT_DIR)/data
GENOMES_DIR := $(DATA_DIR)/Genomes
RESULTS_DIR := $(ROOT_DIR)/results

# TEs
DEV_H4_UNCLEAN_TEs := $(RESULTS_DIR)/Pan_Annotation/H4_NewNames.fa.mod.EDTA.TEanno.gff3
DEV_DN_UNCLEAN_TEs := $(RESULTS_DIR)/Pan_Annotation/DN_NewNames.fa.mod.EDTA.TEanno.gff3
DEV_RR_UNCLEAN_TEs := $(RESULTS_DIR)/Pan_Annotation/RR_NewNames.fa.mod.EDTA.TEanno.gff3

DEV_FVI_UNCLEAN_TEs := $(RESULTS_DIR)/Pan_Annotation/FVI_NewNames.fa.mod.EDTA.TEanno.gff3
DEV_FNI_UNCLEAN_TEs := $(RESULTS_DIR)/Pan_Annotation/FNI_NewNames.fa.mod.EDTA.TEanno.gff3
DEV_FII_UNCLEAN_TEs := $(RESULTS_DIR)/Pan_Annotation/FII_NewNames.fa.mod.EDTA.TEanno.gff3


# NOTE TODO these file paths no longer work
DEV_H4_UNCLEAN_GENES := $(DATA_DIR)/Genomes/H4/H4_GeneAnnotation.gff
DEV_DN_UNCLEAN_GENES := $(DATA_DIR)/Genomes/Del_Norte/DN_GeneAnnotation.gff
DEV_RR_UNCLEAN_GENES := $(DATA_DIR)/Genomes/Royal_Royce/RR_GeneAnnotation.gff

DEV_FVI_UNCLEAN_GENES := $(DATA_DIR)/Genomes/F_viridis/FVI_GeneAnnotation.gff
DEV_FNI_UNCLEAN_GENES := $(DATA_DIR)/Genomes/F_nipponica/FNI_GeneAnnotation.gff
DEV_FII_UNCLEAN_GENES := $(DATA_DIR)/Genomes/F_iinumae/FII_GeneAnnotation.gff

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
	mkdir -p $(RESULTS_DIR)/cleaned_annotations/
	python $(ROOT_DIR)/src/import_strawberry_EDTA.py $(DEV_H4_UNCLEAN_TEs) $(RESULTS_DIR)/cleaned_annotations/ H4
	python $(ROOT_DIR)/src/import_strawberry_EDTA.py $(DEV_DN_UNCLEAN_TEs) $(RESULTS_DIR)/cleaned_annotations/ DN
	python $(ROOT_DIR)/src/import_strawberry_EDTA.py $(DEV_RR_UNCLEAN_TEs) $(RESULTS_DIR)/cleaned_annotations/ RR
	python $(ROOT_DIR)/src/import_strawberry_EDTA.py $(DEV_FVI_UNCLEAN_TEs) $(RESULTS_DIR)/cleaned_annotations/ FVI
	python $(ROOT_DIR)/src/import_strawberry_EDTA.py $(DEV_FNI_UNCLEAN_TEs) $(RESULTS_DIR)/cleaned_annotations/ FNI
	python $(ROOT_DIR)/src/import_strawberry_EDTA.py $(DEV_FII_UNCLEAN_TEs) $(RESULTS_DIR)/cleaned_annotations/ FII

filter_genes:
	@echo Filtering strawberry genes into appropriate format for TE Density
	mkdir -p $(RESULTS_DIR)/cleaned_annotations/
	python $(ROOT_DIR)/src/import_strawberry_gene_anno.py $(DEV_H4_UNCLEAN_GENES) $(RESULTS_DIR)/cleaned_annotations/ H4
	python $(ROOT_DIR)/src/import_strawberry_gene_anno.py $(DEV_DN_UNCLEAN_GENES) $(RESULTS_DIR)/cleaned_annotations/ DN
	python $(ROOT_DIR)/src/import_strawberry_gene_anno.py $(DEV_RR_UNCLEAN_GENES) $(RESULTS_DIR)/cleaned_annotations/ RR
	python $(ROOT_DIR)/src/import_strawberry_gene_anno.py $(DEV_FVI_UNCLEAN_GENES) $(RESULTS_DIR)/cleaned_annotations/ FVI
	python $(ROOT_DIR)/src/import_strawberry_gene_anno.py $(DEV_FNI_UNCLEAN_GENES) $(RESULTS_DIR)/cleaned_annotations/ FNI
	python $(ROOT_DIR)/src/import_strawberry_gene_anno.py $(DEV_FII_UNCLEAN_GENES) $(RESULTS_DIR)/cleaned_annotations/ FII


#-------------------------------------------------------------------#
# Orthology Analysis
#
.PHONY: filter_RR_DN_syntelogs
filter_RR_DN_syntelogs:
	@echo Filtering Royal Royce and Del Norte SynMap results...
	mkdir -p $(DATA_DIR)/orthologs/filtered
	python $(ROOT_DIR)/src/orthologs/syntelogs.py $(DATA_DIR)/orthologs/RR_DN_SynMap.txt DN $(DATA_DIR)/orthologs/filtered/Cleaned_RR_DN_Syntelogs.tsv

.PHONY: filter_RR_H4_syntelogs
filter_RR_H4_syntelogs:
	@echo Filtering Royal Royce and H4 SynMap results
	mkdir -p $(DATA_DIR)/orthologs/filtered
	python $(ROOT_DIR)/src/orthologs/syntelogs.py $(DATA_DIR)/orthologs/RR_H4_SynMap.txt H4 $(DATA_DIR)/orthologs/filtered/Cleaned_RR_H4_Syntelogs.tsv

# Filter the BLAST results for RR and H4
.PHONY: filter_RR_H4_homologs
filter_RR_H4_homologs:
	@echo TODO
	python $(ROOT_DIR)/src/orthologs/reformat_RR_H4_BLAST_results.py $(DATA_DIR)/orthologs/RR_H4.blast $(DATA_DIR)/orthologs/filtered/RR_H4_BLAST_renamed.txt

# Replace the names for the Del Norte BLAST results and filter the results
.PHONY: filter_RR_DN_homologs
filter_RR_DN_homologs:
	@echo Replacing Del Norte BLAST names
	@echo TODO
	python $(ROOT_DIR)/src/orthologs/replace_and_reformat_DN_RR_BLAST_results.py $(DATA_DIR)/orthologs/RR_DN.blast $(DATA_DIR)/orthologs/DN_salt.translation $(DATA_DIR)/orthologs/filtered/RR_DN_BLAST_renamed.txt

.PHONY: create_pan_orthology_table
create_pan_orthology_table:
	mkdir -p $(RESULTS_DIR)/orthologs
	python $(ROOT_DIR)/src/orthologs/pan_orthology_table.py $(DATA_DIR)/orthologs/filtered/Cleaned_RR_H4_Syntelogs.tsv $(DATA_DIR)/orthologs/filtered/RR_H4_BLAST_renamed.txt $(DATA_DIR)/orthologs/filtered/Cleaned_RR_DN_Syntelogs.tsv $(DATA_DIR)/orthologs/filtered/RR_DN_BLAST_renamed.txt $(DATA_DIR)/orthologs/H4-At_Orthologs_March2022.tsv $(RESULTS_DIR)/orthologs/



#-------------------------------------------------------------------#
# TODO this may need to be changed when I start looking at the other genomes
.PHONY: filter_RR_expression
filter_RR_expression:
	@echo Filtering RoyalRoyce expression data
	python $(ROOT_DIR)/src/expression_data.py $(DATA_DIR)/Genomes/Royal_Royce/RoyalRoyce_Cold_Warm_count_matrix.csv $(RESULTS_DIR)

# NOTE DEPRECATED 1/23/2024
# NOTE this translates the Royal Royce CDS FASTA
.PHONY: translate_CDS_to_protein
translate_CDS_to_protein:
	mkdir -p $(DATA_DIR)/Royal_Royce/RR_Proteins
	@echo Converting Royal Royce files
	python $(ROOT_DIR)/src/orthologs/translate_cds_fasta_to_protein.py $(DATA_DIR)/Genomes/Royal_Royce/RR_CDS.fa $(DATA_DIR)/orthologs/RR_CDS_as_Proteins.fa
	@echo Converting H4 files
	python $(ROOT_DIR)/src/orthologs/translate_cds_fasta_to_protein.py $(DATA_DIR)/Genomes/H4/H4_CDS.fa $(DATA_DIR)/orthologs/H4_CDS_as_Proteins.fa

.PHONY: format_protein_database
format_protein_database:
	# TODO this is an sbatch script

# NOTE DEPRECATED
# TODO clean up the paths
.PHONY: clean_RR_tpases
clean_RR_tpases:
	python $(ROOT_DIR)/doc/LINE/clean_CDS.py \
		/mnt/research/edgerpat_lab/Scotty/Strawberry_Domestication/data/Genomes/Royal_Royce/Royal_Royce_CDS.fasta \
		/mnt/research/edgerpat_lab/Scotty/Strawberry_Domestication/doc/LINE/NingTpases_CDS_results.txt \
		/mnt/research/edgerpat_lab/Scotty/Strawberry_Domestication/data/Genomes/Royal_Royce/Royal_Royce_CDS_NewNames.fasta \
		/mnt/research/edgerpat_lab/Scotty/Strawberry_Domestication/results/Fixed_Fastas/ \
