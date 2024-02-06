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

#-------------------------------------------------------------------#
# Filter the genes for TE Density
# Define the file paths for the uncleaned gene annotations
H4_UNCLEAN_GENES := $(DATA_DIR)/Genomes/H4/H4_GeneAnnotation.gff
DN_UNCLEAN_GENES := $(DATA_DIR)/Genomes/Del_Norte/DN_GeneAnnotation.gff
RR_UNCLEAN_GENES := $(DATA_DIR)/Genomes/Royal_Royce/RR_GeneAnnotation.gff
FVI_UNCLEAN_GENES := $(DATA_DIR)/Genomes/F_viridis/FVI_GeneAnnotation.gff
FNI_UNCLEAN_GENES := $(DATA_DIR)/Genomes/F_nipponica/FNI_GeneAnnotation.gff
FII_UNCLEAN_GENES := $(DATA_DIR)/Genomes/F_iinumae/FII_GeneAnnotation.gff

# Define the file paths for the cleaned gene annotations
H4_CLEAN_GENES := $(RESULTS_DIR)/cleaned_annotations/Cleaned_H4_GeneAnnotation.gff
DN_CLEAN_GENES := $(RESULTS_DIR)/cleaned_annotations/Cleaned_DN_GeneAnnotation.gff
RR_CLEAN_GENES := $(RESULTS_DIR)/cleaned_annotations/Cleaned_RR_GeneAnnotation.gff
FVI_CLEAN_GENES := $(RESULTS_DIR)/cleaned_annotations/Cleaned_FVI_GeneAnnotation.gff
FNI_CLEAN_GENES := $(RESULTS_DIR)/cleaned_annotations/Cleaned_FNI_GeneAnnotation.gff
FII_CLEAN_GENES := $(RESULTS_DIR)/cleaned_annotations/Cleaned_FII_GeneAnnotation.gff

# Define a target to create the output directory if it doesn't exist
$(RESULTS_DIR)/cleaned_annotations:
	mkdir -p $@

$(H4_CLEAN_GENES): $(H4_UNCLEAN_GENES) | $(RESULTS_DIR)/cleaned_annotations
	@echo Filtering H4 strawberry genes into appropriate format for TE Density
	python $(ROOT_DIR)/src/import_strawberry_gene_anno.py $< H4 $@

$(DN_CLEAN_GENES): $(DN_UNCLEAN_GENES) | $(RESULTS_DIR)/cleaned_annotations
	@echo Filtering DN strawberry genes into appropriate format for TE Density
	python $(ROOT_DIR)/src/import_strawberry_gene_anno.py $< DN $@

$(RR_CLEAN_GENES): $(RR_UNCLEAN_GENES) | $(RESULTS_DIR)/cleaned_annotations
	@echo Filtering RR strawberry genes into appropriate format for TE Density
	python $(ROOT_DIR)/src/import_strawberry_gene_anno.py $< RR $@

$(FII_CLEAN_GENES): $(FII_UNCLEAN_GENES) | $(RESULTS_DIR)/cleaned_annotations
	@echo Filtering FII strawberry genes into appropriate format for TE Density
	python $(ROOT_DIR)/src/import_strawberry_gene_anno.py $< FII $@

$(FNI_CLEAN_GENES): $(FNI_UNCLEAN_GENES) | $(RESULTS_DIR)/cleaned_annotations
	@echo Filtering FNI strawberry genes into appropriate format for TE Density
	python $(ROOT_DIR)/src/import_strawberry_gene_anno.py $< FNI $@

$(FVI_CLEAN_GENES): $(FVI_UNCLEAN_GENES) | $(RESULTS_DIR)/cleaned_annotations
	@echo Filtering FVI strawberry genes into appropriate format for TE Density
	python $(ROOT_DIR)/src/import_strawberry_gene_anno.py $< FVI $@

# Define a phony target to filter all gene annotations for convenience
.PHONY: filter_all_genes
filter_all_genes: $(H4_CLEAN_GENES) $(DN_CLEAN_GENES) $(RR_CLEAN_GENES) $(FII_CLEAN_GENES) $(FNI_CLEAN_GENES) $(FVI_CLEAN_GENES)

# Define a phony target to run a singular gene filtering for convenience
# FUTURE if you want more individually type the phonies out
.PHONY: filter_H4_genes
filter_H4_genes: $(H4_CLEAN_GENES)

#-------------------------------------------------------------------#
# Orthology Analysis

# Define the file paths for the cleaned syntelogs files
CLEANED_RR_DN_SYNTELOGS := $(DATA_DIR)/orthologs/filtered/Cleaned_RR_DN_Syntelogs.tsv
CLEANED_RR_H4_SYNTELOGS := $(DATA_DIR)/orthologs/filtered/Cleaned_RR_H4_Syntelogs.tsv

# Define the file paths for the renamed BLAST results
RR_H4_BLAST_RENAMED := $(DATA_DIR)/orthologs/filtered/RR_H4_BLAST_renamed.txt
RR_DN_BLAST_RENAMED := $(DATA_DIR)/orthologs/filtered/RR_DN_BLAST_renamed.txt

# Define the file paths for the master orthology table that is the final product
STRAWBERRY_ORTHOLOG_TABLE := $(RESULTS_DIR)/orthologs/Strawberry_Arabidopsis_Ortholog_Table.tsv

# Define a target to create the directory if it doesn't exist
$(DATA_DIR)/orthologs/filtered $(RESULTS_DIR)/orthologs:
	mkdir -p $@

# Define the target for generating cleaned RR_DN_Syntelogs file
$(CLEANED_RR_DN_SYNTELOGS): $(DATA_DIR)/orthologs/RR_DN_SynMap.txt | $(DATA_DIR)/orthologs/filtered
	@echo Filtering Royal Royce and Del Norte SynMap results...
	python $(ROOT_DIR)/src/orthologs/syntelogs.py $< DN $@

# Define the target for generating cleaned RR_H4_Syntelogs file
$(CLEANED_RR_H4_SYNTELOGS): $(DATA_DIR)/orthologs/RR_H4_SynMap.txt | $(DATA_DIR)/orthologs/filtered
	@echo "Filtering Royal Royce and H4 SynMap results..."
	python $(ROOT_DIR)/src/orthologs/syntelogs.py $< H4 $@

# Define the target for renaming RR_DN BLAST results
# NOTE there is a gene renaming step in this script
$(RR_DN_BLAST_RENAMED): $(DATA_DIR)/orthologs/RR_DN.blast $(DATA_DIR)/orthologs/DN_salt.translation | $(DATA_DIR)/orthologs/filtered
	@echo "Reformatting the Royal Royce and Del Norte BLAST results..."
	python $(ROOT_DIR)/src/orthologs/replace_and_reformat_DN_RR_BLAST_results.py $^ $@

# Define the target for renaming RR_H4 BLAST results
$(RR_H4_BLAST_RENAMED): $(DATA_DIR)/orthologs/RR_H4.blast | $(DATA_DIR)/orthologs/filtered
	@echo "Reformatting the Royal Royce and H4 BLAST results..."
	python $(ROOT_DIR)/src/orthologs/reformat_RR_H4_BLAST_results.py $< $@

# Define the target for creating the master strawberry-arabidopsis ortholog table
$(STRAWBERRY_ORTHOLOG_TABLE): $(CLEANED_RR_H4_SYNTELOGS) $(RR_H4_BLAST_RENAMED) $(CLEANED_RR_DN_SYNTELOGS) $(RR_DN_BLAST_RENAMED) $(DATA_DIR)/orthologs/H4-At_Orthologs_March2022.tsv | $(RESULTS_DIR)/orthologs
	@echo "Creating the master strawberry-arabidopsis ortholog table..."
	python $(ROOT_DIR)/src/orthologs/pan_orthology_table.py $^ $(RESULTS_DIR)/orthologs/ $@

# Define a phony target to run the orthology analysis for convenience
.PHONY: orthology_analysis
orthology_analysis: $(STRAWBERRY_ORTHOLOG_TABLE)

# FUTURE TODO add a PHONY target to run the syntelogs and BLAST renaming steps if those need to be run individually


#-------------------------------------------------------------------#
# TE Density analysis of Syntelogs
# Define the file paths for the density data
DN_DENSITY_DIR := $(DATA_DIR)/density/DN
RR_DENSITY_DIR := $(DATA_DIR)/density/RR

# Define the file paths for the GeneData
DN_GENE_DATA := $(DATA_DIR)/density/DN
RR_GENE_DATA := $(DATA_DIR)/density/RR

# Define the file path for the syntelog density table
SYNTELOG_W_DENSITY_TABLE := $(RESULTS_DIR)/density_analysis/Orthologs_w_Density.tsv

# Define a target to create the directory if it doesn't exist
$(RESULTS_DIR)/density_analysis:
	mkdir -p $@

# TODO add the paths for the other genomes
$(SYNTELOG_W_DENSITY_TABLE): $(STRAWBERRY_ORTHOLOG_TABLE) $(DN_CLEAN_GENES) $(RR_CLEAN_GENES) $(DN_DENSITY_DIR) $(RR_DENSITY_DIR) $(RESULTS_DIR) | $(RESULTS_DIR)/density_analysis
	echo TODO
	python $(ROOT_DIR)/src/syntelog_differences/parse_density_data.py $^

# Define a phony target to create the syntelog density table for convenience
.PHONY: create_sytelog_density_table
create_sytelog_density_table: $(SYNTELOG_W_DENSITY_TABLE)

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
