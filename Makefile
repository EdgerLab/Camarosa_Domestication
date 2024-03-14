# scripts for development
# __file__ Makefile
# __author__ Scott Teresi
SHELL=/bin/bash

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
# Prepare the gene annotation inputs for TE Density

# Define the file paths for the uncleaned gene annotations
H4_UNCLEAN_GENES := $(DATA_DIR)/Genomes/H4/H4_GeneAnnotation.gff
DN_UNCLEAN_GENES := $(DATA_DIR)/Genomes/Del_Norte/DN_GeneAnnotation.gff
RR_UNCLEAN_GENES := $(DATA_DIR)/Genomes/Royal_Royce/RR_GeneAnnotation.gff
FVI_UNCLEAN_GENES := $(DATA_DIR)/Genomes/F_viridis/FVI_GeneAnnotation.gff
FNI_UNCLEAN_GENES := $(DATA_DIR)/Genomes/F_nipponica/FNI_GeneAnnotation.gff
FII_UNCLEAN_GENES := $(DATA_DIR)/Genomes/F_iinumae/FII_GeneAnnotation.gff

# Define the file paths for the cleaned gene annotations
H4_CLEAN_GENES := $(RESULTS_DIR)/cleaned_annotations/Cleaned_H4_GeneAnnotation.tsv
DN_CLEAN_GENES := $(RESULTS_DIR)/cleaned_annotations/Cleaned_DN_GeneAnnotation.tsv
RR_CLEAN_GENES := $(RESULTS_DIR)/cleaned_annotations/Cleaned_RR_GeneAnnotation.tsv
FVI_CLEAN_GENES := $(RESULTS_DIR)/cleaned_annotations/Cleaned_FVI_GeneAnnotation.tsv
FNI_CLEAN_GENES := $(RESULTS_DIR)/cleaned_annotations/Cleaned_FNI_GeneAnnotation.tsv
FII_CLEAN_GENES := $(RESULTS_DIR)/cleaned_annotations/Cleaned_FII_GeneAnnotation.tsv

# Define a target to create the output directory if it doesn't exist
$(RESULTS_DIR)/cleaned_annotations:
	mkdir -p $@
#-----------------------------------#

# Define a phony target to filter all gene annotations for convenience
.PHONY: filter_all_genes
filter_all_genes: $(H4_CLEAN_GENES) $(DN_CLEAN_GENES) $(RR_CLEAN_GENES) $(FII_CLEAN_GENES) $(FNI_CLEAN_GENES) $(FVI_CLEAN_GENES)

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

# Define a phony target to run a singular gene filtering for convenience
# FUTURE if you want more individually type the phonies out
# TODO DEV REMOVE FUTURE
.PHONY: filter_H4_genes
filter_H4_genes: $(DN_CLEAN_GENES)


# TODO CLEANUP
#-----------------------------------------------------

UNCLEAN_GO_FILE := $(DATA_DIR)/GO/ATH_GOSLIM.txt
CLEANED_GO_FILE := $(RESULTS_DIR)/go_analysis/GO_ID_w_Term.tsv
GO_OUT_DIR := $(RESULTS_DIR)/go_analysis
$(GO_OUT_DIR):
	mkdir -p $@

# Define a target to clean the GO file
.PHONY: filter_go_slim
filter_go_slim: $(CLEANED_GO_FILE) $(TOP_GO_REFERENCE_FILE)

$(CLEANED_GO_FILE) $(TOP_GO_REFERENCE_FILE): $(UNCLEAN_GO_FILE) | $(GO_OUT_DIR)
	python $(ROOT_DIR)/src/go_analysis/generate_gene_w_GO_term.py $< $(GO_OUT_DIR)

.PHONY: clear_go_output
clear_go_output:
	rm -f $(CLEANED_GO_FILE)
	rm -f $(TOP_GO_REFERENCE_FILE)

#-------------------------------------------------------------------#
# Orthology Analysis:
	# 1. Filter the syntelogs that were generated from SynMap on CoGe
	# 2. Run the BLAST scripts so that we have a dataset supplemental to the CoGe results
	# 3. Filter the BLAST results
	# 4. Generate the master orthology table

# Define the file paths for the regular and cleaned syntelogs files
RR_DN_SYNTELOG_REGULAR := $(DATA_DIR)/orthologs/RR_DN_SynMap.txt
RR_H4_SYNTELOG_REGULAR := $(DATA_DIR)/orthologs/RR_H4_SynMap.txt
CLEANED_RR_DN_SYNTELOGS := $(DATA_DIR)/orthologs/filtered/Cleaned_RR_DN_Syntelogs.tsv
CLEANED_RR_H4_SYNTELOGS := $(DATA_DIR)/orthologs/filtered/Cleaned_RR_H4_Syntelogs.tsv

# Define the file paths for the renamed BLAST results
RR_H4_BLAST_REGULAR := $(DATA_DIR)/orthologs/RR_H4.blast
RR_DN_BLAST_REGULAR := $(DATA_DIR)/orthologs/RR_DN.blast
RR_H4_BLAST_RENAMED := $(DATA_DIR)/orthologs/filtered/RR_H4_BLAST_renamed.tsv
RR_DN_BLAST_RENAMED := $(DATA_DIR)/orthologs/filtered/RR_DN_BLAST_renamed.tsv

# Define the file paths for the master orthology table that is the final product
STRAWBERRY_ORTHOLOG_TABLE := $(RESULTS_DIR)/orthologs/Strawberry_Arabidopsis_Ortholog_Table.tsv

# Define a target to create the directory if it doesn't exist
$(DATA_DIR)/orthologs/filtered $(RESULTS_DIR)/orthologs:
	mkdir -p $@
#-----------------------------------#

.PHONY: filter_RR_DN_syntelogs
filter_RR_DN_syntelogs: $(CLEANED_RR_DN_SYNTELOGS)

.PHONY: filter_RR_H4_syntelogs
filter_RR_H4_syntelogs: $(CLEANED_RR_H4_SYNTELOGS)

$(CLEANED_RR_DN_SYNTELOGS): $(RR_DN_SYNTELOG_REGULAR) | $(DATA_DIR)/orthologs/filtered
	python $(ROOT_DIR)/src/orthologs/syntelogs.py $< DN $@

$(CLEANED_RR_H4_SYNTELOGS): $(RR_H4_SYNTELOG_REGULAR) | $(DATA_DIR)/orthologs/filtered
	python $(ROOT_DIR)/src/orthologs/syntelogs.py $< H4 $@

# Define a target to generate the vanilla BLAST results, must be run on cluster.
.PHONY: generate_BLAST
generate_BLAST: $(RR_H4_BLAST_REGULAR) $(RR_DN_BLAST_REGULAR)

$(RR_H4_BLAST_REGULAR):
	sbatch $(ROOT_DIR)/src/orthologs/rr_h4_blastall.sb

$(RR_DN_BLAST_REGULAR):
	sbatch $(ROOT_DIR)/src/orthologs/rr_dn_blastall.sb

# Define a target to rename the BLAST results
.PHONY: rename_RR_H4_BLAST
rename_RR_H4_BLAST: $(RR_H4_BLAST_RENAMED)

.PHONY: rename_RR_DN_BLAST
rename_RR_DN_BLAST: $(RR_DN_BLAST_RENAMED)

$(RR_H4_BLAST_RENAMED): $(RR_H4_BLAST_REGULAR) $(RR_CLEAN_GENES) $(H4_CLEAN_GENES) | $(DATA_DIR)/orthologs/filtered
	python $(ROOT_DIR)/src/orthologs/reformat_RR_H4_BLAST_results.py $^ $@

# NOTE there is a gene renaming step in this particular script
$(RR_DN_BLAST_RENAMED): $(RR_DN_BLAST_REGULAR) $(DATA_DIR)/orthologs/DN_salt.translation $(RR_CLEAN_GENES) $(DN_CLEAN_GENES) | $(DATA_DIR)/orthologs/filtered
	python $(ROOT_DIR)/src/orthologs/replace_and_reformat_DN_RR_BLAST_results.py $^ $@


# Define a target to generate the ortholog table from the BLAST and SynMap results
.PHONY: generate_ortholog_table
generate_ortholog_table: $(STRAWBERRY_ORTHOLOG_TABLE)

# TODO
# DEFINE A better target for the H4-At_Orthologs_March2022.tsv file
# Define the target for creating the master strawberry-arabidopsis ortholog table
$(STRAWBERRY_ORTHOLOG_TABLE): $(CLEANED_RR_H4_SYNTELOGS) $(RR_H4_BLAST_RENAMED) $(CLEANED_RR_DN_SYNTELOGS) $(RR_DN_BLAST_RENAMED) $(DATA_DIR)/orthologs/H4-At_Orthologs_March2022.tsv $(CLEANED_GO_FILE) $(DN_CLEAN_GENES) $(RR_CLEAN_GENES) $(H4_CLEAN_GENES) $(DATA_DIR)/orthologs/filtered
	python $(ROOT_DIR)/src/orthologs/pan_orthology_table.py $^ $@

#-------------------------------------------------------------------#
# TE Density analysis of Syntelogs
# Define the file paths for the density data
DN_DENSITY_DIR := $(RESULTS_DIR)/density/DN
RR_DENSITY_DIR := $(RESULTS_DIR)/density/RR
H4_DENSITY_DIR := $(RESULTS_DIR)/density/H4

# Define the file path for the syntelog density table
# TODO this is actually never made because the code got refactored into multiple tables

# Define a target to create the directory if it doesn't exist
$(RESULTS_DIR)/density_analysis:
	mkdir -p $@
#-----------------------------------#

# Define a phony target to create the syntelog density table for convenience
# NOTE this script takes a while....
# FUTURE efforts could be made to take the loop out of the parse_density_data.py script and make it so that
# the script can be run in parallel
.PHONY: create_sytelog_density_tables
create_sytelog_density_tables: $(STRAWBERRY_ORTHOLOG_TABLE) $(DN_CLEAN_GENES) $(RR_CLEAN_GENES) $(H4_CLEAN_GENES) $(DN_DENSITY_DIR) $(RR_DENSITY_DIR) $(H4_DENSITY_DIR) $(RESULTS_DIR)/density_analysis | $(RESULTS_DIR)/density_analysis
	python $(ROOT_DIR)/src/syntelog_differences/parse_density_data.py $^

#.PHONY: test_sytelog_density_tables



# TODO add the mkdir as a dependency not in the command
.PHONY: graph_syntelog_plots
graph_syntelog_plots:
	mkdir -p $(RESULTS_DIR)/density_analysis/figures
	ls $(RESULTS_DIR)/density_analysis/*.tsv | parallel python $(ROOT_DIR)/src/syntelog_differences/bargraphs.py {} $(RESULTS_DIR)/density_analysis/figures

# TODO remove in future, DEV testing
.PHONY: test_syntelog_plots
test_syntelog_plots:
	python $(ROOT_DIR)/src/syntelog_differences/bargraphs.py $(RESULTS_DIR)/density_analysis/DN_minus_RR_Copia_1000_Downstream.tsv $(RESULTS_DIR)/density_analysis/figures

#-------------------------------------------------------------------#
# Locate super TE-dense genes
# GO Analysis: Filter the Arabidopsis GO SLIM file
# GO Analysis: Run TopGO on the TE-dense genes
# Generate an UpSet plot for the GO enrichment results, showing the shared GO terms

TOP_GO_REFERENCE_FILE := $(RESULTS_DIR)/go_analysis/ArabidopsisGene_w_GO.tsv

CUTOFF_TABLES_DIR := $(RESULTS_DIR)/density_analysis/cutoff_tables
GO_ENRICHMENT_DIR := $(GO_OUT_DIR)/enrichment
GO_UPSET_PLOT_DIR := $(RESULTS_DIR)/go_analysis/enrichment/upset_plots

# Define a target to create the directory if it doesn't exist
$(CUTOFF_TABLES_DIR):
	mkdir -p $@

$(GO_ENRICHMENT_DIR):
	mkdir -p $@

$(GO_UPSET_PLOT_DIR):
	mkdir -p $@
#-----------------------------------#

# Define a target to create tables of the top X% of TE-dense genes, and the DN - RR difference
.PHONY: find_super_dense_genes
find_super_dense_genes: | $(CUTOFF_TABLES_DIR)
	ls $(RESULTS_DIR)/density_analysis/*.tsv | parallel python $(ROOT_DIR)/src/go_analysis/find_abnormal_genes.py {} 95 5 $(CUTOFF_TABLES_DIR)

# TODO remove in future, DEV testing
.PHONY: test_cutoff
test_cutoff:
	python $(ROOT_DIR)/src/go_analysis/find_abnormal_genes.py $(RESULTS_DIR)/density_analysis/DN_minus_RR_Copia_1000_Downstream.tsv 99 1 $(CUTOFF_TABLES_DIR)



# Define a target to run TopGO on the super dense gene output files
.PHONY: run_topgo
run_topgo: $(TOP_GO_REFERENCE_FILE) | $(GO_ENRICHMENT_DIR)
	rm $(GO_ENRICHMENT_DIR)/*
	ls $(CUTOFF_TABLES_DIR)/*.tsv | parallel Rscript $(ROOT_DIR)/src/go_analysis/TopGO.R {} $(TOP_GO_REFERENCE_FILE) $(GO_ENRICHMENT_DIR)

# to look at output 
# grep -nr "stress" . --exclude="*Tc1*" | cut -f1,20-24
#
#

# Define a target to generate an UpSet plot for the GO enrichment results
# TODO verify one more time that the files are being provided in the correct order
.PHONY: upset_plot
upset_plot: | $(GO_UPSET_PLOT_DIR)
	find $(GO_ENRICHMENT_DIR)/ -maxdepth 1 -type f | sort | xargs -n2 sh -c 'python $(ROOT_DIR)/src/go_analysis/upset_plot.py $$2 $$3 $$1' sh $(GO_UPSET_PLOT_DIR)


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


