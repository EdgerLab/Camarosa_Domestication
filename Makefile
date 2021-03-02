# scripts for development
# __file__ Makefile
# __author__ Scott Teresi

ROOT_DIR:=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))
DATA_DIR := $(ROOT_DIR)/../Strawberry_Data
DOMESTICATION := $(ROOT_DIR)/../Domestication_Data

process_all: ## process the all genomes and generate the ortholog table
	python $(ROOT_DIR)/scripts/AT_Fragaria/generate_pairs.py $(DATA_DIR)/SynMap_Output/SynMapOutput_AT_UnmaskedCam.txt $(DATA_DIR)/BLAST_Data/AT_Unmasked_Camarosa.blast $(DATA_DIR)/SynMap_Output/SynMapOutput_MaskedVesca_UnmaskedCam.txt.txt $(DATA_DIR)/SynMap_Output/SynMapOutput_UnmaskedCam_MaskedDelNorte.txt

process_comparisons_w_s:
	python $(ROOT_DIR)/scripts/compare_density.py $(DOMESTICATION)/Strawberry_AT_Ortholog_Table_with_Supplementation.tsv

process_comparisons_no_s:
	python $(ROOT_DIR)/scripts/compare_density.py $(DOMESTICATION)/Strawberry_AT_Ortholog_Table_without_Supplementation.tsv
