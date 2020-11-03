# scripts for development
# __file__ Makefile
# __author__ Scott Teresi

ROOT_DIR:=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))
DEV_DATA := $(ROOT_DIR)/../Strawberry_Data/CoGe_Blast_Data

DEV_CAM_SYN := $(DEV_DATA)/AT_Camarosa_Synmap_Raw.txt
DEV_CAM_NAME:= "Camarosa"

DEV_H4_SYN := $(DEV_DATA)/AT_H4_Synmap_Raw.txt
DEV_H4_NAME:= "H4"

# TODO change to Del Norte
# DEV_IINUMAE_SYN := $(DEV_DATA)/AT_Iinumae_Synmap_Raw.txt
# DEV_IINUMAE_NAME:= "Iinumae"


process_cam: ## process the camarosa genome
	python $(ROOT_DIR)/scripts/AT_Fragaria/generate_pairs.py ../Strawberry_Data/SynMap_Output/SynMapOutput_AT_UnmaskedCam.txt ../Strawberry_Data/BLAST_Data/AT_Unmasked_Camarosa.blast ../Strawberry_Data/SynMap_Output/SynMapOutput_MaskedVesca_UnmaskedCam.txt.txt

