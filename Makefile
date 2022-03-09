# scripts for development
# __file__ Makefile
# __author__ Scott Teresi

ROOT_DIR:=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))
DATA_DIR := $(ROOT_DIR)/data
GENOMES_DIR := $(DATA_DIR)/Genomes
RESULTS_DIR := $(ROOT_DIR)/results

setup:
	mkdir -p requirements doc results src data

# NOTE: Go through the directories of the strawberry genomes, 
# and create a CDS FASTA file with sequence IDs that
# conform to the sequence ID length requirements of EDTA. 
#
# NOTE must be done with the python environment activated
.ONESHELL:
.SILENT:
.SHELLFLAGS := -euo pipefail -c
generate_cds_and_fix_names:
	for dir in $$(find $(GENOMES_DIR) -mindepth 1 -type d)
	do
	gff_file=
	fasta_file=
	echo 
	#echo $$dir
		# Set the fasta file and gff file arguments for gffread
		for a_file in $$(find -L $$dir -mindepth 1 -type f)
			do
				if [ "$${a_file: -4}" == ".gff" ]
				then
					gff_file=$$a_file
				fi
				if [ "$${a_file: -3}" == ".fa" ]
				then
					fasta_file=$$a_file
				fi	
			done
		if [ -z $${gff_file} ]  && [ -z $${fasta_file} ]
			then
				echo "FAILURE"
			else
				# Define genome name from dir
				genome_name=$$(basename $$dir)
				echo "Genome Name: $$genome_name"
				# Define output dir from dir and path
				output_location=$$(realpath "$(RESULTS_DIR)/$$genome_name")
				mkdir -p $$output_location
				echo "Output Location: $$output_location"
				# Define output (CDS fasta) filename from path
				output_filename="$$output_location/$${genome_name}_CDS.fasta"
				echo "Output CDS Filename: $$output_filename"
				# Begin running GFFRead
				/mnt/research/edgerpat_lab/Scotty/gffread/gffread -x $$output_filename -g $$fasta_file $$gff_file
				# Fix gene names in that output. They are too long, and EDTA won't like it. Makes a new file with the New_Names suffix
				python $(ROOT_DIR)/src/fix_cds_names.py $$output_filename $$genome_name $$output_location
				# Compress the CDS file with the original names, to consere space
				gzip $$output_filename
		fi
		# Move to next dir, create whitespace
	done
