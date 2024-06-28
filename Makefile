# __file__ Makefile
# __author__ Scott Teresi
# Purpose: Makefile for the strawberry domestication project.
# - This Makefile is organized into several sections, each of which is responsible for a different part of the analysis.
# - Sections are denoted by `--------------`
# - At the top of each section is a definition of various paths, to sort input and output files.
#-------------------------------------------------------------------#
# Main

setup:
	mkdir -p requirements doc results src data

SHELL=/bin/bash
ROOT_DIR:=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))
DATA_DIR := $(ROOT_DIR)/data
GENOMES_DIR := $(DATA_DIR)/Genomes
RESULTS_DIR := $(ROOT_DIR)/results
CODE_DIR := $(ROOT_DIR)/src

#-------------------------------------------------------------------#
# Prepare the gene and TE annotation inputs for TE Density

CLEANED_ANNOTATION_DIR := $(RESULTS_DIR)/cleaned_annotations
$(CLEANED_ANNOTATION_DIR):
	mkdir -p $@

# Define the file paths for the uncleaned gene annotations
H4_UNCLEAN_GENES := $(DATA_DIR)/Genomes/H4/H4_GeneAnnotation.gff
DN_UNCLEAN_GENES := $(DATA_DIR)/Genomes/Del_Norte/DN_GeneAnnotation.gff
RR_UNCLEAN_GENES := $(DATA_DIR)/Genomes/Royal_Royce/RR_GeneAnnotation.gff

# Define the file paths for the cleaned gene annotations
H4_CLEAN_GENES := $(CLEANED_ANNOTATION_DIR)/Cleaned_H4_GeneAnnotation.tsv
DN_CLEAN_GENES := $(CLEANED_ANNOTATION_DIR)/Cleaned_DN_GeneAnnotation.tsv
RR_CLEAN_GENES := $(CLEANED_ANNOTATION_DIR)/Cleaned_RR_GeneAnnotation.tsv

# Define the file paths for the uncleaned TE annotations
H4_UNCLEAN_TEs := $(RESULTS_DIR)/Pan_Annotation/H4_NewNames.fa.mod.EDTA.TEanno.gff3
DN_UNCLEAN_TEs := $(RESULTS_DIR)/Pan_Annotation/DN_NewNames.fa.mod.EDTA.TEanno.gff3
RR_UNCLEAN_TEs := $(RESULTS_DIR)/Pan_Annotation/RR_NewNames.fa.mod.EDTA.TEanno.gff3

# Define the file paths for the cleaned TE annotations
H4_CLEAN_TEs := $(CLEANED_ANNOTATION_DIR)/Cleaned_H4_NewNames.fa.mod.EDTA.TEanno.gff3
DN_CLEAN_TEs := $(CLEANED_ANNOTATION_DIR)/Cleaned_DN_NewNames.fa.mod.EDTA.TEanno.gff3
RR_CLEAN_TEs := $(CLEANED_ANNOTATION_DIR)/Cleaned_RR_NewNames.fa.mod.EDTA.TEanno.gff3

$(H4_CLEAN_GENES): $(H4_UNCLEAN_GENES) | $(CLEANED_ANNOTATION_DIR)
	@echo Filtering H4 strawberry genes into appropriate format for TE Density
	python $(CODE_DIR)/import_strawberry_gene_anno.py $< H4 $@

$(DN_CLEAN_GENES): $(DN_UNCLEAN_GENES) | $(CLEANED_ANNOTATION_DIR)
	@echo Filtering DN strawberry genes into appropriate format for TE Density
	python $(CODE_DIR)/import_strawberry_gene_anno.py $< DN $@

$(RR_CLEAN_GENES): $(RR_UNCLEAN_GENES) | $(CLEANED_ANNOTATION_DIR)
	@echo Filtering RR strawberry genes into appropriate format for TE Density
	python $(CODE_DIR)/import_strawberry_gene_anno.py $< RR $@

$(H4_CLEAN_TEs): $(H4_UNCLEAN_TEs) | $(CLEANED_ANNOTATION_DIR)
	@echo Filtering H4 strawberry TEs into appropriate format for TE Density
	python $(CODE_DIR)/import_strawberry_EDTA.py $< H4 $@

$(DN_CLEAN_TEs): $(DN_UNCLEAN_TEs) | $(CLEANED_ANNOTATION_DIR)
	@echo Filtering DN strawberry TEs into appropriate format for TE Density
	python $(CODE_DIR)/import_strawberry_EDTA.py $< DN $@

$(RR_CLEAN_TEs): $(RR_UNCLEAN_TEs) | $(CLEANED_ANNOTATION_DIR)
	@echo Filtering RR strawberry TEs into appropriate format for TE Density
	python $(CODE_DIR)/import_strawberry_EDTA.py $< RR $@

# Define a phony target to filter all gene annotations for convenience
.PHONY: filter_all_genes
filter_all_genes: $(H4_CLEAN_GENES) $(DN_CLEAN_GENES) $(RR_CLEAN_GENES)

# Define a phony target to filter all TE annotations for convenience
.PHONY: filter_all_transposons
filter_all_transposons: $(H4_CLEAN_TEs) $(DN_CLEAN_TEs) $(RR_CLEAN_TEs)

#-------------------------------------------------------------------#
# Prepare the syntelog and BLAST results to generate the ortholog table

# Define the file paths for the unclean syntelog files
RR_DN_SYNTELOG_REGULAR := $(DATA_DIR)/orthologs/RR_DN_SynMap.txt
RR_H4_SYNTELOG_REGULAR := $(DATA_DIR)/orthologs/RR_H4_SynMap.txt
# Define the file paths for the unclean BLAST results
RR_H4_BLAST_REGULAR := $(DATA_DIR)/orthologs/RR_H4.blast
RR_DN_BLAST_REGULAR := $(DATA_DIR)/orthologs/RR_DN.blast

# Define the file paths for the clean syntelog files
CLEANED_RR_DN_SYNTELOGS := $(DATA_DIR)/orthologs/filtered/Cleaned_RR_DN_Syntelogs.tsv
CLEANED_RR_H4_SYNTELOGS := $(DATA_DIR)/orthologs/filtered/Cleaned_RR_H4_Syntelogs.tsv
# Define the file paths for the cleaned BLAST results
RR_H4_BLAST_RENAMED := $(DATA_DIR)/orthologs/filtered/RR_H4_BLAST_renamed.tsv
RR_DN_BLAST_RENAMED := $(DATA_DIR)/orthologs/filtered/RR_DN_BLAST_renamed.tsv

# Define the file paths for processing the GO terms
UNCLEAN_GO_FILE := $(DATA_DIR)/GO/ATH_GOSLIM.txt
CLEANED_GO_FILE := $(RESULTS_DIR)/go_analysis/GO_ID_w_Term.tsv
TOP_GO_REFERENCE_FILE := $(RESULTS_DIR)/go_analysis/ArabidopsisGene_w_GO.tsv
GO_OUT_DIR := $(RESULTS_DIR)/go_analysis

# Define the file paths for the master orthology table --- this is the final product
STRAWBERRY_ORTHOLOG_TABLE := $(RESULTS_DIR)/orthologs/Strawberry_Arabidopsis_Ortholog_Table.tsv

# Define a target to create the directory if it doesn't exist
$(DATA_DIR)/orthologs/filtered $(RESULTS_DIR)/orthologs:
	mkdir -p $@

$(GO_OUT_DIR):
	mkdir -p $@


.PHONY: filter_RR_DN_syntelogs
filter_RR_DN_syntelogs: $(CLEANED_RR_DN_SYNTELOGS)

$(CLEANED_RR_DN_SYNTELOGS): $(RR_DN_SYNTELOG_REGULAR) | $(DATA_DIR)/orthologs/filtered
	python $(CODE_DIR)/orthologs/syntelogs.py $< DN $@

.PHONY: filter_RR_H4_syntelogs
filter_RR_H4_syntelogs: $(CLEANED_RR_H4_SYNTELOGS)

$(CLEANED_RR_H4_SYNTELOGS): $(RR_H4_SYNTELOG_REGULAR) | $(DATA_DIR)/orthologs/filtered
	python $(CODE_DIR)/orthologs/syntelogs.py $< H4 $@


$(RR_H4_BLAST_REGULAR):
	@echo  BLASTING RR and H4, must be run on the cluster
	sbatch $(CODE_DIR)/orthologs/rr_h4_blastall.sb

$(RR_DN_BLAST_REGULAR):
	@echo  BLASTING RR and DN, must be run on the cluster
	sbatch $(CODE_DIR)/orthologs/rr_dn_blastall.sb

# Define a target to generate the vanilla BLAST results, must be run on cluster.
.PHONY: generate_BLAST
generate_BLAST: $(RR_H4_BLAST_REGULAR) $(RR_DN_BLAST_REGULAR)


$(RR_H4_BLAST_RENAMED): $(RR_H4_BLAST_REGULAR) $(RR_CLEAN_GENES) $(H4_CLEAN_GENES) | $(DATA_DIR)/orthologs/filtered
	python $(CODE_DIR)/orthologs/reformat_RR_H4_BLAST_results.py $^ $@

# NOTE there is a gene renaming step in this particular script
$(RR_DN_BLAST_RENAMED): $(RR_DN_BLAST_REGULAR) $(DATA_DIR)/orthologs/DN_salt.translation $(RR_CLEAN_GENES) $(DN_CLEAN_GENES) | $(DATA_DIR)/orthologs/filtered
	python $(CODE_DIR)/orthologs/replace_and_reformat_DN_RR_BLAST_results.py $^ $@

# Define a target to rename the BLAST results
.PHONY: rename_RR_H4_BLAST
rename_RR_H4_BLAST: $(RR_H4_BLAST_RENAMED)

.PHONY: rename_RR_DN_BLAST
rename_RR_DN_BLAST: $(RR_DN_BLAST_RENAMED)

# Define a target to clean the GO file
.PHONY: filter_go_slim
filter_go_slim: $(CLEANED_GO_FILE) $(TOP_GO_REFERENCE_FILE)

$(CLEANED_GO_FILE) $(TOP_GO_REFERENCE_FILE): $(UNCLEAN_GO_FILE) | $(GO_OUT_DIR)
	python $(CODE_DIR)/go_analysis/generate_gene_w_GO_term.py $< $(GO_OUT_DIR)

.PHONY: clear_go_output
clear_go_output:
	rm -f $(CLEANED_GO_FILE)
	rm -f $(TOP_GO_REFERENCE_FILE)

# Define a target to generate the ortholog table from the BLAST and SynMap results
.PHONY: generate_ortholog_table
generate_ortholog_table: $(STRAWBERRY_ORTHOLOG_TABLE)

# TODO DEFINE a better target for the H4-At_Orthologs_March2022.tsv file
# Define the target for creating the master strawberry-arabidopsis ortholog table
$(STRAWBERRY_ORTHOLOG_TABLE): $(CLEANED_RR_H4_SYNTELOGS) $(RR_H4_BLAST_RENAMED) $(CLEANED_RR_DN_SYNTELOGS) $(RR_DN_BLAST_RENAMED) $(DATA_DIR)/orthologs/H4-At_Orthologs_March2022.tsv $(CLEANED_GO_FILE) $(DN_CLEAN_GENES) $(RR_CLEAN_GENES) $(H4_CLEAN_GENES) $(DATA_DIR)/orthologs/filtered
	python $(CODE_DIR)/orthologs/pan_orthology_table.py $^ $@

#-------------------------------------------------------------------#
# TE Density analysis of Syntelogs
# Define the file paths for the density data
DN_DENSITY_DIR := $(RESULTS_DIR)/density/DN/output_data
RR_DENSITY_DIR := $(RESULTS_DIR)/density/RR/output_data
H4_DENSITY_DIR := $(RESULTS_DIR)/density/H4/output_data
DENSITY_TABLE_DIR := $(RESULTS_DIR)/density_analysis/tables
SYNTELOG_PLOT_DIR := $(RESULTS_DIR)/density_analysis/figures

# Define a target to create the directory if it doesn't exist
$(DENSITY_TABLE_DIR):
	mkdir -p $@

$(SYNTELOG_PLOT_DIR):
	mkdir -p $@

# Define a phony target to create the syntelog density table for convenience
# NOTE this script takes a while....
.PHONY: generate_density_tables
generate_density_tables: $(STRAWBERRY_ORTHOLOG_TABLE) $(DN_CLEAN_GENES) $(RR_CLEAN_GENES) $(H4_CLEAN_GENES) $(DN_DENSITY_DIR) $(RR_DENSITY_DIR) $(H4_DENSITY_DIR) $(DENSITY_TABLE_DIR) | $(DENSITY_TABLE_DIR)
	python $(CODE_DIR)/syntelog_differences/parse_density_data.py $^

.PHONY: generate_syntelog_plots
generate_syntelog_plots: | $(SYNTELOG_PLOT_DIR)
	find $(DENSITY_TABLE_DIR) -type f -name '*minus*' | parallel python $(CODE_DIR)/syntelog_differences/bargraphs.py {} $(SYNTELOG_PLOT_DIR)

# TODO remove in future, DEV testing
.PHONY: test_syntelog_plots
test_syntelog_plots:
	python $(CODE_DIR)/syntelog_differences/bargraphs.py $(DENSITY_TABLE_DIR)/DN_minus_RR_Copia_1000_Downstream.tsv $(SYNTELOG_PLOT_DIR)



#-------------------------------------------------------------------#
# Define the file paths for dotplot output directory
DOTPLOT_OUT_DIR := $(RESULTS_DIR)/dotplot
DOTPLOT_CONFIG_FILE := $(CODE_DIR)/dotplot/dotplot_config_parameters.ini

# Define a target to create the directory if it doesn't exist
$(DOTPLOT_OUT_DIR):
	mkdir -p $@

.PHONY: generate_RR_dotplot
generate_RR_dotplot: $(RR_DENSITY_DIR) $(RR_CLEAN_GENES) $(DOTPLOT_OUT_DIR) $(STRAWBERRY_ORTHOLOG_TABLE) $(DOTPLOT_CONFIG_FILE) | $(DOTPLOT_OUT_DIR)
	python $(CODE_DIR)/dotplot/generate_dotplots.py $^ 'RR_(.*?).h5' 'RR' 

.PHONY: generate_DN_dotplot
generate_DN_dotplot: $(DN_DENSITY_DIR) $(DN_CLEAN_GENES) $(DOTPLOT_OUT_DIR) $(STRAWBERRY_ORTHOLOG_TABLE) $(DOTPLOT_CONFIG_FILE) | $(DOTPLOT_OUT_DIR)
	python $(CODE_DIR)/dotplot/generate_dotplots.py $^ 'DN_(.*?).h5' 'DN' 

#-------------------------------------------------------------------#
# Get the AED scores for the genes
# Define the file paths for the uncleaned gene annotations
H4_AED_SCORE := $(RESULTS_DIR)/AED/H4_AED.tsv
DN_AED_SCORE := $(RESULTS_DIR)/AED/DN_AED.tsv
RR_AED_SCORE := $(RESULTS_DIR)/AED/RR_AED.tsv
AED_SCORE_DIR := $(RESULTS_DIR)/AED

$(AED_SCORE_DIR):
	mkdir -p $@

# Define a phony target to generate all AED tables for convenience
.PHONY: generate_AED_score_tables
generate_AED_score_tables: $(H4_AED_SCORE) $(DN_AED_SCORE) $(RR_AED_SCORE)

$(H4_AED_SCORE): $(H4_UNCLEAN_GENES) | $(AED_SCORE_DIR)
	@echo Filtering H4 strawberry genes into appropriate format for TE Density
	python $(CODE_DIR)/extract_AED_score.py $< H4 $(AED_SCORE_DIR)/H4_AED_distribution.png $@

$(DN_AED_SCORE): $(DN_UNCLEAN_GENES) | $(AED_SCORE_DIR)
	@echo Filtering DN strawberry genes into appropriate format for TE Density
	python $(CODE_DIR)/extract_AED_score.py $< DN $(AED_SCORE_DIR)/DN_AED_distribution.png $@

$(RR_AED_SCORE): $(RR_UNCLEAN_GENES) | $(AED_SCORE_DIR)
	@echo Filtering RR strawberry genes into appropriate format for TE Density
	python $(CODE_DIR)/extract_AED_score.py $< RR $(AED_SCORE_DIR)/RR_AED_distribution.png $@

#-------------------------------------------------------------------#
# Locate super TE-dense genes

SUPER_DENSE_CUTOFF_TABLE_DIR := $(RESULTS_DIR)/density_analysis/cutoff_tables
GO_ENRICHMENT_DIR := $(GO_OUT_DIR)/enrichment
GO_UPSET_PLOT_DIR := $(RESULTS_DIR)/go_analysis/enrichment/upset_plots
GO_UPSET_NONSYNTENIC_PLOT_DIR := $(RESULTS_DIR)/go_analysis/enrichment/upset_plots/nonsyntenic
GO_UPSET_SYNTENIC_PLOT_DIR := $(RESULTS_DIR)/go_analysis/enrichment/upset_plots/syntenic
INTERMEDIATE_OUTPUT_DIFF_UNIQ := $(GO_ENRICHMENT_DIR)/difference_intermediate.tsv
FINAL_OUTPUT_DIFF_UNIQ := $(GO_ENRICHMENT_DIR)/Total_Diff_1K_unique_GO_terms.tsv
GO_UNIQ_ENRICHMENT_DIR := $(GO_ENRICHMENT_DIR)/unique_GO_terms

# NOTE this was made manually with a vim macro
GO_NONSYNTENIC_UPSET_KEY := $(CODE_DIR)/go_analysis/go_nonsyntenic_upset_key.tsv
GO_SYNTENIC_UPSET_KEY := $(CODE_DIR)/go_analysis/go_syntenic_upset_key.tsv

# Define a target to create the directory if it doesn't exist
$(SUPER_DENSE_CUTOFF_TABLE_DIR):
	mkdir -p $@

$(GO_ENRICHMENT_DIR):
	mkdir -p $@

$(GO_UNIQ_ENRICHMENT_DIR):
	mkdir -p $@

$(GO_UPSET_NONSYNTENIC_PLOT_DIR):
	mkdir -p $@

$(GO_UPSET_SYNTENIC_PLOT_DIR):
	mkdir -p $@

$(GO_UPSET_PLOT_DIR):
	mkdir -p $@
#-----------------------------------#

# NOTE takes a long time to run
# NOTE be careful interpreting the ortholog survival rate plots, as at large window sizes,
# the 95th percentile can have a low TE value that you need to beat, this especially matters for the
# uncommon TE superfamilies. For example the TE value to beat for the 95th percentile for the 10B window for 
# CACTA might be ~0.3 which isn't a lot of TE
# Total TE category might be most informative and easiest to interpret because of high cutoff value
# NOTE unused because I don't want all the different TE types
.PHONY: generate_super_dense_gene_tables_single_genome
generate_super_dense_gene_tables_single_genome: $(STRAWBERRY_ORTHOLOG_TABLE) $(SUPER_DENSE_CUTOFF_TABLE_DIR)
	# Lol I made the manifest file with Vim in like 1 minute, I'm not going to automate that
	mkdir -p $(SUPER_DENSE_CUTOFF_TABLE_DIR)/no_Arabidopsis
	mkdir -p $(SUPER_DENSE_CUTOFF_TABLE_DIR)/ortholog_analysis
	cat $(CODE_DIR)/go_analysis/cutoff_single_genome_manifest.tsv | parallel -a - -C '\t' python $(CODE_DIR)/go_analysis/find_abnormal_genes.py $(DENSITY_TABLE_DIR)/{1} $(RESULTS_DIR)/AED/{2} 95 5 $^

prefixes := RR DN H4
types := Total_TE_Density_5000_Upstream TIR_5000_Upstream LTR_5000_Upstream
# Loop over a specific set
super_dense_single_genome_specific: $(STRAWBERRY_ORTHOLOG_TABLE) $(SUPER_DENSE_CUTOFF_TABLE_DIR)
	mkdir -p $(SUPER_DENSE_CUTOFF_TABLE_DIR)/no_Arabidopsis
	mkdir -p $(SUPER_DENSE_CUTOFF_TABLE_DIR)/ortholog_analysis
	$(foreach prefix,$(prefixes), \
		$(foreach type,$(types), \
			python $(CODE_DIR)/go_analysis/find_abnormal_genes.py \
				$(DENSITY_TABLE_DIR)/$(prefix)_$(type).tsv \
				$(RESULTS_DIR)/AED/$(prefix)_AED.tsv \
				95 \
				5 \
				$^ ; \
		) \
	)


# NOTE unused because I don't want all the different TE types,
# Use super_dense_two_genome_specific instead
.PHONY: deprecated_generate_super_dense_gene_tables_differing_syntelogs
deprecated_generate_super_dense_gene_tables_differing_syntelogs: $(SUPER_DENSE_CUTOFF_TABLE_DIR)
	# Lol I made the manifest file with Vim in like 1 minute, I'm not going to automate that
	mkdir -p $(SUPER_DENSE_CUTOFF_TABLE_DIR)/no_Arabidopsis
	mkdir -p $(SUPER_DENSE_CUTOFF_TABLE_DIR)/ortholog_analysis
	cat $(CODE_DIR)/go_analysis/cutoff_differing_syntelogs_manifest.tsv | parallel -a - -C '\t' python $(CODE_DIR)/go_analysis/find_differing_syntelogs.py $(DENSITY_TABLE_DIR)/{1} $(DN_AED_SCORE) $(RR_AED_SCORE) $^

double_prefixes := DN_minus_RR
.PHONY: super_dense_two_genome_specific
super_dense_two_genome_specific: $(DN_AED_SCORE) $(RR_AED_SCORE) $(SUPER_DENSE_CUTOFF_TABLE_DIR)
	mkdir -p $(SUPER_DENSE_CUTOFF_TABLE_DIR)/no_Arabidopsis
	mkdir -p $(SUPER_DENSE_CUTOFF_TABLE_DIR)/ortholog_analysis
	$(foreach prefix,$(double_prefixes), \
		$(foreach type,$(types), \
			python $(CODE_DIR)/go_analysis/find_differing_syntelogs.py \
	       			$(DENSITY_TABLE_DIR)/$(prefix)_$(type).tsv $^ ; \
		) \
	)

# Define a target to run TopGO on the super dense gene output files
.PHONY: generate_go_enrichments
generate_go_enrichments: $(TOP_GO_REFERENCE_FILE) | $(GO_ENRICHMENT_DIR)
	find $(SUPER_DENSE_CUTOFF_TABLE_DIR)/ -maxdepth 1 -type f -name '*.tsv' | parallel Rscript $(CODE_DIR)/go_analysis/TopGO.R {} $(TOP_GO_REFERENCE_FILE) $(GO_ENRICHMENT_DIR)


# Creates text files that are more legible
# TODO this is untested
.ONESHELL: extract_unique_enriched_GO_terms
.PHONY: extract_unique_enriched_GO_terms
extract_unique_enriched_GO_terms: | $(GO_UNIQ_ENRICHMENT_DIR)
	for file in $(GO_ENRICHMENT_DIR)/*.tsv; do \
		cut -d$$'\t' -f 1,2 "$$file" | sort | uniq | sed \$$d > $(GO_UNIQ_ENRICHMENT_DIR)/Unique_$$(basename $$file); \
	done

# TODO this is untested
.ONESHELL: extract_diff_total_unique_enriched_GO_terms
.PHONY: extract_diff_total_unique_enriched_GO_terms
extract_diff_total_unique_enriched_GO_terms:
	cut -d$$'\t' -f 1,2 $(GO_ENRICHMENT_DIR)/Overrepresented_Difference_Total_TE_Density_1000_Upstream_Biased_Towards_DN_95_density_percentile.tsv | sort | uniq | sed \$$d  > $(FINAL_OUTPUT_DIFF_UNIQ)
	#rm -f $(INTERMEDIATE_OUTPUT_DIFF_UNIQ)
	#cut -d$$'\t' -f 1,2 $(GO_ENRICHMENT_DIR)/Overrepresented_Difference_Total_TE_Density_1000_Upstream_Biased_Towards_DN_95_density_percentile.tsv >> $(INTERMEDIATE_OUTPUT_DIFF_UNIQ)
	#cut -d$$'\t' -f 1,2 $(INTERMEDIATE_OUTPUT_DIFF_UNIQ) | sort | uniq > $(FINAL_OUTPUT_DIFF_UNIQ)



.PHONY: generate_dn_rr_nonsyntenic_upset_plot
generate_dn_rr_nonsyntenic_upset_plot: | $(GO_UPSET_NONSYNTENIC_PLOT_DIR)
	cat $(GO_NONSYNTENIC_UPSET_KEY) | while IFS=$$'\t' read -r dn_data rr_data; do \
		python $(CODE_DIR)/go_analysis/dn_vs_rr_upset.py $(GO_ENRICHMENT_DIR)/$$dn_data $(GO_ENRICHMENT_DIR)/$$rr_data $(GO_UPSET_NONSYNTENIC_PLOT_DIR); \
	done

.PHONY: generate_dn_rr_syntenic_upset_plot
generate_dn_rr_syntenic_upset_plot: | $(GO_UPSET_SYNTENIC_PLOT_DIR)
	cat $(GO_SYNTENIC_UPSET_KEY) | while IFS=$$'\t' read -r dn_data rr_data; do \
		python $(CODE_DIR)/go_analysis/dn_vs_rr_upset.py $(GO_ENRICHMENT_DIR)/$$dn_data $(GO_ENRICHMENT_DIR)/$$rr_data $(GO_UPSET_SYNTENIC_PLOT_DIR) --syntelog; \
	done
	
	
.PHONY: test_nonsyntenic_upset
test_nonsyntenic_upset: | $(GO_UPSET_NONSYNTENIC_PLOT_DIR)
	python $(CODE_DIR)/go_analysis/dn_vs_rr_upset.py $(GO_ENRICHMENT_DIR)/Overrepresented_DN_Total_TE_Density_5000_Upstream_Upper_95_density_percentile.tsv $(GO_ENRICHMENT_DIR)/Overrepresented_RR_Total_TE_Density_5000_Upstream_Upper_95_density_percentile.tsv $(GO_UPSET_NONSYNTENIC_PLOT_DIR)

# grep -nr "flower development" Overrepresented_RR_Total_TE_Density_5000_Upstream_Upper_95_density_percentile.tsv | awk -F'\t' '{print $1, $2, $8}'
# 

.PHONY: test_syntenic
test_syntenic: | $(GO_UPSET_SYNTENIC_PLOT_DIR)
	python $(CODE_DIR)/go_analysis/dn_vs_rr_upset.py $(GO_ENRICHMENT_DIR)/Overrepresented_Difference_Total_TE_Density_5000_Upstream_Syntelogs_Biased_Towards_DN.tsv $(GO_ENRICHMENT_DIR)/Overrepresented_Difference_Total_TE_Density_5000_Upstream_Syntelogs_Biased_Towards_RR.tsv $(GO_UPSET_SYNTENIC_PLOT_DIR) --syntelog

# Define a target to generate an UpSet plot for the GO enrichment results
# TODO verify one more time that the files are being provided in the correct order
# NOTE THIS DOES NOT WORK FOR THE 'DIFFERENCE' FILES
# TODO this is untested
# # TODO MAKE SURE THE ORDER OF INPUTS IS RIGHT
.PHONY: generate_all_upset_plot
generate_all_upset_plot: | $(GO_UPSET_PLOT_DIR)
	find $(GO_ENRICHMENT_DIR)/ -maxdepth 1 -type f -not -name "*Difference*" | sort | xargs -n3 sh -c 'python $(CODE_DIR)/go_analysis/upset_plot.py $$1 $$2 $$3 $(GO_UPSET_PLOT_DIR)' sh

# NOTE deprecated
.PHONY: test_upset
test_upset:
	mkdir -p $(GO_UPSET_PLOT_DIR)
	python $(CODE_DIR)/go_analysis/upset_plot.py \
		$(GO_ENRICHMENT_DIR)/Overrepresented_DN_LTR_5000_Upstream_Upper_95_density_percentile.tsv \
		$(GO_ENRICHMENT_DIR)/Overrepresented_RR_LTR_5000_Upstream_Upper_95_density_percentile.tsv \
		$(GO_ENRICHMENT_DIR)/Overrepresented_H4_LTR_5000_Upstream_Upper_95_density_percentile.tsv \
		$(GO_UPSET_PLOT_DIR)
#-------------------------------------------------------------------#
# TODO get an outputfilepath variable name
#
.PHONY: generate_gene_distance_plots
generate_gene_distance_plots:
	python $(CODE_DIR)/gene_distances/gene_distances.py $(DN_CLEAN_GENES) DN $(RESULTS_DIR)/gene_distances
	python $(CODE_DIR)/gene_distances/gene_distances.py $(RR_CLEAN_GENES) RR $(RESULTS_DIR)/gene_distances

#-------------------------------------------------------------------#
# Single Copy Ortholog Analaysis
# Do genes that are single copy orthologs have lower TE densities than other genes?
# What are the GO terms associated with single copy ortholog genes?

SCO_TABLE := $(DATA_DIR)/orthologs/single_copy_orthologs/sd01.tsv
SCO_GO_OUTPUT := $(GO_ENRICHMENT_DIR)/Overrepresented_SCO_GO_terms.tsv

.PHONY: single_copy_orthologs
single_copy_orthologs: $(STRAWBERRY_ORTHOLOG_TABLE) $(SCO_TABLE) $(RESULTS_DIR)
	python src/orthologs/single_copy_orthologs.py $(DENSITY_TABLE_DIR)/H4_Total_TE_Density_5000_Upstream.tsv $(DENSITY_TABLE_DIR)/RR_Total_TE_Density_5000_Upstream.tsv $^ 


# Define a target to run TopGO on the super dense gene output files
$(SCO_GO_OUTPUT): $(SCO_TABLE) $(TOP_GO_REFERENCE_FILE) 
	@echo "Running TopGO on single copy orthologs"
	Rscript $(CODE_DIR)/go_analysis/sco_go.R $(SCO_TABLE) $(TOP_GO_REFERENCE_FILE) $@

.PHONY: generate_sco_go_enrichments
generate_sco_go_enrichments: $(SCO_GO_OUTPUT)

#-------------------------------------------------------------------#
# Intersect the selection/domestication sweeps with my list of TE-dense genes,
# these genes are the output of TopGO.

# Define a file path for the domestication sweeps data from the publication
SWEEP_TABLE := $(DATA_DIR)/domestication_sweeps.csv
INTERSECTION_MANIFEST := $(CODE_DIR)/go_analysis/intersect_sweeps_manifest.tsv


.PHONY: intersect_sweeps_with_TE_dense_genes
intersect_sweeps_with_TE_dense_genes: 
	mkdir -p $(GO_ENRICHMENT_DIR)/sweep_intersections
	cat $(INTERSECTION_MANIFEST) | parallel -a - -C '\t' python $(CODE_DIR)/intersect_sweeps_w_enriched_terms.py $(GO_ENRICHMENT_DIR)/{1} ${GO_ENRICHMENT_DIR}/sweep_intersections/{2} $(SWEEP_TABLE) 

.PHONY: dev_intersect_sweeps_with_TE_dense_genes
dev_intersect_sweeps_with_TE_dense_genes: 
	python $(CODE_DIR)/intersect_sweeps_w_enriched_terms.py $(SWEEP_TABLE) $(GO_ENRICHMENT_DIR)/Overrepresented_RR_Total_TE_Density_5000_Upstream_Upper_95_density_percentile.tsv test_intersect.tsv
	#python $(CODE_DIR)/intersect_sweeps_w_enriched_terms.py $(SWEEP_TABLE) $(GO_ENRICHMENT_DIR)/Overrepresented_Difference_Total_TE_Density_5000_Upstream_Syntelogs_Biased_Towards_RR.tsv


#-------------------------------------------------------------------#
# TODO unfishied
# TODO this may need to be changed when I start looking at the other genomes
.PHONY: filter_RR_expression
filter_RR_expression:
	@echo Filtering RoyalRoyce expression data
	python $(CODE_DIR)/expression_data.py $(DATA_DIR)/Genomes/Royal_Royce/RoyalRoyce_Cold_Warm_count_matrix.csv $(RESULTS_DIR)

# NOTE DEPRECATED 1/23/2024
# NOTE this translates the Royal Royce CDS FASTA
.PHONY: translate_CDS_to_protein
translate_CDS_to_protein:
	mkdir -p $(DATA_DIR)/Royal_Royce/RR_Proteins
	@echo Converting Royal Royce files
	python $(CODE_DIR)/orthologs/translate_cds_fasta_to_protein.py $(DATA_DIR)/Genomes/Royal_Royce/RR_CDS.fa $(DATA_DIR)/orthologs/RR_CDS_as_Proteins.fa
	@echo Converting H4 files
	python $(CODE_DIR)/orthologs/translate_cds_fasta_to_protein.py $(DATA_DIR)/Genomes/H4/H4_CDS.fa $(DATA_DIR)/orthologs/H4_CDS_as_Proteins.fa

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


