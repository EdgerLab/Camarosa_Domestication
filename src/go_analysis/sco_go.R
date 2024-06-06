# __author__ Scott Teresi
# Purpose:
	# Perform a GO enrichment on a list of Arabidopsis genes pre-identified as
	# being single-copy orthologs.

# Context:
	# TODO

suppressPackageStartupMessages(library(topGO))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))

# NOTE not sourcing the other TopGO script, just copying the function because I don't want to run stuff from the main environment.
# Apparently you can't just source functions, oh well, would have been better to not write stuff in main, but I'm not going to change it now.
# NOTE future, just write an Rscript that only has functions, don't mix main and functins if you ever want to re-use those functions
# source('/home/scott/Documents/Uni/Research/Projects/Strawberry_Domestication/src/go_analysis/TopGO.R')
#---------------------------------------------------------------------#
# Define input args and do simple set-up
# Check to see you have the args
args= commandArgs(trailingOnly=TRUE)

# TODO CHANGE THIS NAME
sco_arabidopsis_gene_table = read.delim(file=args[1], header=TRUE, sep='\t')


# Define your mapping (gene universe), MAGIC input arg 2
geneID2GO = readMappings(file=args[2], sep='\t') # NB global
master_genes = names(geneID2GO)  # NB global

# Must be tsv
outfile = args[3]

my_interesting_genes = sco_arabidopsis_gene_table$gene_id

# NOTE this fuction is copied from TopGO.R
run_topgo = function(master_genes, geneID2GO, my_interesting_genes, ontology_group){
	geneList = factor(as.integer(master_genes %in% my_interesting_genes))
	names(geneList) = master_genes # NB give it names

	GOdata = new('topGOdata',
		     ontology = ontology_group,
		     allGenes = geneList,
		     annot = annFUN.gene2GO,
		     gene2GO = geneID2GO)

	# NOTE get the number of genes in this interesting list
	# NOTE how many genes am I losing when I perform this step?
	# Do some genes not have the 'MF' ontology or something?
	# print(GOdata)  # NB shows the number of genes lost

	sig_genes = sigGenes(GOdata)  # the significant genes
	num_sig_genes = numSigGenes(GOdata) # NB the no. of signifcant genes,
	# not necessarily the nodes in the graph
	num_nodes = length(usedGO(object=GOdata))  # NB all nodes in the graph,
	# not necessarily the significant nodes
	# FUTURE would this be useful information to store?

	resultFisher_weight = runTest(GOdata, algorithm='weight01', statistic='fisher')
	#resultFisher_classic = runTest(GOdata, algorithm='classic', statistic='fisher')

	# NB from Nolan:
		# algorithm="classic" WILL NOT take GO hierarchy into account
	    	# ^ The limitation of using "classic" is that all genes
       		# annotated to a GO terms will be automatically annotated
       		# to its parents as well, therefore a GO term might look
       		# enriched just because its children are enriched
	    	# algorithm="weight01" WILL take GO hierarchy into account
		# These p-values have not been corrected for multiple testing

	# NB from Pat:
       		# After discussion with Pat, since I already have a cutoff applied
       		# to my modules, and the GO hieracrhy is useful information, I don't want
       		# to penalize the stats further, so I am going with classic and won't do
       		# an FDR.

	# list the top significant GO terms
	# This is for the weight option
    	allRes = GenTable(GOdata, classicFisher=resultFisher_weight, 
			  orderBy='classicFisher',
			  ranksOf='classicFisher',
			  topNodes=num_nodes,
			  numChar=1000) # NB avoid truncation of string
	# This is for classic option
    	# allRes = GenTable(GOdata, classicFisher=resultFisher_classic, 
	# 		  orderBy='classicFisher',
	# 		  ranksOf='classicFisher',
	# 		  topNodes=num_nodes,
	# 		  numChar=1000) # NB avoid truncation of string
			      
			      
	# MAGIC p-val threshold
       	pval_threshold = 0.05
	
	# NB apply the p-val threshold
	allRes_significant = allRes[which(allRes$classicFisher < pval_threshold),]


	# NB add a string to a new column that says 'yes' if we have
       	# more significant than expected
	allRes_significant$Overrepresented = ifelse(allRes_significant$Significant > allRes_significant$Expected, 'Y', 'N')

	return(allRes_significant)

}	

# Run the function
muh_table = run_topgo(master_genes, geneID2GO, my_interesting_genes, 'BP')
# Print to console where I am writing the file
print(paste('Writing to:', outfile))
write.table(muh_table, file=file.path(outfile), sep='\t', quote=FALSE, row.names=FALSE)
