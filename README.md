# Contact Information:
| Role          | Name          | GitHub                                                  | Email              |
|---------------|---------------|---------------------------------------------------------|--------------------|
| Project Lead: | Scott Teresi  | [Personal GitHub](https://github.com/sjteresi) | <teresisc@msu.edu> |
| PI:           | Patrick Edger | [Lab GitHub](https://github.com/EdgerLab)               | <edgerpat@msu.edu> |

# Genome Table and Files
| **Genome**   | **Regular Fasta**                                                                                         | **CDS FASTA**                                                                                                                          | **Gene Annotation**                                                                                                          | **TE Density** |
|--------------|-----------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------|----------------|
| Del Norte    | /mnt/research/strawberry/Reference_Genomes/Octoploids/Del_Norte/DelNorte_ragtag.scaffolds.fasta            | (TODO) Make via GFFRead: /mnt/research/edgerpat_lab/Scotty/Strawberry_Domestication/results/Del_Norte/Del_Norte_CDS_NewNames.fa.gz     | /mnt/research/strawberry/Reference_Genomes/Octoploids/Del_Norte/DelNorte-salt_maker_annotation.renamed.gff                   | TODO           |
| Royal Royce  | /mnt/research/strawberry/Reference_Genomes/Octoploids/Royal_Royce/farr1.fa                                | (TODO) Make via GFFRead: /mnt/research/edgerpat_lab/Scotty/Strawberry_Domestication/results/Royal_Royce/Royal_Royce_CDS_NewNames.fa.gz | /mnt/research/strawberry/Reference_Genomes/Octoploids/Royal_Royce/farr1.gene_models.gff                                      | TODO           |
| Hawaii 4     | /mnt/research/strawberry/Reference_Genomes/Diploids/Fragaria_vesca/Hawaii-4/H4_FASTA.fa                   | TODO) Make via GFFRead: /mnt/research/edgerpat_lab/Scotty/Strawberry_Domestication/results/H4/H4_CDS_NewNames.fa.gz                    | /mnt/research/strawberry/Reference_Genomes/Diploids/Fragaria_vesca/Hawaii-4/Fragaria_vesca_v4.0.a1.transcripts.gff3          | TODO           |
| Vesca 502    | /mnt/research/strawberry/Reference_Genomes/Diploids/Fragaria_vesca/502/502.ragtag.scaffolds.renamed.fasta | TODO) Make via GFFRead: /mnt/research/edgerpat_lab/Scotty/Strawberry_Domestication/results/Vesca_502/Vesca_502_CDS_NewNames.fa.gz      | /mnt/research/strawberry/Reference_Genomes/Diploids/Fragaria_vesca/502/502_maker_annotation.gff.v1.2.agat.gff                | Not creating   |
| Vesca 562    | /mnt/research/strawberry/Reference_Genomes/Diploids/Fragaria_vesca/562/562.ragtag.scaffolds.renamed.fasta | (TODO) Make via GFFRead: /mnt/research/edgerpat_lab/Scotty/Strawberry_Domestication/results/Vesca_562/Vesca_562_CDS_NewNames.fa.gz     | /mnt/research/strawberry/Reference_Genomes/Diploids/Fragaria_vesca/562/562_maker_annotation.gff.v1.2.agat.gff                | Not creating   |
| Vesca 2339   | /mnt/research/strawberry/Reference_Genomes/Diploids/Fragaria_vesca/1008/Renamed_2339/2339.final.fasta     | (TODO) Make via GFFRead: /mnt/research/edgerpat_lab/Scotty/Strawberry_Domestication/results/Vesca_2339/Vesca_2339_CDS_NewNames.fa.gz   | /mnt/research/strawberry/Reference_Genomes/Diploids/Fragaria_vesca/1008/Renamed_2339/2339_maker_annotation.gff.v1.2.agat.gff | Not creating   |
| F. iinumae   | /mnt/research/strawberry/Reference_Genomes/Diploids/F_iinumae/FII.fasta                                   | Pre-existing: /mnt/research/strawberry/Reference_Genomes/Diploids/F_iinumae/FII.cds                                                    | Not using                                                                                                                    | Not creating   |
| F. nipponica | /mnt/research/strawberry/Reference_Genomes/Diploids/F_nipponica/nipponica.ragtag.scaffolds.fasta          | Not creating or using                                                                                                                  | Not available                                                                                                                | Not creating   |
| F. viridis   | /mnt/research/strawberry/Reference_Genomes/Diploids/F_viridis/viridis.ragtag.scaffolds.fasta              | Not creating or using                                                                                                                  | Not available                                                                                                                | Not creating   |



# New
H4
https://genomevolution.org/r/1r9uy

DN
https://genomevolution.org/r/1r9uk


# EDTA Install:
```
- module purge
- module load Conda/3
- conda activate EDTA
- Install from YML file
- conda install bedtools
- conda install samtools
```

# EDTA Running:
I was having issues running EDTA v2.2.1 on RR (Royal Royce) genome.
It was crashing part way through because it was encountering a TE type not in some sort of internal EDTA whitelist.
However, v2.2.1 had a better panEDTA script, so it was desirable to use the latest version.
I was also having trouble with my annotation quality, because I did not originally provide CDS of each genome to EDTA.
This resulted in my genomes having a lot of genes annotated as TEs.
So my objective was to improve the annotation quality and make RR work with EDTA v2.2.1.

The solution I centered on was running RR through EDTA v2.1.1 to get the `RepeatMasker.out` file, and then run EDTA v2.2.1 on RR with the RepeatMasker as an additional argument.
This worked, and I was able to avoid the whitelist error.
I ran DN and H4 from the get-go with EDTA v2.2.1 (with CDSs) and then ran panEDTA on everything.

# Ortholog Analysis:
Generate the orthologs for the strawberries, and generate an ortholog table

	# 1. Filter the syntelogs that were generated from SynMap on CoGe
	# 2. Run the BLAST scripts so that we have a dataset supplemental to the CoGe results
	# 3. Filter the BLAST results
	# 4. Generate the master orthology table

