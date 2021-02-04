# Goal:
Develop requisite data structures and code for comparing TE Density data between strawberries.

# Contact Information:
| Role          | Name          | GitHub                                                  | Email              |
|---------------|---------------|---------------------------------------------------------|--------------------|
| Project Lead: | Scott Teresi  | [Personal GitHub](https://github.com/huckleberry-hound) | <teresisc@msu.edu> |
| PI:           | Patrick Edger | [Lab GitHub](https://github.com/EdgerLab)               | <edgerpat@msu.edu> |

# Genomes Utilized:
| Regular CoGe ID                   | Masked ID                                    | Peptide Source | Shorthand Name |
|-----------------------------------|----------------------------------------------|----------------|-------------|
| Arabidopsis thaliana Col-0 (id 1) | CNS PL.20 Masked repeats 50X (v10, id 16746) | CoGe Derived | AT |
| Fragaria x ananass subsp. ananassa (id 41115) | unmasked (v1, id 35789) | CoGe Derived | Camarosa |
| Fragaria vesca subsp. vesca (vv4, id 57743) | NCBI WindowMasker Hard | Not Used | H4 |
| Fragaria chiloensis subsp. chiloensis (Del Norte) (v1, id60430) | Hard Masked | Not Used | Del Norte | 

# Identification of Syntelogs and Homologs:
First I need to create a list of the syntelogs and homologs between the strawberry species and with *Arabidopsis* so that I may later collect gene function information. The output of this portion of the project will primarily yield a `.tsv` file that provides something similar to the following:

| Arabidopsis Gene | Camarosa Gene | H4 Gene      | Del Norte Gene |
|------------------|---------------|--------------|----------------|
| ATG_1            | Dummy_gene_1  | Dummy_gene_1 | Dummy_gene_1   |
| ATG_2            | Dummy_gene_2  | Dummy_gene_2 | Dummy_gene_2   |


## Workflow:
* Generate a list of syntelogs and homologs from Camarosa to *Arabidopsis*.
* Generate a list of syntelogs for each non-Camarosa strawberry genome to Camarosa.
* Generate a table containing the above information in an easily searchable format for later use in conjunction with TE Density data.


## Running SynMap:
This section describes the methods to run [SynMap](https://genomevolution.org/CoGe/SynMap.pl) on CoGe. I ran SynMap with mostly [default options](https://genomevolution.org/wiki/index.php/SynMap), I did change one option: under *Merge Syntenic Blocks* I set it to `Quota Align Merge`. Here is the [link](https://genomevolution.org/r/1gl8p) for *Arabidopsis vs Fragaria x ananassa*. Here is the [link](https://genomevolution.org/r/1gla3) for the analysis between the H4 and Camarosa. And finally, here is the link to regenerate the analysis for Del Norte and Camarosa

## Running BLAST:
We are doing this step to identify homologs that may have been missed using a synteny-based approach. Genes that could have been missed by the synteny search include single-gene transpositions (and others). We are going to use a BLAST database of protein predictions. 

Caveat: I had to run the code in `scripts/AT_Fragaria/fix_fasta.py` for the Camarosa protein fasta because of the initial format of the FASTA file. The initial format was not crashing during BLAST, but it produced an erroneous output file that did not have the correct Camarosa gene names; this was likely due to incorrect splitting of the sequence ID in the original FASTA. My `fix_fasta.py` script merely renames the sequence IDs and is a little hard-coded to meet this goal. My guess is that it did not like the spaces in the sequence ID string and just took the first item, resulting in each "gene name" becoming "Fragaria".

In regards as to actually running BLAST, refer to the script at `scripts/AT_Fragaria/blastall.sb` for more information. First we generate a BLAST database and prepare it for protein indices. Then we can run the BLAST algorithm on the two sequence files, this may take awhile. For notes on the options for `blastall`, please refer to the [documentation](https://www.ncbi.nlm.nih.gov/Class/BLAST/blastallopts.txt). Caveat: When I was test running the script on a **masked** fasta file for Camarosa there were a few hundred lines of warnings (representing genes) in the output file stating that the Karlin-Altschul parameters could not be calculated due to an odd query sequence, I determined this to be the result of using a masked fasta file (some sequences were partially if not fully represented by N), and this was causing the warning.
