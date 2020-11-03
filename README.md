# Contact Information:
| Role          | Name          | GitHub                                                  | Email              |
|---------------|---------------|---------------------------------------------------------|--------------------|
| Project Lead: | Scott Teresi  | [Personal GitHub](https://github.com/huckleberry-hound) | <teresisc@msu.edu> |
| PI:           | Patrick Edger | [Lab GitHub](https://github.com/EdgerLab)               | <edgerpat@msu.edu> |


# Identification of Syntenic and Homologous Genes:

## Design Plan:
* Generate a list of syntelogs and homologs from Camarosa to *Arabidopsis*.
* Generate a list of syntelogs for each non-Camarosa strawberry genome to Camarosa
* Generate a table containing the above information in an easily searchable format

## Genomes:
| Regular CoGe ID                   | Masked ID                                    | Peptide Source |
|-----------------------------------|----------------------------------------------|----------------|
| Arabidopsis thaliana Col-0 (id 1) | CNS PL.20 Masked repeats 50X (v10, id 16746) | CoGe Derived |
| Fragaria x ananass subsp. ananassa (id 41115) | masked (v1, id 35784) | CoGe Derived |
| Fragaria x ananass subsp. ananassa (id 41115) | unmasked (v1, id 35789) | CoGe Derived |
| Fragaria vesca subsp vesca (vv4, id 57743) | NCBI WindowMasker Hard | CoGe Derived [^1] |
| Fragaria iinumae (id 54118) | Unmasked | CoGe Derived |

## Running SynMap:
This section describes the methods to run [SynMap](https://genomevolution.org/CoGe/SynMap.pl) on CoGe. I ran SynMap with mostly [default options](https://genomevolution.org/wiki/index.php/SynMap), I did change one option: under *Merge Syntenic Blocks* I set it to `Quota Align Merge`. Here is the [link](https://genomevolution.org/r/1fo40) to regenerate the analysis on CoGe for *Arabidopsis vs Fragaria x ananassa*. But now I just re-did the Synmap with the unmasked version and here is the [link](https://genomevolution.org/r/1gl8p). Here is the [link](https://genomevolution.org/r/1gla3) for the analysis between the masked vesca and unmasked Cam.

### NOTE may be able to remove this portion
 here is the [link](https://genomevolution.org/r/1fo4c) for *Arabidopsis vs Fragaria vesca*, and here is the [link](https://genomevolution.org/r/1fvm3) to regernerate the analysis of *Arabidopsis vs Fragaria iinumae* (link may be restricted). For the sake of simplicity, the Arabidopsis genome was always the first in the comparison order.

## Running BLAST:
We are doing this step to identify homologs that may have been missed using a synteny-based approach. Genes that could have been missed by the synteny search include single-gene transpositions (and others). We are going to use a BLAST database of protein predictions. 

I had to run the code in `scripts/AT_Fragaria/fix_fasta.py` for the Camarosa protein fasta because of the initial format of the FASTA file. The initial format was not crashing during BLAST, but it produced an erroneous output file that did not have the correct Camarosa gene names; this was likely due to incorrect splitting of the sequence ID in the original FASTA. My `fix_fasta.py` script merely renames the sequence IDs and is a little hard-coded to meet this goal. My guess is that it did not like the spaces in the sequence ID string and just took the first item, resulting in each "gene name" becoming "Fragaria".

Please refer to the script at `scripts/AT_Fragaria/blastall.sb` for more information. First we generate a BLAST database and prepare it for protein indices. Then we can run the BLAST algorithm on the two sequence files, this may take awhile. For notes on the options for `blastall`, please refer to the [documentation](https://www.ncbi.nlm.nih.gov/Class/BLAST/blastallopts.txt). Caveat, when running `scripts/AT_Fragaria/blastall.sb` for Camarosa there were a few hundred lines of warnings (representing genes) in the output file stating that the Karlin-Altschul parameters could not be calculated due to an odd query sequence, I determined this to be the result of using a masked fasta file, and proceeded to continue with the analysis.


# Footnotes
[^1]: Here is the [publication](http://gigadb.org/dataset/100372) for more data.
