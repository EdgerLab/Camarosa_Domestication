# Contact Information:
| Role          | Name          | GitHub                                                  | Email              |
|---------------|---------------|---------------------------------------------------------|--------------------|
| Project Lead: | Scott Teresi  | [Personal GitHub](https://github.com/huckleberry-hound) | <teresisc@msu.edu> |
| PI:           | Patrick Edger | [Lab GitHub](https://github.com/EdgerLab)               | <edgerpat@msu.edu> |

# Identification of Syntenic Genes:

## Genome:
| Regular CoGe ID                   | Masked ID                                    |
|-----------------------------------|----------------------------------------------|
| Arabidopsis thaliana Col-0 (id 1) | CNS PL.20 Masked repeats 50X (v10, id 16746) |
| Fragaria x ananass subsp. ananassa (id 41115) | masked (v1, id 35784) |
| Fragaria vesva subsp vesca (id 34820) | NCBI WindowMasker Hard (vv4, id 57743) |

## Running SynMap:
This section describes the methods to run [SynMap](https://genomevolution.org/CoGe/SynMap.pl) on CoGe. I ran SynMap with mostly [default options](https://genomevolution.org/wiki/index.php/SynMap), I did change one option: under *Merge Syntenic Blocks* I set it to `Quota Align Merge`. Here is the [link](https://genomevolution.org/r/1fo40) to regenerate the analysis on CoGe for *Arabidopsis vs Fragaria x ananassa* and here is the [link](https://genomevolution.org/r/1fo4c) for *Arabidopsis vs Fragaria vesca*.

## Running BLAST:
We are doing this step to identify homologs that may have been missed using a synteny-based approach. Genes that could have been missed by the synteny search include single-gene transpositions (and others). We are going to use a BLAST database of protein predictions.

Please refer to the script at `scripts/At_BB/blastall.sb` for more information. First we generate a BLAST database and prepare it for protein indices. Then we can run the BLAST algorithm on the two sequence files, this may take awhile. For notes on the options for `blastall`, please refer to the [documentation](https://www.ncbi.nlm.nih.gov/Class/BLAST/blastallopts.txt).
