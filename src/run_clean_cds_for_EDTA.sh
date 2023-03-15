#--------------------------------------------------------
# Define base paths
OUT_DATA_DIR=/mnt/research/edgerpat_lab/Scotty/Strawberry_Domestication/data/Genomes/Fixed_Fastas

#--------------------------------------------------------
# Make directories for fixed FASTA files
DN_OUT_DIR=$OUT_DATA_DIR/Del_Norte
RR_OUT_DIR=$OUT_DATA_DIR/Royal_Royce
H4_OUT_DIR=$OUT_DATA_DIR/H4
FNI_OUT_DIR=$OUT_DATA_DIR/F_nipponica
FVI_OUT_DIR=$OUT_DATA_DIR/F_viridis
FII_OUT_DIR=$OUT_DATA_DIR/F_iinumae
mkdir -p $DN_OUT_DIR $RR_OUT_DIR $H4_OUT_DIR $FNI_OUT_DIR $FVI_OUT_DIR $FII_OUT_DIR
#--------------------------------------------------------
# FUTURE: perhaps there is a way to integrate this with the "cleaning" of the CDS files prior to running BLAST.
# FUTURE: maybe things can be integrated and I can use a config file? There is a lot of repeated code here...

# Get filepaths for each genome
DN="DN"
DN_CDS_FASTA=/mnt/research/edgerpat_lab/Scotty/Strawberry_Domestication/data/Genomes/Fixed_Fastas/Del_Norte/DN_CDS_NewNames.fa
BLACKLISTED_DN=/mnt/research/edgerpat_lab/Scotty/Strawberry_Domestication/data/Genomes/Fixed_Fastas/Del_Norte/DN_CDS_sans_blacklist.fa
DN_BLAST=/mnt/research/edgerpat_lab/Scotty/Strawberry_Domestication/doc/LINE/results/Del_Norte_NingTpases_RESULTS.txt

RR="RR"
RR_CDS_FASTA=/mnt/research/edgerpat_lab/Scotty/Strawberry_Domestication/data/Genomes/Fixed_Fastas/Royal_Royce/RR_CDS_NewNames.fa
BLACKLISTED_RR=/mnt/research/edgerpat_lab/Scotty/Strawberry_Domestication/data/Genomes/Fixed_Fastas/Royal_Royce/RR_CDS_sans_blacklist.fa
RR_BLAST=/mnt/research/edgerpat_lab/Scotty/Strawberry_Domestication/doc/LINE/results/Royal_Royce_NingTpases_RESULTS.txt

H4="H4"
H4_CDS_FASTA=/mnt/research/edgerpat_lab/Scotty/Strawberry_Domestication/data/Genomes/Fixed_Fastas/H4/H4_CDS_NewNames.fa
BLACKLISTED_H4=/mnt/research/edgerpat_lab/Scotty/Strawberry_Domestication/data/Genomes/Fixed_Fastas/H4/H4_CDS_sans_blacklist.fa
H4_BLAST=/mnt/research/edgerpat_lab/Scotty/Strawberry_Domestication/doc/LINE/results/H4_NingTpases_RESULTS.txt

FII="FII"
FII_CDS_FASTA=/mnt/research/edgerpat_lab/Scotty/Strawberry_Domestication/data/Genomes/Fixed_Fastas/F_iinumae/FII_CDS_NewNames.fa
BLACKLISTED_FII=/mnt/research/edgerpat_lab/Scotty/Strawberry_Domestication/data/Genomes/Fixed_Fastas/F_iinumae/FII_CDS_sans_blacklist.fa
FII_BLAST=/mnt/research/edgerpat_lab/Scotty/Strawberry_Domestication/doc/LINE/results/Iinumae_NingTpases_RESULTS.txt

FVI="FVI"
FVI_CDS_FASTA=/mnt/research/edgerpat_lab/Scotty/Strawberry_Domestication/data/Genomes/Fixed_Fastas/F_viridis/FVI_CDS_NewNames.fa
BLACKLISTED_FVI=/mnt/research/edgerpat_lab/Scotty/Strawberry_Domestication/data/Genomes/Fixed_Fastas/F_viridis/FVI_CDS_sans_blacklist.fa
FVI_BLAST=/mnt/research/edgerpat_lab/Scotty/Strawberry_Domestication/doc/LINE/results/Viridis_NingTpases_RESULTS.txt

FNI="FNI"
FNI_CDS_FASTA=/mnt/research/edgerpat_lab/Scotty/Strawberry_Domestication/data/Genomes/Fixed_Fastas/F_nipponica/FNI_CDS_NewNames.fa
BLACKLISTED_FNI=/mnt/research/edgerpat_lab/Scotty/Strawberry_Domestication/data/Genomes/Fixed_Fastas/F_nipponica/FNI_CDS_sans_blacklist.fa
FNI_BLAST=/mnt/research/edgerpat_lab/Scotty/Strawberry_Domestication/doc/LINE/results/Nipponica_NingTpases_RESULTS.txt

#--------------------------------------------------------
# Load modules for Python work (FASTA renaming)
module purge
module load GCCcore/11.3.0 Python/3.10.4 && source /mnt/research/edgerpat_lab/Scotty/venvs/S_Domestication/bin/activate

#--------------------------------------------------------
# Prepare the CDS FASTAs
SCRIPT=/mnt/research/edgerpat_lab/Scotty/Strawberry_Domestication/src/remove_hits_from_CDS.py



# NOTE, done
# Del Norte
python $SCRIPT $DN_CDS_FASTA $DN_BLAST $BLACKLISTED_DN $DN $DN_OUT_DIR

# Royal Royce
python $SCRIPT $RR_CDS_FASTA $RR_BLAST $BLACKLISTED_RR $RR $RR_OUT_DIR

# H4
python $SCRIPT $H4_CDS_FASTA $H4_BLAST $BLACKLISTED_H4 $H4 $H4_OUT_DIR 

# FII
python $SCRIPT $FII_CDS_FASTA $FII_BLAST $BLACKLISTED_FII $FII $FII_OUT_DIR 

# FVI
python $SCRIPT $FVI_CDS_FASTA $FVI_BLAST $BLACKLISTED_FVI $FVI $FVI_OUT_DIR 

# FNI
python $SCRIPT $FNI_CDS_FASTA $FNI_BLAST $BLACKLISTED_FNI $FNI $FNI_OUT_DIR 
