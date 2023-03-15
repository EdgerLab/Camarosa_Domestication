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
# Get filepaths for each genome
DN="DN"
DN_CDS_FASTA=/mnt/research/edgerpat_lab/Scotty/Strawberry_Domestication/data/Genomes/Del_Norte/DN_CDS.fa
DN_NEW_CDS_FASTA=$DN_OUT_DIR/DN_CDS_NewNames.fa

RR="RR"
RR_CDS_FASTA=/mnt/research/edgerpat_lab/Scotty/Strawberry_Domestication/data/Genomes/Royal_Royce/RR_CDS.fa
RR_NEW_CDS_FASTA=$RR_OUT_DIR/RR_CDS_NewNames.fa

H4="H4"
H4_CDS_FASTA=/mnt/research/edgerpat_lab/Scotty/Strawberry_Domestication/data/Genomes/H4/H4_CDS.fa
H4_NEW_CDS_FASTA=$H4_OUT_DIR/H4_CDS_NewNames.fa

FII="FII"
FII_CDS_FASTA=/mnt/research/edgerpat_lab/Scotty/Strawberry_Domestication/data/Genomes/F_iinumae/FII_CDS.fa
FII_NEW_CDS_FASTA=$FII_OUT_DIR/FII_CDS_NewNames.fa

FVI="FVI"
FVI_CDS_FASTA=/mnt/research/edgerpat_lab/Scotty/Strawberry_Domestication/data/Genomes/F_viridis/FVI_CDS.fa
FVI_NEW_CDS_FASTA=$FVI_OUT_DIR/FVI_CDS_NewNames.fa

FNI="FNI"
FNI_CDS_FASTA=/mnt/research/edgerpat_lab/Scotty/Strawberry_Domestication/data/Genomes/F_nipponica/FNI_CDS.fa
FNI_NEW_CDS_FASTA=$FNI_OUT_DIR/FNI_CDS_NewNames.fa

#--------------------------------------------------------
# Load modules for Python work (FASTA renaming)
module purge
module load GCCcore/11.3.0 Python/3.10.4 && source /mnt/research/edgerpat_lab/Scotty/venvs/S_Domestication/bin/activate

#--------------------------------------------------------
# Prepare the CDS FASTAs
SCRIPT=/mnt/research/edgerpat_lab/Scotty/Strawberry_Domestication/src/trim_CDS_names_for_blast.py

# Del Norte
python $SCRIPT $DN_CDS_FASTA $DN_NEW_CDS_FASTA $DN
# Royal Royce
python $SCRIPT $RR_CDS_FASTA $RR_NEW_CDS_FASTA $RR
# H4
python $SCRIPT $H4_CDS_FASTA $H4_NEW_CDS_FASTA $H4
# FII
python $SCRIPT $FII_CDS_FASTA $FII_NEW_CDS_FASTA $FII
# FVI
python $SCRIPT $FVI_CDS_FASTA $FVI_NEW_CDS_FASTA $FVI
# FNI
python $SCRIPT $FNI_CDS_FASTA $FNI_NEW_CDS_FASTA $FNI
