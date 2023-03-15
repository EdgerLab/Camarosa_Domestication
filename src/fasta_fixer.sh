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
# TODO think about putting everything in a congif file because this script is mega ugly.
# TODO this could be refactored with the Fix CDS
# FUTURE could just use a different function for each.

#--------------------------------------------------------
# PREPARE THE FASTA files

# NOTE, Del Norte needs a reformatted regular FASTA
DN="DN"
DN_REG_FASTA=/mnt/research/edgerpat_lab/Scotty/Strawberry_Domestication/data/Genomes/Del_Norte/DN.fa
DN_NEW_REG_FASTA=$DN_OUT_DIR/DN_NewNames.fa

# Viridis needs a reformatted regular FASTA
FVI="FVI"
FVI_REG_FASTA=/mnt/research/edgerpat_lab/Scotty/Strawberry_Domestication/data/Genomes/F_viridis/FVI.fa
FVI_NEW_REG_FASTA=$FVI_OUT_DIR/FVI_NewNames.fa

# Nipponica needs a reformatted regular FASTA
FNI="FNI"
FNI_REG_FASTA=/mnt/research/edgerpat_lab/Scotty/Strawberry_Domestication/data/Genomes/F_nipponica/FNI.fa
FNI_NEW_REG_FASTA=$FNI_OUT_DIR/FNI_NewNames.fa

# These genomes are OK, I am going to make a symbolic link to the file for the purposes of having similar file suffixes and locations
# NOTE, FII has a regular FASTA that is already properly formatted
FII='FII'
FII_REG_FASTA=/mnt/research/edgerpat_lab/Scotty/Strawberry_Domestication/data/Genomes/F_iinumae/FII.fa
FII_NEW_REG_FASTA=$FII_OUT_DIR/FII_NewNames.fa
# NOTE, RR has a regular FASTA that is already properly formatted
RR="RR"
RR_REG_FASTA=/mnt/research/edgerpat_lab/Scotty/Strawberry_Domestication/data/Genomes/Royal_Royce/RR.fa
RR_NEW_REG_FASTA=$RR_OUT_DIR/RR_NewNames.fa
# NOTE, H4 has a regular FASTA that is already properly formatted
H4="H4"
H4_REG_FASTA=/mnt/research/edgerpat_lab/Scotty/Strawberry_Domestication/data/Genomes/H4/H4.fa
H4_NEW_REG_FASTA=$H4_OUT_DIR/H4_NewNames.fa


#--------------------------------------------------------
# Load modules for Python work (FASTA renaming)
module purge
module load GCCcore/11.3.0 Python/3.10.4 && source /mnt/research/edgerpat_lab/Scotty/venvs/S_Domestication/bin/activate

#--------------------------------------------------------
# Preparing Regular FASTAs
# This must be done because the names are too long
python /mnt/research/edgerpat_lab/Scotty/Strawberry_Domestication/src/fix_fasta_names.py $DN_REG_FASTA $DN $DN_OUT_DIR $DN_NEW_REG_FASTA
python /mnt/research/edgerpat_lab/Scotty/Strawberry_Domestication/src/fix_fasta_names.py $FVI_REG_FASTA $FVI $FVI_OUT_DIR $FVI_NEW_REG_FASTA
python /mnt/research/edgerpat_lab/Scotty/Strawberry_Domestication/src/fix_fasta_names.py $FNI_REG_FASTA $FNI $FNI_OUT_DIR $FNI_NEW_REG_FASTA

echo "Make symbolic links for the FASTA files that are already in proper format."
ln -svf $FII_REG_FASTA $FII_NEW_REG_FASTA
ln -svf $RR_REG_FASTA $RR_NEW_REG_FASTA
ln -svf $H4_REG_FASTA $H4_NEW_REG_FASTA

echo "Gzip the original files that were re-made with edited names, to save disk space"
gzip $DN_REG_FASTA $FVI_REG_FASTA $FNI_REG_FASTA
