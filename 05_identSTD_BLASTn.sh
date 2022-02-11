####
#### Standard STD identification
####

#cd XXXX

DBPATH=seqdata_STD/FishSTD_MiFish
QUERYPATH=03_OTUClusteringOut/OTU_seqs.fa
OUTPUT_DIR=05_identSTD_BLASTnOut
EVALUE_SET=1e-50

mkdir ${OUTPUT_DIR}
blastn -db ${DBPATH} -query ${QUERYPATH} -evalue ${EVALUE_SET} -outfmt 6 -out ${OUTPUT_DIR}/STDseqOut.txt

