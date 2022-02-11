####
#### Manual demultiplex from iSeq data
#### 2021.11.29: RIR-031 data
####

#---------------------------------------------------#
# REQUIRED: seqkit (https://bioinf.shenwei.me/seqkit)
# REQUIRED: cutadapt (https://github.com/marcelm/cutadapt)

WORKING_DIR="/XXXX/01_Demultiplex"
SEQ_DIR="../seqdata_raw/3_rep_vol_test"
SAMPLE_DATA="../../../sampledata/RIR-031/RIR-031_IndexInfo.csv"
OUTPUT_DIR="01_Demultiplex_3_rep_vol_testOut"
cd ${WORKING_DIR}
mkdir ${OUTPUT_DIR}


# ------------------------------------------------------------------------------------- #
# Step 1. Extract barcode region of fastq (the first 1-8 bp region)
# ------------------------------------------------------------------------------------- #
mkdir ${OUTPUT_DIR}/01_Out
cd ${SEQ_DIR}
for file in *.fastq.gz; do
  seqkit subseq -r 1:8 ${file} | gzip -c > ../../01_Demultiplex/${OUTPUT_DIR}/01_Out/${file%.fastq.gz}_index.fastq.gz
done


# ------------------------------------------------------------------------------------- #
# Step 2. Extract IDs of index sequences
# ------------------------------------------------------------------------------------- #
cd ../../01_Demultiplex/${OUTPUT_DIR}/01_Out
mkdir ../02_Out

# Processing sequence reads
while read row; do
  sample_id=$(echo ${row} | cut -d , -f 5)
  index=$(echo ${row} | cut -d , -f 10)
  seqdata_file=$(echo ${row} | cut -d , -f 3)
  seqkit grep -srip ^$index ${seqdata_file}_R1_001_index.fastq.gz -o ../02_Out/${sample_id}_R1_ID.fastq.gz
  seqkit grep -srip ^$index ${seqdata_file}_R2_001_index.fastq.gz -o ../02_Out/${sample_id}_R2_ID.fastq.gz
done < ${SAMPLE_DATA}


# ------------------------------------------------------------------------------------- #
# Step 2. Get common IDs of R1 and R2 reads
# ------------------------------------------------------------------------------------- #
cd ../02_Out
mkdir ../03_Out
### manually remove non-paired files
## ID sequences of R1 reads that match R1_ID
for file in *_R1_ID.fastq.gz; do
  seqkit grep -f <(seqkit seq -ni ${file}) ${file%_R1_ID.fastq.gz}_R2_ID.fastq.gz | seqkit seq -ni > ../03_Out/${file%_R1_ID.fastq.gz}_common_ID.txt
done


# ------------------------------------------------------------------------------------- #
# Step 3. Get sequences based on the common IDs
# (Trim index sequences here)
# ------------------------------------------------------------------------------------- #
cd ../03_Out
mkdir ../04_Out

# Demultiplex all fastq files
while read row; do
  sample_id=$(echo ${row} | cut -d , -f 5)
  seqdata_file=$(echo ${row} | cut -d , -f 3)
  seqkit grep -f ${sample_id}_common_ID.txt ../../${SEQ_DIR}/${seqdata_file}_R1_001.fastq.gz | seqkit subseq -r 9:-1 | gzip -c > ../04_Out/${sample_id}_R1.fastq.gz
  seqkit grep -f ${sample_id}_common_ID.txt ../../${SEQ_DIR}/${seqdata_file}_R2_001.fastq.gz | seqkit subseq -r 9:-1 | gzip -c > ../04_Out/${sample_id}_R2.fastq.gz
done < ${SAMPLE_DATA}


# ------------------------------------------------------------------------------------- #
# Step 4. Trim primers
# ------------------------------------------------------------------------------------- #
cd ../04_Out/
mkdir ../05_Out
for file in *_R1.fastq.gz; do
  # Single Primer removal (MiFish primers)
  cutadapt -j 4 \
  -g GTCGGTAAAACTCGTGCCAGC \
  -a CAAACTGGGATTAGATACCCCACTATG \
  -G CATAGTGGGGTATCTAATCCCAGTTTG \
  -A GCTGGCACGAGTTTTACCGAC \
  -n 2 \
  -o ../05_Out/${file%_R1.fastq.gz}_trimmed_R1.fastq.gz \
  -p ../05_Out/${file%_R1.fastq.gz}_trimmed_R2.fastq.gz \
  ${file} \
  ${file%_R1.fastq.gz}_R2.fastq.gz
done

#-g Forward \
#-a Reverse-RevComp \
#-G Reverse \
#-A Forward-RevComp \

# ------------------------------------------------------------------------------------- #
# Step 5. Summarize the number of sequence reads
# ------------------------------------------------------------------------------------- #
cd ..
mkdir stats
seqkit stats -a  ../${SEQ_DIR}/*.fastq.gz > stats/00_OriginalSeqdata.txt
seqkit stats -a  ../${SEQ_DIR}/*_R1_001.fastq.gz > stats/00_OriginalSeqdata_R1_only.txt
seqkit stats -a  01_Out/*.fastq.gz > stats/01_BarcodeSeq.txt
seqkit stats -a  02_Out/*.fastq.gz > stats/02_IndexSeq.txt
seqkit stats -a  04_Out/*.fastq.gz > stats/04_Demultiplexed.txt
seqkit stats -a  04_Out/*_R1.fastq.gz > stats/04_Demultiplexed_R1_only.txt
seqkit stats -a  05_Out/*.fastq.gz > stats/05_PrimerTrimmed.txt
seqkit stats -a  05_Out/*_R1.fastq.gz > stats/05_PrimerTrimmed_R1_only.txt

## Move primer-trimmed sequences
mv 05_Out/*.fastq.gz ../../seqdata_demultiplexed/3_rep_vol_test

## Delete temporal files
rm -r 01_Out
rm -r 02_Out
rm -r 03_Out
rm -r 04_Out
rm -r 05_Out
