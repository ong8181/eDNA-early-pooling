####
#### Manual demultiplex from iSeq data
#### 2021.11.29: RIR-027 data
####

#---------------------------------------------------#
# REQUIRED: seqkit (https://bioinf.shenwei.me/seqkit)
# REQUIRED: cutadapt (https://github.com/marcelm/cutadapt)

WORKING_DIR="/XXXX/01_Demultiplex"
SEQ_DIR="../seqdata_raw/1_1st_2nd_indexing"
SAMPLE_DATA="../../../sampledata/RIR-027/RIR-027_IndexInfo.csv"
FILENAME_DATA="../../sampledata/RIR-027/RIR-027_file_rename.csv"
OUTPUT_DIR="01_Demultiplex_1_1st_2nd_indexingOut"
cd ${WORKING_DIR}
mkdir ${OUTPUT_DIR}


# ------------------------------------------------------------------------------------- #
# Step 1. Extract barcode region of fastq (the first 1-8 bp region)
# ------------------------------------------------------------------------------------- #
mkdir ${OUTPUT_DIR}/01_Out
cd ${SEQ_DIR}
# Manual demultiplexing only for Undetermined_*.fastq.gz
for file in Undetermined_*.fastq.gz; do
seqkit subseq -r 1:8 ${file} | gzip -c > ../../01_Demultiplex/${OUTPUT_DIR}/01_Out/${file%.fastq.gz}_index.fastq.gz
done


# ------------------------------------------------------------------------------------- #
# Step 2. Extract IDs of index sequences
# ------------------------------------------------------------------------------------- #
cd ../../01_Demultiplex/${OUTPUT_DIR}/01_Out
mkdir ../02_Out
### Forward read processing
while read row; do
  sample_name=$(echo ${row} | cut -d , -f 5)
  index=$(echo ${row} | cut -d , -f 10)
  seqkit grep -srip ^$index Undetermined_S0_L001_R1_001_index.fastq.gz -o ../02_Out/${sample_name}_R1_ID.fastq.gz
  seqkit grep -srip ^$index Undetermined_S0_L001_R2_001_index.fastq.gz -o ../02_Out/${sample_name}_R2_ID.fastq.gz
done < ${SAMPLE_DATA}


# ------------------------------------------------------------------------------------- #
# Step 2. Get common IDs of R1 and R2 reads
# ------------------------------------------------------------------------------------- #
cd ../02_Out
mkdir ../03_Out
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
for file in *_common_ID.txt; do
seqkit grep -f ${file} ../../${SEQ_DIR}/Undetermined_S0_L001_R1_001.fastq.gz | seqkit subseq -r 9:-1 | gzip -c > ../04_Out/${file%_common_ID.txt}_R1.fastq.gz
seqkit grep -f ${file} ../../${SEQ_DIR}/Undetermined_S0_L001_R2_001.fastq.gz | seqkit subseq -r 9:-1 | gzip -c > ../04_Out/${file%_common_ID.txt}_R2.fastq.gz
done


# ------------------------------------------------------------------------------------- #
# Step 4. Rename fastq files from Treatment 2 & 3
# ------------------------------------------------------------------------------------- #
cd ../../${SEQ_DIR}
# Use filename list
while read row; do
  sample_name=$(echo ${row} | cut -d , -f 1)
  read_id=$(echo ${row} | cut -d , -f 2)
  file_name=$(echo ${row} | cut -d , -f 3)
  cp ${file_name} ../../01_Demultiplex/01_Demultiplex_1_1st_2nd_indexingOut/04_Out/${sample_name}_${read_id}.fastq.gz
done < ${FILENAME_DATA}


# ------------------------------------------------------------------------------------- #
# Step 5. Trim primers
# ------------------------------------------------------------------------------------- #
cd ../../01_Demultiplex/01_Demultiplex_1_1st_2nd_indexingOut/04_Out
mkdir ../05_Out
for file in *_R1.fastq.gz; do
  # Single Primer removal
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
mv 05_Out/*.fastq.gz ../../seqdata_demultiplexed/1_1st_2nd_indexing

## Delete temporal files
rm -r 01_Out
rm -r 02_Out
rm -r 03_Out
rm -r 04_Out
rm -r 05_Out
