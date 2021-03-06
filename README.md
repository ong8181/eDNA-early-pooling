# Analysis codes for Ushio et al. (2022) "An efficient early-pooling protocol for environmental DNA metabarcoding" _Environmental DNA_
[![DOI](https://zenodo.org/badge/458195380.svg)](https://zenodo.org/badge/latestdoi/458195380)

This repository contains analysis codes to reproduce the results in Ushio et al. (2022) _Environmental DNA_ https://doi.org/10.1002/edn3.337. A preprint version is also available in
_bioRxiv_ https://doi.org/10.1101/2022.02.15.480497.

:heavy_exclamation_mark: If you want to perform the analyses from the demultiplexing step using the original FASTQ files, you should download FASTQ files from DDBJ (see the section "Downloading sequence data" below). Note that the FASTQ file in DDBJ are renamed to follow DDBJ instructions. The FASTQ files should be re-renamed by executing the commands in the section "Downloading sequence data".

Alternatively, please contact the corresponding author (Masayuki Ushio; ong8181@gmail.com).

# License
See LICENSE.

# Analysis flow
Results of the analysis are stored in `xxxxOut` folers.


## Step 1. Sequence data preprocessing
- `01_Demultiplex/`: Demultiplex sequence data using `seqkit`<br>
Codes in this folder were used to rename and demultiplex the original FASTQ files (generated by Basespace functionality), and are not required to reproduce the results. The codes are placed here to show how the original FASTQ files were preprocessed.

## Step 2. Sequence data processing (`dada2` and `DECIPHER`)
- `02_DADA2/`: Create an ASV-sample matrix using `dada2` <br>
- `03_OTUClustering.R`: Cluster ASVs into OTUs using `DECIPHER`. `DECIPHER` version should be <= 2.22.0.<br>


## Step 3. Taxa assignments (`Claident`)
- `04_TaxaAssignment.sh`: Assign taxa using `Claident`<br>
- `05_identSTD_BLASTn.sh`: Identify standard DNA sequences<br>
- `06_TaxaSTDcombine.R`: Compile taxa information<br>

## Step 4. Statistical analysis (`phyloseq`) and visualization (`ggplot2`)
- `07_CompilePhyloseq.R`: Prepare phyloseq objects<br>
- `08_Exp1_1st2nd/`: Analyze Experiment I<br>
- `09_Exp2_exosap/`: Analyze Experiment II<br>
- `10_Exp3_repvol/`: Analyze Experiment III<br>
- `11_SameIndexTest/`: Test index-specific biases<br>
- `12_NagahamaNMDS/`: Test protocol-specific biases<br>
- `13_RareOTU/`: Examine how rare OTUs are detected<br>


## Step 5. Format figures (ggplot2)
- `FigCode/`: Format figures<br>
- `FigCode/Fig_Combine.Rmd`: Combine figures<br>


# Package versions
see text files in `00_SessionInfo/`


# Downloading sequence data
Demultiplexed FASTQ files are available in DDBJ DRA (accession number = DRA013399). If you want to reproduce the analysis, place the downloaded FASTQ files in `seqdata_demultiplexed` folder using the following commands. You can re-run the analysis from DADA2 analysis.

:heavy_exclamation_mark: Each sample folder usually contains two FASTQ files (paired-end files). However, four sample folders (DRX333539, DRX333586, DRX333596, and DRX333675) contain three FASTQ files. The third FASTQ file contains sequences of which paired-sequences were not found, which can simply be discarded when you want to reproduce the analyses.

```
# Go to the repository folder
cd eDNA-early-pooling

# Create folders to place demultiplexed FASTQ files for the three experiments
mkdir seqdata_demultiplexed/1_1st_2nd_indexing
mkdir seqdata_demultiplexed/2_exosap
mkdir seqdata_demultiplexed/3_rep_vol_test

# Create temporal folders
mkdir seqdata_raw/temp
mkdir seqdata_raw/temp/xml
cd seqdata_raw/temp

# Download data (downlaod only *.bz2 and *.xml files)
wget -r --no-parent -A "*.bz2","*.xml" "ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/DRA013/DRA013399/"

# If you need to specify the proxy:
#wget -r --no-parent -A "*.bz2","*.xml" -e FTP_PROXY=proxy.xxx.xx:xxxx "ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/DRA013/DRA013399/"

# Move files
mv ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/DRA013/DRA013399/*/*.fastq.bz2 ./
mv ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/DRA013/DRA013399/*.xml xml

# Delete temporal folder
rm -r ftp.ddbj.nig.ac.jp

# Convert bz2 to gz format
for f in *.bz2; do
  bzcat "$f" | gzip -c >"${f%.bz2}.gz"
  rm $f
done

# Move files to each experiment folder
## Experiment I
count=0
while read row || [ -n "${row}" ]; do
  if ((count > 0)); then
    file_name=$(echo ${row} | cut -d , -f 14)
    echo "Moving ${file_name} files..."
    mv ${file_name}* ../../seqdata_demultiplexed/1_1st_2nd_indexing
  fi
  count=`expr $count + 1`
done < ../../sampledata/DRA_metadata/ushio-0017_Exp1_1st_2nd_indexing.csv

## Experiment II
count=0
while read row || [ -n "${row}" ]; do
  if ((count > 0)); then
    file_name=$(echo ${row} | cut -d , -f 14)
    echo "Moving ${file_name} files..."
    mv ${file_name}* ../../seqdata_demultiplexed/2_exosap
  fi
  count=`expr $count + 1`
done < ../../sampledata/DRA_metadata/ushio-0017_Exp2_exosap.csv

## Experiment III
count=0
while read row || [ -n "${row}" ]; do
  if ((count > 0)); then
    file_name=$(echo ${row} | cut -d , -f 14)
    echo "Moving ${file_name} files..."
    mv ${file_name}* ../../seqdata_demultiplexed/3_rep_vol_test
  fi
  count=`expr $count + 1`
done < ../../sampledata/DRA_metadata/ushio-0017_Exp3_rep_vol_test.csv

# Delete the temporal folder
mv xml ..
cd ..
rm -r temp

# Change the names of the FASTQ files to the original ones
## Experiment I
cd ../seqdata_demultiplexed/1_1st_2nd_indexing
count=0
while read row || [ -n "${row}" ]; do
  if ((count > 0)); then
    file_name1=$(echo ${row} | cut -d , -f 14)
    file_name2=$(echo ${row} | cut -d , -f 17 | cut -c 2-22)
    echo "Renaming ${file_name1} files..."
    mv ${file_name1}_1.fastq.gz ${file_name2}_R1.fastq.gz
    mv ${file_name1}_2.fastq.gz ${file_name2}_R2.fastq.gz
  fi
  count=`expr $count + 1`
done < ../../sampledata/DRA_metadata/ushio-0017_Exp1_1st_2nd_indexing.csv

## Experiment II
cd ../2_exosap
count=0
while read row || [ -n "${row}" ]; do
  if ((count > 0)); then
    file_name1=$(echo ${row} | cut -d , -f 14)
    file_name2=$(echo ${row} | cut -d , -f 17 | cut -c 2-22)
    echo "Renaming ${file_name1} files..."
    mv ${file_name1}_1.fastq.gz ${file_name2}_R1.fastq.gz
    mv ${file_name1}_2.fastq.gz ${file_name2}_R2.fastq.gz
  fi
  count=`expr $count + 1`
done < ../../sampledata/DRA_metadata/ushio-0017_Exp2_exosap.csv

## Experiment III
cd ../3_rep_vol_test
count=0
while read row || [ -n "${row}" ]; do
  if ((count > 0)); then
    file_name1=$(echo ${row} | cut -d , -f 14)
    file_name2=$(echo ${row} | cut -d , -f 17 | cut -c 2-22)
    echo "Renaming ${file_name1} files..."
    mv ${file_name1}_1.fastq.gz ${file_name2}_R1.fastq.gz
    mv ${file_name1}_2.fastq.gz ${file_name2}_R2.fastq.gz
  fi
  count=`expr $count + 1`
done < ../../sampledata/DRA_metadata/ushio-0017_Exp3_rep_vol_test.csv
```

Now you can move to `02_DADA2` folder and start the analysis in your computer!!
