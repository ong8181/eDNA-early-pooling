####
#### DADA2 analysis of fastq files
#### Paired-end (iSeq 150PE)
#### 2021.11.29 Ushio
#### R 4.1.2
####

#--------------- DADA2 processing ---------------#
# Set random seeds (for reproduction)
ran.seed <- 1234
set.seed(ran.seed)
dir.create("../00_SessionInfo")

# Load library and functions
library(dada2); packageVersion("dada2") #1.22.0, 2021.11.16
library(ShortRead); packageVersion("ShortRead") #1.52.0, 2021.11.16
library(tidyverse); packageVersion("tidyverse") #1.3.1, 2021.8.4
source("../functions_R/F01_HelperFunctions.R")

# Create output directory
od <- basename(rstudioapi::getSourceEditorContext()$path)
(output_folder <- paste0(str_sub(od, end = -3), "Out")); rm(od)
dir.create(output_folder)

# Load sequence reads
path <- "../seqdata_demultiplexed/2_exosap"
fnFs <- sort(list.files(path, pattern="R1.fastq.gz", full.names = T)) # Forward read files
fnRs <- sort(list.files(path, pattern="R2.fastq.gz", full.names = T)) # Reverse read files
# Get sample names, assuming files named as so: SAMPLENAME_XXX.fastq
sample_names <- sapply(strsplit(basename(fnFs), "_trimmed"), `[`, 1)

# Visualize quality
#plotQualityProfile(fnFs[1])
#plotQualityProfile(fnRs[1])

# ------------------------ Primer removal check ---------------------------- #
# Identify primers
FWD <- "GTCGGTAAAACTCGTGCCAGC" # MiFish-F
REV <- "CATAGTGGGGTATCTAATCCCAGTTTG" # MiFish-R
FWD_orients <- AllOrients(FWD)
REV_orients <- AllOrients(REV)

# Identify primers
seq_id <- 1
rbind(FWD.ForwardReads = sapply(FWD_orients, PrimerHits, fn = fnFs[[seq_id]]),
      FWD.ReverseReads = sapply(FWD_orients, PrimerHits, fn = fnRs[[seq_id]]), 
      REV.ForwardReads = sapply(REV_orients, PrimerHits, fn = fnFs[[seq_id]]), 
      REV.ReverseReads = sapply(REV_orients, PrimerHits, fn = fnRs[[seq_id]]))

# Performing filtering and trimming
filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample_names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample_names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen=c(110, 110), # Truncate the end of reads
                     minLen = 100, # Remove unexpectedly short reads
                     maxN = 0, maxEE = c(2,2), truncQ = 2, rm.phix = T,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
out
#plotQualityProfile(filtFs[1])
#plotQualityProfile(filtRs[1])

# Exclude 0 seq samples, rename filtFs and filtRs
if(length(sample_names[out[,2]<1 | out[,1]<1]) > 0){
  filtFs <- file.path(filt_path, paste0(sample_names[out[,2]>0 & out[,1]>0], "_F_filt.fastq.gz"))
  filtRs <- file.path(filt_path, paste0(sample_names[out[,2]>0 & out[,1]>0], "_R_filt.fastq.gz"))
}

# Learn the error rates
min_nbases <- 110 * sum(out[,2]) # Use a small number of bases to speed up the analysis
errF <- learnErrors(filtFs, multithread = TRUE, randomize = TRUE, MAX_CONSIST = 20, nbases = min_nbases)
errR <- learnErrors(filtRs, multithread = TRUE, randomize = TRUE, MAX_CONSIST = 20, nbases = min_nbases)

# Visualize errors
ggsave(sprintf("%s/errF.pdf", output_folder), plotErrors(errF, nominalQ = T), width = 10, height = 10)
ggsave(sprintf("%s/errR.pdf", output_folder), plotErrors(errR, nominalQ = T), width = 10, height = 10)

# Dereplicatin
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample_names[out[,2]>0 & out[,1]>0]
names(derepRs) <- sample_names[out[,2]>0 & out[,1]>0]

# Sample inference
dadaFs <- dada(derepFs, err=errF, multithread = TRUE, pool = TRUE)
dadaRs <- dada(derepRs, err=errR, multithread = TRUE, pool = TRUE)

# Merging paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE, minOverlap = 20)

# Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab) # 61 samples with 265 ASVs (pooling)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
# Cutting unexpected length sequences
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(150,190)]
table(nchar(getSequences(seqtab2)))

# Remove chimeras
seqtab_nochim <- removeBimeraDenovo(seqtab2, method = "consensus", multithread = TRUE, verbose = TRUE)
dim(seqtab_nochim)
sum(seqtab_nochim)/sum(seqtab2)

# Track reads thourhg the pipeline
out2 <- out[out[,2]>0 & out[,1]>0,]
getN <- function(x) sum(getUniques(x))
track <- cbind(out2, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab2), rowSums(seqtab_nochim),  rowSums(seqtab_nochim)/out2[,1])
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "tabled2", "nonchim", "prop(last/first)")
rownames(track) <- sample_names[out[,2]>0 & out[,1]>0]
head(track)

# Taxa output for claident tax assignment
seqs <- colnames(seqtab_nochim)
seqs_out <- as.matrix(c(rbind(sprintf(">ASV%06d", 1:length(seqs)), seqs)), ncol = 1)
seqtab_nochim2 <- seqtab_nochim
colnames(seqtab_nochim2) <- sprintf("ASV%06d", 1:length(seqs)) # Temporal object for saving

# Save outputs
saveRDS(seqtab_nochim, paste0(output_folder, "/seqtab_nochim.obj"))
write.csv(seqs, paste0(output_folder, "/seq_only.csv"), row.names = colnames(seqtab_nochim2))
write.table(seqs_out, paste0(output_folder, "/ASV_seqs.fa"), col.names = FALSE, row.names = FALSE, quote = FALSE)
write.csv(seqtab_nochim2, paste0(output_folder, "/seqtab_nochim.csv"))
write.csv(track, paste0(output_folder, "/track.csv"))

# Save workspace
rm(seqtab_nochim2)
rm(derepFs)
rm(derepRs)
rm(dadaFs)
rm(dadaRs)
save(list = ls(all.names = TRUE),
     file = paste0(output_folder, "/", output_folder, ".RData"))

# Save session info
writeLines(capture.output(sessionInfo()),
           paste0("../00_SessionInfo/", output_folder, "_", substr(Sys.time(), 1, 10), ".txt"))
