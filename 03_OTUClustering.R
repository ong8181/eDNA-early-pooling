####
#### Clustering DADA2 results
#### 2021.11.29 Ushio (R4.1.2)
#### R 4.1.2
####

# Load library and functions
library(tidyverse); packageVersion("tidyverse") #1.3.1, 2021.11.11
library(phyloseq); packageVersion("phyloseq") #1.38.0, 2021.11.11
library(ShortRead); packageVersion("ShortRead") #1.52.0, 2021.11.11
library(dada2); packageVersion("dada2") #1.22.0, 2021.11.11
library(cowplot); packageVersion("cowplot") # 1.1.1, 2021.6.12
theme_set(theme_cowplot())
# Packages that are required but not loaded:
# library(DECIPHER)
# library(Biostrings)

# Set random seeds (for reproduction)
ran.seed <- 1234
set.seed(ran.seed)
od <- basename(rstudioapi::getSourceEditorContext()$path)
(output_folder <- paste0(str_sub(od, end = -3), "Out")); rm(od)
dir.create(output_folder)


# <-----------------------------------------------------> #
#                       Load data
# <-----------------------------------------------------> #
# Load seqtab_nochim data
seqtab_nochim1 <- readRDS("02_DADA2/02_DADA2_1_1st_2nd_indexingOut/seqtab_nochim.obj")
seqtab_nochim2 <- readRDS("02_DADA2/02_DADA2_2_exosapOut/seqtab_nochim.obj")
seqtab_nochim3 <- readRDS("02_DADA2/02_DADA2_3_rep_vol_testOut/seqtab_nochim.obj")

# Load sample data
sample_sheet <- read.csv("sampledata/eDNAseq_SampleSheet_All.csv")
sample_sheet1 <- read.csv("sampledata/RIR-027/RIR-027_SampleSheet_All.csv")
sample_sheet2 <- read.csv("sampledata/RIR-036/RIR-036_SampleSheet_All.csv")
sample_sheet3 <- read.csv("sampledata/RIR-031/RIR-031_SampleSheet_All.csv")

# Check data size
dim(sample_sheet)
dim(seqtab_nochim1); dim(seqtab_nochim2); dim(seqtab_nochim3)
sum(seqtab_nochim1); sum(seqtab_nochim2); sum(seqtab_nochim3)


# <-----------------------------------------------------> #
# Compile data
# <-----------------------------------------------------> #
# Extract sequences
seqs1 <- colnames(seqtab_nochim1)
seqs2 <- colnames(seqtab_nochim2)
seqs3 <- colnames(seqtab_nochim3)
seqs_all <- c(seqs1, seqs2, seqs3)
seqs1_ncol <- 1:length(seqs1)
seqs2_ncol <- (max(seqs1_ncol)+1):(max(seqs1_ncol)+length(seqs2))
seqs3_ncol <- (max(seqs2_ncol)+1):(max(seqs2_ncol)+length(seqs3))

# Prepare seqtab_nochim_all
seqtab_nochim_all <- matrix(0, ncol = length(seqs_all), nrow = nrow(sample_sheet))
rownames(seqtab_nochim_all) <- rownames(sample_sheet) <- sample_sheet$Sample_Name2
colnames(seqtab_nochim_all) <- seqs_all

# Fill empty samples with 0
seqtab_nochim_all[rownames(seqtab_nochim1),seqs1_ncol] <- seqtab_nochim1
seqtab_nochim_all[rownames(seqtab_nochim2),seqs2_ncol] <- seqtab_nochim2
seqtab_nochim_all[rownames(seqtab_nochim3),seqs3_ncol] <- seqtab_nochim3


# <-----------------------------------------------------> #
#       Merge identical sequence
# <-----------------------------------------------------> #
seqtab_nochim_all <- seqtab_nochim_all %>% t %>%
  rowsum(colnames(seqtab_nochim_all)) %>% t
dim(seqtab_nochim_all)
seqs_all2 <- colnames(seqtab_nochim_all)


# <-----------------------------------------------------> #
#                  Clustering by DECIPHER
# <-----------------------------------------------------> #
# Preparation
n_cores <- parallel::detectCores() # set to number of cpus/processors to use for the clustering
dna <- Biostrings::DNAStringSet(seqs_all2) # "seqs" is an object from DADA2 output
## Find clusters of ASVs to form the new OTUs
aln <- DECIPHER::AlignSeqs(dna, processors = n_cores)
aln_dist <- DECIPHER::DistanceMatrix(aln, processors = n_cores)
# use `cutoff = 0.03` for a 97% OTU 
clusters <- DECIPHER::IdClusters(aln_dist, method = "complete",
                                 processors = n_cores, type = "clusters",
                                 cutoff = 0.03, showPlot = F)
colnames(clusters) <- "cluster_id"
clusters$cluster_id <- factor(clusters$cluster_id, levels = unique(clusters$cluster_id))


# <-----------------------------------------------------> #
#             Make OTU-based phyloseq objects
# <-----------------------------------------------------> #
# Import DADA2 ASV output to phyloseq
ps <- phyloseq(otu_table(seqtab_nochim_all, taxa_are_rows = FALSE))
# Merge taxa in "ps" using cluster information
ps_otu <- speedyseq::merge_taxa_vec(ps, group = clusters$cluster_id)
otu_seqs <- colnames(otu_table(ps_otu))

# Quick taxa assignment
otu_only <- data.frame(taxa_id = sprintf("OTU%05d", 1:length(otu_seqs)),
                       seq = otu_seqs)
write.csv(otu_only, sprintf("%s/otu_only.csv", output_folder), row.names = T)

# Check corresponence between ASV sequences and OTU representative sequences
clusters$seq_sum <- colSums(seqtab_nochim_all)
clusters$asv_seq <- seqs_all2
asv_1st <- clusters %>% group_by(cluster_id) %>% summarize(otu_seqs = asv_seq[[1]])
clusters$otu_seq <- unlist(asv_1st[match(clusters$cluster_id, asv_1st$cluster_id), "otu_seqs"])
write.csv(clusters, sprintf("%s/cluster_summary.csv", output_folder), row.names = F)

# Save OTU table
otu_out <- as.matrix(c(rbind(sprintf(">OTU%05d", 1:length(otu_seqs)), otu_seqs)), ncol=1)
write.table(otu_out, sprintf("%s/OTU_seqs.fa", output_folder), col.names = FALSE, row.names = FALSE, quote = FALSE)

otu_mat <- as.data.frame(otu_table(ps_otu))
colnames(otu_mat) <- sprintf("OTU%05d", 1:length(otu_seqs))
write.csv(otu_mat, sprintf("%s/otu_table.csv", output_folder), row.names = T)


# <-----------------------------------------------------> #
#                     Save workspace
# <-----------------------------------------------------> #
# Save workspace
save(list = ls(all.names = TRUE),
     file = paste0(output_folder, "/", output_folder, ".RData"))

# Save session info
writeLines(capture.output(sessionInfo()),
           paste0("00_SessionInfo/", output_folder, "_", substr(Sys.time(), 1, 10), ".txt"))
