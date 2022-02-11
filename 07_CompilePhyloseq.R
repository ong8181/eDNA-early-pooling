####
#### No.7 Compile phyloseq objects
#### 2022.01.06 Ushio
#### R 4.1.2
####

# Set random seeds (for reproduction)
ran.seed <- 1234
set.seed(ran.seed)

# Load library and functions
library(tidyverse); packageVersion("tidyverse") # 1.3.1, 2021.10.16
library(phyloseq); packageVersion("phyloseq") # 1.38.0, 2021.11.18
library(cowplot); packageVersion("cowplot") # 1.1.1, 2021.6.13
library(ggsci); packageVersion("ggsci") # 2.9, 2021.6.13
library(RColorBrewer); packageVersion("RColorBrewer") # 1.1.2, 2021.6.13
theme_set(theme_cowplot())
source("functions_R/F02_HelperFunctions.R") # Helper function for visualization

# Generate output folder
od <- basename(rstudioapi::getSourceEditorContext()$path)
(output_folder <- paste0(str_sub(od, end = -3), "Out")); rm(od)
dir.create(output_folder)

# Load sample data
seqtab_data <- read.csv("03_OTUClusteringOut/otu_table.csv", row.names = 1)
sample_sheet <- read.csv("sampledata/eDNAseq_SampleSheet_All.csv")
tax_sheet <- read.csv("06_TaxaSTDcombineOut/claident_tax_revise.csv")
tax_seq <- read.csv("03_OTUClusteringOut/otu_only.csv", row.names = 1)
std_names <- c("STD_MiFish02", "STD_MiFish04", "STD_MiFish05", "STD_MiFish08", "STD_MiFish09",
               "STD_Anguilla", "STD_cyprinus", "STD_takifugu", "STD_lateolabrax", "STD_Engraulis")

# ----------------------------------------------- #
#         Compile data
# ----------------------------------------------- #
# Check structure
dim(seqtab_data)
dim(sample_sheet)
dim(tax_sheet); dim(tax_seq)

# Adjust row- and col-names
rownames(sample_sheet) <- rownames(seqtab_data) <- sample_sheet$Sample_Name2
tax_sheet$seq <- tax_seq[,2]
tax_sheet$seqlen <- nchar(tax_seq[,2])
colnames(seqtab_data) <- rownames(tax_sheet) <- paste0("OTU", str_sub(tax_seq[,1], start = 6, end = 8))
sample_sheet$site <- factor(sample_sheet$site, levels = c("Sea_Nagahama", "Sea_Otomi", "River_Seta", "STD_Mix"))


# ----------------------------------------------- #
#         Import to phyloseq
# ----------------------------------------------- #
# Import to phyloseq
ps_all0 <- ps_all <- phyloseq(otu_table(seqtab_data, taxa_are_rows = FALSE),
                              sample_data(sample_sheet),
                              tax_table(as.matrix(tax_sheet)))

# Merge STD sequences into one OTU
for (i in 1:length(std_names)) {
  merged_taxa_i <- taxa_names(subset_taxa(ps_all, species == std_names[i]))
  if(length(merged_taxa_i) > 1) {
    ps_all <- merge_taxa(ps_all, merged_taxa_i, archetype = merged_taxa_i[1])
    tax_table(ps_all)[merged_taxa_i[1],] <- tax_table(ps_all0)[merged_taxa_i[1],]
  }
}

# Re-load tax_table with habitat information
tax_sheet_w_habitat <- read.csv(sprintf("%s/tax_table_all_RM.csv", output_folder), row.names = 1)
rownames(tax_sheet_w_habitat) <- paste0("OTU", str_sub(rownames(tax_sheet_w_habitat), start = -3))
## Check tax_sheet contents
dim(tax_sheet_w_habitat); dim(tax_table(ps_all))
all(rownames(tax_sheet_w_habitat) == rownames(tax_table(ps_all)))

# Replace the original tax_sheet with tax_sheet_w_habitat
tax_table(ps_all) <- tax_table(as.matrix(tax_sheet_w_habitat))

# Extract fish species
ps_fish <- subset_taxa(ps_all, class == "Actinopteri")


# ----------------------------------------------- #
# Divide data into each experiment
# ----------------------------------------------- #
ps_exp1 <- subset_samples(ps_fish, Experiment_ID == "E01") %>% prune_taxa(taxa_sums(.) > 0, .)
ps_exp2 <- subset_samples(ps_fish, Experiment_ID == "E02") %>% prune_taxa(taxa_sums(.) > 0, .)
ps_exp3 <- subset_samples(ps_fish, Experiment_ID == "E03") %>% prune_taxa(taxa_sums(.) > 0, .)


# ----------------------------------------------- #
#         Save data
# ----------------------------------------------- #
# Re-output data
write.csv(otu_table(ps_all), sprintf("%s/otu_table_all.csv", output_folder))
write.csv(sample_data(ps_all), sprintf("%s/sample_data_all.csv", output_folder))
write.csv(as.data.frame(tax_table(ps_all)), sprintf("%s/tax_table_all.csv", output_folder))

# Save separated phyloseq objects
saveRDS(ps_all, sprintf("%s/ps_all.obj", output_folder))
saveRDS(ps_fish, sprintf("%s/ps_fish_all.obj", output_folder))
saveRDS(ps_exp1, sprintf("%s/ps_exp1.obj", output_folder))
saveRDS(ps_exp2, sprintf("%s/ps_exp2.obj", output_folder))
saveRDS(ps_exp3, sprintf("%s/ps_exp3.obj", output_folder))

# Save session info
writeLines(capture.output(sessionInfo()),
           paste0("00_SessionInfo/", output_folder, "_", substr(Sys.time(), 1, 10), ".txt"))


