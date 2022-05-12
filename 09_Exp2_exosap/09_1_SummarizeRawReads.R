####
#### No.9.1 Summarize Experiment II
#### R 4.1.2
####

# Set working directory
if(basename(getwd()) != "09_Exp2_exosap") setwd("09_Exp2_exosap")

# Set random seeds (for reproduction)
ran.seed <- 1234
set.seed(ran.seed)

# Load library and functions
library(tidyverse); packageVersion("tidyverse") # 1.3.1, 2021.10.16
library(phyloseq); packageVersion("phyloseq") # 1.38.0, 2021.11.18
library(cowplot); packageVersion("cowplot") # 1.1.1, 2021.6.13
library(ggsci); packageVersion("ggsci") # 2.9, 2021.6.13
library(RColorBrewer); packageVersion("RColorBrewer") # 1.1.3, 2021.6.13
theme_set(theme_cowplot())
source("../functions_R/F02_HelperFunctions.R") # Helper function for visualization

# Generate output folder
od <- basename(rstudioapi::getSourceEditorContext()$path)
(output_folder <- paste0(str_sub(od, end = -3), "Out")); rm(od)
dir.create(output_folder)


# ----------------------------------------------- #
#    Load phyloseq object
# ----------------------------------------------- #
ps_all <- readRDS("../07_CompilePhyloseqOut/ps_all.obj")
ps_exp2 <- readRDS("../07_CompilePhyloseqOut/ps_exp2.obj")


# ----------------------------------------------- #
#         Visualize pattern: Reads
# ----------------------------------------------- #
get_palette <- colorRampPalette(brewer.pal(8, "Paired"))
# Visualize
target_rank <- "family"
ps_rename <- taxa_name_summarize(ps_exp2, target_rank, top_taxa_n = 10)
ps_m1 <- speedyseq::psmelt(ps_rename)
ps_m2 <- stats::aggregate(ps_m1$Abundance, by=list(ps_m1$Sample, ps_m1$family), "sum") # Summed up to make phylum sum
ps_m3 <- stats::aggregate(ps_m1$Abundance, by=list(ps_m1$Sample, ps_m1$rep_tax), "sum") # Summed up to make phylum sum
colnames(ps_m2) <- c("sample", target_rank, "abundance")
colnames(ps_m3) <- c("sample", "rep_tax", "abundance")
# Figures
f1 <- ggplot(ps_m2, aes(x = sample, y = abundance, group = family, fill = family)) +
  geom_bar(stat = "identity", colour = NA) +
  theme(axis.text.x = element_text(angle = -90, hjust = 1, vjust = 0.5, size = 6)) + 
  scale_fill_manual(values = get_palette(29)) +
  xlab(NULL) + ylab("Sequence reads") +
  NULL
f2 <- ggplot(ps_m3, aes(x = sample, y = abundance, group = rep_tax, fill = rep_tax)) +
  geom_bar(stat = "identity", colour = NA) + theme(axis.text.x = element_text(angle = -90, hjust = 1, vjust = 0.5, size = 6)) + 
  xlab(NULL) + ylab("Sequence reads") +
  scale_fill_manual(values = get_palette(11)) +
  NULL


# ----------------------------------------------- #
#     Normalize reads against "sample" well
# ----------------------------------------------- #
## Normalization using total reads of the sample well
ps_exp2_cat01 <- ps_exp2 %>% subset_samples(test_category == "Exo_005min_4C") %>%
  transform_sample_counts(fun = function(x) x/max(sample_sums(.)))
ps_exp2_cat02 <- ps_exp2 %>% subset_samples(test_category == "Exo_030min_4C") %>%
  transform_sample_counts(fun = function(x) x/max(sample_sums(.)))
ps_exp2_cat03 <- ps_exp2 %>% subset_samples(test_category == "Exo_120min_4C") %>%
  transform_sample_counts(fun = function(x) x/max(sample_sums(.)))
ps_exp2_cat04 <- ps_exp2 %>% subset_samples(test_category == "NoExo_005min_4C") %>%
  transform_sample_counts(fun = function(x) x/max(sample_sums(.)))
ps_exp2_cat05 <- ps_exp2 %>% subset_samples(test_category == "NoExo_030min_4C") %>%
  transform_sample_counts(fun = function(x) x/max(sample_sums(.)))
ps_exp2_cat06 <- ps_exp2 %>% subset_samples(test_category == "NoExo_120min_4C") %>%
  transform_sample_counts(fun = function(x) x/max(sample_sums(.)))
ps_exp2_cat07 <- ps_exp2 %>% subset_samples(test_category == "Exo_005min_RT") %>%
  transform_sample_counts(fun = function(x) x/max(sample_sums(.)))
ps_exp2_cat08 <- ps_exp2 %>% subset_samples(test_category == "Exo_030min_RT") %>%
  transform_sample_counts(fun = function(x) x/max(sample_sums(.)))
ps_exp2_cat09 <- ps_exp2 %>% subset_samples(test_category == "Exo_120min_RT") %>%
  transform_sample_counts(fun = function(x) x/max(sample_sums(.)))
ps_exp2_cat10 <- ps_exp2 %>% subset_samples(test_category == "NoExo_005min_RT") %>%
  transform_sample_counts(fun = function(x) x/max(sample_sums(.)))
ps_exp2_cat11 <- ps_exp2 %>% subset_samples(test_category == "NoExo_030min_RT") %>%
  transform_sample_counts(fun = function(x) x/max(sample_sums(.)))
ps_exp2_cat12 <- ps_exp2 %>% subset_samples(test_category == "NoExo_120min_RT") %>%
  transform_sample_counts(fun = function(x) x/max(sample_sums(.)))

# Re-merge phyloseq objects
ps_exp2_norm <- merge_phyloseq(ps_exp2_cat01, ps_exp2_cat02, ps_exp2_cat03, ps_exp2_cat04,
                               ps_exp2_cat05, ps_exp2_cat06, ps_exp2_cat07, ps_exp2_cat08,
                               ps_exp2_cat09, ps_exp2_cat10, ps_exp2_cat11, ps_exp2_cat12)

## Normalization using total reads of each OTU of the sample well
column_wise_normalization <- function(ps_object){
  ps_subset <- ps_object
  new_otu_table <- matrix(NaN, ncol = dim(otu_table(ps_subset))[2], nrow = dim(otu_table(ps_subset))[1])
  for(i in 1:dim(otu_table(ps_subset))[2]) {
    otu_read_sample_i <- as.numeric(otu_table(ps_subset)[1,i])
    new_otu_table[,i] <- otu_table(ps_subset)[,i]/otu_read_sample_i
  }
  colnames(new_otu_table) <- taxa_names(ps_subset)
  rownames(new_otu_table) <- sample_names(ps_subset)
  otu_table(ps_subset) <- otu_table(new_otu_table, taxa_are_rows = FALSE)
  return(ps_subset)
}

ps_exp2_cat21 <- column_wise_normalization(subset_samples(ps_exp2, test_category == "Exo_005min_4C"))
ps_exp2_cat22 <- column_wise_normalization(subset_samples(ps_exp2, test_category == "Exo_030min_4C"))
ps_exp2_cat23 <- column_wise_normalization(subset_samples(ps_exp2, test_category == "Exo_120min_4C"))
ps_exp2_cat24 <- column_wise_normalization(subset_samples(ps_exp2, test_category == "NoExo_005min_4C"))
ps_exp2_cat25 <- column_wise_normalization(subset_samples(ps_exp2, test_category == "NoExo_030min_4C"))
ps_exp2_cat26 <- column_wise_normalization(subset_samples(ps_exp2, test_category == "NoExo_120min_4C"))
ps_exp2_cat27 <- column_wise_normalization(subset_samples(ps_exp2, test_category == "Exo_005min_RT"))
ps_exp2_cat28 <- column_wise_normalization(subset_samples(ps_exp2, test_category == "Exo_030min_RT"))
ps_exp2_cat29 <- column_wise_normalization(subset_samples(ps_exp2, test_category == "Exo_120min_RT"))
ps_exp2_cat30 <- column_wise_normalization(subset_samples(ps_exp2, test_category == "NoExo_005min_RT"))
ps_exp2_cat31 <- column_wise_normalization(subset_samples(ps_exp2, test_category == "NoExo_030min_RT"))
ps_exp2_cat32 <- column_wise_normalization(subset_samples(ps_exp2, test_category == "NoExo_120min_RT"))

# Re-merge phyloseq objects
ps_exp2_norm2 <- merge_phyloseq(ps_exp2_cat21, ps_exp2_cat22, ps_exp2_cat23, ps_exp2_cat24,
                                ps_exp2_cat25, ps_exp2_cat26, ps_exp2_cat27, ps_exp2_cat28,
                                ps_exp2_cat29, ps_exp2_cat30, ps_exp2_cat31, ps_exp2_cat32)


# ----------------------------------------------- #
#         Save data
# ----------------------------------------------- #
# Re-output data
write.csv(otu_table(ps_exp2), sprintf("%s/otu_table.csv", output_folder))
write.csv(sample_data(ps_exp2), sprintf("%s/sample_data.csv", output_folder))
write.csv(as.data.frame(tax_table(ps_exp2)), sprintf("%s/tax_table.csv", output_folder))
saveRDS(ps_exp2_norm, sprintf("%s/ps_exp2_norm.obj", output_folder))
saveRDS(ps_exp2_norm2, sprintf("%s/ps_exp2_norm2.obj", output_folder))

# Save session info
writeLines(capture.output(sessionInfo()),
           paste0("../00_SessionInfo/", output_folder, "_", substr(Sys.time(), 1, 10), ".txt"))

