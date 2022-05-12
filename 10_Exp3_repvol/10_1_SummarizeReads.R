####
#### No.10.1 Summarize Experiment III: Replicate and volume
#### 2022.05.12 revision for Environmental DNA
#### R 4.1.2
####

# Set working directory
if(basename(getwd()) != "10_Exp3_repvol") setwd("10_Exp3_repvol")

# Set random seeds (for reproduction)
ran.seed <- 1234
set.seed(ran.seed)

# Load library and functions
library(tidyverse); packageVersion("tidyverse") # 1.3.1, 2021.10.16
library(phyloseq); packageVersion("phyloseq") # 1.38.0, 2021.11.18
library(cowplot); packageVersion("cowplot") # 1.1.1, 2021.6.13
library(ggsci); packageVersion("ggsci") # 2.9, 2021.6.13
library(RColorBrewer); packageVersion("RColorBrewer") # 1.1.3, 2022.5.12
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
ps_exp3 <- readRDS("../07_CompilePhyloseqOut/ps_exp3.obj")


# ----------------------------------------------- #
#    Subtract NC from ps_exp3
# ----------------------------------------------- #
# Prepare NC OTU table
ps_exp3_nc <- ps_exp3 %>% subset_samples(sample_nc == "nc")
## Extract cat_comb for all data
group_nc_df <- ps_exp3 %>% sample_data %>% as_tibble %>%
  select(sample_nc, rep_method)
## Extract NC reads
nc_reads <- ps_exp3_nc %>% otu_table %>% as.matrix %>% 
  .[match(group_nc_df$rep_method, ps_exp3_nc %>% sample_data %>% as_tibble %>% pull(rep_method)),]
## Subtract NC reads from sample reads
ps_exp3_flt1 <- ps_exp3
otu_table(ps_exp3_flt1) <- ((ps_exp3 %>% otu_table %>% as.matrix) - nc_reads) %>% 
  replace(., . < 0, 0) %>% otu_table(., taxa_are_rows = FALSE)
taxa_sums(ps_exp3_flt1)
ps_exp3_flt1 <- ps_exp3_flt1 %>% prune_taxa(taxa_sums(.) > 0, .)


# ----------------------------------------------- #
#         Visualize pattern: Reads
# ----------------------------------------------- #
get_palette <- colorRampPalette(brewer.pal(8, "Paired"))
# Visualize
target_rank <- "family"
ps_rename <- taxa_name_summarize(ps_exp3_flt1, target_rank, top_taxa_n = 10)
ps_m1 <- speedyseq::psmelt(ps_rename)
ps_m2 <- stats::aggregate(ps_m1$Abundance, by=list(ps_m1$Sample, ps_m1$family), "sum") # Summed up to make phylum sum
ps_m3 <- stats::aggregate(ps_m1$Abundance, by=list(ps_m1$Sample, ps_m1$rep_tax), "sum") # Summed up to make phylum sum
colnames(ps_m2) <- c("sample", target_rank, "abundance")
colnames(ps_m3) <- c("sample", "rep_tax", "abundance")
# Figures
f1 <- ggplot(ps_m2, aes(x = sample, y = abundance, group = family, fill = family)) +
  geom_bar(stat = "identity", colour = NA) +
  theme(axis.text.x = element_text(angle = -90, hjust = 1, vjust = 0.5, size = 6)) + 
  scale_fill_manual(values = get_palette(45)) +
  #scale_fill_igv() +
  xlab(NULL) + ylab("Sequence reads") +
  NULL
f2 <- ggplot(ps_m3, aes(x = sample, y = abundance, group = rep_tax, fill = rep_tax)) +
  geom_bar(stat = "identity", colour = NA) + theme(axis.text.x = element_text(angle = -90, hjust = 1, vjust = 0.5, size = 6)) + 
  xlab(NULL) + ylab("Sequence reads") +
  scale_fill_manual(values = get_palette(11)) +
  NULL

f3_1 <- subset_samples(ps_exp3_flt1, sample_nc == "sample" & site == "Sea_Nagahama") %>%
  subset_taxa(species != "STDseqs") %>%
  plot_richness(measures = "Observed", x = "replicate") +
  theme(axis.text.x = element_text(size = 6)) +
  facet_wrap(~ test_name + test_category, ncol = 4) +
  ylab("No. of OTUs") + panel_border() +
  xlab(NULL) + ylim(0,40) +
  NULL
sp_data <- f3_1$data[,c("site", "test_name", "test_category", "value")]
sp_data %>%
  group_by(site, test_name, test_category) %>%
  summarize(mean_div = mean(value),
            sd_div = sd(value))

f3_2 <- subset_samples(ps_exp3_flt1, sample_nc == "sample" & site == "STD_Mix") %>%
  subset_taxa(genus == "STDseqs") %>%
  plot_richness(measures = "Observed", x = "replicate") +
  theme(axis.text.x = element_text(size = 6)) +
  facet_wrap(~ test_name + test_category, ncol = 4) +
  ylab("No. of OTUs") + panel_border() +
  scale_y_continuous(breaks = seq(4,11), limits = c(4,11)) +
  geom_hline(yintercept = 10, linetype = 2) +
  xlab(NULL) +
  NULL
sp_data <- f3_2$data[,c("site", "test_name", "test_category", "value")]
sp_data %>%
  group_by(site, test_name, test_category) %>%
  summarize(mean_div = mean(value),
            sd_div = sd(value))


# ----------------------------------------------- #
#  Visualize pattern: Additional visualization
# ----------------------------------------------- #
# Extract top taxa > 100 reads
top_taxa <- taxa_names(ps_exp3_flt1)[taxa_sums(ps_exp3) > 100]
ps_top <- prune_taxa(top_taxa, ps_exp3_flt1)
f4 <- plot_richness(ps_top, measures = "Observed") + xlab(NULL) + theme(axis.text.x = element_text(size = 6))
f5 <- ggplot(ps_m1, aes(x = replicate, y = Abundance, fill = rep_tax)) +
  geom_bar(stat = "identity", colour = NA) +
  theme(axis.text.x = element_text(angle = -90, hjust = 1, vjust = 0.5, size = 6)) + 
  facet_wrap(~ site + test_name + test_category) +
  scale_fill_brewer(NULL, palette = "Paired") +
  xlab(NULL) + ylab("Sequence reads") + panel_border()

## Select or remove STD sequences
ps_s1 <- ps_m1 %>% filter(family == "STDseqs") 
ps_s2 <- ps_m1 %>% filter(family != "STDseqs")
f6 <- ggplot(ps_s1, aes(x = replicate, y = Abundance, fill = species)) +
  geom_bar(stat = "identity", colour = NA) +
  theme(axis.text.x = element_text(angle = -90, hjust = 1, vjust = 0.5, size = 6)) + 
  facet_wrap(~ site + test_name + test_category) +
  scale_fill_manual(values = get_palette(12)) +
  xlab(NULL) + ylab("Sequence reads") + panel_border()
f7 <- ggplot(ps_s2, aes(x = replicate, y = Abundance, fill = family)) +
  geom_bar(stat = "identity", colour = NA) +
  theme(axis.text.x = element_text(angle = -90, hjust = 1, vjust = 0.5, size = 6)) + 
  facet_wrap(~ site + test_name + test_category) +
  scale_fill_manual(values = get_palette(length(unique(ps_s2$family)))) +
  xlab(NULL) + ylab("Sequence reads") + panel_border() + theme(legend.position = "none")


# ----------------------------------------------- #
#         Relative abundance data
# ----------------------------------------------- #
# Check relative abundance
ps_rename2 <- taxa_name_summarize(ps_exp3_flt1, "family", top_taxa_n = 10)
ps_rel <- transform_sample_counts(ps_rename2, function(x) x/sum(x))
otu_table(ps_rel)[is.na(otu_table(ps_rel))] <- 0
ps_rel <- prune_taxa(taxa_sums(ps_rel) > 0, ps_rel)
ps_rel2 <- filter_taxa(ps_rel, function(x) mean(x) > 1e-02, TRUE)
ps_r1 <- speedyseq::psmelt(ps_rel)
ps_r2 <- stats::aggregate(ps_r1$Abundance, by=list(ps_r1$Sample, ps_r1$genus), "sum") # Summed up to make phylum sum
ps_r3 <- stats::aggregate(ps_r1$Abundance, by=list(ps_r1$Sample, ps_r1$rep_tax), "sum") # Summed up to make phylum sum
colnames(ps_r2) <- c("sample", "genus", "abundance")
colnames(ps_r3) <- c("sample", "rep_tax", "abundance")

# Figures
r1 <- ggplot(ps_r1, aes(x = replicate, y = Abundance, fill = rep_tax)) +
  geom_bar(stat = "identity", colour = NA) +
  theme(axis.text.x = element_text(angle = -90, hjust = 1, vjust = 0.5, size = 6)) + 
  facet_wrap(~ site + test_name + test_category) +
  scale_fill_brewer(NULL, palette = "Paired") +
  xlab(NULL) + ylab("Relative abundance") + panel_border()
r2 <- ps_r1 %>% filter(sample_nc == "sample") %>%
  ggplot(aes(x = replicate, y = Abundance, fill = rep_tax)) +
  geom_bar(stat = "identity", colour = NA) +
  theme(axis.text.x = element_text(angle = -90, hjust = 1, vjust = 0.5, size = 6)) + 
  facet_wrap(~ site + test_name + test_category) +
  scale_fill_brewer("genus", palette = "Paired") +
  xlab(NULL) + ylab("Relative abundance") + panel_border()
r3 <- ps_r1 %>% filter(family == "STDseqs") %>%
  ggplot(aes(x = replicate, y = Abundance, fill = species)) +
  geom_bar(stat = "identity", colour = NA) +
  theme(axis.text.x = element_text(angle = -90, hjust = 1, vjust = 0.5, size = 6)) + 
  facet_wrap(~ site + test_name + test_category) +
  scale_fill_manual(values = get_palette(12)) +
  xlab(NULL) + ylab("Relative abundance") + panel_border()
r4 <- ps_r1 %>% filter(family != "STDseqs") %>%
  ggplot(aes(x = replicate, y = Abundance, fill = family)) +
  geom_bar(stat = "identity", colour = NA) +
  theme(axis.text.x = element_text(angle = -90, hjust = 1, vjust = 0.5, size = 6)) + 
  facet_wrap(~ site + test_name + test_category) +
  scale_fill_manual(values = get_palette(length(unique(ps_r1$family)))) +
  xlab(NULL) + ylab("Sequence reads") + panel_border() + theme(legend.position = "none")


# ----------------------------------------------- #
#         Save data
# ----------------------------------------------- #
# Re-output data
write.csv(otu_table(ps_exp3), sprintf("%s/otu_table.csv", output_folder))
write.csv(sample_data(ps_exp3), sprintf("%s/sample_data.csv", output_folder))
write.csv(as.data.frame(tax_table(ps_exp3)), sprintf("%s/tax_table.csv", output_folder))
saveRDS(ps_exp3_flt1, sprintf("%s/ps_exp3_flt1.obj", output_folder))

# Save session info
writeLines(capture.output(sessionInfo()),
           paste0("../00_SessionInfo/", output_folder, "_", substr(Sys.time(), 1, 10), ".txt"))

