####
#### Rarefying sequence depth and compare diversity
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
ps_exp3_flt1 <- readRDS("10_1_SummarizeReadsOut/ps_exp3_flt1.obj")


# ----------------------------------------------- #
#    Rarefying sequence depth
# ----------------------------------------------- #
# Split samples and rarefy to the minimum reads in each site
## Exclude E03_T001_S001 as it contains only 304 sequences
ps_exp3_site1 <- prune_samples(sample_names(ps_exp3_flt1) != "E03_T001_S001", ps_exp3_flt1) %>%
  subset_samples(sample_nc == "sample" & site == "Sea_Nagahama") %>%
  subset_taxa(habitat != "freshwater" & habitat != "std") %>%
  rarefy_even_depth(rngseed = ran.seed, replace = FALSE, trimOTUs = FALSE)
ps_exp3_site2 <- prune_samples(sample_names(ps_exp3_flt1) != "E03_T001_S001", ps_exp3_flt1) %>%
  subset_samples(sample_nc == "sample" & site == "STD_Mix") %>%
  subset_taxa(habitat == "std") %>%
  rarefy_even_depth(rngseed = ran.seed, replace = FALSE, trimOTUs = FALSE)
sample_sums(ps_exp3_site1) # 15394 reads/sample
sample_sums(ps_exp3_site2) # 3025 reads/sample

# Replace OTUs of which relative abundance is < 0.01% with 0
ps_exp3_site1 <- transform_sample_counts(ps_exp3_site1, function(x){replace(x, x < ceiling(mean(sample_sums(ps_exp3_site1))*0.0001), 0)}) 
ps_exp3_site2 <- transform_sample_counts(ps_exp3_site2, function(x){replace(x, x < ceiling(mean(sample_sums(ps_exp3_site2))*0.0001), 0)}) 
mean(sample_sums(ps_exp3_site1)) # 15393.44 reads for each sample
mean(sample_sums(ps_exp3_site2)) # 3025 reads for each sample

# Re-merge phyloseq objects
ps_exp3_even <- merge_phyloseq(ps_exp3_site1, ps_exp3_site2) %>%
  prune_taxa(taxa_sums(.) > 0, .)
write.csv(as.data.frame(tax_table(ps_exp3_even)), sprintf("%s/ps_exp1_tax_table.csv", output_folder))


# ----------------------------------------------- #
#   Visualize pattern
# ----------------------------------------------- #
get_palette <- colorRampPalette(brewer.pal(8, "Paired"))
target_rank <- "family"
ps_rename <- taxa_name_summarize(ps_exp3_even, target_rank, top_taxa_n = 10)
# Replace rep_tax names of STD
std_names <- c("STD_MiFish02", "STD_MiFish04", "STD_MiFish05", "STD_MiFish08", "STD_MiFish09",
               "STD_Anguilla", "STD_cyprinus", "STD_Engraulis", "STD_lateolabrax", "STD_takifugu")
new_tax_table <- data.frame(tax_table(ps_rename))
for(i in 1:length(std_names)) {
  new_tax_table$rep_tax[new_tax_table$species == std_names[i]] <- std_names[i]
}
## Adjust names
non_std_names <- sort(unique(new_tax_table$rep_tax)[!unique(new_tax_table$rep_tax) %in% std_names])
new_tax_table$rep_tax <- factor(new_tax_table$rep_tax, levels = c(non_std_names, std_names))
tax_table(ps_rename) <- tax_table(as.matrix(new_tax_table))

# Melt phyloseq object
ps_m1 <- speedyseq::psmelt(ps_rename)
ps_m2 <- stats::aggregate(ps_m1$Abundance, by=list(ps_m1$Sample, ps_m1$family), "sum") # Summed up to make phylum sum
ps_m3 <- stats::aggregate(ps_m1$Abundance, by=list(ps_m1$Sample, ps_m1$rep_tax), "sum") # Summed up to make phylum sum
colnames(ps_m2) <- c("sample", target_rank, "abundance")
colnames(ps_m3) <- c("sample", "rep_tax", "abundance")
## Adjust name orders
ps_m1$rep_tax <- factor(ps_m1$rep_tax, levels = c(non_std_names, std_names))
ps_m3$rep_tax <- factor(ps_m3$rep_tax, levels = c(non_std_names, std_names))

# Figures
f1 <- ggplot(ps_m2, aes(x = sample, y = abundance, group = family, fill = family)) +
  geom_bar(stat = "identity", colour = NA) +
  theme(axis.text.x = element_text(angle = -90, hjust = 1, vjust = 0.5, size = 6)) + 
  scale_fill_manual(values = get_palette(45)) +
  xlab(NULL) + ylab("Sequence reads") +
  NULL
f2 <- ggplot(ps_m3, aes(x = sample, y = abundance, group = rep_tax, fill = rep_tax)) +
  geom_bar(stat = "identity", colour = NA) + theme(axis.text.x = element_text(angle = -90, hjust = 1, vjust = 0.5, size = 6)) + 
  xlab(NULL) + ylab("Sequence reads") +
  scale_fill_manual(values = get_palette(20)) +
  NULL


# ----------------------------------------------- #
#    Visualize pattern: Reads
# ----------------------------------------------- #
f3 <- ggplot(ps_m1, aes(x = replicate, y = Abundance, fill = rep_tax)) +
    geom_bar(stat = "identity", colour = NA) +
    theme(axis.text.x = element_text(angle = -90, hjust = 1, vjust = 0.5, size = 6)) + 
    facet_wrap(~ site + test_name + test_category) +
    scale_fill_manual(values = get_palette(20)) +
    xlab(NULL) + ylab("Sequence reads") + panel_border()

## Select or remove STD sequences
ps_s1 <- ps_m1 %>% filter(family == "STDseqs")
ps_s2 <- ps_m1 %>% filter(family != "STDseqs")
f4 <- ggplot(ps_s1, aes(x = replicate, y = Abundance, fill = species)) +
  geom_bar(stat = "identity", colour = NA) +
  theme(axis.text.x = element_text(angle = -90, hjust = 1, vjust = 0.5, size = 6)) + 
  facet_wrap(~ site + test_name + test_category) +
  scale_fill_manual(values = get_palette(20)) +
  xlab(NULL) + ylab("Sequence reads") + panel_border()
f5 <- ggplot(ps_s2, aes(x = replicate, y = Abundance, fill = family)) +
  geom_bar(stat = "identity", colour = NA) +
  theme(axis.text.x = element_text(angle = -90, hjust = 1, vjust = 0.5, size = 6)) + 
  facet_wrap(~ site + test_name + test_category) +
  scale_fill_manual(values = get_palette(length(unique(ps_s2$family)))) +
  xlab(NULL) + ylab("Sequence reads") + panel_border()


ps_rename2 <- taxa_name_summarize(ps_exp3_even %>% transform_sample_counts(function(x) x/sum(x)), target_rank, top_taxa_n = 10)
new_tax_table2 <- data.frame(tax_table(ps_rename2))
for(i in 1:length(std_names)) {
  new_tax_table2$rep_tax[new_tax_table2$species == std_names[i]] <- std_names[i]
}
## Adjust names
non_std_names2 <- sort(unique(new_tax_table2$rep_tax)[!unique(new_tax_table2$rep_tax) %in% std_names])
new_tax_table2$rep_tax <- factor(new_tax_table2$rep_tax, levels = c(non_std_names2, std_names))
tax_table(ps_rename2) <- tax_table(as.matrix(new_tax_table2))
# Melt phyloseq object
ps_f1 <- speedyseq::psmelt(ps_rename2)
ps_f1$rep_tax <- factor(ps_f1$rep_tax, levels = c(non_std_names2, std_names))
f6 <- ggplot(ps_f1, aes(x = replicate, y = Abundance, fill = rep_tax)) +
  geom_bar(stat = "identity", colour = NA) +
  theme(axis.text.x = element_text(angle = -90, hjust = 1, vjust = 0.5, size = 6)) + 
  facet_wrap(~ site + test_name + test_category) +
  scale_fill_manual(values = get_palette(20)) +
  xlab(NULL) + ylab("Sequence reads") + panel_border()


# ----------------------------------------------- #
#    Visualize pattern: Diversity
# ----------------------------------------------- #
f3_1 <- subset_samples(ps_exp3_even, sample_nc == "sample" & site == "Sea_Nagahama") %>%
  subset_taxa(species != "STDseqs") %>%
  plot_richness(measures = "Observed", x = "replicate") +
  theme(axis.text.x = element_text(size = 6)) +
  facet_wrap(~ test_name + test_category, ncol = 4) +
  ylab("No. of OTUs") + panel_border() +
  xlab(NULL) + ylim(0,40) +
  NULL
f3_1_data <- f3_1$data[,c("site", "test_name", "test_category", "value")]
f3_1_data %>%
  group_by(site, test_name, test_category) %>%
  summarize(mean_div = mean(value),
            sd_div = sd(value))

f3_2 <- subset_samples(ps_exp3_even, sample_nc == "sample" & site == "STD_Mix") %>%
  subset_taxa(genus == "STDseqs") %>%
  plot_richness(measures = "Observed", x = "replicate") +
  theme(axis.text.x = element_text(size = 6)) +
  facet_wrap(~ test_name + test_category, ncol = 4) +
  ylab("No. of OTUs") + panel_border() +
  scale_y_continuous(breaks = seq(4,11), limits = c(4,11)) +
  geom_hline(yintercept = 10, linetype = 2) +
  xlab(NULL) +
  NULL
f3_2_data <- f3_2$data[,c("site", "test_name", "test_category", "value")]
f3_2_data %>%
  group_by(site, test_name, test_category) %>%
  summarize(mean_div = mean(value),
            sd_div = sd(value))


# ----------------------------------------------- #
#    Visualize pattern: Diversity
# ----------------------------------------------- #
# Jitter plot
## Nagahama samples
f3_1_data$reaction_scale <- c(rep(1, 9), rep(2, 10), rep(4, 10), rep(8, 10))
f3_1_data$reaction_scale_fac <- factor(c(rep(1, 9), rep(2, 10), rep(4, 10), rep(8, 10)), levels = c(1,2,4,8))
d1_1 <- f3_1_data %>% ggplot(aes(x = log2(reaction_scale) + 1, y = value, fill = test_name)) + 
  geom_boxplot(aes(x = reaction_scale_fac, y = value), outlier.size = 0, outlier.shape = NA, alpha = 0.5) +
  scale_fill_manual(values = c("red3", "royalblue")) + 
  geom_smooth(color = "black", fill = "gray") +
  geom_jitter(width = 0.1, height = 0) +
  ylim(0, 40) +
  xlab("Reaction scale (no. of replicates or volume of DNA template)") +
  panel_border() +
  facet_wrap(. ~ test_name) +
  ylab("No. of detected fish species eDNA") +
  ggtitle("Effects of PCR reaction volume and replicates") +
  NULL
d1_2 <- f3_1_data %>% ggplot(aes(x = reaction_scale_fac, y = value, fill = test_name)) + 
  geom_boxplot(outlier.size = 0, outlier.shape = NA, alpha = 0.5) +
  scale_fill_manual(values = c("red3", "royalblue")) + 
  geom_jitter(width = 0.1, height = 0) +
  ylim(0, 40) +
  xlab("Reaction scale (no. of replicates or volume of DNA template)") +
  panel_border() +
  facet_wrap(. ~ test_name) +
  ylab("No. of detected fish species eDNA") +
  ggtitle("Effects of PCR reaction volume and replicates") +
  NULL

## STD_Mix samples
f3_2_data$reaction_scale <- c(rep(1, 10), rep(2, 10), rep(4, 10), rep(8, 10))
f3_2_data$reaction_scale_fac <- factor(c(rep(1, 10), rep(2, 10), rep(4, 10), rep(8, 10)), levels = c(1,2,4,8))
d2_1 <- f3_2_data %>% ggplot(aes(x = log2(reaction_scale) + 1, y = value, fill = test_name)) + 
  geom_boxplot(aes(x = reaction_scale_fac, y = value), outlier.size = 0, outlier.shape = NA, alpha = 0.5) +
  scale_fill_manual(values = c("red3", "royalblue")) + 
  geom_smooth(color = "black", fill = "gray") +
  geom_jitter(width = 0.1, height = 0) +
  ylim(0, 12) +
  geom_hline(yintercept = 10, linetype = 2) +
  xlab("Reaction scale (no. of replicates or volume of DNA template)") +
  panel_border() +
  facet_wrap(. ~ test_name) +
  ylab("No. of detected fish species eDNA") +
  ggtitle("Effects of PCR reaction volume and replicates") +
  NULL
d2_2 <- f3_2_data %>% ggplot(aes(x = reaction_scale_fac, y = value, fill = test_name)) + 
  geom_boxplot(outlier.size = 0, outlier.shape = NA, alpha = 0.5) +
  scale_fill_manual(values = c("red3", "royalblue")) + 
  geom_jitter(width = 0.1, height = 0) +
  ylim(0, 12) +
  geom_hline(yintercept = 10, linetype = 2) +
  xlab("Reaction scale (no. of replicates or volume of DNA template)") +
  panel_border() +
  facet_wrap(. ~ test_name) +
  ylab("No. of detected fish species eDNA") +
  ggtitle("Effects of PCR reaction volume and replicates") +
  NULL

## Brief statistics
seaN_glm <- glm(value ~ reaction_scale * test_name, data = f3_1_data, family = poisson(link = "log"))
stdM_glm <- glm(value ~ reaction_scale * test_name, data = f3_2_data, family = poisson(link = "log"))
summary(seaN_glm)
summary(stdM_glm)

seaN_glm1 <- glm(value ~ test_name, family = poisson(link = "log"),
                 data = f3_1_data %>% filter(reaction_scale_fac == 1))
seaN_glm2 <- glm(value ~ test_name, family = poisson(link = "log"),
                 data = f3_1_data %>% filter(reaction_scale_fac == 2))
seaN_glm4 <- glm(value ~ test_name, family = poisson(link = "log"),
                 data = f3_1_data %>% filter(reaction_scale_fac == 4))
seaN_glm8 <- glm(value ~ test_name, family = poisson(link = "log"),
                 data = f3_1_data %>% filter(reaction_scale_fac == 8))
summary(seaN_glm1); summary(seaN_glm2); summary(seaN_glm4); summary(seaN_glm8)

stdM_glm1 <- glm(value ~ test_name, family = poisson(link = "log"),
                 data = f3_2_data %>% filter(reaction_scale_fac == 1))
stdM_glm2 <- glm(value ~ test_name, family = poisson(link = "log"),
                 data = f3_2_data %>% filter(reaction_scale_fac == 2))
stdM_glm4 <- glm(value ~ test_name, family = poisson(link = "log"),
                 data = f3_2_data %>% filter(reaction_scale_fac == 4))
stdM_glm8 <- glm(value ~ test_name, family = poisson(link = "log"),
                 data = f3_2_data %>% filter(reaction_scale_fac == 8))
summary(stdM_glm1); summary(stdM_glm2); summary(stdM_glm4); summary(stdM_glm8)

# Multiplecomparison
#library(multcomp); packageVersion("multcomp") # 1.4.17, 2021.11.24
#multcomp::glht(seaN_glm, linfct = mcp(test_name = "Tukey"))


# ----------------------------------------------- #
#         Save data
# ----------------------------------------------- #
# Save figure objects
#dir.create("../FigCode"); dir.create("../FigCode/00_RawFigs")
fig_dir <- "../FigCode/00_RawFigs/"
saveRDS(list(f3_1, f3_2, f6), paste0(fig_dir, "10_2_Fig_Exp3_RepVolDiveristy.obj"))
saveRDS(list(d1_1, d1_2, d2_1, d2_2), paste0(fig_dir, "10_2_Fig_Exp3_DivBoxplot.obj"))

# Re-output data
write.csv(otu_table(ps_exp3_even), sprintf("%s/otu_table_exp3.csv", output_folder))
write.csv(sample_data(ps_exp3_even), sprintf("%s/sample_data_exp3.csv", output_folder))
write.csv(as.data.frame(tax_table(ps_exp3_even)), sprintf("%s/tax_table_exp3.csv", output_folder))
saveRDS(ps_exp3_even, sprintf("%s/ps_exp3_even.obj", output_folder))

# Save session info
writeLines(capture.output(sessionInfo()),
           paste0("../00_SessionInfo/", output_folder, "_", substr(Sys.time(), 1, 10), ".txt"))

