####
#### No.8.1 Summarize Experiment I
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
source("../functions_R/F02_HelperFunctions.R") # Helper function for visualization

# Generate output folder
od <- basename(rstudioapi::getSourceEditorContext()$path)
(output_folder <- paste0(str_sub(od, end = -3), "Out")); rm(od)
dir.create(output_folder)


# ----------------------------------------------- #
#    Load phyloseq object
# ----------------------------------------------- #
ps_all <- readRDS("../07_CompilePhyloseqOut/ps_all.obj")
ps_fish <- readRDS("../07_CompilePhyloseqOut/ps_fish_all.obj")
ps_exp1 <- readRDS("../07_CompilePhyloseqOut/ps_exp1.obj")


# ----------------------------------------------- #
#    Visualize pattern: Reads
# ----------------------------------------------- #
get_palette <- colorRampPalette(brewer.pal(8, "Paired"))

target_rank <- "family"
ps_rename <- taxa_name_summarize(ps_exp1, target_rank, top_taxa_n = 10)
ps_m1 <- speedyseq::psmelt(ps_rename)
ps_m2 <- stats::aggregate(ps_m1$Abundance, by=list(ps_m1$Sample, ps_m1$family), "sum") # Summed up to make phylum sum
ps_m3 <- stats::aggregate(ps_m1$Abundance, by=list(ps_m1$Sample, ps_m1$rep_tax), "sum") # Summed up to make phylum sum
colnames(ps_m2) <- c("sample", target_rank, "abundance")
colnames(ps_m3) <- c("sample", "rep_tax", "abundance")
# Figures
(f1 <- ggplot(ps_m2, aes(x = sample, y = abundance, group = family, fill = family)) +
    geom_bar(stat = "identity", colour = NA) +
    theme(axis.text.x = element_text(angle = -90, hjust = 1, vjust = 0.5, size = 6)) + 
    scale_fill_manual(values = get_palette(45)) +
    xlab(NULL) + ylab("Sequence reads") +
    NULL)
(f2 <- ggplot(ps_m3, aes(x = sample, y = abundance, group = rep_tax, fill = rep_tax)) +
    geom_bar(stat = "identity", colour = NA) + theme(axis.text.x = element_text(angle = -90, hjust = 1, vjust = 0.5, size = 6)) + 
    xlab(NULL) + ylab("Sequence reads") +
    scale_fill_manual(values = get_palette(11)) +
    NULL)

# ----------------------------------------------- #
#    Visualize pattern: Reads
# ----------------------------------------------- #
(f3 <- ggplot(ps_m1, aes(x = replicate, y = Abundance, fill = rep_tax)) +
    geom_bar(stat = "identity", colour = NA) +
    theme(axis.text.x = element_text(angle = -90, hjust = 1, vjust = 0.5, size = 6)) + 
    facet_wrap(~ site + index_method + enzyme, ncol = 3) +
    scale_fill_brewer("family", palette = "Paired") +
    xlab(NULL) + ylab("Sequence reads") + panel_border())

## Select or remove STD sequences
ps_s1 <- ps_m1 %>% filter(family == "STDseqs") 
ps_s2 <- ps_m1 %>% filter(family != "STDseqs")
(f4 <- ggplot(ps_s1, aes(x = replicate, y = Abundance, fill = species)) +
    geom_bar(stat = "identity", colour = NA) +
    theme(axis.text.x = element_text(angle = -90, hjust = 1, vjust = 0.5, size = 6)) + 
    facet_wrap(~ site + index_method + enzyme, ncol = 3) +
    scale_fill_manual(values = get_palette(12)) +
    xlab(NULL) + ylab("Sequence reads") + panel_border())


# ----------------------------------------------- #
#    Visualize pattern: Diversity
# ----------------------------------------------- #
# Exclude apparent cross-contamination (e.g., STD reads from samples)
ps_exp1_natural <- subset_samples(ps_exp1, site != "STD_Mix") %>% subset_taxa(genus != "STDseqs")
ps_exp1_std <- subset_samples(ps_exp1, site == "STD_Mix") %>% subset_taxa(genus == "STDseqs")
ps_exp1_decontam <- merge_phyloseq(ps_exp1_natural, ps_exp1_std)

f5 <- plot_richness(subset_samples(ps_exp1_decontam, sample_nc == "sample"), measures = "Observed", x = "enzyme") +
  theme(axis.text.x = element_text(size = 6)) +
  facet_wrap(~ site + index_method, scales = "free_x",  ncol = 2) +
  scale_x_discrete(expand = c(0.8,0.5)) +
  ylab("No. of OTUs") +
  xlab(NULL) +
  NULL
f5$layers <- f5$layers[-1]
f5 <- f5 + geom_boxplot(width = 0.5, outlier.shape = NULL, outlier.size = 0, outlier.colour = "white") + 
  geom_jitter(size = 1.5, width = 0.1) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
sp_data <- f5$data[,c("site", "enzyme", "index_method", "value")]
sp_data %>%
  group_by(site, enzyme, index_method) %>%
  summarize(mean_div = mean(value),
            sd_div = sd(value))


# ----------------------------------------------- #
#         Relative abundance data
# ----------------------------------------------- #
# Check relative abundance
ps_rename2 <- taxa_name_summarize(ps_exp1, "family", top_taxa_n = 10)
ps_rel <- transform_sample_counts(ps_rename2, function(x) x/sum(x))
otu_table(ps_rel)[is.na(otu_table(ps_rel))] <- 0
ps_rel <- prune_taxa(taxa_sums(ps_rel) > 0, ps_rel)
ps_rel2 <- filter_taxa(ps_rel, function(x) mean(x) > 1e-02, TRUE)
ps_r1 <- speedyseq::psmelt(ps_rel)
ps_r2 <- stats::aggregate(ps_r1$Abundance, by=list(ps_r1$Sample, ps_r1$family), "sum") # Summed up to make phylum sum
ps_r3 <- stats::aggregate(ps_r1$Abundance, by=list(ps_r1$Sample, ps_r1$rep_tax), "sum") # Summed up to make phylum sum
colnames(ps_r2) <- c("sample", target_rank, "abundance")
colnames(ps_r3) <- c("sample", "rep_tax", "abundance")

# Figures
(r1 <- ggplot(ps_r1, aes(x = replicate, y = Abundance, fill = rep_tax)) +
    geom_bar(stat = "identity", colour = NA) +
    theme(axis.text.x = element_text(angle = -90, hjust = 1, vjust = 0.5, size = 6)) + 
    facet_wrap(~ site + index_method + enzyme, ncol = 3) +
    scale_fill_brewer("family", palette = "Paired") +
    xlab(NULL) + ylab("Relative abundance") + panel_border())
(r2 <- ps_r1 %>% filter(sample_nc == "sample") %>%
    ggplot(aes(x = replicate, y = Abundance, fill = rep_tax)) +
    geom_bar(stat = "identity", colour = NA) +
    theme(axis.text.x = element_text(angle = -90, hjust = 1, vjust = 0.5, size = 6)) + 
    facet_wrap(~ site + index_method + enzyme, ncol = 3) +
    scale_fill_brewer("family", palette = "Paired") +
    xlab(NULL) + ylab("Relative abundance") + panel_border())
(r3 <- ps_r1 %>% filter(family == "STDseqs") %>%
    ggplot(aes(x = replicate, y = Abundance, fill = species)) +
    geom_bar(stat = "identity", colour = NA) +
    theme(axis.text.x = element_text(angle = -90, hjust = 1, vjust = 0.5, size = 6)) + 
    facet_wrap(~ site + index_method + enzyme, ncol = 3) +
    scale_fill_manual(values = get_palette(12)) +
    xlab(NULL) + ylab("Relative abundance") + panel_border())
(r4 <- ps_r1 %>% filter(query != "STDseqs") %>%
    ggplot(aes(x = replicate, y = Abundance, fill = family)) +
    geom_bar(stat = "identity", colour = NA) +
    theme(axis.text.x = element_text(angle = -90, hjust = 1, vjust = 0.5, size = 6)) + 
    facet_wrap(~ site + index_method + enzyme, ncol = 3) +
    scale_fill_manual(values = get_palette(length(unique(ps_r1$family)))) +
    xlab(NULL) + ylab("Sequence reads") + panel_border() + theme(legend.position = "none"))


# ----------------------------------------------- #
#         Save figures
# ----------------------------------------------- #
#pdf(file = sprintf("%s/Summary.pdf", output_folder), width = 18, height = 6)
#plot_grid(f2, plot_richness(ps_exp1, measures = "Observed"),
#          ncol = 2, rel_widths = c(1,0.7))
#dev.off()
#ggsave(file = sprintf("%s/SummaryReads.pdf", output_folder), plot = f3, width = 12, height = 12)
#ggsave(file = sprintf("%s/SummaryDiversity.pdf", output_folder), plot = f5, width = 12, height = 12)
#ggsave(file = sprintf("%s/STDReads.pdf", output_folder), plot = f4, width = 12, height = 12)
#ggsave(file = sprintf("%s/SummaryRelative.pdf", output_folder), plot = r1, width = 12, height = 12)
#ggsave(file = sprintf("%s/SummaryRelative_noNC.pdf", output_folder), plot = r2, width = 12, height = 12)
#ggsave(file = sprintf("%s/STDRelative.pdf", output_folder), plot = r3, width = 12, height = 12)


# ----------------------------------------------- #
#         Save data
# ----------------------------------------------- #
# Re-output data
write.csv(otu_table(ps_exp1), sprintf("%s/otu_table_exp1.csv", output_folder))
write.csv(sample_data(ps_exp1), sprintf("%s/sample_data_exp1.csv", output_folder))
write.csv(as.data.frame(tax_table(ps_exp1)), sprintf("%s/tax_table_exp1.csv", output_folder))

# Save session info
writeLines(capture.output(sessionInfo()),
           paste0("../00_SessionInfo/", output_folder, "_", substr(Sys.time(), 1, 10), ".txt"))


# ----------------------------------------------- #
#   Save figures for publication
# ----------------------------------------------- #
dir.create("../FigCode")
dir.create("../FigCode/00_RawFigs")
fig_dir <- "../FigCode/00_RawFigs/"
# Relabeling the facet strips
saveRDS(f3, paste0(fig_dir, "8_1_Fig_Exp1_SummaryReads.obj"))

