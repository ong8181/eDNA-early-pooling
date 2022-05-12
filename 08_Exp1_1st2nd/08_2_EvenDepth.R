####
#### Rarefying sequence depth and compare composition and diversity
#### 2022.05.12 revision for Environmental DNA
#### R 4.1.2
####

# Set working directory
if(basename(getwd()) != "08_Exp1_1st2nd") setwd("08_Exp1_1st2nd")

# Set random seeds (for reproduction)
ran.seed <- 1234
set.seed(ran.seed)

# Load library and functions
library(tidyverse); packageVersion("tidyverse") # 1.3.1, 2021.10.16
library(phyloseq); packageVersion("phyloseq") # 1.38.0, 2021.11.18
library(cowplot); packageVersion("cowplot") # 1.1.1, 2021.6.13
library(ggsci); packageVersion("ggsci") # 2.9, 2021.6.13
library(RColorBrewer); packageVersion("RColorBrewer") # 1.1.3, 2022.5.5
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
## Use NC subtracted OTU table
ps_exp1_flt1 <- readRDS("08_1_SummarizeRawReadsOut/ps_exp1_flt1.obj")


# ----------------------------------------------- #
#    Rarefying sequence depth
# ----------------------------------------------- #
# Split samples and rarefy to the minimum reads in each site
ps_exp1_site1 <- subset_samples(ps_exp1_flt1, sample_nc == "sample" & site == "Sea_Nagahama") %>%
  subset_taxa(habitat != "freshwater") %>% subset_taxa(habitat != "std") %>%
  rarefy_even_depth(rngseed = ran.seed, replace = FALSE, trimOTUs = FALSE)
ps_exp1_site2 <- subset_samples(ps_exp1_flt1, sample_nc == "sample" & site == "Sea_Otomi") %>%
  subset_taxa(habitat != "freshwater") %>% subset_taxa(habitat != "std") %>%
  rarefy_even_depth(rngseed = ran.seed, replace = FALSE, trimOTUs = FALSE)
ps_exp1_site3 <- subset_samples(ps_exp1_flt1, sample_nc == "sample" & site == "River_Seta") %>%
  subset_taxa(habitat != "marine") %>% subset_taxa(habitat != "std") %>%
  rarefy_even_depth(rngseed = ran.seed, replace = FALSE, trimOTUs = FALSE)
ps_exp1_site4 <- subset_samples(ps_exp1_flt1, sample_nc == "sample" & site == "STD_Mix") %>%
  subset_taxa(habitat == "std") %>%
  rarefy_even_depth(rngseed = ran.seed, replace = FALSE, trimOTUs = FALSE)
sample_sums(ps_exp1_site1) # 20750 reads for each sample
sample_sums(ps_exp1_site2) # 8351 reads for each sample
sample_sums(ps_exp1_site3) # 10148 reads for each sample
sample_sums(ps_exp1_site4) # 9574 reads for each sample

# Replace OTUs of which relative abundance is < 0.01% with 0
ps_exp1_site1 <- transform_sample_counts(ps_exp1_site1, function(x){replace(x, x < ceiling(mean(sample_sums(ps_exp1_site1))*0.0001), 0)}) 
ps_exp1_site2 <- transform_sample_counts(ps_exp1_site2, function(x){replace(x, x < ceiling(mean(sample_sums(ps_exp1_site2))*0.0001), 0)}) 
ps_exp1_site3 <- transform_sample_counts(ps_exp1_site3, function(x){replace(x, x < ceiling(mean(sample_sums(ps_exp1_site3))*0.0001), 0)}) 
ps_exp1_site4 <- transform_sample_counts(ps_exp1_site4, function(x){replace(x, x < ceiling(mean(sample_sums(ps_exp1_site4))*0.0001), 0)}) 
mean(sample_sums(ps_exp1_site1)) # 20747.67 reads for each sample
mean(sample_sums(ps_exp1_site2)) # 8351 reads for each sample
mean(sample_sums(ps_exp1_site3)) # 10147.87 reads for each sample
mean(sample_sums(ps_exp1_site4)) # 9574 reads for each sample

# Re-merge phyloseq objects
ps_exp1_even <- merge_phyloseq(ps_exp1_site1, ps_exp1_site2, ps_exp1_site3, ps_exp1_site4) %>%
  prune_taxa(taxa_sums(.) > 0, .)
write.csv(as.data.frame(tax_table(ps_exp1_even)), sprintf("%s/ps_exp1_tax_table.csv", output_folder))


# ----------------------------------------------- #
#   Visualize pattern
# ----------------------------------------------- #
get_palette <- colorRampPalette(brewer.pal(8, "Paired"))
target_rank <- "family"
ps_rename <- taxa_name_summarize(ps_exp1_even, target_rank, top_taxa_n = 10)
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
  xlab(NULL) + ylab("Sequence reads") +
  NULL
f2 <- ggplot(ps_m3, aes(x = sample, y = abundance, group = rep_tax, fill = rep_tax)) +
  geom_bar(stat = "identity", colour = NA) + theme(axis.text.x = element_text(angle = -90, hjust = 1, vjust = 0.5, size = 6)) + 
  xlab(NULL) + ylab("Sequence reads") +
  scale_fill_manual(values = get_palette(11)) +
  NULL


# ----------------------------------------------- #
#    Visualize pattern: Reads
# ----------------------------------------------- #
f3 <- ggplot(ps_m1, aes(x = replicate, y = Abundance, fill = rep_tax)) +
  geom_bar(stat = "identity", colour = NA) +
  theme(axis.text.x = element_text(angle = -90, hjust = 1, vjust = 0.5, size = 6)) + 
  facet_wrap(~ site + index_method + enzyme, ncol = 3) +
  scale_fill_brewer("family", palette = "Paired") +
  xlab(NULL) + ylab("Sequence reads") + panel_border()

## Select or remove STD sequences
ps_s1 <- ps_m1 %>% filter(family == "STDseqs") 
ps_s2 <- ps_m1 %>% filter(family != "STDseqs")
f4 <- ggplot(ps_s1, aes(x = replicate, y = Abundance, fill = species)) +
  geom_bar(stat = "identity", colour = NA) +
  theme(axis.text.x = element_text(angle = -90, hjust = 1, vjust = 0.5, size = 6)) + 
  facet_wrap(~ site + index_method + enzyme, ncol = 3) +
  scale_fill_manual(values = get_palette(12)) +
  xlab(NULL) + ylab("Sequence reads") + panel_border()
f5 <- ggplot(ps_s2, aes(x = replicate, y = Abundance, fill = family)) +
  geom_bar(stat = "identity", colour = NA) +
  theme(axis.text.x = element_text(angle = -90, hjust = 1, vjust = 0.5, size = 6)) + 
  facet_wrap(~ site + index_method + enzyme, ncol = 3) +
  scale_fill_manual(values = get_palette(length(unique(ps_s2$family)))) +
  xlab(NULL) + ylab("Sequence reads") + panel_border()


# ----------------------------------------------- #
#    Visualize pattern: Diversity
# ----------------------------------------------- #
# Exclude apparent cross-contamination (e.g., STD reads from samples)
ps_exp1_natural <- subset_samples(ps_exp1_even, site != "STD_Mix") %>% subset_taxa(genus != "STDseqs")
ps_exp1_std <- subset_samples(ps_exp1_even, site == "STD_Mix") %>% subset_taxa(genus == "STDseqs")
ps_exp1_decontam <- merge_phyloseq(ps_exp1_natural, ps_exp1_std)

f6 <- plot_richness(subset_samples(ps_exp1_decontam, sample_nc == "sample"), measures = "Observed", x = "test_category") +
  theme(axis.text.x = element_text(size = 8, angle = -90, hjust = 0)) +
  facet_wrap(~ site, scales = "free_x",  ncol = 2) +
  #scale_x_discrete(expand = c(0.8,0.5)) +
  ylab("No. of OTUs") +
  xlab(NULL) +
  NULL
f6$layers <- f6$layers[-1]
f6 <- f6 + geom_boxplot(width = 0.8, outlier.shape = NULL, outlier.size = 0, outlier.colour = "white") + 
  geom_jitter(size = 1.5, width = 0.1) + ylim(0, 35) + panel_border() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5)) +
  NULL
sp_data <- f6$data[,c("site", "enzyme", "index_method", "value")]
sp_data %>%
  group_by(site, enzyme, index_method) %>%
  summarize(mean_div = mean(value),
            sd_div = sd(value))

## GLM assuming Poisson distributions
sp_data$method <- paste0(sp_data$index_method, "_", sp_data$enzyme)
### Sea_Nagahama
m1_1 <- sp_data %>% filter(site == "Sea_Nagahama") %>%
  glm(value ~ method, data = ., family = poisson(link = "log"))
m1_2 <- sp_data %>% filter(site == "Sea_Nagahama") %>%
  glm(value ~ 1, data = ., family = poisson(link = "log"))
anova(m1_2, m1_1, test="Chisq")
### Sea_Otomi
m2_1 <- sp_data %>% filter(site == "Sea_Otomi") %>%
  glm(value ~ method, data = ., family = poisson(link = "log"))
m2_2 <- sp_data %>% filter(site == "Sea_Otomi") %>%
  glm(value ~ 1, data = ., family = poisson(link = "log"))
anova(m2_2, m2_1, test="Chisq")
### River_Seta
m3_1 <- sp_data %>% filter(site == "River_Seta") %>%
  glm(value ~ method, data = ., family = poisson(link = "log"))
m3_2 <- sp_data %>% filter(site == "River_Seta") %>%
  glm(value ~ 1, data = ., family = poisson(link = "log"))
anova(m3_2, m3_1, test="Chisq")
### River_Seta
m4_1 <- sp_data %>% filter(site == "STD_Mix") %>%
  glm(value ~ method, data = ., family = poisson(link = "log"))
m4_2 <- sp_data %>% filter(site == "STD_Mix") %>%
  glm(value ~ 1, data = ., family = poisson(link = "log"))
anova(m4_2, m4_1, test="Chisq")


# ----------------------------------------------- #
#         Relative abundance data
# ----------------------------------------------- #
ps_rename2 <- taxa_name_summarize(ps_exp1_even, "family", top_taxa_n = 10)
# Replace rep_tax names of STD
std_names <- c("STD_MiFish02", "STD_MiFish04", "STD_MiFish05", "STD_MiFish08", "STD_MiFish09",
               "STD_Anguilla", "STD_cyprinus", "STD_Engraulis", "STD_lateolabrax", "STD_takifugu")
new_tax_table <- data.frame(tax_table(ps_rename2))
for(i in 1:length(std_names)) {
  new_tax_table$rep_tax[new_tax_table$species == std_names[i]] <- std_names[i]
}
## Adjust names
non_std_names <- sort(unique(new_tax_table$rep_tax)[!unique(new_tax_table$rep_tax) %in% std_names])
new_tax_table$rep_tax <- factor(new_tax_table$rep_tax, levels = c(non_std_names, std_names))
tax_table(ps_rename2) <- tax_table(as.matrix(new_tax_table))

# Visualize
ps_rel <- transform_sample_counts(ps_rename2, function(x) x/sum(x))
otu_table(ps_rel)[is.na(otu_table(ps_rel))] <- 0
ps_rel <- prune_taxa(taxa_sums(ps_rel) > 0, ps_rel)
ps_rel2 <- filter_taxa(ps_rel, function(x) mean(x) > 1e-02, TRUE)
ps_r1 <- speedyseq::psmelt(ps_rel)
ps_r2 <- stats::aggregate(ps_r1$Abundance, by=list(ps_r1$Sample, ps_r1$family), "sum") # Summed up to make phylum sum
ps_r3 <- stats::aggregate(ps_r1$Abundance, by=list(ps_r1$Sample, ps_r1$rep_tax), "sum") # Summed up to make phylum sum
colnames(ps_r2) <- c("sample", target_rank, "abundance")
colnames(ps_r3) <- c("sample", "rep_tax", "abundance")
## Adjust name orders
ps_r1$rep_tax <- factor(ps_r1$rep_tax, levels = c(non_std_names, std_names))
ps_r3$rep_tax <- factor(ps_r3$rep_tax, levels = c(non_std_names, std_names))

# Figures
r1 <- ggplot(ps_r1, aes(x = replicate, y = Abundance, fill = rep_tax)) +
  geom_bar(stat = "identity", colour = NA) +
  theme(axis.text.x = element_text(angle = -90, hjust = 1, vjust = 0.5, size = 6)) + 
  facet_wrap(~ site + index_method + enzyme, ncol = 3) +
  scale_fill_igv() +
  xlab(NULL) + ylab("Relative abundance") + panel_border()
r2 <- ps_r1 %>% filter(family == "STDseqs") %>%
  ggplot(aes(x = replicate, y = Abundance, fill = species)) +
  geom_bar(stat = "identity", colour = NA) +
  theme(axis.text.x = element_text(angle = -90, hjust = 1, vjust = 0.5, size = 6)) + 
  facet_wrap(~ site + index_method + enzyme, ncol = 3) +
  scale_fill_manual(values = get_palette(12)) +
  xlab(NULL) + ylab("Relative abundance") + panel_border()
r3 <- ps_r1 %>% filter(genus != "STDseqs") %>%
  ggplot(aes(x = replicate, y = Abundance, fill = family)) +
  geom_bar(stat = "identity", colour = NA) +
  theme(axis.text.x = element_text(angle = -90, hjust = 1, vjust = 0.5, size = 6)) + 
  facet_wrap(~ site + index_method + enzyme, ncol = 3) +
  scale_fill_manual(values = get_palette(length(unique(ps_r1$family)))) +
  xlab(NULL) + ylab("Sequence reads") + panel_border() + theme(legend.position = "none")


# ----------------------------------------------- #
#         Save data
# ----------------------------------------------- #
# Re-output data
write.csv(otu_table(ps_exp1_even), sprintf("%s/otu_table_exp1.csv", output_folder))
write.csv(sample_data(ps_exp1_even), sprintf("%s/sample_data_exp1.csv", output_folder))
write.csv(as.data.frame(tax_table(ps_exp1_even)), sprintf("%s/tax_table_exp1.csv", output_folder))
saveRDS(ps_exp1_even, sprintf("%s/ps_exp1_even.obj", output_folder))

# Save session info
writeLines(capture.output(sessionInfo()),
           paste0("../00_SessionInfo/", output_folder, "_", substr(Sys.time(), 1, 10), ".txt"))

# ----------------------------------------------- #
#   Save figures for publication
# ----------------------------------------------- #
# Save figure objects
#dir.create("../FigCode"); dir.create("../FigCode/00_RawFigs")
fig_dir <- "../FigCode/00_RawFigs/"
# Save figures for publication
saveRDS(list(f6, r1), paste0(fig_dir, "8_2_Fig_Exp1_CompositionDiveristy.obj"))

