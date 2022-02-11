####
#### No.13.1 Summarizing rare OTUs
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
get_palette <- colorRampPalette(brewer.pal(8, "Paired"))

# Generate output folder
od <- basename(rstudioapi::getSourceEditorContext()$path)
(output_folder <- paste0(str_sub(od, end = -3), "Out")); rm(od)
dir.create(output_folder)


# ----------------------------------------------- #
#    Load phyloseq object
# ----------------------------------------------- #
ps_exp3_even <- readRDS("../10_Exp3_repvol/10_2_EvenDepthOut/ps_exp3_even.obj")
sample_sums(ps_exp3_even)

# ----------------------------------------------- #
#    Extract treatment-specific taxa
# ----------------------------------------------- #
# Replication test
ps_exp3_1rep <- ps_exp3_even %>% subset_samples(test_category == "1rep") %>%
  subset_samples(site == "Sea_Nagahama") %>%
  prune_taxa(taxa_sums(.) > 0, .)
ps_exp3_2rep <- ps_exp3_even %>% subset_samples(test_category == "2rep") %>%
  subset_samples(site == "Sea_Nagahama") %>%
  prune_taxa(taxa_sums(.) > 0, .)
ps_exp3_4rep <- ps_exp3_even %>% subset_samples(test_category == "4rep") %>%
  subset_samples(site == "Sea_Nagahama") %>%
  prune_taxa(taxa_sums(.) > 0, .)
ps_exp3_8rep <- ps_exp3_even %>% subset_samples(test_category == "8rep") %>%
  subset_samples(site == "Sea_Nagahama") %>%
  prune_taxa(taxa_sums(.) > 0, .)
# Volume test
ps_exp3_1vol <- ps_exp3_even %>% subset_samples(test_category == "1ul") %>%
  subset_samples(site == "Sea_Nagahama") %>%
  prune_taxa(taxa_sums(.) > 0, .)
ps_exp3_2vol <- ps_exp3_even %>% subset_samples(test_category == "2ul") %>%
  subset_samples(site == "Sea_Nagahama") %>%
  prune_taxa(taxa_sums(.) > 0, .)
ps_exp3_4vol <- ps_exp3_even %>% subset_samples(test_category == "4ul") %>%
  subset_samples(site == "Sea_Nagahama") %>%
  prune_taxa(taxa_sums(.) > 0, .)
ps_exp3_8vol <- ps_exp3_even %>% subset_samples(test_category == "8ul") %>%
  subset_samples(site == "Sea_Nagahama") %>%
  prune_taxa(taxa_sums(.) > 0, .)

# Identify treatment-specific taxa
## Replication test
taxa_newin_1rep <- taxa_names(ps_exp3_1rep)
taxa_newin_2rep <- taxa_names(ps_exp3_2rep)[!taxa_names(ps_exp3_2rep) %in% taxa_names(ps_exp3_1rep)]
taxa_newin_4rep <- taxa_names(ps_exp3_4rep)[!taxa_names(ps_exp3_4rep) %in% unique(c(taxa_names(ps_exp3_2rep), taxa_names(ps_exp3_1rep)))]
taxa_newin_8rep <- taxa_names(ps_exp3_8rep)[!taxa_names(ps_exp3_8rep) %in% unique(c(taxa_names(ps_exp3_4rep), taxa_names(ps_exp3_2rep), taxa_names(ps_exp3_1rep)))]
length(unique(c(taxa_names(ps_exp3_1rep), taxa_names(ps_exp3_2rep), taxa_names(ps_exp3_4rep), taxa_names(ps_exp3_8rep))))
length(unique(c(taxa_newin_1rep, taxa_newin_2rep, taxa_newin_4rep, taxa_newin_8rep)))
sum(taxa_names(ps_exp3_1rep) %in% taxa_names(ps_exp3_8rep))
taxa_newin_2rep %in% taxa_names(ps_exp3_8rep)
taxa_newin_1rep %in% taxa_names(ps_exp3_8rep)
## Taxa detected in all replication treatments
taxa_in_allrep <- taxa_names(ps_exp3_1rep) %>%
  intersect(taxa_names(ps_exp3_2rep)) %>%
  intersect(taxa_names(ps_exp3_4rep)) %>%
  intersect(taxa_names(ps_exp3_8rep)) %>% sort()

## Volume test
taxa_newin_1vol <- taxa_names(ps_exp3_1vol)
taxa_newin_2vol <- taxa_names(ps_exp3_2vol)[!taxa_names(ps_exp3_2vol) %in% taxa_names(ps_exp3_1vol)]
taxa_newin_4vol <- taxa_names(ps_exp3_4vol)[!taxa_names(ps_exp3_4vol) %in% unique(c(taxa_names(ps_exp3_2vol), taxa_names(ps_exp3_1vol)))]
taxa_newin_8vol <- taxa_names(ps_exp3_8vol)[!taxa_names(ps_exp3_8vol) %in% unique(c(taxa_names(ps_exp3_4vol), taxa_names(ps_exp3_2vol), taxa_names(ps_exp3_1vol)))]
length(unique(c(taxa_names(ps_exp3_1vol), taxa_names(ps_exp3_2vol), taxa_names(ps_exp3_4vol), taxa_names(ps_exp3_8vol))))
length(unique(c(taxa_newin_1vol, taxa_newin_2vol, taxa_newin_4vol, taxa_newin_8vol)))
sum(taxa_names(ps_exp3_1vol) %in% taxa_names(ps_exp3_8vol))
taxa_newin_2vol %in% taxa_names(ps_exp3_8vol)
taxa_newin_1vol %in% taxa_names(ps_exp3_8vol)
## Taxa detected in all replication treatments
taxa_in_allvol <- taxa_names(ps_exp3_1vol) %>%
  intersect(taxa_names(ps_exp3_2vol)) %>%
  intersect(taxa_names(ps_exp3_4vol)) %>%
  intersect(taxa_names(ps_exp3_8vol)) %>% sort()


# ----------------------------------------------- #
# Assgin detected treatment categories
# ----------------------------------------------- #
# Replication test
tax_table_rep <- merge_phyloseq(ps_exp3_1rep, ps_exp3_2rep, ps_exp3_4rep, ps_exp3_8rep) %>%
  tax_table %>% data.frame %>% mutate(detected_tr = NaN) %>% mutate(detected_tr2 = NaN)
## Detected treatment 1
tax_table_rep$detected_tr[rownames(tax_table_rep) %in% taxa_newin_1rep] <- "detected_in_1_rep"
tax_table_rep$detected_tr[rownames(tax_table_rep) %in% taxa_newin_2rep] <- "detected_in_2_rep"
tax_table_rep$detected_tr[rownames(tax_table_rep) %in% taxa_newin_4rep] <- "detected_in_4_rep"
tax_table_rep$detected_tr[rownames(tax_table_rep) %in% taxa_newin_8rep] <- "detected_in_8_rep"
## Detected treatment 1
tax_table_rep$detected_tr2[rownames(tax_table_rep) %in% taxa_in_allrep] <- "detected_in_all_rep"
tax_table_rep$detected_tr2[!rownames(tax_table_rep) %in% taxa_in_allrep] <- "detected_in_some_rep"
## Re-input tax_table
ps_exp3_rep_all <- merge_phyloseq(ps_exp3_1rep, ps_exp3_2rep, ps_exp3_4rep, ps_exp3_8rep) %>%
  transform_sample_counts(function(x) x/sum(x))
tax_table(ps_exp3_rep_all) <- tax_table(as.matrix(tax_table_rep))

# Volume test
tax_table_vol <- merge_phyloseq(ps_exp3_1vol, ps_exp3_2vol, ps_exp3_4vol, ps_exp3_8vol) %>%
  tax_table %>% data.frame %>% mutate(detected_tr = NaN) %>% mutate(detected_tr2 = NaN)
## Detected treatment 1
tax_table_vol$detected_tr[rownames(tax_table_vol) %in% taxa_newin_1vol] <- "detected_in_1ul"
tax_table_vol$detected_tr[rownames(tax_table_vol) %in% taxa_newin_2vol] <- "detected_in_2ul"
tax_table_vol$detected_tr[rownames(tax_table_vol) %in% taxa_newin_4vol] <- "detected_in_4ul"
tax_table_vol$detected_tr[rownames(tax_table_vol) %in% taxa_newin_8vol] <- "detected_in_8ul"
## Detected treatment 1
tax_table_vol$detected_tr2[rownames(tax_table_vol) %in% taxa_in_allvol] <- "detected_in_all_vol"
tax_table_vol$detected_tr2[!rownames(tax_table_vol) %in% taxa_in_allvol] <- "detected_in_some_vol"
## Re-input tax_table
ps_exp3_vol_all <- merge_phyloseq(ps_exp3_1vol, ps_exp3_2vol, ps_exp3_4vol, ps_exp3_8vol) %>%
  transform_sample_counts(function(x) x/sum(x))
tax_table(ps_exp3_vol_all) <- tax_table(as.matrix(tax_table_vol))


# ----------------------------------------------- #
# List up rare OTU names
# ----------------------------------------------- #
# Calculate shared OTU ratio
shared_otu_ratio <- function(ps_object1, ps_object2) {
  n_shared_otu <- length(intersect(taxa_names(ps_object1), taxa_names(ps_object2)))
  mean_otu <- mean(length(taxa_names(ps_object1)), length(taxa_names(ps_object2)))
  res_vec <- c(n_shared_otu, mean_otu, round(n_shared_otu/mean_otu, 5))
  names(res_vec) <- c("n_shared_otu", "mean_n_otu", "shared_otu_ratio")
  return(res_vec)
}
shared_otu_ratio(ps_exp3_rep_all, ps_exp3_vol_all)
shared_otu_ratio(ps_exp3_1rep, ps_exp3_1vol)
shared_otu_ratio(ps_exp3_2rep, ps_exp3_2vol)
shared_otu_ratio(ps_exp3_4rep, ps_exp3_4vol)
shared_otu_ratio(ps_exp3_8rep, ps_exp3_8vol)

## Check unique OTUs of "rep" against "vol" treatment
length(setdiff(taxa_names(ps_exp3_1rep), taxa_names(ps_exp3_8vol)))/length(taxa_names(ps_exp3_1rep)) # 90% detected in 8vol
length(setdiff(taxa_names(ps_exp3_2rep), taxa_names(ps_exp3_8vol)))/length(taxa_names(ps_exp3_2rep)) # 83% detected in 8vol
length(setdiff(taxa_names(ps_exp3_4rep), taxa_names(ps_exp3_8vol)))/length(taxa_names(ps_exp3_4rep)) # 72% detected in 8vol
length(setdiff(taxa_names(ps_exp3_8rep), taxa_names(ps_exp3_8vol)))/length(taxa_names(ps_exp3_8rep)) # 71% detected in 8vol
## Check unique OTUs of "vol" against "rep" treatment
length(setdiff(taxa_names(ps_exp3_1vol), taxa_names(ps_exp3_8rep)))/length(taxa_names(ps_exp3_1vol)) # 96% detected in 8rep
length(setdiff(taxa_names(ps_exp3_2vol), taxa_names(ps_exp3_8rep)))/length(taxa_names(ps_exp3_2vol)) # 85% detected in 8rep
length(setdiff(taxa_names(ps_exp3_4vol), taxa_names(ps_exp3_8rep)))/length(taxa_names(ps_exp3_4vol)) # 87% detected in 8rep
length(setdiff(taxa_names(ps_exp3_8vol), taxa_names(ps_exp3_8rep)))/length(taxa_names(ps_exp3_8vol)) # 76% detected in 8rep

# Check OTU characteristics
taxa_rep_vol <- intersect(taxa_names(ps_exp3_8rep), taxa_names(ps_exp3_8vol))
taxa_only_rep <- setdiff(taxa_names(ps_exp3_8rep), taxa_names(ps_exp3_8vol))
taxa_only_vol <- setdiff(taxa_names(ps_exp3_8vol), taxa_names(ps_exp3_8rep))
ps_exp3_tax_all <- data.frame(tax_table(ps_exp3_even))
ps_exp3_tax_all$total_reads <- taxa_sums(ps_exp3_even)
ps_exp3_tax_all[taxa_rep_vol,c("family","genus","species","total_reads","jp_name_memo")]
ps_exp3_tax_all[taxa_only_rep,c("family","genus","species","total_reads","jp_name_memo")]
ps_exp3_tax_all[taxa_only_vol,c("family","genus","species","total_reads","jp_name_memo")]

# Venn diagram
venn_rep_data <- list(rep1 = taxa_names(ps_exp3_1rep),
                      rep2 = taxa_names(ps_exp3_2rep),
                      rep4 = taxa_names(ps_exp3_4rep),
                      rep8 = taxa_names(ps_exp3_8rep))
venn_vol_data <- list(vol1 = taxa_names(ps_exp3_1vol),
                      vol2 = taxa_names(ps_exp3_2vol),
                      vol4 = taxa_names(ps_exp3_4vol),
                      vol8 = taxa_names(ps_exp3_8vol))
venn_vol_rep <- list(rep8 = taxa_names(ps_exp3_8rep),
                     rep4 = taxa_names(ps_exp3_4rep),
                     vol4 = taxa_names(ps_exp3_4vol),
                     vol8 = taxa_names(ps_exp3_8vol))
v1 <- ggVennDiagram::ggVennDiagram(venn_rep_data)
v2 <- ggVennDiagram::ggVennDiagram(venn_vol_data)
v3 <- ggVennDiagram::ggVennDiagram(venn_vol_rep)


# ----------------------------------------------- #
# Visualize pattern: Which OTUs can be detected in each treatment?
# ----------------------------------------------- #
# Replication test
otu_table(ps_exp3_rep_all); sample_sums(ps_exp3_rep_all)
ps_rep_m1 <- speedyseq::psmelt(ps_exp3_rep_all)
g1 <- ps_rep_m1 %>%
  filter(Abundance > 0) %>%
  mutate(Abundance = replace(Abundance, Abundance < 0.005, 0.005)) %>%
  ggplot(aes(x = detected_tr, y = Abundance*100, color = test_category)) + 
  geom_jitter(width = 0.2, height = 0) +
  scale_y_log10(limits = c(0.5, 100)) +
  scale_color_startrek() + theme(axis.text.x = element_text(angle = 45, vjust = 0.8, hjust = 0.8)) +
  ggtitle("Replication test") +
  xlab(NULL) + ylab("Relative abundance (%)") +
  NULL
g2 <- ps_rep_m1 %>%
  filter(Abundance > 0) %>%
  mutate(Abundance = replace(Abundance, Abundance < 0.005, 0.005)) %>%
  ggplot(aes(x = detected_tr2, y = Abundance*100, color = test_category)) + 
  geom_jitter(width = 0.2, height = 0) +
  scale_y_log10(limits = c(0.5, 100)) +
  scale_color_startrek() + theme(axis.text.x = element_text(angle = 45, vjust = 0.8, hjust = 0.8)) +
  ggtitle("Replication test") +
  xlab(NULL) + ylab("Relative abundance (%)") +
  NULL

# Volume test
otu_table(ps_exp3_vol_all); sample_sums(ps_exp3_vol_all)
ps_vol_m1 <- speedyseq::psmelt(ps_exp3_vol_all)
g3 <- ps_vol_m1 %>%
  filter(Abundance > 0) %>%
  mutate(Abundance = replace(Abundance, Abundance < 0.005, 0.005)) %>%
  ggplot(aes(x = detected_tr, y = Abundance * 100, color = test_category)) + 
  geom_jitter(width = 0.2, height = 0) +
  scale_y_log10(limits = c(0.5, 100)) +
  scale_color_startrek() + theme(axis.text.x = element_text(angle = 45, vjust = 0.8, hjust = 0.8)) +
  ggtitle("Volume test") +
  xlab(NULL) + ylab("Relative abundance (%)") +
  NULL
g4 <- ps_vol_m1 %>%
  filter(Abundance > 0) %>%
  mutate(Abundance = replace(Abundance, Abundance < 0.005, 0.005)) %>%
  ggplot(aes(x = detected_tr2, y = Abundance * 100, color = test_category)) + 
  geom_jitter(width = 0.2, height = 0) +
  scale_y_log10(limits = c(0.5, 100)) +
  scale_color_startrek() + theme(axis.text.x = element_text(angle = 45, vjust = 0.8, hjust = 0.8)) +
  ggtitle("Volume test") +
  xlab(NULL) + ylab("Relative abundance (%)") +
  NULL


# ----------------------------------------------- #
#         Save data
# ----------------------------------------------- #
# Save figures
ggsave(file = sprintf("%s/RareOTUs_SeaNagahama_Rep.pdf", output_folder),
       plot = plot_grid(g1, g2, rel_widths = c(1,0.7), align = "hv", ncol = 2), width = 12, height = 6)
ggsave(file = sprintf("%s/RareOTUs_SeaNagahama_Vol.pdf", output_folder),
       plot = plot_grid(g3, g4, rel_widths = c(1,0.7), align = "hv", ncol = 2), width = 12, height = 6)

# Save figure objects
fig_dir <- "../FigCode/00_RawFigs/"
saveRDS(list(g1, g2, g3, g4), paste0(fig_dir, "13_1_Fig_RareOTUs.obj"))

# Save session info
writeLines(capture.output(sessionInfo()),
           paste0("../00_SessionInfo/", output_folder, "_", substr(Sys.time(), 1, 10), ".txt"))

