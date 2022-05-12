####
#### No.11.2 Visualize index-dependence
#### 2022.05.12 revision for Environmental DNA
#### R 4.1.2
####

# Set working directory
if(basename(getwd()) != "11_SameIndexTest") setwd("11_SameIndexTest")

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
get_palette <- colorRampPalette(brewer.pal(8, "Paired"))

# Generate output folder
od <- basename(rstudioapi::getSourceEditorContext()$path)
(output_folder <- paste0(str_sub(od, end = -3), "Out")); rm(od)
dir.create(output_folder)


# ----------------------------------------------- #
#    Load phyloseq object
# ----------------------------------------------- #
ps_exp1_even <- readRDS("../08_Exp1_1st2nd/08_2_EvenDepthOut/ps_exp1_even.obj")
ps_exp3_even <- readRDS("../10_Exp3_repvol/10_2_EvenDepthOut/ps_exp3_even.obj")

# Edit I5_Index_ID2 and create Index_name
sample_data(ps_exp1_even)[,"Index_name"] <- 
  sample_data(ps_exp1_even)[,"I5_Index_ID2"] %>%
  unlist %>% str_extract(regex("IPg...."))
sample_data(ps_exp3_even)[,"Index_name"] <- 
  sample_data(ps_exp3_even)[,"I5_Index_ID2"] %>%
  unlist %>% str_extract(regex("IPg...."))

# ----------------------------------------------- #
#   Visualize pattern (Experiment 3)
# ----------------------------------------------- #
target_rank <- "family"
ps_exp3_even_stdM <- ps_exp3_even %>%
  subset_samples(site == "STD_Mix") %>%
  subset_taxa(habitat == "std") %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  transform_sample_counts(function(x) x/sum(x))
sample_sums(ps_exp3_even_stdM)

ps_rename2 <- taxa_name_summarize(ps_exp3_even_stdM, target_rank, top_taxa_n = 10)
ps_m2 <- speedyseq::psmelt(ps_rename2)
# Figures
f1 <- ggplot(ps_m2, aes(x = test_category, y = Abundance, group = rep_tax, fill = rep_tax)) +
    geom_bar(stat = "identity", colour = NA) +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5, size = 8)) + 
    xlab("Time") + ylab("Sequence reads") +
    facet_wrap(~ test_name + Index_name, ncol = 5) +
    scale_fill_manual(values = get_palette(11)) +
    ggtitle("Data from Experiment 3") +
    NULL

# Extract top taxa and visualize the relative abundance for each index
top_rank <- 1:3
top_taxa <- names(sort(taxa_sums(ps_exp3_even_stdM), decreasing = TRUE))[top_rank]
ps_exp3_top_melt <- ps_exp3_even_stdM %>% prune_taxa(top_taxa, .) %>%
  speedyseq::psmelt()
f2 <- ggplot(ps_exp3_top_melt, aes(x = Index_name, y = Abundance)) +
  geom_boxplot(outlier.shape = NA, outlier.size = 0) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.7) +
  facet_wrap(.~ OTU, nrow = 3) +
  panel_border() +
  ylim(0, 0.45) +
  xlab(NULL) + ylab("Relative abundance") +
  ggtitle("Experiment 3") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5, size = 10)) + 
  NULL

# GLM test # No effects of tag
ps_exp3_melt2 <- ps_exp3_even %>%
  subset_samples(site == "STD_Mix") %>%
  subset_taxa(habitat == "std") %>%
  prune_taxa(taxa_sums(.) > 0, .) %>% prune_taxa(top_taxa, .) %>%
  speedyseq::psmelt()
summary(aov(glm(cbind(Abundance, sample_sums(ps_exp3_even)[ps_exp3_melt2$Sample]) ~ OTU + Index_name,
                data = ps_exp3_melt2,
    family = binomial(link = "logit"))))


# ----------------------------------------------- #
#  Visualize pattern: Combine all results for each tag
# ----------------------------------------------- #
# Data from experiment 1
ps_exp1_stdM <- ps_exp1_even %>%
  subset_samples(site == "STD_Mix") %>%
  prune_taxa(taxa_sums(.) > 0, .)
sample_sums(ps_exp1_stdM)

ps_exp1_top_melt <- ps_exp1_stdM %>%
  subset_samples(enzyme == "Platinum") %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  prune_taxa(top_taxa, .) %>%
  speedyseq::psmelt()
ps_exp1_top_melt$Index_name[is.na(ps_exp1_top_melt$Index_name)] <-
  ps_exp1_top_melt$Index_name[is.na(ps_exp1_top_melt$Index_name)] <-
  "2nd PCR indexing"
f3 <- ggplot(ps_exp1_top_melt, aes(x = Index_name, y = Abundance)) +
  geom_boxplot(outlier.shape = NA, outlier.size = 0) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.7) +
  facet_wrap(.~ OTU, nrow = 3) +
  panel_border() +
  ylim(0, 0.45) +
  xlab(NULL) +
  ylab("Relative abundance") +
  ggtitle("Experiment 1") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5, size = 10)) + 
  NULL
mean_2nd_indexing <- f3$data %>% group_by(Index_name, OTU) %>%
  summarize(mean_prop = mean(Abundance)) %>% .[1:3,] %>%
  data.frame()

# GLM test # No effects of tag
ps_exp1_top_melt2 <- ps_exp1_stdM %>%
  subset_samples(enzyme == "Platinum") %>%
  prune_taxa(top_taxa, .) %>% 
  speedyseq::psmelt()
summary(aov(glm(cbind(Abundance,
                      sample_sums(ps_exp1_stdM)[ps_exp1_top_melt2$Sample]) ~ OTU + Index_name,
                data = ps_exp1_top_melt2,
                family = binomial(link = "logit"))))


# ----------------------------------------------- #
#         Save data
# ----------------------------------------------- #
# Save figures
prop_all <- plot_grid(f3 + geom_hline(yintercept = mean_2nd_indexing$mean_prop, linetype = 2, col = "red3"),
          f2 + geom_hline(yintercept = mean_2nd_indexing$mean_prop, linetype = 2, col = "red3"),
          align = "hv", axis = "lrtb", nrow = 1, rel_widths = c(1,1.5))
ggsave(file = sprintf("%s/TagDependence.pdf", output_folder),
       plot = prop_all, width = 14, height = 8)

# Save session info
writeLines(capture.output(sessionInfo()),
           paste0("../00_SessionInfo/", output_folder, "_", substr(Sys.time(), 1, 10), ".txt"))

