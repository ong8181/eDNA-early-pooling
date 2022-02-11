####
#### No.11.1 Visualize index-specific biases
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
ps_all <- readRDS("../07_CompilePhyloseqOut/ps_all.obj")
ps_exp1_even <- readRDS("../08_Exp1_1st2nd/08_2_EvenDepthOut/ps_exp1_even.obj")
ps_exp2 <- readRDS("../07_CompilePhyloseqOut/ps_exp2.obj")
ps_exp3_even <- readRDS("../10_Exp3_repvol/10_2_EvenDepthOut/ps_exp3_even.obj")

# Compile px_exp2
ps_exp2_even <- ps_exp2 %>% subset_samples(sample_nc == "sample") %>%
  subset_taxa(habitat != "freshwater" & habitat != "std") %>%
  prune_taxa(taxa_sums(.) > 0, .) %>% 
  rarefy_even_depth(rngseed = ran.seed, replace = FALSE, trimOTUs = FALSE)
sample_sums(ps_exp2_even)


# ----------------------------------------------- #
#   Visualize pattern (Experiment 2)
# ----------------------------------------------- #
# Visualize
target_rank <- "family"
ps_rename <- taxa_name_summarize(ps_exp2_even, target_rank, top_taxa_n = 10)
ps_m1 <- speedyseq::psmelt(ps_rename)
ps_m1$time <- factor(ps_m1$time)
# Figures
f1 <- ggplot(ps_m1, aes(x = Sample, y = Abundance, group = rep_tax, fill = rep_tax)) +
    geom_bar(stat = "identity", colour = NA) + theme(axis.text.x = element_text(angle = -90, hjust = 1, vjust = 0.5, size = 6)) + 
    xlab(NULL) + ylab("Sequence reads") +
    scale_fill_manual(values = get_palette(11)) +
    NULL
f2 <- ggplot(ps_m1, aes(x = time, y = Abundance, group = rep_tax, fill = rep_tax)) +
    geom_bar(stat = "identity", colour = NA) + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5, size = 8)) + 
    xlab("Time") + ylab("Sequence reads") +
    facet_wrap(. ~ purification_after_1st_pcr + temperature, ncol = 2) +
    scale_fill_manual(values = get_palette(11)) +
    NULL


# Visualize the relative abundance of top OTUs
top_taxa <- names(sort(taxa_sums(ps_exp2_even), decreasing = TRUE))[c(1,2,4)]
ps_exp2_top_melt <- ps_exp2_even %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  prune_taxa(top_taxa, .) %>%
  speedyseq::psmelt()
(f2_1 <- ggplot(ps_exp2_top_melt, aes(x = I5_Index_ID2, y = Abundance)) +
    geom_boxplot(outlier.shape = NA, outlier.size = 0) +
    geom_jitter(width = 0.2, height = 0, alpha = 0.7) +
    facet_wrap(.~ OTU, nrow = 3) +
    ggtitle("Experiment 2") +
    panel_border() +
    ylim(0,0.8) +
    xlab(NULL) + ylab("Relative abundance") +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5, size = 10)) + 
    NULL)


# ----------------------------------------------- #
#   Visualize pattern (Experiment 3)
# ----------------------------------------------- #
ps_exp3_even_seaN <- ps_exp3_even %>%
  subset_samples(site == "Sea_Nagahama") %>%
  subset_taxa(habitat != "std" & habitat != "freshwater") %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  transform_sample_counts(function(x) x/sum(x))
sample_sums(ps_exp3_even_seaN)

ps_rename2 <- taxa_name_summarize(ps_exp3_even_seaN, target_rank, top_taxa_n = 10)
ps_m2 <- speedyseq::psmelt(ps_rename2)
# Figures
f3 <- ggplot(ps_m2, aes(x = test_category, y = Abundance, group = rep_tax, fill = rep_tax)) +
    geom_bar(stat = "identity", colour = NA) +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5, size = 8)) + 
    xlab("Time") + ylab("Sequence reads") +
    facet_wrap(~ test_name + I5_Index_ID2, ncol = 5) +
    scale_fill_manual(values = get_palette(11)) +
    ggtitle("Data from Experiment 3") +
    NULL

 
# Extract top taxa and visualize the relative abundance for each index
top_rank <- 1:3
top_taxa <- names(sort(taxa_sums(ps_exp3_even_seaN), decreasing = TRUE))[top_rank]
ps_exp3_top_melt <- ps_exp3_even_seaN %>% prune_taxa(top_taxa, .) %>%
  speedyseq::psmelt()
f4 <- ggplot(ps_exp3_top_melt, aes(x = I5_Index_ID2, y = Abundance)) +
    geom_boxplot(outlier.shape = NA, outlier.size = 0) +
    geom_jitter(width = 0.2, height = 0, alpha = 0.7) +
    facet_wrap(.~ OTU, nrow = 3) +
    panel_border() +
    ylim(0,0.8) +
    xlab(NULL) + ylab("Relative abundance") +
    ggtitle("Experiment 3") +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5, size = 10)) + 
    NULL

# GLM test # No effects of tag
ps_exp3_melt2 <- ps_exp3_even %>%
  subset_samples(site == "Sea_Nagahama") %>%
  subset_taxa(habitat != "std" & habitat != "freshwater") %>%
  prune_taxa(taxa_sums(.) > 0, .) %>% prune_taxa(top_taxa, .) %>%
  speedyseq::psmelt()
summary(aov(glm(cbind(Abundance, sample_sums(ps_exp3_even)[1]) ~ OTU + I5_Index_ID2,
                data = ps_exp3_melt2,
    family = binomial(link = "logit"))))

summary(aov(glm(cbind(Abundance, sample_sums(ps_exp3_even)[1]) ~ I5_Index_ID2,
                data = ps_exp3_melt2 %>% filter(OTU == "OTU002"),
                family = binomial(link = "logit"))))
summary(aov(glm(cbind(Abundance, sample_sums(ps_exp3_even)[1]) ~ I5_Index_ID2,
                data = ps_exp3_melt2 %>% filter(OTU == "OTU007"),
                family = binomial(link = "logit"))))
summary(aov(glm(cbind(Abundance, sample_sums(ps_exp3_even)[1]) ~ I5_Index_ID2,
                data = ps_exp3_melt2 %>% filter(OTU == "OTU056"),
                family = binomial(link = "logit"))))


# ----------------------------------------------- #
#  Visualize pattern: Combine all results for each tag
# ----------------------------------------------- #
# Data from experiment 1
ps_exp1_seaN <- ps_exp1_even %>%
  subset_samples(site == "Sea_Nagahama") %>%
  prune_taxa(taxa_sums(.) > 0, .)
sample_sums(ps_exp1_seaN)

ps_exp1_top_melt <- ps_exp1_seaN %>%
  subset_samples(enzyme == "Platinum") %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  prune_taxa(top_taxa, .) %>%
  speedyseq::psmelt()
ps_exp1_top_melt$I5_Index_ID2[is.na(ps_exp1_top_melt$I5_Index_ID2)] <-
  ps_exp1_top_melt$I5_Index_ID2[is.na(ps_exp1_top_melt$I5_Index_ID2)] <-
  "2nd PCR indexing"
(f5 <- ggplot(ps_exp1_top_melt, aes(x = I5_Index_ID2, y = Abundance)) +
    geom_boxplot(outlier.shape = NA, outlier.size = 0) +
    geom_jitter(width = 0.2, height = 0, alpha = 0.7) +
    facet_wrap(.~ OTU, nrow = 3) +
    panel_border() +
    ylim(0,0.8) +
    xlab(NULL) +
    ylab("Relative abundance") +
    ggtitle("Experiment 1") +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5, size = 10)) + 
    NULL)
mean_2nd_indexing <- f5$data %>% group_by(I5_Index_ID2, OTU) %>%
  summarize(mean_prop = mean(Abundance)) %>% .[1:3,] %>%
  data.frame()

# GLM test # No effects of tag
rarefied_counts1 <- min(sample_sums(ps_exp1_seaN %>%
  subset_samples(enzyme == "Platinum") %>%
  prune_taxa(top_taxa, .)))
ps_exp1_top_melt2 <- ps_exp1_seaN %>%
  subset_samples(enzyme == "Platinum") %>%
  prune_taxa(top_taxa, .) %>% 
  speedyseq::psmelt()
summary(aov(glm(cbind(Abundance, rarefied_counts1) ~ OTU + I5_Index_ID2,
                data = ps_exp1_top_melt2,
                family = binomial(link = "logit"))))


# ----------------------------------------------- #
#         Save data
# ----------------------------------------------- #
# Save figures
f5_2 <- f5 + geom_hline(yintercept = mean_2nd_indexing$mean_prop, linetype = 2, col = "red3")
f2_2 <- f2_1 + geom_hline(yintercept = mean_2nd_indexing$mean_prop, linetype = 2, col = "red3")
f4_2 <- f4 + geom_hline(yintercept = mean_2nd_indexing$mean_prop, linetype = 2, col = "red3")
prop_all <- plot_grid(f5_2, f2_2, f4_2,
          align = "hv", axis = "lrtb", nrow = 1, rel_widths = c(3,1.3,5))
ggsave(file = sprintf("%s/Exp2_ComComposition.pdf", output_folder),
       plot = f1, width = 8, height = 8)
ggsave(file = sprintf("%s/Exp2_ComComposition2.pdf", output_folder),
       plot = f2, width = 8, height = 8)
ggsave(file = sprintf("%s/Exp3_ComComposition2.pdf", output_folder),
       plot = f3, width = 8, height = 8)
ggsave(file = sprintf("%s/IndexDependence.pdf", output_folder),
       plot = prop_all, width = 16, height = 8)

# Save figure objects
fig_dir <- "../FigCode/00_RawFigs/"
saveRDS(list(f1, f2, f3), paste0(fig_dir, "11_1_Fig_SameIndexComp_SeaN.obj"))
saveRDS(list(f5_2, f2_2, f4_2), paste0(fig_dir, "11_1_Fig_SameIndexBoxplot.obj"))

# Save session info
writeLines(capture.output(sessionInfo()),
           paste0("../00_SessionInfo/", output_folder, "_", substr(Sys.time(), 1, 10), ".txt"))
