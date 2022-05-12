####
#### No.12.1 Nagahama NMDS all
#### 2022.05.12 revision for Environmental DNA
#### R 4.1.2
####

# Set working directory
if(basename(getwd()) != "12_NagahamaNMDS") setwd("12_NagahamaNMDS")

# Set random seeds (for reproduction)
ran.seed <- 1234
set.seed(ran.seed)

# Load library and functions
library(tidyverse); packageVersion("tidyverse") # 1.3.1, 2021.10.16
library(phyloseq); packageVersion("phyloseq") # 1.38.0, 2021.11.18
library(cowplot); packageVersion("cowplot") # 1.1.1, 2021.6.13
library(ggsci); packageVersion("ggsci") # 2.9, 2021.6.13
library(RColorBrewer); packageVersion("RColorBrewer") # 1.1.3, 2022.5.5
library(ggrepel); packageVersion("ggrepel") # 0.9.1, 2021.11.29
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
ps_exp2 <- readRDS("../07_CompilePhyloseqOut/ps_exp2.obj")
ps_exp3_even <- readRDS("../10_Exp3_repvol/10_2_EvenDepthOut/ps_exp3_even.obj")

# Compile px_exp2
ps_exp2_even <- ps_exp2 %>% subset_samples(sample_nc == "sample") %>%
  subset_taxa(habitat != "freshwater" & habitat != "std") %>%
  prune_taxa(taxa_sums(.) > 0, .) %>% 
  rarefy_even_depth(rngseed = ran.seed, replace = FALSE, trimOTUs = FALSE)
sample_sums(ps_exp2_even)
## Remove OTUs < 0.01% in each sample
ps_exp2_even <- transform_sample_counts(ps_exp2_even, function(x){replace(x, x < ceiling(mean(sample_sums(ps_exp2_even))*0.0001), 0)}) 
mean(sample_sums(ps_exp2_even)) # 95321.33 reads for each sample


# ----------------------------------------------- #
#   Compile all Sea_Nagahama samples and rarefy
# ----------------------------------------------- #
ps_seaN_even <- merge_phyloseq(ps_exp1_even, ps_exp2_even, ps_exp3_even) %>%
  subset_samples((site == "River_Seta" | site == "Sea_Nagahama" | site == "Sea_Otomi") & sample_nc == "sample") %>%
  subset_taxa(habitat != "freshwater" & habitat != "std") %>%
  rarefy_even_depth(rngseed = ran.seed, replace = FALSE, trimOTUs = FALSE)
sample_sums(ps_seaN_even)


# ----------------------------------------------- #
#   Visualize pattern (Experiment 2)
# ----------------------------------------------- #
# Visualize
target_rank <- "family"
ps_seaN_new_df <- sample_data(ps_seaN_even) %>% data.frame %>%
  mutate(test_category2 = test_category)
ps_seaN_new_df$test_category2[ps_seaN_new_df$test_name == "ExoSAP_test" & ps_seaN_new_df$purification_after_1st_pcr == "ExoSAP"] <- "Purified_by_ExoSAP"
ps_seaN_new_df$test_category2[ps_seaN_new_df$test_name == "ExoSAP_test" & ps_seaN_new_df$purification_after_1st_pcr == "none"] <- "Not_purified"
sample_data(ps_seaN_even) <- sample_data(ps_seaN_new_df)
ps_rename <- taxa_name_summarize(ps_seaN_even %>% subset_samples(site == "Sea_Nagahama"),
                                 target_rank, top_taxa_n = 10)
ps_m1 <- speedyseq::psmelt(ps_rename)
# Figures
f1 <- ggplot(ps_m1, aes(x = Sample, y = Abundance, group = rep_tax, fill = rep_tax)) +
  geom_bar(stat = "identity", colour = NA) + theme(axis.text.x = element_text(angle = -90, hjust = 1, vjust = 0.5, size = 6)) + 
  xlab(NULL) + ylab("Sequence reads") +
  scale_fill_manual(name = "Family", values = get_palette(11)) +
  #facet_wrap(. ~ Experiment_ID, scale = "free_x") +
  facet_grid(. ~ Experiment_ID + test_category2, scale = "free", space = "free") +
  NULL


# ----------------------------------------------- #
#   NMDS
# ----------------------------------------------- #
# Nonmetric multidimensional scaling
set.seed(ran.seed); ps_seaN_bray <- ordinate(ps_seaN_even, "NMDS", "bray")
set.seed(ran.seed)
(n1 <- plot_ordination(ps_seaN_even, ps_seaN_bray,
                       shape = "purification_after_1st_pcr",
                       color = "test_category") +
    stat_ellipse(geom = "polygon", alpha = 0.3, aes(fill = site), color = "gray", type = "t") +
    geom_point(size = 4) +
    scale_color_igv(name = "Test category") +
    scale_shape_manual(name = "Post 1st PCR purification", values = c(16, 17)) +
    #scale_color_manual(values = get_palette(30)) +
    scale_fill_igv(name = "Site") +
    xlab("Axis 1") + ylab("Axis 2") +
    geom_text_repel(aes(label = test_category), max.overlaps = 30) +
    ggtitle(NULL) +
    NULL)

set.seed(ran.seed)
n2 <- plot_ordination(ps_seaN_even, ps_seaN_bray,
                      shape = "site",
                      color = "test_category2")
n2_df <- n2$data
n2_df$test_category2 <- paste0(n2_df$index_method, "_", n2_df$purification_after_1st_pcr)
n2$data <- n2_df
n2 <- n2 + stat_ellipse(geom = "polygon", alpha = 0.3, aes(fill = site), color = "gray", type = "t") +
  geom_point(size = 3, aes(color = test_category2)) +
  scale_color_igv(name = "Protocol") +
  scale_fill_igv(name = "Site") +
  scale_shape_manual(name = "Site", values = c(16, 15, 17)) +
  #geom_text_repel(aes(label = index_method), max.overlaps = 20) +
  xlab("Axis 1") + ylab("Axis 2") +
  ggtitle(NULL) +
  NULL


# ----------------------------------------------- #
#         Save data
# ----------------------------------------------- #
# Save figure objects
#dir.create("../FigCode"); dir.create("../FigCode/00_RawFigs")
fig_dir <- "../FigCode/00_RawFigs/"
saveRDS(list(f1, n1, n2), paste0(fig_dir, "12_1_Fig_NagahamaNMDS_all.obj"))

# Re-output data
write.csv(otu_table(ps_seaN_even), sprintf("%s/otu_table_seaN_even.csv", output_folder))
write.csv(sample_data(ps_seaN_even), sprintf("%s/sample_data_seaN_even.csv", output_folder))
write.csv(as.data.frame(tax_table(ps_seaN_even)), sprintf("%s/tax_table_seaN_even.csv", output_folder))

# Save session info
writeLines(capture.output(sessionInfo()),
           paste0("../00_SessionInfo/", output_folder, "_", substr(Sys.time(), 1, 10), ".txt"))

