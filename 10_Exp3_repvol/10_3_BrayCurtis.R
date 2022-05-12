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
library(RColorBrewer); packageVersion("RColorBrewer") # 1.1.3, 2021.6.13
theme_set(theme_cowplot())
get_palette <- colorRampPalette(brewer.pal(8, "Paired"))

# Generate output folder
od <- basename(rstudioapi::getSourceEditorContext()$path)
(output_folder <- paste0(str_sub(od, end = -3), "Out")); rm(od)
dir.create(output_folder)


# ----------------------------------------------- #
#    Load phyloseq object
# ----------------------------------------------- #
ps_exp3_even <- readRDS("10_2_EvenDepthOut/ps_exp3_even.obj")

## Divide the phyloseq object into each site x experiment
ps_exp3_seaN_rep <- ps_exp3_even %>%
  subset_samples(site == "Sea_Nagahama" & test_name == "Replicate_test") %>%
  prune_taxa(taxa_sums(.) > 0, .)
ps_exp3_seaN_vol <- ps_exp3_even %>%
  subset_samples(site == "Sea_Nagahama" & test_name == "Volume_test") %>%
  prune_taxa(taxa_sums(.) > 0, .)
ps_exp3_stdM_rep <- ps_exp3_even %>%
  subset_samples(site == "STD_Mix" & test_name == "Replicate_test") %>%
  prune_taxa(taxa_sums(.) > 0, .)
ps_exp3_stdM_vol <- ps_exp3_even %>%
  subset_samples(site == "STD_Mix" & test_name == "Volume_test") %>%
  prune_taxa(taxa_sums(.) > 0, .)


# ----------------------------------------------- #
#   Calculate Bray-Curtis dissimilarity
# ----------------------------------------------- #
# Bary-curtis dissimilarity
bray1 <- vegan::vegdist(as.matrix(otu_table(ps_exp3_seaN_rep)), method = "bray") %>% as.matrix()
bray2 <- vegan::vegdist(as.matrix(otu_table(ps_exp3_seaN_vol)), method = "bray") %>% as.matrix()
bray3 <- vegan::vegdist(as.matrix(otu_table(ps_exp3_stdM_rep)), method = "bray") %>% as.matrix()
bray4 <- vegan::vegdist(as.matrix(otu_table(ps_exp3_stdM_vol)), method = "bray") %>% as.matrix()

# Seq_Nagahama, Replication test
seaN_rep_scale1 <- bray1[  1:4,  1:4][lower.tri(bray1[  1:4,  1:4])] %>% c(., rep(NA, 4))
seaN_rep_scale2 <- bray1[  5:9,  5:9][lower.tri(bray1[  5:9,  5:9])]
seaN_rep_scale4 <- bray1[10:14,10:14][lower.tri(bray1[10:14,10:14])]
seaN_rep_scale8 <- bray1[15:19,15:19][lower.tri(bray1[15:19,15:19])]

# Seq_Nagahama, Volume test
seaN_vol_scale1 <- bray2[  1:5,  1:5][lower.tri(bray2[  1:5,  1:5])]
seaN_vol_scale2 <- bray2[ 6:10, 6:10][lower.tri(bray2[ 6:10, 6:10])]
seaN_vol_scale4 <- bray2[11:15,11:15][lower.tri(bray2[11:15,11:15])]
seaN_vol_scale8 <- bray2[16:20,16:20][lower.tri(bray2[16:20,16:20])]

# Seq_Nagahama, Replication test
stdM_rep_scale1 <- bray3[  1:5,  1:5][lower.tri(bray3[  1:5,  1:5])]
stdM_rep_scale2 <- bray3[ 6:10, 6:10][lower.tri(bray3[ 6:10, 6:10])]
stdM_rep_scale4 <- bray3[11:15,11:15][lower.tri(bray3[11:15,11:15])]
stdM_rep_scale8 <- bray3[16:20,16:20][lower.tri(bray3[16:20,16:20])]

# Seq_Nagahama, Replication test
stdM_vol_scale1 <- bray4[  1:5,  1:5][lower.tri(bray4[  1:5,  1:5])]
stdM_vol_scale2 <- bray4[ 6:10, 6:10][lower.tri(bray4[ 6:10, 6:10])]
stdM_vol_scale4 <- bray4[11:15,11:15][lower.tri(bray4[11:15,11:15])]
stdM_vol_scale8 <- bray4[16:20,16:20][lower.tri(bray4[16:20,16:20])]


# ----------------------------------------------- #
#         Generate a new data.frame
# ----------------------------------------------- #
# Within method variations
bray_curtis_values <- c(seaN_rep_scale1, seaN_rep_scale2, seaN_rep_scale4, seaN_rep_scale8,
                        seaN_vol_scale1, seaN_vol_scale2, seaN_vol_scale4, seaN_vol_scale8,
                        stdM_rep_scale1, stdM_rep_scale2, stdM_rep_scale4, stdM_rep_scale8,
                        stdM_vol_scale1, stdM_vol_scale2, stdM_vol_scale4, stdM_vol_scale8)
bray_df <- data.frame(data_id = sprintf("ID%03d", 1:(40*4)),
                      rep_id = rep(sprintf("Rep%02d", 1:10), 16),
                      site = c(rep("Sea_Nagahama", 80), rep("STD_Mix", 80)),
                      test_name = rep(c(rep("Replicate_test",40), rep("Volume_test",40)), 2),
                      reaction_scale = rep(c(rep("1_rep",10), rep("2_rep",10), rep("4_rep",10), rep("8_rep",10),
                                             rep("1_vol",10), rep("2_vol",10), rep("4_vol",10), rep("8_vol",10)), 2)
                      ) %>%
  mutate(treatment = sprintf("%s-%s-%s", site, test_name, reaction_scale)) %>%
  mutate(bray_curtis = bray_curtis_values) %>% na.omit()


# ----------------------------------------------- #
#         Visualize pattern
# ----------------------------------------------- #
(g1 <- bray_df %>% ggplot(aes(x = reaction_scale, y = bray_curtis)) +
    geom_boxplot(outlier.color = NA, outlier.size = 0) +
    geom_jitter(height = 0, width = 0.2) +
    facet_wrap(~ site + test_name, ncol = 2, scales = "free_x") +
    panel_border() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    xlab(NULL) + ylab("Bray-Curtis dissimilarity") +
    ylim(0, 1) +
    NULL)


# ----------------------------------------------- #
#   Custom bootstrap test
# ----------------------------------------------- #
bray_boot <- function(bray_df_object, n_iter = 1000, ran.seed = 1234) {
  # Convert to matrix
  bray_df_object <- as.matrix(bray_df_object)
  n_rows <- dim(bray_df_object)[1]
  # Set ID values
  if (n_rows < 20) {
    id_list <- list(df1 = 1:4, df2 = 5:9, df3 = 10:14, df4 = 15:19)
    id_method <- c(6,10,10,10)
  } else {
    id_list <- list(df1 = 1:5, df2 = 6:10, df3 = 11:15, df4 = 16:20)
    id_method <- c(10,10,10,10)
  }
  # Prepare output object
  explained_boot_all <- c(NULL)
  
  # Bootstrap loop
  for(i in 1:n_iter) {
    # Generate random index
    # Labels "1st PCR indexing" or "2nd PCR indexing" were randomly shuffled
    # And the resultant difference between bray-curtis values was compared with the original one
    set.seed(ran.seed+i); (index_id_boot <- sample(n_rows, n_rows, replace = FALSE))
    # Shuffle the distance matrix
    bray_df_boot <- bray_df_object[index_id_boot, index_id_boot]
    # Extract (pseudo) treatment
    bray_df1 <- bray_df_boot[id_list$df1,id_list$df1]
    bray_df2 <- bray_df_boot[id_list$df2,id_list$df2]
    bray_df3 <- bray_df_boot[id_list$df3,id_list$df3]
    bray_df4 <- bray_df_boot[id_list$df4,id_list$df4]
    bray_vals <- c(bray_df1[lower.tri(bray_df1)], bray_df2[lower.tri(bray_df2)],
                   bray_df3[lower.tri(bray_df3)], bray_df4[lower.tri(bray_df4)])
    bray_df_boot2 <- data.frame(method = c(rep("method1", id_method[1]), rep("method2", id_method[2]),
                                           rep("method3", id_method[3]), rep("method4", id_method[4])),
                                values = bray_vals) %>% tibble()
    # Calculate mean residuals
    bray_boot_summary <- summary(aov(bray_df_boot2$values ~ bray_df_boot2$method))
    explained_boot <- bray_boot_summary[[1]]$`Sum Sq`[1]
    explained_boot_all <- c(explained_boot_all, explained_boot)
  }
  
  # Get original values
  bray_o_df1 <- bray_df_object[id_list$df1,id_list$df1]
  bray_o_df2 <- bray_df_object[id_list$df2,id_list$df2]
  bray_o_df3 <- bray_df_object[id_list$df3,id_list$df3]
  bray_o_df4 <- bray_df_object[id_list$df4,id_list$df4]
  bray_o_vals <- c(bray_o_df1[lower.tri(bray_o_df1)], bray_o_df2[lower.tri(bray_o_df2)],
                   bray_o_df3[lower.tri(bray_o_df3)], bray_o_df4[lower.tri(bray_o_df4)])
  bray_o_df <- data.frame(method = c(rep("method1", id_method[1]), rep("method2", id_method[2]),
                                     rep("method3", id_method[3]), rep("method4", id_method[4])),
                          values = bray_o_vals) %>% tibble()
  bray_o_summary <- summary(aov(bray_o_df$values ~ bray_o_df$method))
  explained_sq <- bray_o_summary[[1]]$`Sum Sq`[1]
  # Calculate mean difference
  quant_explained <- quantile(explained_boot_all, probs = c(0.9, 0.95, 0.975, 0.99, 0.995, 0.999))
  boot_p_val <- sum(explained_boot_all - explained_sq > 0)/n_iter
  res_all <- c(quant_explained, explained_sq, boot_p_val)
  names(res_all)[(length(res_all)-1):length(res_all)] <- c("observed", "p_val")
  return(res_all)
}

# Test Bary-Curtis dissimilarity values
## Reaction_scale influences Bray-Curtis dissimilarities?
(seaN_rep <- bray_boot(as.dist(bray1))) # P = 0.039
(seaN_vol <- bray_boot(as.dist(bray2))) # P = 0.043
(stdM_rep <- bray_boot(as.dist(bray3))) # P = 0.003
(stdM_vol <- bray_boot(as.dist(bray4))) # P = 0.006


# ----------------------------------------------- #
#         Save figures
# ----------------------------------------------- #
ggsave(file = sprintf("%s/BrayCurtis_repvol.pdf", output_folder), plot = g1, width = 6, height = 6)


# ----------------------------------------------- #
#         Save data
# ----------------------------------------------- #
# Save figure objects
#dir.create("../FigCode"); dir.create("../FigCode/00_RawFigs")
fig_dir <- "../FigCode/00_RawFigs/"
saveRDS(g1, paste0(fig_dir, "10_3_Fig_Exp3_RepVolBrayCurtis.obj"))

# Save session info
writeLines(capture.output(sessionInfo()),
           paste0("../00_SessionInfo/", output_folder, "_", substr(Sys.time(), 1, 10), ".txt"))

