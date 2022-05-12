####
#### Community difference
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
library(RColorBrewer); packageVersion("RColorBrewer") # 1.1.3, 2022.5.5
theme_set(theme_cowplot())

# Generate output folder
od <- basename(rstudioapi::getSourceEditorContext()$path)
(output_folder <- paste0(str_sub(od, end = -3), "Out")); rm(od)
dir.create(output_folder)

# ----------------------------------------------- #
#    Load phyloseq object
# ----------------------------------------------- #
ps_all <- readRDS("../07_CompilePhyloseqOut/ps_all.obj")
ps_exp1_even <- readRDS("08_2_EvenDepthOut/ps_exp1_even.obj")


# ----------------------------------------------- #
#         Separate by site
# ----------------------------------------------- #
## Sea_Nagahama
ps_fish1 <- subset_samples(ps_exp1_even, site == "Sea_Nagahama") %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  transform_sample_counts(function(x) x/sum(x))
sample_data(ps_fish1) <- sample_data(ps_fish1) %>% data.frame %>%
  mutate(libprep = paste0(.$index_method, "_", .$enzyme)) %>% sample_data()
## Sea_Otomi
ps_fish2 <- subset_samples(ps_exp1_even, site == "Sea_Otomi") %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  transform_sample_counts(function(x) x/sum(x))
sample_data(ps_fish2) <- sample_data(ps_fish2) %>% data.frame %>%
  mutate(libprep = paste0(.$index_method, "_", .$enzyme)) %>% sample_data()
## River_Seta
ps_fish3 <- subset_samples(ps_exp1_even, site == "River_Seta") %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  transform_sample_counts(function(x) x/sum(x))
sample_data(ps_fish3) <- sample_data(ps_fish3) %>% data.frame %>%
  mutate(libprep = paste0(.$index_method, "_", .$enzyme)) %>% sample_data()
## STD_Mix
ps_fish4 <- subset_samples(ps_exp1_even, site == "STD_Mix") %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  transform_sample_counts(function(x) x/sum(x))
sample_data(ps_fish4) <- sample_data(ps_fish4) %>% data.frame %>%
  mutate(libprep = paste0(.$index_method, "_", .$enzyme)) %>% sample_data()
sample_sums(ps_fish1); sample_sums(ps_fish2)
sample_sums(ps_fish3); sample_sums(ps_fish4)


# ----------------------------------------------- #
#         Calculating community difference
# ----------------------------------------------- #
# Bary-curtis dissimilarity
bray1 <- vegan::vegdist(as.matrix(otu_table(ps_fish1)), method = "bray") %>% as.matrix()
bray2 <- vegan::vegdist(as.matrix(otu_table(ps_fish2)), method = "bray") %>% as.matrix()
bray3 <- vegan::vegdist(as.matrix(otu_table(ps_fish3)), method = "bray") %>% as.matrix()
bray4 <- vegan::vegdist(as.matrix(otu_table(ps_fish4)), method = "bray") %>% as.matrix()

# Extract ID for stats
plat_id <- c(1:10)
index2nd_id <- c(6:15)
newold_id <- c(1:5,11:15)

# Define grid - category
seaN_index1st_plat <- bray1[  1:5,  1:5][lower.tri(bray1[  1:5,  1:5])]
seaN_index2nd_plat <- bray1[ 6:10, 6:10][lower.tri(bray1[ 6:10, 6:10])]
seaN_index2nd_kapa <- bray1[11:15,11:15][lower.tri(bray1[11:15,11:15])]
plat_only1 <- as.dist(bray1[plat_id, plat_id])
id2nd_only1 <- as.dist(bray1[index2nd_id, index2nd_id])
new_old1 <- as.dist(bray1[newold_id,newold_id])

seaO_index1st_plat <- bray2[  1:5,  1:5][lower.tri(bray2[  1:5,  1:5])]
seaO_index2nd_plat <- bray2[ 6:10, 6:10][lower.tri(bray2[ 6:10, 6:10])]
seaO_index2nd_kapa <- bray2[11:15,11:15][lower.tri(bray2[11:15,11:15])]
plat_only2 <- as.dist(bray2[plat_id, plat_id])
id2nd_only2 <- as.dist(bray2[index2nd_id, index2nd_id])
new_old2 <- as.dist(bray2[newold_id,newold_id])

rivS_index1st_plat <- bray3[  1:5,  1:5][lower.tri(bray3[  1:5,  1:5])]
rivS_index2nd_plat <- bray3[ 6:10, 6:10][lower.tri(bray3[ 6:10, 6:10])]
rivS_index2nd_kapa <- bray3[11:15,11:15][lower.tri(bray3[11:15,11:15])]
plat_only3 <- as.dist(bray3[plat_id, plat_id])
id2nd_only3 <- as.dist(bray3[index2nd_id, index2nd_id])
new_old3 <- as.dist(bray3[newold_id,newold_id])

stdM_index1st_plat <- bray4[  1:5,  1:5][lower.tri(bray4[  1:5,  1:5])]
stdM_index2nd_plat <- bray4[ 6:10, 6:10][lower.tri(bray4[ 6:10, 6:10])]
stdM_index2nd_kapa <- bray4[11:15,11:15][lower.tri(bray4[11:15,11:15])]
plat_only4 <- as.dist(bray4[plat_id, plat_id])
id2nd_only4 <- as.dist(bray4[index2nd_id, index2nd_id])
new_old4 <- as.dist(bray4[newold_id,newold_id])

# ----------------------------------------------- #
#         Generate a new data.frame
# ----------------------------------------------- #
# Within method variations
bray_df <- data.frame(data_id = sprintf("ID%03d", 1:120),
                      rep_id = rep(sprintf("Rep%02d", 1:10), 12),
                      site = c(rep("Sea_Nagahama", 30), rep("Sea_Otomi", 30), rep("River_Seta", 30), rep("STD_Mix", 30)),
                      index_method = rep(c(rep("1stPCR_indexing", 10), rep("2ndPCR_indexing", 20)), 4),
                      enzyme = rep(c(rep("Platinum", 20), rep("KAPA", 10)), 4))
bray_df$libprep <- paste0(bray_df$index_method, "_", bray_df$enzyme)
bray_df$bray_curtis <- c(seaN_index1st_plat, seaN_index2nd_plat, seaN_index2nd_kapa,
                         seaO_index1st_plat, seaO_index2nd_plat, seaO_index2nd_kapa,
                         rivS_index1st_plat, rivS_index2nd_plat, rivS_index2nd_kapa,
                         stdM_index1st_plat, stdM_index2nd_plat, stdM_index2nd_kapa)


# ----------------------------------------------- #
#         Visualize pattern
# ----------------------------------------------- #
g1 <- bray_df %>% ggplot(aes(x = libprep, y = bray_curtis)) +
  geom_boxplot(outlier.color = NA, outlier.size = 0) +
  geom_jitter(height = 0, width = 0.2) +
  facet_wrap(~ site, ncol = 2) +
  panel_border() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  xlab(NULL) + ylab("Bray-Curtis dissimilarity") +
  ylim(0, 1) +
  NULL


# ----------------------------------------------- #
#    Statistics (GLM, ANOSIM)
# ----------------------------------------------- #
# Brief statistics
libprep_id_all <- c(rep("1stPCR_indexing_Platinum", 5),
                rep("2ndPCR_indexing_Platinum", 5),
                rep("2ndPCR_indexing_KAPA", 5))
plat_id_all <- c(rep("1stPCR_indexing_Platinum", 5),
                    rep("2ndPCR_indexing_Platinum", 5))
index2nd_id_all <- c(rep("2ndPCR_indexing_Platinum", 5),
                    rep("2ndPCR_indexing_KAPA", 5))
## 1st PCR indexing v.s. 2nd PCR indexing (Platinum)
set.seed(ran.seed); vegan::anosim(bray1, libprep_id_all) # Sea_Nagahama, P = 0.019
set.seed(ran.seed); vegan::anosim(bray2, libprep_id_all) # Sea_Otomi, P = 0.632
set.seed(ran.seed); vegan::anosim(bray3, libprep_id_all) # River_Seta, P = 0.981
set.seed(ran.seed); vegan::anosim(bray4, libprep_id_all) # STD_Mix, P = 0.221
## Platinum v.s. KAPA (2nd PCR indexing)
set.seed(ran.seed); vegan::anosim(id2nd_only1, index2nd_id_all) # Sea_Nagahama, P = 0.065
set.seed(ran.seed); vegan::anosim(id2nd_only2, index2nd_id_all) # Sea_Otomi, P = 0.484
set.seed(ran.seed); vegan::anosim(id2nd_only3, index2nd_id_all) # River_Seta, P = 0.974
set.seed(ran.seed); vegan::anosim(id2nd_only4, index2nd_id_all) # STD_Mix, P = 0.698
## 1st PCR indexing v.s. 2nd PCR indexing (Platinum)
set.seed(ran.seed); vegan::anosim(plat_only1, plat_id_all) # Sea_Nagahama, P = 0.440
set.seed(ran.seed); vegan::anosim(plat_only2, plat_id_all) # Sea_Otomi, P = 0.757
set.seed(ran.seed); vegan::anosim(plat_only3, plat_id_all) # River_Seta, P = 0.935
set.seed(ran.seed); vegan::anosim(plat_only4, plat_id_all) # STD_Mix, P = 0.109

## Difference in Bray-Curtis dissimilarity between 1st PCR indexing v.s. 2nd PCR indexing
seaN_tr <- subset(bray_df, site == "Sea_Nagahama") %>% mutate(log_bray = -log(bray_curtis))
seaO_tr <- subset(bray_df, site == "Sea_Otomi") %>% mutate(log_bray = -log(bray_curtis))
rivS_tr <- subset(bray_df, site == "River_Seta") %>% mutate(log_bray = -log(bray_curtis))
stdM_tr <- subset(bray_df, site == "STD_Mix") %>% mutate(log_bray = -log(bray_curtis))

# Perform GLM using gamma distribution (log link)
summary(glm(log_bray ~ index_method + enzyme, data = seaN_plat_tr, family = Gamma(link = "log")))
## Sea_Nagahama, 2ndPCR_indexing, P = 0.00906, *
summary(glm(log_bray ~ index_method + enzyme, data = seaO_tr, family = Gamma(link = "log")))
## Sea_Otomi, 2ndPCR_indexing, P = 0.4566, NS
summary(glm(log_bray ~ index_method + enzyme, data = rivS_tr, family = Gamma(link = "log")))
## River_Seta, 2ndPCR_indexing, P = 0.3831, NS
summary(glm(log_bray ~ index_method + enzyme, data = stdM_tr, family = Gamma(link = "log")))
## STD_Mix, 2ndPCR_indexing, P < 0.0001, ***

# ----------------------------------------------- #
#   Custom bootstrap test
# ----------------------------------------------- #
bray_boot <- function(bray_df_object, n_iter = 1000, ran.seed = 1234) {
  # Convert to matrix
  bray_df_object <- as.matrix(bray_df_object)
  # Prepare output object
  diff_all <- c(NULL)
  for(i in 1:n_iter) {
    # Generate random index
    # Labels "1st PCR indexing" or "2nd PCR indexing" were randomly shuffled
    # And the resultant difference between bray-curtis values was compared with the original one
    set.seed(ran.seed+i); (index_id_boot <- sample(10, 10, replace = FALSE))
    # Shuffle the distance matrix
    bray_df_boot <- bray_df_object[index_id_boot, index_id_boot]
    # Extract (pseudo) treatment
    bray_1st_df <- bray_df_boot[1:5,1:5]
    bray_2nd_df <- bray_df_boot[6:10,6:10]
    bray_vals <- c(bray_1st_df[lower.tri(bray_1st_df)], bray_2nd_df[lower.tri(bray_2nd_df)])
    bray_df_boot2 <- data.frame(method = c(rep("method1", 10), rep("method2", 10)),
                                values = bray_vals) %>% tibble()
    # Calculate mean difference
    diff_tmp <- diff(as.numeric(unlist(bray_df_boot2 %>% group_by(method) %>% summarize(mean_bray = mean(values)))[3:4]))
    diff_all <- c(diff_all, abs(diff_tmp))
  }
  # Get original values
  bray_1st_df0 <- bray_df_object[1:5,1:5]
  bray_2nd_df0 <- bray_df_object[6:10,6:10]
  bray_vals0 <- c(bray_1st_df0[lower.tri(bray_1st_df0)], bray_2nd_df0[lower.tri(bray_2nd_df0)])
  bray_df_boot0 <- data.frame(method = c(rep("method1", 10), rep("method2", 10)),
                              values = bray_vals0) %>% tibble()
  # Calculate mean difference
  obs_diff <- abs(diff(as.numeric(unlist(bray_df_boot0 %>% group_by(method) %>% summarize(mean_bray = mean(values)))[3:4])))
  quant_diff <- quantile(diff_all, probs = c(0.9, 0.95, 0.975, 0.99, 0.995, 0.999))
  boot_p_val <- sum(diff_all - obs_diff > 0)/n_iter
  res_diff <- c(quant_diff, obs_diff, boot_p_val)
  names(res_diff)[(length(res_diff)-1):length(res_diff)] <- c("observed", "p_val")
  return(res_diff)
}

# Test Bary-Curtis dissimilarity values
## 1st_indexing_Platinum v.s. 2nd_indexing_KAPA
(seaN_method_ci <- bray_boot(new_old1)) # NS, P = 0.291
(seaO_method_ci <- bray_boot(new_old2)) # NS, P = 0.811
(rivS_method_ci <- bray_boot(new_old3)) # NS, P = 0.607
(stdM_method_ci <- bray_boot(new_old4)) # NS, P = 0.065 -
## 2nd_indexing, KAPA v.s. Platinum
(seaN_id2nd_ci <- bray_boot(id2nd_only1)) # NS, P = 0.945
(seaO_id2nd_ci <- bray_boot(id2nd_only2)) # NS, P = 0.946
(rivS_id2nd_ci <- bray_boot(id2nd_only3)) # NS, P = 0.846
(stdM_id2nd_ci <- bray_boot(id2nd_only4)) # NS, P = 0.065 -
## Platinum, 1st_indexing v.s. 2nd_indexing
(seaN_plat_ci <- bray_boot(plat_only1)) # NS, P = 0.308
(seaO_plat_ci <- bray_boot(plat_only2)) # NS, P = 0.759
(rivS_plat_ci <- bray_boot(plat_only3)) # NS, P = 0.603
(stdM_plat_ci <- bray_boot(plat_only4)) # Significant, P < 0.001


# ----------------------------------------------- #
#  Relative abundance of the most dominant taxon
# ----------------------------------------------- #
top_rank_list <- list()
top_taxa_list <- list()
for(i in 1:5){
  # Set top rank
  top_rank <- i
  top_sp1 <- names(sort(taxa_sums(ps_fish1), decreasing = T))[top_rank]
  top_sp2 <- names(sort(taxa_sums(ps_fish2), decreasing = T))[top_rank]
  top_sp3 <- names(sort(taxa_sums(ps_fish3), decreasing = T))[top_rank]
  top_sp4 <- names(sort(taxa_sums(ps_fish4), decreasing = T))[top_rank]
  
  ## Generate a new data.frame
  new_df1 <- cbind(sample_data(ps_fish1), c(otu_table(ps_fish1)[,top_sp1]))
  new_df2 <- cbind(sample_data(ps_fish2), c(otu_table(ps_fish2)[,top_sp2]))
  new_df3 <- cbind(sample_data(ps_fish3), c(otu_table(ps_fish3)[,top_sp3]))
  new_df4 <- cbind(sample_data(ps_fish4), c(otu_table(ps_fish4)[,top_sp4]))
  colnames(new_df1)[dim(new_df1)[2]] <-
    colnames(new_df2)[dim(new_df1)[2]] <-
    colnames(new_df3)[dim(new_df1)[2]] <-
    colnames(new_df4)[dim(new_df1)[2]] <- "top_prop"
  rel_df <- rbind(new_df1[,c("Sample_ID", "replicate", "site", "index_method", "enzyme", "libprep", "top_prop")],
                  new_df2[,c("Sample_ID", "replicate", "site", "index_method", "enzyme", "libprep", "top_prop")],
                  new_df3[,c("Sample_ID", "replicate", "site", "index_method", "enzyme", "libprep", "top_prop")],
                  new_df4[,c("Sample_ID", "replicate", "site", "index_method", "enzyme", "libprep", "top_prop")])
  
  ## Save figure
  (g4 <- rel_df %>% ggplot(aes(x = libprep, y = top_prop)) +
      geom_boxplot(outlier.color = NA, outlier.size = 0) +
      geom_jitter(height = 0, width = 0.2) +
      facet_wrap(~ site, ncol = 4) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
      xlab(NULL) + ylab("Proportion of the top taxon") +
      panel_border() +
      NULL)
  ggsave(file = sprintf("%s/TopAbundance_%02d.pdf", output_folder, top_rank),
         plot = g4, width = 8, height = 8)
  top_info1 <- c(top_sp1, tax_table(ps_exp1_even)[top_sp1,c("family","genus")])
  top_info2 <- c(top_sp2, tax_table(ps_exp1_even)[top_sp2,c("family","genus")])
  top_info3 <- c(top_sp3, tax_table(ps_exp1_even)[top_sp3,c("family","genus")])
  top_info4 <- c(top_sp4, tax_table(ps_exp1_even)[top_sp4,c("genus","species")])
  name_df <- data.frame(query = c(top_info1[1], top_info2[1], top_info3[1], top_info4[1]),
                        info1 = c(top_info1[2], top_info2[2], top_info3[2], top_info4[2]),
                        info2 = c(top_info1[3], top_info2[3], top_info3[3], top_info4[3]))
  # Collect ggplot objects
  top_rank_list <- c(top_rank_list, list(g4))
  top_taxa_list <- c(top_taxa_list, list(name_df))
}


# ----------------------------------------------- #
#         Save figures
# ----------------------------------------------- #
# Save figure objects
#dir.create("../FigCode"); dir.create("../FigCode/00_RawFigs")
fig_dir <- "../FigCode/00_RawFigs/"
saveRDS(g1, paste0(fig_dir, "8_4_BrayCurtis.obj"))
saveRDS(top_rank_list, paste0(fig_dir, "8_4_TopProportion.obj"))
saveRDS(top_taxa_list, paste0(fig_dir, "8_4_TopPropTaxaInfo.obj"))

# Save session info
writeLines(capture.output(sessionInfo()),
           paste0("../00_SessionInfo/", output_folder, "_", substr(Sys.time(), 1, 10), ".txt"))

