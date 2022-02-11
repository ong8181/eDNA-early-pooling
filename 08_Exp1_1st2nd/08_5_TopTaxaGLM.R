####
#### No.8.4 Analyzing top taxa
#### R 4.1.2
####

# Set random seeds (for reproduction)
ran.seed <- 1234
set.seed(ran.seed)

# Load library and functions
library(tidyverse); packageVersion("tidyverse") # 1.3.1, 2021.10.16
library(phyloseq); packageVersion("phyloseq") # 1.38.0, 2021.11.18
library(cowplot); packageVersion("cowplot") # 1.1.1, 2021.6.13
theme_set(theme_cowplot())

# Generate output folder
od <- basename(rstudioapi::getSourceEditorContext()$path)
(output_folder <- paste0(str_sub(od, end = -3), "Out")); rm(od)
dir.create(output_folder)


# ----------------------------------------------- #
# Load objects
# ----------------------------------------------- #
ps_exp1_even <- readRDS("08_2_EvenDepthOut/ps_exp1_even.obj")
top_rank_list <- readRDS("../FigCode/00_RawFigs/8_4_TopProportion.obj")


# ----------------------------------------------- #
# Compile data
# ----------------------------------------------- #
## Sea_Nagahama
ps_fish1 <- subset_samples(ps_exp1_even, site == "Sea_Nagahama") %>% prune_taxa(taxa_sums(.) > 0, .)
## Sea_Otomi
ps_fish2 <- subset_samples(ps_exp1_even, site == "Sea_Otomi") %>% prune_taxa(taxa_sums(.) > 0, .)
## River_Seta
ps_fish3 <- subset_samples(ps_exp1_even, site == "River_Seta") %>% prune_taxa(taxa_sums(.) > 0, .)
## STD_Mix
ps_fish4 <- subset_samples(ps_exp1_even, site == "STD_Mix") %>% prune_taxa(taxa_sums(.) > 0, .)
sample_sums(ps_fish1); sample_sums(ps_fish2); sample_sums(ps_fish3); sample_sums(ps_fish4)

# Define function
add_reads <- function (list_df) {
  list_df$reads <- NA
  list_df$reads[list_df$site == "Sea_Nagahama"] <- list_df$top_prop[list_df$site == "Sea_Nagahama"] * unique(sample_sums(ps_fish1))
  list_df$reads[list_df$site == "Sea_Otomi"] <- list_df$top_prop[list_df$site == "Sea_Otomi"] * unique(sample_sums(ps_fish2))
  list_df$reads[list_df$site == "River_Seta"] <- list_df$top_prop[list_df$site == "River_Seta"] * unique(sample_sums(ps_fish3))
  list_df$reads[list_df$site == "STD_Mix"] <- list_df$top_prop[list_df$site == "STD_Mix"] * unique(sample_sums(ps_fish4))
  return(list_df)
}

# Add sequence reads information
top_list1 <- add_reads(top_rank_list[[1]]$data)
top_list2 <- add_reads(top_rank_list[[2]]$data)
top_list3 <- add_reads(top_rank_list[[3]]$data)
top_list4 <- add_reads(top_rank_list[[4]]$data)
top_list5 <- add_reads(top_rank_list[[5]]$data)
# Top 1
top_list1_seaN <- top_list1 %>% filter(site == "Sea_Nagahama") %>% droplevels()
top_list1_seaO <- top_list1 %>% filter(site == "Sea_Otomi") %>% droplevels()
top_list1_rivS <- top_list1 %>% filter(site == "River_Seta") %>% droplevels()
top_list1_stdM <- top_list1 %>% filter(site == "STD_Mix") %>% droplevels()
# Top 2
top_list2_seaN <- top_list2 %>% filter(site == "Sea_Nagahama") %>% droplevels()
top_list2_seaO <- top_list2 %>% filter(site == "Sea_Otomi") %>% droplevels()
top_list2_rivS <- top_list2 %>% filter(site == "River_Seta") %>% droplevels()
top_list2_stdM <- top_list2 %>% filter(site == "STD_Mix") %>% droplevels()
# Top 3
top_list3_seaN <- top_list3 %>% filter(site == "Sea_Nagahama") %>% droplevels()
top_list3_seaO <- top_list3 %>% filter(site == "Sea_Otomi") %>% droplevels()
top_list3_rivS <- top_list3 %>% filter(site == "River_Seta") %>% droplevels()
top_list3_stdM <- top_list3 %>% filter(site == "STD_Mix") %>% droplevels()


# ----------------------------------------------- #
# GLM
# ----------------------------------------------- #
## Top 1
TukeyHSD(aov(reads ~ libprep, data = top_list1_seaN)) # N.S.
TukeyHSD(aov(reads ~ libprep, data = top_list1_seaO)) # N.S.
TukeyHSD(aov(reads ~ libprep, data = top_list1_rivS)) # N.S.
TukeyHSD(aov(reads ~ libprep, data = top_list1_stdM)) # N.S.
## Top 2
TukeyHSD(aov(reads ~ libprep, data = top_list2_seaN)) # 2ndKAPA - 1stPlat: P < 0.01 (0.0026), 2ndPlat - 1stPlat; P = 0.080
TukeyHSD(aov(reads ~ libprep, data = top_list2_seaO)) # N.S.
TukeyHSD(aov(reads ~ libprep, data = top_list2_rivS)) # N.S.
TukeyHSD(aov(reads ~ libprep, data = top_list2_stdM)) # N.S.
## Top 3
TukeyHSD(aov(reads ~ libprep, data = top_list3_seaN)) # N.S.
TukeyHSD(aov(reads ~ libprep, data = top_list3_seaO)) # N.S.
TukeyHSD(aov(reads ~ libprep, data = top_list3_rivS)) # N.S.
TukeyHSD(aov(reads ~ libprep, data = top_list3_stdM)) # N.S.


# ----------------------------------------------- #
#  Save results
# ----------------------------------------------- #
# Save session info
writeLines(capture.output(sessionInfo()),
           paste0("../00_SessionInfo/", output_folder, "_", substr(Sys.time(), 1, 10), ".txt"))

