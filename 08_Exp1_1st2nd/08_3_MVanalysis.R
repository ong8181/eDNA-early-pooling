####
#### Multivariate analysis
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
#         Split by site
# ----------------------------------------------- #
ps_fish1 <- subset_samples(ps_exp1_even, site == "Sea_Nagahama") %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  transform_sample_counts(function(x) x/sum(x))
ps_fish2 <- subset_samples(ps_exp1_even, site == "Sea_Otomi") %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  transform_sample_counts(function(x) x/sum(x))
ps_fish3 <- subset_samples(ps_exp1_even, site == "River_Seta") %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  transform_sample_counts(function(x) x/sum(x))
ps_fish4 <- subset_samples(ps_exp1_even, site == "STD_Mix") %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  transform_sample_counts(function(x) x/sum(x))
sample_sums(ps_fish1); sample_sums(ps_fish2)
sample_sums(ps_fish3); sample_sums(ps_fish4)


# ----------------------------------------------- #
#         NMDS
# ----------------------------------------------- #
# Nonmetric multidimensional scaling
set.seed(ran.seed); ps_bray1 <- ordinate(ps_fish1, "NMDS", "bray")
set.seed(ran.seed); ps_bray2 <- ordinate(ps_fish2, "NMDS", "bray")
set.seed(ran.seed); ps_bray3 <- ordinate(ps_fish3, "NMDS", "bray")
set.seed(ran.seed); ps_bray4 <- ordinate(ps_fish4, "NMDS", "bray", trymax = 100)

n1 <- plot_ordination(ps_fish1, ps_bray1, shape = "enzyme", color = "index_method") +
  stat_ellipse(geom = "polygon", alpha = 0.1, aes(fill=index_method)) +
  geom_point(size = 2) + scale_color_startrek() + scale_fill_startrek() +
  xlab("Axis 1") + ylab("Axis 2") + ggtitle(NULL) + facet_wrap(~ site)
n2 <- plot_ordination(ps_fish2, ps_bray2, shape = "enzyme", color = "index_method") +
  stat_ellipse(geom = "polygon", alpha = 0.1, aes(fill=index_method)) +
  geom_point(size = 2) + scale_color_startrek() + scale_fill_startrek() +
  xlab("Axis 1") + ylab("Axis 2") + ggtitle(NULL) + facet_wrap(~ site)
n3 <- plot_ordination(ps_fish3, ps_bray3, shape = "enzyme", color = "index_method") +
  stat_ellipse(geom = "polygon", alpha = 0.1, aes(fill=index_method)) +
  geom_point(size = 2) + scale_color_startrek() + scale_fill_startrek() +
  xlab("Axis 1") + ylab("Axis 2") + ggtitle(NULL) + facet_wrap(~ site)
n4 <- plot_ordination(ps_fish4, ps_bray4, shape = "enzyme", color = "index_method") +
  stat_ellipse(geom = "polygon", alpha = 0.1, aes(fill=index_method)) +
  geom_point(size = 2) + scale_color_startrek() + scale_fill_startrek() +
  xlab("Axis 1") + ylab("Axis 2") + ggtitle(NULL) + facet_wrap(~ site)
n_all <- plot_grid(n1, n2, n3, n4, ncol = 2, align = "hv")


# ----------------------------------------------- #
#         Save figures
# ----------------------------------------------- #
# Output figures
ggsave(file = sprintf("%s/NMDS_all.pdf", output_folder),
       plot = n_all, width = 12, height = 8)

# Save figure objects
fig_dir <- "../FigCode/00_RawFigs/"
saveRDS(list(n1, n2, n3, n4), paste0(fig_dir, "8_3_NMDS_all.obj"))

# Save session info
writeLines(capture.output(sessionInfo()),
           paste0("../00_SessionInfo/", output_folder, "_", substr(Sys.time(), 1, 10), ".txt"))

