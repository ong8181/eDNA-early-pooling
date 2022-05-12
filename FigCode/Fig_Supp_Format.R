####
#### Format figures for eDNA-seq
#### 2022.05.12 revision for Environmental DNA
#### R 4.1.2
####

# Set working directory
if(basename(getwd()) != "FigCode") setwd("FigCode")

# Load libraries
## For general
library(tidyverse); packageVersion("tidyverse") # 1.3.1, 2021.8.25
library(lubridate); packageVersion("lubridate") # 1.8.0, 2021.11.10
## For image
library(ggimage); packageVersion("ggimage") # 0.3.0, 2021.12.8
library(magick); packageVersion("magick") # 2.7.3, 2021.12.8
## For ggplot
library(cowplot); packageVersion("cowplot") # 1.1.1, 2021.3.1
library(ggrepel); packageVersion("ggrepel") # 0.9.1, 2021.11.29
library(ggsci); packageVersion("ggsci") # 2.9, 2021.8.26
theme_set(theme_cowplot())

# Prepare output folders
fig_dir <- "00_RawFigs/"
fig_dir_out <- "00_FormattedFigs/"
# Generate output directory
dir.create(fig_dir_out)

# Prepare color palette
library(RColorBrewer); packageVersion("RColorBrewer") # 1.1.3, 2022.5.5
get_palette <- colorRampPalette(brewer.pal(8, "Paired"))


# <---------------------------------------------> #
#  Load figure objects
# <---------------------------------------------> #
# phyloseq object to check taxa
ps_fish <- readRDS("../07_CompilePhyloseqOut/ps_fish_all.obj")
taxa1 <- phyloseq::tax_table(ps_fish)
# Additional analyses
fig_11_sameidbox <- readRDS(paste0(fig_dir, "11_1_Fig_SameIndexBoxplot.obj"))
fig_11_sameidcomp <- readRDS(paste0(fig_dir, "11_1_Fig_SameIndexComp_SeaN.obj"))
fig_12_nagahama <- readRDS(paste0(fig_dir, "12_1_Fig_NagahamaNMDS_all.obj"))
fig_14_contamrate <- readRDS(paste0(fig_dir, "14_1_Fig_ContamRate_Exp2_Fxp3.obj"))


# <---------------------------------------------> #
#  Same index test
# <---------------------------------------------> #
## Relabel
taxa1["OTU002"]; taxa1["OTU007"]; taxa1["OTU056"]
## Prepare annotation
id_text <- data.frame(x = 1.5, y = 0.75, lab = "Tag effect: N.S.")
## Same index clustering
fig_same_index <- plot_grid(fig_11_sameidbox[[1]] + ggtitle("Experiment I"),
                            fig_11_sameidbox[[2]] + ggtitle("Experiment II"),
                            fig_11_sameidbox[[3]] + ggtitle("Experiment III") +
                              geom_text(data = id_text, aes(x = x, y = y, label = lab)),
                            align = "hv", axis = "lrtb",
                            nrow = 1, rel_widths = c(3,1,5))
## Statistics
df3 <- fig_11_sameidbox[[3]]$data
summary(aov(Abundance ~ I7_Index_ID2, data = df3 %>% filter(OTU == "OTU002")))
summary(aov(Abundance ~ I7_Index_ID2, data = df3 %>% filter(OTU == "OTU007")))
summary(aov(Abundance ~ I7_Index_ID2, data = df3 %>% filter(OTU == "OTU056")))


# <---------------------------------------------> #
#  Sea_Nagahama test
# <---------------------------------------------> #
## Relabel
df12 <- fig_12_nagahama[[1]]$data
### Experiment
df12$Experiment_ID[df12$Experiment_ID == "E01"] <- "Exp. I"
df12$Experiment_ID[df12$Experiment_ID == "E02"] <- "Exp. II"
df12$Experiment_ID[df12$Experiment_ID == "E03"] <- "Exp. III"
### ExoSAP-IT
cat2 <- unique(df12$test_category2)
df12$test_category2[df12$test_category2 == cat2[1]] <- "2-µl"
df12$test_category2[df12$test_category2 == cat2[2]] <- "1-µl"
df12$test_category2[df12$test_category2 == cat2[3]] <- "w/ exonuclease"
df12$test_category2[df12$test_category2 == cat2[4]] <- "1st-tagging, Platinum"
df12$test_category2[df12$test_category2 == cat2[5]] <- "4-µl"
df12$test_category2[df12$test_category2 == cat2[6]] <- "8-µl"
df12$test_category2[df12$test_category2 == cat2[7]] <- "1-rep."
df12$test_category2[df12$test_category2 == cat2[8]] <- "2-rep."
df12$test_category2[df12$test_category2 == cat2[9]] <- "4-rep."
df12$test_category2[df12$test_category2 == cat2[10]] <- "w/o exonuclease"
df12$test_category2[df12$test_category2 == cat2[11]] <- "2nd-indexing, KAPA"
df12$test_category2[df12$test_category2 == cat2[12]] <- "2nd-indexing, Platinum"
df12$test_category2[df12$test_category2 == cat2[13]] <- "8-rep."
### Re-import df12
fig_12_nagahama[[1]]$data <- df12

# Re-label fig_12_nagahama[[2]]
df12_2 <- fig_12_nagahama[[2]]$data
cat12_2 <- unique(df12_2$test_category)
df12_2$test_category[df12_2$test_category == cat12_2[1]] <- "1st_tagging_Platinum"
cat12_3 <- unique(df12_2$test_category)
df12_2$test_category <- factor(df12_2$test_category, levels = cat12_3[1:23])
fig_12_nagahama[[2]]$data <- df12_2

fig_nagahama <- plot_grid(fig_12_nagahama[[1]] + scale_fill_igv(name = "Family"),
                          fig_12_nagahama[[2]],
                          ncol = 1, labels = letters[1:2])


# <---------------------------------------------> #
# Combined to the same index figure
# <---------------------------------------------> #
df13 <- fig_12_nagahama[[3]]$data
cat3 <- unique(df13$test_category2)
df13$test_category2[df13$test_category2 == cat3[1]] <- "1st PCR tagging, w/o exonuclease"
df13$test_category2[df13$test_category2 == cat3[2]] <- "2nd PCR indexing, w/o exonuclease"
df13$test_category2[df13$test_category2 == cat3[3]] <- "1st PCR tagging, w exonuclease"
fig_12_nagahama[[3]]$data <- df13

fig_sameid_nagahama <- plot_grid(fig_same_index,
                                 fig_12_nagahama[[3]],
                                 ncol = 1, labels = letters[1:2],
                                 rel_heights = c(1, 0.8))


# <---------------------------------------------> #
# Save figures
# <---------------------------------------------> #
# Main figures
## Figure 5
ggsave(file = sprintf("%s/Figure_05.pdf", fig_dir_out),
       plot = fig_sameid_nagahama, width = 16, height = 14)


# Supplementary figures
## Figure S5
ggsave(file = sprintf("%s/Figure_S05.pdf", fig_dir_out),
       plot = fig_nagahama, width = 28, height = 22)


