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


# --------------------------------------------------------- #
#  Load figure objects
# --------------------------------------------------------- #
# Experiment Design
fig_exp_design <- ggdraw() + draw_image(image_read("00_ExpDesign/eDNAseq_ExpDesign.jpg"))
# Experiment I
fig_exp1_summary <- readRDS(paste0(fig_dir, "8_1_Fig_Exp1_SummaryReads.obj"))
fig_exp1_compdiv <- readRDS(paste0(fig_dir, "8_2_Fig_Exp1_CompositionDiveristy.obj"))
fig_exp1_nmds <- readRDS(paste0(fig_dir, "8_3_NMDS_all.obj"))
fig_exp1_bray <- readRDS(paste0(fig_dir, "8_4_BrayCurtis.obj"))
fig_exp1_topprop <- readRDS(paste0(fig_dir, "8_4_TopProportion.obj"))
tab_exp1_topprop <- readRDS(paste0(fig_dir, "8_4_TopPropTaxaInfo.obj"))
## For relabeling
var_char1 <- c("1stPCR_indexing", "2ndPCR_indexing")
## Label
discrete_label <- c("2nd PCR\nindexing\nKAPA",
                    "2nd PCR\nindexing\nPlatinum",
                    "1st PCR\ntagging\nPlatinum")
site_level <- c("Sea_Nagahama", "Sea_Otomi", "River_Seta", "STD_Mix")


# --------------------------------------------------------- #
#   NMDS
# --------------------------------------------------------- #
### Relabeling
fig_exp1_nmds[[1]] <- fig_exp1_nmds[[1]] + labs(color = "Method", shape = "Enzyme", fill = "Method")
fig_exp1_nmds[[1]]$data$index_method[fig_exp1_nmds[[1]]$data$index_method == var_char1[1]] <- "1st PCR tagging"
fig_exp1_nmds[[1]]$data$index_method[fig_exp1_nmds[[1]]$data$index_method == var_char1[2]] <- "2nd PCR indexing"
nmds_legend <- get_legend(fig_exp1_nmds[[1]])
fig_exp1_nmds_all <- plot_grid(fig_exp1_nmds[[1]] + theme(legend.position = "none"),
                               fig_exp1_nmds[[2]] + theme(legend.position = "none"),
                               nmds_legend,
                               fig_exp1_nmds[[3]] + theme(legend.position = "none"),
                               fig_exp1_nmds[[4]] + theme(legend.position = "none"),
                               nmds_legend,
                               align = "hv", axis = "lrbt", rel_widths = c(1,1,0.4,1,1,0.4),
                               ncol = 3,
                               labels = c("a","b",NA,"c","d",NA))


# --------------------------------------------------------- #
#   Barplot and diversity
# --------------------------------------------------------- #
## Relabeling
fig_exp1_compdiv[[1]]$data$test_category <- factor(fig_exp1_compdiv[[1]]$data$test_category,
                                                   levels = c("2nd_indexing_KAPA", "2nd_indexing_Platinum", "1st_indexing_Platinum"))
fig_exp1_compdiv[[1]]$data$site <- factor(fig_exp1_compdiv[[1]]$data$site, levels = site_level)
fig_exp1_compdiv[[2]]$data$index_method[fig_exp1_compdiv[[2]]$data$index_method == var_char1[1]] <- "1st PCR tagging"
fig_exp1_compdiv[[2]]$data$index_method[fig_exp1_compdiv[[2]]$data$index_method == var_char1[2]] <- "2nd PCR indexing"
fig_exp1_compdiv[[2]]$data$site <- factor(fig_exp1_compdiv[[2]]$data$site, levels = c("Sea_Nagahama", "Sea_Otomi", "River_Seta", "STD_Mix"))
fig_exp1_compdiv[[2]]$data$index_method <- factor(fig_exp1_compdiv[[2]]$data$index_method, levels = c("2nd PCR indexing", "1st PCR tagging"))
fig_exp1_compdiv[[2]]$data$enzyme <- factor(fig_exp1_compdiv[[2]]$data$enzyme, levels = c("KAPA", "Platinum"))
## Prepare N.S. text for diversity plot
text_div <- data.frame(x = 2.3, y = 32, lab = "N.S.")

## Composition and diversity
fig_exp1_bar_div <- plot_grid(fig_exp1_compdiv[[1]] + # For OTU richness, "methods" is not significant
                                scale_x_discrete(labels = discrete_label) +
                                theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
                                facet_wrap(~ site, ncol = 4) + ylim(0,32) +
                                geom_text(data = text_div, aes(x = x, y = y, label = lab)) +
                                ylab("No. of fish OTUs"),
                              fig_exp1_compdiv[[2]] + scale_fill_igv(name = "Family \nor STD name"),
                              ncol = 1, labels = letters[1:2],
                              align = "hv", axis = "lrbt",
                              rel_heights = c(0.4,1))


# --------------------------------------------------------- #
#  Bray-Curtis
# --------------------------------------------------------- #
### Relabel
fig_exp1_bray$data$site <- factor(fig_exp1_bray$data$site, levels = c("Sea_Nagahama", "Sea_Otomi", "River_Seta", "STD_Mix"))
fig_exp1_bray$data$libprep <- factor(fig_exp1_bray$data$libprep,
                                     levels = c("2ndPCR_indexing_KAPA",
                                                "2ndPCR_indexing_Platinum",
                                                "1stPCR_indexing_Platinum"))
### Add significance labels
text_bray1 <- data.frame(x = 2, y = 1, lab = c("N.S."), site = factor("Sea_Nagahama", levels = levels(fig_exp1_bray$data$site)))
text_bray2 <- data.frame(x = 2, y = 1, lab = c("N.S."), site = factor("Sea_Otomi", levels = levels(fig_exp1_bray$data$site)))
text_bray3 <- data.frame(x = 2, y = 1, lab = c("N.S."), site = factor("River_Seta", levels = levels(fig_exp1_bray$data$site)))
text_bray4 <- data.frame(x = c(1.5, 2.0, 2.5),
                         y = c(0.55, 0.65, 0.75),
                         lab = c("italic(P) == 0.065", "italic(P) == 0.065", "italic(P) < 0.001"),
                         site = factor("STD_Mix", levels = levels(fig_exp1_bray$data$site)))
line_bray4 <- data.frame(x = c(1, 1, 2), xend = c(2, 3, 3),
                         y = c(0.52, 0.62, 0.72), yend = c(0.52, 0.62, 0.72),
                         site = factor("STD_Mix", levels = levels(fig_exp1_bray$data$site)))

fig_exp1_bray2 <- fig_exp1_bray +
  scale_x_discrete(labels = discrete_label) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  facet_wrap(~ site, ncol = 4) +
  geom_text(data = text_bray1, aes(x = x, y = y, label = lab), size = 4) +
  geom_text(data = text_bray2, aes(x = x, y = y, label = lab), size = 4) +
  geom_text(data = text_bray3, aes(x = x, y = y, label = lab), size = 4) +
  geom_text(data = text_bray4, aes(x = x, y = y, label = lab), parse = T, size = 4) +
  geom_segment(data = line_bray4, aes(x = x, xend = xend, y = y, yend = yend)) +
  NULL


# --------------------------------------------------------- #
#  Combine figures from Experiment I
# --------------------------------------------------------- #
fig_exp1_main <- plot_grid(fig_exp1_nmds_all,
                           fig_exp1_bray2,
                           ncol = 1, rel_heights = c(1,0.6),
                           labels = c(NA, "e"))
                           

# --------------------------------------------------------- #
#  Supplementary Figures
# --------------------------------------------------------- #
## Exp1 sequence reads
fig_exp1_summary$data$index_method[fig_exp1_summary$data$index_method == var_char1[1]] <- "1st PCR tagging"
fig_exp1_summary$data$index_method[fig_exp1_summary$data$index_method == var_char1[2]] <- "2nd PCR indexing"
fig_exp1_summary$data$index_method <- factor(fig_exp1_summary$data$index_method, levels = c("2nd PCR indexing", "1st PCR tagging"))
fig_exp1_summary$data$enzyme <- factor(fig_exp1_summary$data$enzyme, levels = c("KAPA", "Platinum"))
fig_exp1_summary <- fig_exp1_summary + scale_fill_igv(name = "Family")

## Top proportion test
fig_exp1_topprop[[1]]$data$site <- factor(fig_exp1_topprop[[1]]$data$site, levels = c("Sea_Nagahama", "Sea_Otomi", "River_Seta", "STD_Mix"))
fig_exp1_topprop[[2]]$data$site <- factor(fig_exp1_topprop[[2]]$data$site, levels = c("Sea_Nagahama", "Sea_Otomi", "River_Seta", "STD_Mix"))
fig_exp1_topprop[[3]]$data$site <- factor(fig_exp1_topprop[[3]]$data$site, levels = c("Sea_Nagahama", "Sea_Otomi", "River_Seta", "STD_Mix"))
### Change variable order
fig_exp1_topprop[[1]]$data$libprep <- factor(fig_exp1_topprop[[1]]$data$libprep, levels = c("2ndPCR_indexing_KAPA", "2ndPCR_indexing_Platinum", "1stPCR_indexing_Platinum"))
fig_exp1_topprop[[2]]$data$libprep <- factor(fig_exp1_topprop[[2]]$data$libprep, levels = c("2ndPCR_indexing_KAPA", "2ndPCR_indexing_Platinum", "1stPCR_indexing_Platinum"))
fig_exp1_topprop[[3]]$data$libprep <- factor(fig_exp1_topprop[[3]]$data$libprep, levels = c("2ndPCR_indexing_KAPA", "2ndPCR_indexing_Platinum", "1stPCR_indexing_Platinum"))

### Top 1
tab_exp1_topprop[[1]]
tax_prop1 <- data.frame(x = rep(2,4), y = rep(1,4),
                         lab = c("italic(Acanthopagrus)", "italic(Hexagrammos)",
                                 "Cyprinidae", "STD_MiFish04"),
                         site = factor(site_level))
tax_prop1_ns <- data.frame(x = rep(2,4), y = rep(0.92,4),
                        lab = c("N.S."), site = factor(site_level))
fig_exp1_topprop1 <- fig_exp1_topprop[[1]] +
  scale_x_discrete(labels = discrete_label) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  geom_text(data = tax_prop1, aes(x = x, y = y, label = lab), parse = T, size = 4) +
  geom_text(data = tax_prop1_ns, aes(x = x, y = y, label = lab), parse = T, size = 4) +
  ylim(0, 1) +
  NULL
### Top 2
tab_exp1_topprop[[2]]
tax_prop2 <- data.frame(x = rep(2,4), y = rep(0.5,4),
                         lab = c("italic(Parablennius)", "italic(Saurida)",
                                 "Cyprinidae", "STD_MiFish09"),
                         site = factor(site_level))
tax_prop2_ns <- data.frame(x = rep(2,3), y = rep(0.46,3),
                           lab = c("N.S."), site = factor(site_level[2:4]))
text_prop2 <- data.frame(x = c(2, 2.5), y = c(0.37, 0.30),
                         lab = c("italic(P) < 0.01", "italic(P) == 0.068"),
                         site = factor("Sea_Nagahama", levels = levels(fig_exp1_topprop[[2]]$data$site)))
line_prop <- data.frame(x = c(1, 2), xend = c(3, 3), y = c(0.35, 0.28), yend = c(0.35, 0.28),
                         site = factor("Sea_Nagahama", levels = levels(fig_exp1_topprop[[2]]$data$site)))
fig_exp1_topprop2 <- fig_exp1_topprop[[2]] +
  scale_x_discrete(labels = discrete_label) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  geom_text(data = tax_prop2, aes(x = x, y = y, label = lab), parse = T, size = 4) +
  geom_text(data = tax_prop2_ns, aes(x = x, y = y, label = lab), parse = T, size = 4) +
  geom_text(data = text_prop2, aes(x = x, y = y, label = lab), parse = T, size = 4) +
  geom_segment(data = line_prop, aes(x = x, xend = xend, y = y, yend = yend)) +
  ylim(0, 0.5) +
  NULL
### Top 3
tab_exp1_topprop[[3]]
tax_prop3 <- data.frame(x = rep(2,4), y = rep(0.35,4),
                        lab = c("italic(Takifugu)", "italic(Entomacrodus)",
                                "italic(Lepomis)", "STD_lateolabrax"),
                        site = factor(site_level))
tax_prop3_ns <- data.frame(x = rep(2,4), y = rep(0.32,4),
                           lab = c("N.S."), site = factor(site_level))
fig_exp1_topprop3 <- fig_exp1_topprop[[3]] +
  scale_x_discrete(labels = discrete_label) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  geom_text(data = tax_prop3, aes(x = x, y = y, label = lab), parse = T, size = 4) +
  geom_text(data = tax_prop3_ns, aes(x = x, y = y, label = lab), parse = T, size = 4) +
  ylim(0,0.35) +
  NULL

## Relative abundance of top-rank taxa
fig_exp1_top_grid <- plot_grid(fig_exp1_topprop1 +
                                 theme(axis.text.x = element_blank()) +
                                 ggtitle("The most abundant taxon"),
                               fig_exp1_topprop2 +
                                 theme(axis.text.x = element_blank()) +
                                 ggtitle("Second most-abundant taxon"),
                               fig_exp1_topprop3 +
                                 ggtitle("Third most-abundant taxon"),
                               ncol = 1, rel_heights = c(1,1,1.1), labels = letters[1:3])


# --------------------------------------------------------- #
#  Save figures
# --------------------------------------------------------- #
# Main figures
## Figure 1
ggsave(file = sprintf("%s/Figure_01.pdf", fig_dir_out),
       plot = fig_exp_design, width = 12, height = 14)
## Figure 2
ggsave(file = sprintf("%s/Figure_02.pdf", fig_dir_out),
       plot = fig_exp1_main, width = 12, height = 14)


# Supplementary figures
## Figure S1
ggsave(file = sprintf("%s/Figure_S01.pdf", fig_dir_out),
       plot = fig_exp1_summary, width = 12, height = 12)
## Figure S2
ggsave(file = sprintf("%s/Figure_S02.pdf", fig_dir_out),
       plot = fig_exp1_bar_div, width = 12, height = 16)
## Figure S3
ggsave(file = sprintf("%s/Figure_S03.pdf", fig_dir_out),
       plot = fig_exp1_top_grid, width = 10, height = 14)


