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


# <---------------------------------------------> #
#  Load figure objects
# <---------------------------------------------> #
# Experiment III
fig_exp3_divraw <- readRDS(paste0(fig_dir, "10_2_Fig_Exp3_RepVolDiveristy.obj"))
fig_exp3_divbox <- readRDS(paste0(fig_dir, "10_2_Fig_Exp3_DivBoxplot.obj"))
# Additional analyses
fig_exp3_bray <- readRDS(paste0(fig_dir, "10_3_Fig_Exp3_RepVolBrayCurtis.obj"))
# Additional analyses
fig_13_rareotu <- readRDS(paste0(fig_dir, "13_1_Fig_RareOTUs.obj"))


# <---------------------------------------------> #
# Experiment 3
# <---------------------------------------------> #
## Sea_Nagahama
test_cat1 <- c(rep("1-rep.",4), rep("1-µl",5), rep("2-rep.",5),rep("2-µl",5),
               rep("4-rep.",5), rep("4-µl",5), rep("8-rep.",5),rep("8-µl",5))
fig_exp3_divbox[[2]]$data$test_name[fig_exp3_divbox[[2]]$data$test_name == "Replicate_test"] <- "Replicate test"
fig_exp3_divbox[[2]]$data$test_name[fig_exp3_divbox[[2]]$data$test_name == "Volume_test"] <- "Volume test"
fig_exp3_divbox[[2]]$data$reaction_scale_fac <- test_cat1
fig_exp3_divbox[[2]]$layers[[1]] <- fig_exp3_divbox[[2]]$layers[[2]] <- NULL
### Format figure
g1 <- fig_exp3_divbox[[2]] +
  geom_boxplot(width = 0.5, alpha = 0.5, outlier.shape = NA, fill = NA) +
  geom_jitter(width = 0.1, height = 0, alpha = 0.65) +
  facet_wrap(. ~ test_name, scales = "free_x") +
  labs(title = "Sea_Nagahama",
       subtitle = bquote("Reaction scale:"~italic(P)~"< 0.0001,"~"Test category:"~"N.S.")) +
  ylab("No. of fish OTUs") + xlab("Treatment") +
  theme(legend.position = "none") + ylim(0,35)
### Statistics
summary(glm(value ~ test_name * reaction_scale, data = g1$data, family = poisson(link = "log")))
summary(glm(value ~ test_category, data = g1$data %>% filter(reaction_scale == 1), family = poisson(link = "log")))
summary(glm(value ~ test_category, data = g1$data %>% filter(reaction_scale == 2), family = poisson(link = "log")))
summary(glm(value ~ test_category, data = g1$data %>% filter(reaction_scale == 4), family = poisson(link = "log")))
summary(glm(value ~ test_category, data = g1$data %>% filter(reaction_scale == 8), family = poisson(link = "log")))

## STD_Mix
test_cat2 <- c(rep("1-rep.",5), rep("1-µl",5), rep("2-rep.",5),rep("2-µl",5),
               rep("4-rep.",5), rep("4-µl",5), rep("8-rep.",5),rep("8-µl",5))
fig_exp3_divbox[[4]]$data$test_name[fig_exp3_divbox[[4]]$data$test_name == "Replicate_test"] <- "Replicate test"
fig_exp3_divbox[[4]]$data$test_name[fig_exp3_divbox[[4]]$data$test_name == "Volume_test"] <- "Volume test"
fig_exp3_divbox[[4]]$data$reaction_scale_fac <- test_cat2
fig_exp3_divbox[[4]]$layers[[1]] <- fig_exp3_divbox[[4]]$layers[[2]] <- NULL
### Format figure
g2 <- fig_exp3_divbox[[4]] +
  geom_boxplot(width = 0.5, alpha = 0.5, outlier.shape = NA, fill = NA) +
  geom_jitter(width = 0.1, height = 0, alpha = 0.65) +
  facet_wrap(. ~ test_name, scales = "free_x") +
  labs(title = "STD_Mix",
       subtitle = bquote("Reaction scale:"~italic(P)~"= 0.103,"~"Test category:"~"N.S.")) +
  ylab("No. of standard DNA OTUs") + xlab("Treatment") +
  theme(legend.position = "none") + scale_y_continuous(limits = c(0,11), breaks = seq(0,12,2))
### Statistics
summary(glm(value ~ test_name * reaction_scale, data = g2$data, family = poisson(link = "log")))
summary(glm(value ~ test_category, data = g2$data %>% filter(reaction_scale == 1), family = poisson(link = "log")))
summary(glm(value ~ test_category, data = g2$data %>% filter(reaction_scale == 2), family = poisson(link = "log")))
summary(glm(value ~ test_category, data = g2$data %>% filter(reaction_scale == 4), family = poisson(link = "log")))
summary(glm(value ~ test_category, data = g2$data %>% filter(reaction_scale == 8), family = poisson(link = "log")))


### Bray-Curtis dissimilarity
test_cat3 <- c(rep("1-rep.",6), rep("2-rep.",10), rep("4-rep.",10), rep("8-rep.",10),
               rep("1-µl", 10), rep("2-µl.", 10), rep("4-µl.", 10), rep("8-µl.",10))
test_cat4 <- c(rep("1-rep.",10), rep("2-rep.",10), rep("4-rep.",10), rep("8-rep.",10),
               rep("1-µl", 10), rep("2-µl.", 10), rep("4-µl.", 10), rep("8-µl.",10))
## Sea_Nagahama
fig_exp3_bray$data$test_name[fig_exp3_bray$data$test_name == "Replicate_test"] <- "Replicate test"
fig_exp3_bray$data$test_name[fig_exp3_bray$data$test_name == "Volume_test"] <- "Volume test"
fig_exp3_bray$data$reaction_scale <- c(test_cat3, test_cat4)
### Prepare text
bray1_p_val1 <- data.frame(x = rep(3.2,2), y = rep(0.95,2),
                           lab = c("Reaction scale:",
                                   "Reaction scale:"),
                           test_name = c("Replicate test", "Volume test"))
bray1_p_val2 <- data.frame(x = rep(3.2,2), y = rep(0.85,2),
                          lab = c("italic(P) == 0.039",
                                  "italic(P) == 0.043"),
                          test_name = c("Replicate test", "Volume test"))
### Format figure
g3 <- fig_exp3_bray$data %>% filter(site == "Sea_Nagahama") %>%
  ggplot(aes(x = reaction_scale, y = bray_curtis)) +
  geom_boxplot(width = 0.5, alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.1, height = 0, alpha = 0.65) +
  facet_wrap(. ~ test_name, scales = "free_x") +
  geom_text(data = bray1_p_val1, aes(x = x, y = y, label = lab), size = 3.5) +
  geom_text(data = bray1_p_val2, aes(x = x, y = y, label = lab), parse = T, size = 3.5) +
  panel_border() + ggtitle("Sea_Nagahama") +
  ylab("Bray-Curtis dissimilarity") + xlab("Treatment") +
  theme(legend.position = "none") + ylim(0,1)

## Std_Mix
### Prepare text
bray2_p_val1 <- data.frame(x = rep(3.2,2), y = rep(0.95,2),
                           lab = c("Reaction scale:",
                                   "Reaction scale:"),
                           test_name = c("Replicate test", "Volume test"))
bray2_p_val2 <- data.frame(x = rep(3.2,2), y = rep(0.85,2),
                           lab = c("italic(P) == 0.003",
                                   "italic(P) == 0.006"),
                           test_name = c("Replicate test", "Volume test"))
### Format figure
g4 <- fig_exp3_bray$data %>% filter(site == "STD_Mix") %>%
  ggplot(aes(x = reaction_scale, y = bray_curtis)) +
  geom_boxplot(width = 0.5, alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.1, height = 0, alpha = 0.65) +
  facet_wrap(. ~ test_name, scales = "free_x") +
  geom_text(data = bray2_p_val1, aes(x = x, y = y, label = lab), size = 3.5) +
  geom_text(data = bray2_p_val2, aes(x = x, y = y, label = lab), parse = T, size = 3.5) +
  panel_border() + ggtitle("STD_Mix") +
  ylab("Bray-Curtis dissimilarity") + xlab("Treatment") +
  theme(legend.position = "none") + ylim(0,1)


### Rare OTUs
g5 <- fig_13_rareotu[[1]] + scale_x_discrete(labels = c("Detected\nOTUs\nin 1-rep.",
                                                        "Newly\ndetected\nOTUs\nin 2-rep.",
                                                        "Newly\ndetected\nOTUs\nin 4-rep.",
                                                        "Newly\ndetected\nOTUs\nin 8-rep.")) +
  scale_color_startrek(name = "Treatment", labels = c("1-rep.", "2-rep.", "4-rep.", "8-rep.")) + 
  ggtitle("Sea_Nagahama, Replication test") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1))
### Volume test
g6 <- fig_13_rareotu[[3]] + scale_x_discrete(labels = c("Detected\nOTUs\nin 1-µl",
                                                        "Newly\ndetected\nOTUs\nin 2-µl",
                                                        "Newly\ndetected\nOTUs\nin 4-µl",
                                                        "Newly\ndetected\nOTUs\nin 8-µl")) +
  scale_color_startrek(name = "Treatment", labels = c("1-µl", "2-µl", "4-µl", "8-µl")) + 
  ggtitle("Sea_Nagahama, Volume test") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1))


# <---------------------------------------------> #
#  Combine main figures
# <---------------------------------------------> #
fig_exp3_divbray <- plot_grid(g1, g3, g2, g4,
                              ncol = 2, align = "hv", axis = "lrtb",
                              byrow = F, rel_widths = c(1, 1),
                              labels = c("a","b","c","d"))
fig_rareotu <- plot_grid(g5,# + theme(legend.position = c(0.7, 0.8)),
                         g6,# + theme(legend.position = c(0.7, 0.8)),
                         ncol = 2, labels = c("e","f"),
                         align = "hv")
fig_exp3_main <- plot_grid(fig_exp3_divbray, fig_rareotu,
                           ncol = 1, rel_heights = c(1,0.5))


# <---------------------------------------------> #
#  Supplementary figure
# <---------------------------------------------> #
## Relabel
df1 <- fig_exp3_divraw[[3]]$data
df1$test_name[df1$test_name == "Replicate_test"] <- "Replicate test"
df1$test_name[df1$test_name == "Volume_test"] <- "Volume test"
## Replicate test names
df1$test_category[df1$test_category == "1rep"] <- "1-rep."
df1$test_category[df1$test_category == "2rep"] <- "2-rep."
df1$test_category[df1$test_category == "4rep"] <- "4-rep."
df1$test_category[df1$test_category == "8rep"] <- "8-rep."
## Volume test names
df1$test_category[df1$test_category == "1ul"] <- "1-µl"
df1$test_category[df1$test_category == "2ul"] <- "2-µl"
df1$test_category[df1$test_category == "4ul"] <- "4-µl"
df1$test_category[df1$test_category == "8ul"] <- "8-µl"

fig_exp3_divraw[[3]]$data <- df1
fig_exp3_divraw[[3]] <- fig_exp3_divraw[[3]] + scale_fill_igv(name = "Family\nor STD name")

# <---------------------------------------------> #
# Save figures
# <---------------------------------------------> #
# Main figures
## Figure 4
ggsave(file = sprintf("%s/Figure_04.pdf", fig_dir_out),
       plot = fig_exp3_main, width = 12, height = 12)

# Supplementary figures
## Figure S4
ggsave(file = sprintf("%s/Figure_S04.pdf", fig_dir_out),
       plot = fig_exp3_divraw[[3]], width = 12, height = 12)

