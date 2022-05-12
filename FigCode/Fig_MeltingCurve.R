####
#### Melting curve figure for eDNA-seq
#### 2022.05.12 revision for Environmental DNA
#### R 4.1.2
####

# Set working directory
if(basename(getwd()) != "FigCode") setwd("FigCode")

# Load libraries
## For general
library(tidyverse); packageVersion("tidyverse") # 1.3.1, 2021.8.25
library(cowplot); packageVersion("cowplot") # 1.1.1, 2021.3.1
library(ggsci); packageVersion("ggsci") # 2.9, 2021.8.26
theme_set(theme_cowplot())

# Prepare output folders
fig_dir <- "00_RawFigs/"
fig_dir_out <- "00_FormattedFigs/"


# --------------------------------------------------------- #
#  Load data
# --------------------------------------------------------- #
d <- read_csv("data_fig/data_melting_curve.csv",
              col_types = cols(temp = col_double(),
                               value = col_double(),
                               site = col_factor(),
                               replicate = col_factor()))

# --------------------------------------------------------- #
#  Visualize
# --------------------------------------------------------- #
g1 <- d %>% 
  ggplot(aes(x = temp, y = value, color = site, group = site:replicate)) +
  geom_line(alpha = 0.75) +
  scale_color_startrek(name = "Site") +
  xlab(expression(paste("Temperature (", degree, "C)"))) +
  ylab("– (d/dT) Fluorescence (465 nm – 510 nm)") +
  theme_gray() +
  NULL

# --------------------------------------------------------- #
#  Save figures
# --------------------------------------------------------- #
cairo_pdf(file = sprintf("%s/Figure_S_MeltingCurve.pdf", fig_dir_out),
          width = 8, height = 4)
g1; dev.off()
