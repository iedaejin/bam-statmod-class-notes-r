# Shared setup used by the class notes (source this first).
# Shared setup for every chapter: source("_common.R") in the first chunk.
knitr::opts_chunk$set(
  echo = TRUE, message = FALSE, warning = FALSE, comment = "#>",
  collapse = FALSE, fig.align = "center", fig.width = 6.5, fig.height = 3.5,
  out.width = "80%", fig.pos = "htbp", cache = FALSE
)
options(digits = 4, width = 72, pillar.sigfig = 4)
suppressPackageStartupMessages(library(ggplot2))
ie_navy <- "#002C8B"; ie_coral <- "#D2463C"; ie_blue <- "#0B6FD0"; ie_green <- "#2E8B57"
ie_pal <- c(ie_navy, ie_coral, ie_blue, ie_green, "#9A6B00", "#6A3D9A")
theme_set(theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", colour = ie_navy),
        axis.title = element_text(colour = ie_navy),
        panel.grid.minor = element_blank(), legend.position = "bottom"))
set.seed(2027)
