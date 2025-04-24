library(tidyverse)
library(here)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

## The following file builds Fig. S2: the taxonomic history of the leopard frog 
## species complex.

##    FILES REQUIRED:
##            rana_taxonomy.txt # Taxonomic history of the group


tax <- read_tsv(here("data", "4_Data_visualization", "data_files_input_into_scripts", "rana_taxonomy.txt"), col_names = TRUE)

plot <-
  tax %>% 
  # filter(Year < 1988) %>% 
  ggplot() +
  annotate("rect", xmin = 1780, xmax = 1940, ymin = 0, ymax = 33,
           alpha = 0.5, fill = "#abb285") +
    annotate("rect", xmin = 1940, xmax = 1960, ymin = 0, ymax = 33,
             alpha = 0.5, fill = "#ccb3d7") +
    annotate("rect", xmin = 1960, xmax = 2000, ymin = 0, ymax = 33,
             alpha = 0.5, fill = "#f39e87") +
    annotate("rect", xmin = 2000, xmax = 2026, ymin = 0, ymax = 33,
             alpha = 0.5, fill = "#88a1c6") +
  geom_line(aes(x = Year, y = To_plot), color = "black") +
  # geom_line(aes(x = Year, y = Cumulative_recognized_today)) +
  ylab("Number of species") +
  scale_x_continuous(limits = c(1780, 2026), breaks = seq(1780, 2026, by = 20), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.title = element_text(face = "bold"))

# Add vertical lines
plot <-
  plot +
  geom_vline(xintercept = 1782, linetype = "dashed", color = "#787878") +
  geom_vline(xintercept = 1920, linetype = "dashed", color = "#787878") +
  geom_vline(xintercept = 1942, linetype = "dashed", color = "#787878") +
  geom_vline(xintercept = 2005, linetype = "dashed", color = "#787878") +
  geom_vline(xintercept = 2024, linetype = "dashed", color = "#787878")

# export plot at 10x6"