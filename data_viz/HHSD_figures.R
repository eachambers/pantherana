library(tidyverse)
library(here)
library(cowplot)
library(maps)
library(mxmaps) # devtools::install_github("diegovalle/mxmaps")
library(usmap)
theme_set(theme_cowplot())

source(here("analysis", "4_hhsd", "hhsd_functions.R"))
source(here("data_viz", "rana_colors.R"))

## The following file makes plots with HHSD results for main set of migpriors (Fig. 4) and
## also all tested migration priors (no figure).


# Only one migprior -------------------------------------------------------

# Colors for plotting
rana_cols <- retrieve_kcols(analysis = "hhsd")

# HHSD provides us with a decision csv with the support for gdi values. I manually
# renamed the decision files such that they're appended with the control file prefixes
# to specify the run.
# files <- list.files(path = here("data", "hhsd"), pattern = "decision.csv")
files <- list.files(path = here("data", "hhsd", "full_runs"), pattern = "decision.csv")

hhsd <- 
  1:length(files) %>% 
  lapply(function(x) {
    all <- read_csv(here("data", "hhsd", "full_runs", files[x]), col_names = TRUE) %>% 
      mutate(filename = paste0(files[x])) %>% 
      # Concatenate nodes to make a comparison col
      unite(c("node 1", "node 2"), col = "comparison", sep = "_", remove = FALSE) %>% 
      # Clean up filename into param settings
      separate(filename, into = c("assembly", "algorithm", "tmp"), sep = "_", remove = TRUE) %>% 
      dplyr::select(-c(`...1`, "tmp"))
    
    if (is.null(all$iteration)) all <- all %>% mutate(iteration = 1)
    
    df_1 <- all %>% 
      dplyr::select("comparison", "node 1", "gdi 1", "assembly", "algorithm", "iteration") %>% 
      rename(node = "node 1", gdi = "gdi 1") %>% 
      unite(c("node", "algorithm"), col = "node_alg", sep = "_", remove = FALSE)
    df_2 <- all %>% 
      dplyr::select("comparison", "node 2", "gdi 2", "assembly", "algorithm", "iteration") %>% 
      rename(node = "node 2", gdi = "gdi 2") %>% 
      unite(c("node", "algorithm"), col = "node_alg", sep = "_", remove = FALSE)
    
    bind_rows(df_1, df_2)
  }) %>% 
  bind_rows()

# TODO add remaining levels in below so they're ordered properly
# hhsd$node = factor(hhsd$node, levels = c("chic", "spec", "berl", "neov"))

# oneset <- hhsd %>% 
#   filter(migprior == "0110")

# Generate plots for each assembly
# p_mxpl <-
#   gdi_plot(dat = oneset %>% filter(assembly == "mxpl"), migpriors = "single") +
#   scale_fill_manual(values = c("chic" = rana_cols %>% filter(pop == "chic") %>% pull(color),
#                                "berl" = rana_cols %>% filter(pop == "berl") %>% pull(color),
#                                "spec" = rana_cols %>% filter(pop == "spec") %>% pull(color),
#                                "neov" = rana_cols %>% filter(pop == "neov") %>% pull(color)))

hhsd$iteration <- as.character(hhsd$iteration)
hhsd$iteration = factor(hhsd$iteration, levels = c("4", "3", "2", "1"))

p_mxpls <- hhsd %>% 
  filter(assembly == "mxpl") %>%
  filter(algorithm == "split") %>%
  ggplot(aes(x = gdi, y = node)) + 
  geom_rect(aes(xmin = 0, xmax = 1, ymin = 0, ymax = 2.5),
            fill = "lightgrey", alpha = 0.5) +
  geom_vline(xintercept = c(0.2, 0.7), linetype = "dashed", color = "grey50", linewidth = 1) +
  # geom_errorbar(aes(xmin = hpd_25, xmax = hpd_975), linewidth = 0.5,
  #               position = position_dodge(0.05)) +
  geom_line() +
  geom_point(fill = "grey", color = "black", pch = 21, size = 5) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), expand = c(0, 0)) +
  theme(axis.title.y = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(linewidth = 1.5),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 15, face = "italic")) +
  facet_grid(rows = vars(iteration), scales = "free_y")
  # facet_grid(rows = vars(iteration), cols = vars(algorithm),
  #            scales = "free_y", switch = "y")

p_mxplm <- hhsd %>% 
  filter(assembly == "mxpl") %>%
  filter(algorithm == "merge") %>% 
  ggplot(aes(x = gdi, y = node)) + 
  geom_rect(aes(xmin = 0, xmax = 1, ymin = 0, ymax = 2.5),
            fill = "lightgrey", alpha = 0.5) +
  geom_vline(xintercept = c(0.2, 0.7), linetype = "dashed", color = "grey50", linewidth = 1) +
  # geom_errorbar(aes(xmin = hpd_25, xmax = hpd_975), linewidth = 0.5,
  #               position = position_dodge(0.05)) +
  geom_line() +
  geom_point(fill = "grey", color = "black", pch = 21, size = 5) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), expand = c(0, 0)) +
  theme(axis.title.y = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(linewidth = 1.5),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 15, face = "italic")) +
  facet_grid(rows = vars(iteration), scales = "free_y")

plot_grid(p_mxpls, p_mxplm, nrow = 1)
ggsave(here("plots", "gdi_single_mxpl.pdf"), width = 6.5, height = 3.75)

########

p_forrs <- hhsd %>% 
  filter(assembly == "forr06K5") %>%
  filter(algorithm == "split") %>%
  ggplot(aes(x = gdi, y = node)) + 
  geom_rect(aes(xmin = 0, xmax = 1, ymin = 0, ymax = 2.5),
            fill = "lightgrey", alpha = 0.5) +
  geom_vline(xintercept = c(0.2, 0.7), linetype = "dashed", color = "grey50", linewidth = 1) +
  # geom_errorbar(aes(xmin = hpd_25, xmax = hpd_975), linewidth = 0.5,
  #               position = position_dodge(0.05)) +
  geom_line() +
  geom_point(fill = "grey", color = "black", pch = 21, size = 5) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), expand = c(0, 0)) +
  theme(axis.title.y = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(linewidth = 1.5),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 15, face = "italic")) +
  facet_grid(rows = vars(iteration), scales = "free_y")
# facet_grid(rows = vars(iteration), cols = vars(algorithm),
#            scales = "free_y", switch = "y")

p_forrm <- hhsd %>% 
  filter(assembly == "forr06K5") %>%
  filter(algorithm == "merge") %>% 
  ggplot(aes(x = gdi, y = node)) + 
  geom_rect(aes(xmin = 0, xmax = 1, ymin = 0, ymax = 2.5),
            fill = "lightgrey", alpha = 0.5) +
  geom_vline(xintercept = c(0.2, 0.7), linetype = "dashed", color = "grey50", linewidth = 1) +
  # geom_errorbar(aes(xmin = hpd_25, xmax = hpd_975), linewidth = 0.5,
  #               position = position_dodge(0.05)) +
  geom_line() +
  geom_point(fill = "grey", color = "black", pch = 21, size = 5) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), expand = c(0, 0)) +
  theme(axis.title.y = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(linewidth = 1.5),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 15, face = "italic")) +
  facet_grid(rows = vars(iteration), scales = "free_y")

plot_grid(p_forrs, p_forrm, nrow = 1)
ggsave(here("plots", "gdi_single_forr.pdf"), width = 6.5, height = 4.5)

########

p_foots <- hhsd %>% 
  filter(assembly == "foothills") %>%
  filter(algorithm == "split") %>%
  ggplot(aes(x = gdi, y = node)) + 
  geom_rect(aes(xmin = 0, xmax = 1, ymin = 0, ymax = 2.5),
            fill = "lightgrey", alpha = 0.5) +
  geom_vline(xintercept = c(0.2, 0.7), linetype = "dashed", color = "grey50", linewidth = 1) +
  # geom_errorbar(aes(xmin = hpd_25, xmax = hpd_975), linewidth = 0.5,
  #               position = position_dodge(0.05)) +
  geom_line() +
  geom_point(fill = "grey", color = "black", pch = 21, size = 5) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), expand = c(0, 0)) +
  theme(axis.title.y = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(linewidth = 1.5),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 15, face = "italic")) +
  facet_grid(rows = vars(iteration), scales = "free_y")
# facet_grid(rows = vars(iteration), cols = vars(algorithm),
#            scales = "free_y", switch = "y")

p_footm <- hhsd %>% 
  filter(assembly == "foothills") %>%
  filter(algorithm == "merge") %>% 
  ggplot(aes(x = gdi, y = node)) + 
  geom_rect(aes(xmin = 0, xmax = 1, ymin = 0, ymax = 2.5),
            fill = "lightgrey", alpha = 0.5) +
  geom_vline(xintercept = c(0.2, 0.7), linetype = "dashed", color = "grey50", linewidth = 1) +
  # geom_errorbar(aes(xmin = hpd_25, xmax = hpd_975), linewidth = 0.5,
  #               position = position_dodge(0.05)) +
  geom_line() +
  geom_point(fill = "grey", color = "black", pch = 21, size = 5) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), expand = c(0, 0)) +
  theme(axis.title.y = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(linewidth = 1.5),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 15, face = "italic")) +
  facet_grid(rows = vars(iteration), scales = "free_y")

plot_grid(p_foots, p_footm, nrow = 1)
ggsave(here("plots", "gdi_single_foot.pdf"), width = 6.5, height = 5.5)


# All migpriors -----------------------------------------------------------

# Generate plots for each
p_mxpl <- gdi_plot(dat = hhsd %>% filter(assembly == "mxpl"), migpriors = "multiple") +
  scale_fill_manual(values = c("chic" = rana_cols %>% filter(pop == "chic") %>% pull(color),
                               "berl" = rana_cols %>% filter(pop == "berl") %>% pull(color),
                               "spec" = rana_cols %>% filter(pop == "spec") %>% pull(color),
                               "neov" = rana_cols %>% filter(pop == "neov") %>% pull(color)))

# plot_grid(p_forr, p_foot, p_mxpl, nrow = 3) # export XXX
# ggsave(here("plots", "gdi_singlemigpriors.pdf"), width = X, height = X)


# GRAVEYARD ---------------------------------------------------------------

hhsd <- 
  1:length(files) %>% 
  lapply(function(x) {
    all <- read_csv(here("data", "hhsd", files[x]), col_names = TRUE) %>% 
      mutate(filename = paste0(files[x])) %>% 
      # Concatenate nodes to make a comparison col
      unite(c("node 1", "node 2"), col = "comparison", sep = "_", remove = FALSE) %>% 
      # Clean up filename into param settings
      separate(filename, into = c("assembly", "algorithm", "migprior", "tmp"), sep = "_", remove = TRUE)
    
    df_1 <- all %>% 
      dplyr::select("comparison", "node 1", "GDI 1", "2.5% HPD...3", "97.5% HPD...4", "assembly", "algorithm", "migprior") %>% 
      rename(node = "node 1", gdi = "GDI 1", hpd_25 = "2.5% HPD...3", hpd_975 = "97.5% HPD...4") %>% 
      unite(c("node", "algorithm"), col = "node_alg", sep = "_", remove = FALSE)
    df_2 <- all %>% 
      dplyr::select("comparison", "node 2", "GDI 2", "2.5% HPD...7", "97.5% HPD...8", "assembly", "algorithm", "migprior") %>% 
      rename(node = "node 2", gdi = "GDI 2", hpd_25 = "2.5% HPD...7", hpd_975 = "97.5% HPD...8") %>% 
      unite(c("node", "algorithm"), col = "node_alg", sep = "_", remove = FALSE)
    
    bind_rows(df_1, df_2)
  }) %>% 
  bind_rows()


mdvals <- c("#613a42", "#a13549", "#c1847e", "#ecbeb6", "#9ab8c9", "#7992a6", "#2073b7")

# theAuk = list(c("#273a66", "#2073b7", "#7992a6", "#9ab8c9", "#d4d8e2", "#ecbeb6", "#c1847e", "#a13549", "#613a42"), type = "diverging"),
# pal <- MVZ_palette("theAuk", 100, type = "continuous")

hcols <- read_tsv(here("data", "color_coding.txt"))

hcols %>% 
  ggplot(aes(x = group, y = gdi, color = gdi)) +
  geom_point(size = 5) +
  geom_label(aes(label = gdi), nudge_x = 0.25) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_gradientn(colors = rev(mdvals))
  # scale_color_gradientn(limits = c(0, 1), colors = pal)

