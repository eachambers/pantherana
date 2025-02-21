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
## also all tested migration priors (Fig. S6).


# Only one migprior -------------------------------------------------------

# Colors for plotting
rana_cols <- retrieve_kcols(analysis = "hhsd")

# HHSD provides us with a decision csv with the support for gdi values. I manually
# renamed the decision files such that they're appended with the control file prefixes
# to specify the run.
files <- list.files(path = here("data", "hhsd"), pattern = "decision.csv")

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

# TODO add remaining levels in below so they're ordered properly
hhsd$node = factor(hhsd$node, levels = c("chic", "spec", "berl", "neov"))

oneset <- hhsd %>% 
  filter(migprior == "0110")

# TODO generate plots for all three assemblies
# Generate plots for each
p_mxpl <-
  gdi_plot(dat = oneset %>% filter(assembly == "mxpl"), migpriors = "single") +
  scale_fill_manual(values = c("chic" = rana_cols %>% filter(pop == "chic") %>% pull(color),
                               "berl" = rana_cols %>% filter(pop == "berl") %>% pull(color),
                               "spec" = rana_cols %>% filter(pop == "spec") %>% pull(color),
                               "neov" = rana_cols %>% filter(pop == "neov") %>% pull(color)))

# plot_grid(p_forr, p_foot, p_mxpl, nrow = 3) # export XXX
# ggsave(here("plots", "gdi_singlemigpriors.pdf"), width = X, height = X)


# All migpriors -----------------------------------------------------------

# TODO generate plots for all three assemblies
# Generate plots for each
p_mxpl <- gdi_plot(dat = hhsd %>% filter(assembly == "mxpl"), migpriors = "multiple") +
  scale_fill_manual(values = c("chic" = rana_cols %>% filter(pop == "chic") %>% pull(color),
                               "berl" = rana_cols %>% filter(pop == "berl") %>% pull(color),
                               "spec" = rana_cols %>% filter(pop == "spec") %>% pull(color),
                               "neov" = rana_cols %>% filter(pop == "neov") %>% pull(color)))

# plot_grid(p_forr, p_foot, p_mxpl, nrow = 3) # export XXX
# ggsave(here("plots", "gdi_singlemigpriors.pdf"), width = X, height = X)
