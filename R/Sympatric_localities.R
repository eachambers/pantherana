## Find sympatric localities

library(tidyverse)
library(here)
library(ggplot2)
library(cowplot)
library(maps)
library(mxmaps) # devtools::install_github("diegovalle/mxmaps")
library(ggrepel)
library(repel) # devtools::install_github("hms-dbmi/repel")
library(usmap)
library(sf)
library(sfheaders)
library(maptools)

theme_set(theme_cowplot())
setwd("~/Box Sync/Rana project/ddRADseq/ALL_RANA/Separate_assemblies/ADMIXTURE/")

source(here("R", "Pop_gen.R"))
source(here("R", "rana_colors.R"))


# Import relevant ADMIXTURE data ------------------------------------------

atlmx <- import_admix_data(path = ".", file = "ATL_MXPL/ATL_MXPL_relaxed_0.25miss_ldp", K_values = 7)
dat <- atlmx[[1]]
centam <- import_admix_data(path = ".", file = "CENTAM/new_CENTAM_relaxed_0.25miss_ldp", K_values = 4)
dat <- centam[[1]]
pacmx <- import_admix_data(path = ".", file = "PACMX/new_PACMX_relaxed_0.25miss_ldp", K_values = 6)
dat <- pacmx[[1]]

forreri <- import_admix_data(path = ".", file = "forreri/forreri_0.25miss_ldp", K_values = 5)


# Process data ------------------------------------------------------------

dataset_name = "PACMX" # "PACMX" / "CENTAM"

# if (dataset_name == "ATL_MXPL") dat <- dat %>% dplyr::filter(Bioinformatics_ID != "IRL57_LCA")
K_value = ncol(dat)-5
kcols <- retrieve_kcols(K_value, dataset = "admixture", dataset_name)
metadata <- read_tsv(paste0(here("data"), "/", dataset_name, "_metadata.txt"), col_names = TRUE)
sites <- read_tsv(paste0(here("data"), "/", dataset_name, "_sites.txt"), col_names = TRUE)

final <- left_join(dat, metadata, by = "Bioinformatics_ID")
final <- left_join(final, sites)

# Sympatric sites are those that have 2 max_Ks
check <-
  final %>%
  group_by(site_name, max_K) %>%
  count()

# Extract sympatric localities only; N.B.: RioelAmate and RioPapagayo are identical
if (dataset_name == "PACMX") symps <- final %>% dplyr::filter(site_name == "PuenteTehuantepec" | site_name == "RiolasTejas" | 
                  site_name == "RioelAmate" | site_name == "RioPapagayo") %>% dplyr::arrange(site_name, max_K) %>% rownames_to_column(var = "new_order")
if (dataset_name == "CENTAM") symps <- final %>% dplyr::filter(site_name == "CerroSPLaLoma") %>% dplyr::arrange(site_name, max_K) %>% rownames_to_column(var = "new_order")

symps$new_order <- factor(symps$new_order, levels = symps$new_order)

symps <-
  symps %>% 
  tidyr::pivot_longer(names_to = "cluster", values_to = "proportion", 
                      -c(Bioinformatics_ID, K_val, order, new_order, max_K, max_K_val, lat, long, site_name, dataset, state))
p <- symps %>% 
  ggplot2::ggplot(aes(x = new_order, y = proportion, fill = cluster)) +
  ggplot2::geom_bar(stat = "identity") + 
  panel_border() + 
  scale_y_continuous(expand = c(0,0)) +
  # scale_x_discrete(expand = c(0,0)) +
  scale_fill_manual(values = kcols) +
  theme_cowplot() %+replace% theme(axis.line = element_line(colour = "black"),
                                   axis.text.x = element_blank(),
                                   axis.text.y = element_blank(),
                                   axis.ticks = element_blank(),
                                   axis.title.x = element_blank(),
                                   axis.title.y = element_blank(),
                                   legend.position = "none",
                                   panel.border = element_rect(fill = NA, colour = "black", linetype = "solid", linewidth = 1.5),
                                   strip.text.y = element_text(size = 30, face = "bold"),
                                   strip.background = element_rect(colour = "white", fill = "white"),
                                   panel.spacing = unit(-0.1, "lines"))

ids <- as.list(unique(symps$Bioinformatics_ID))
p <- p + 
  scale_x_discrete(expand = c(0,0), labels = ids) +
  theme(axis.text.x = element_text(angle = 90, size = 5)) # 10x4 PACMX / 2.3x3 CENTAM


# Guerrero samples for PACMX ----------------------------------------------

data("mxstate.map")
mxstate.map$group <- as.numeric(mxstate.map$group)
gueoax <- mxstate.map %>% dplyr::filter(id == 12 | id == 20)

dataset_name = "PACMX"
K_value = ncol(dat)-5
kcols <- retrieve_kcols(K_value, dataset = "admixture", dataset_name)
metadata <- read_tsv(paste0(here("data"), "/", dataset_name, "_metadata.txt"), col_names = TRUE)
sites <- read_tsv(paste0(here("data"), "/", dataset_name, "_sites.txt"), col_names = TRUE)

final <- left_join(dat, metadata, by = "Bioinformatics_ID")
gro <- left_join(final, sites) %>% 
  dplyr::filter(state == "Guerrero" | state == "Oaxaca") %>% 
  dplyr::arrange(site_name, max_K) %>% 
  rownames_to_column(var = "new_order")
gro$new_order <- factor(gro$new_order, levels = gro$new_order)

# Get observation numbers (optional)
obs <-
  gro %>% 
  # filter(state == "Oaxaca") %>% 
  dplyr::select(site_name) %>% 
  group_by(site_name) %>% 
  summarize(no_inds = n())

gro <-
  gro %>% 
  tidyr::pivot_longer(names_to = "cluster", values_to = "proportion", 
                      -c(Bioinformatics_ID, K_val, order, new_order, max_K, max_K_val, lat, long, site_name, dataset, state))

# Plot structure plot, ordered by site
gro %>% 
  filter(state == "Oaxaca") %>%
  ggplot2::ggplot(aes(x = new_order, y = proportion, fill = cluster)) +
  ggplot2::geom_bar(stat = "identity") + 
  panel_border() + 
  scale_y_continuous(expand = c(0,0)) +
  # scale_x_discrete(expand = c(0,0)) +
  scale_fill_manual(values = kcols) +
  theme_cowplot() %+replace% theme(axis.line = element_line(colour = "black"),
                                   axis.text.x = element_blank(),
                                   axis.text.y = element_blank(),
                                   axis.ticks = element_blank(),
                                   axis.title.x = element_blank(),
                                   axis.title.y = element_blank(),
                                   legend.position = "none",
                                   panel.border = element_rect(fill = NA, colour = "black", linetype = "solid", linewidth = 1.5),
                                   strip.text.y = element_text(size = 30, face = "bold"),
                                   strip.background = element_rect(colour = "white", fill = "white"),
                                   panel.spacing = unit(-0.1, "lines")) # export 12x8

# Make the map with site localities
ggplot() +
  geom_polygon(data = gueoax, aes(x = long, y = lat, group = group), color = "white", fill = "#e6e6e6", linewidth = 0.25) + # color = "#969696" if CENTAM
  theme_map() +
  coord_fixed() +
  geom_point(data = gro, aes(long, lat), size = 3) +
  geom_text(data = gro, aes(x = long, y = lat, label = site_name), color = "black", size = 3) # export 10x8


             
             