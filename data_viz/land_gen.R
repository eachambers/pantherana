library(algatr) # devtools::install_github("TheWangLab/algatr")
library(terra)
library(raster)
library(tidyverse)
library(tidyterra)
library(here)
library(RStoolbox) # devtools::install_github("bleutner/RStoolbox")
library(viridis)
library(cowplot)
library(ggrepel)
library(MVZlibrary)
library(mxmaps)
theme_set(theme_cowplot())

## This code visualizes environmental data used for landscape genomic analyses for the R. forreri complex:
##     (1) Fig. S8A: Map out raster PCA results
##     (2) Fig. S8B: See how bioclimatic variables load onto raster PCs

##    FILES REQUIRED:
##          XXX
##          XXX


# (1) Import data ---------------------------------------------------------

forr_env <- terra::rast(here("data", "3_Analyses", "4_landgen", "forreri_PCenv.tif"))
xlim <- c(xmin(forr_env), xmax(forr_env))
ylim <- c(ymin(forr_env), ymax(forr_env))

# Read in borders from shapefile
border <- sf::st_read(here("data", "2_Data_processing", "data_files_input_into_scripts"), layer = "mx_centam_merged")
sf::st_crs(border) = 4326 # WGS 84

pal <- MVZ_palette("LifeHistories",  type = "discrete")
palette = grDevices::colorRampPalette(pal)(100)

data("mxstate.map")

coords <- read_tsv(here("data", "2_Data_processing", "data_files_input_into_scripts", "forreri_metadata.txt"))


# (2) Figure S8A: raster PCA results --------------------------------------

ggplot() +
  geom_spatraster(data = forr_env[[3]]) +
  scale_fill_gradientn(colors = palette, na.value = NA, name = "PC3") +
  theme_map() +
  geom_sf(data = border, fill = NA, color = "white", size = 0.5) +
  coord_sf(xlim = xlim, ylim = ylim) +
  geom_polygon(data = mxstate.map, aes(x = long, y = lat, group = group), color = "white", fill = NA, linewidth = 0.25) +
  geom_point(data = coords, aes(x = long, y = lat))

ggsave(here("plots", "FigS8A_forrPC3_map.pdf"), width = 10, height = 8, units = "in")


# (3) Figure S8B: bioclimatic vars in each PC -----------------------------

# Retrieve `env_pcs` object from the env_data.R script
l <- loadings(env_pcs$model)
loaddf <- data.frame(matrix(as.numeric(l), attributes(l)$dim, dimnames = attributes(l)$dimnames))
bioclim <- read_tsv(here("data", "2_Data_processing", "data_files_input_into_scripts", "bioclim_vars.txt"), col_names = c("BIO", "Description"))

data <-
  loaddf %>% 
  rownames_to_column(var = "BIO") %>% 
  left_join(., bioclim)

ggplot() +
  geom_hline(yintercept = 0, linewidth = 0.5, color = "gray") +
  geom_vline(xintercept = 0, linewidth = 0.5, color = "gray") +
  geom_segment(data = loaddf, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.3), arrow = arrow(length = unit(0.2, "cm"))) +
  geom_text(data = data, aes(x = Comp.1, y = Comp.3, label = BIO), size = 4) +
  # geom_text_repel(data = data, aes(x = Comp.1, y = Comp.3, label = Description), size = 4) +
  coord_equal() +
  xlab("PC1") +
  ylab("PC3") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 10))
ggsave(here("plots", "FigS8B_forrPC3_loadings.pdf"), width = 2.5, units = "in")
