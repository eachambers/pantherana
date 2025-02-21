library(algatr) # devtools::install_github("TheWangLab/algatr")
library(terra)
library(raster)
library(tidyverse)
library(here)
library(RStoolbox) # devtools::install_github("bleutner/RStoolbox")
library(viridis)
library(cowplot)
library(ggrepel)
library(MVZlibrary)
library(mxmaps)
theme_set(theme_cowplot())

## This code processes environmental data used for landscape genomic analyses for the R. forreri complex:
##     (1) Retrieve metadata
##     (2) Retrieve environmental data
##     (3) Check for collinearity among enviro layers
##     (4) Run raster PCA
##     (5) Export data
##     (6) Fig. S7A: Map out raster PCA results
##     (7) Fig. S7B: See how bioclimatic variables load onto raster PCs

##    FILES GENERATED:
##          forreri_envlayers.tif
##          forreri_envpcs_results.txt


# (1) Retrieve metadata ---------------------------------------------------

samps <- read_tsv(here("data", "forreri_metadata.txt"),
                  col_names = TRUE) %>% 
  dplyr::rename(x = long,
                y = lat) %>% 
  na.omit() # V2797 does not have coords

# names <- samps$Bioinformatics_ID


# (2) Retrieve enviro data ------------------------------------------------

# Retrieve worldclim tiles for coords
env <- algatr::get_worldclim(coords = samps %>% dplyr::select(x, y),
                             res = 0.5,
                             buff = 0.01,
                             save_output = TRUE)

# Take a look at WorldClim tiles that were retrieved
plot(env[[1]], col = magma(100), axes = FALSE)
points(samps %>% dplyr::select(x, y), pch = 19)


# (3) Check for collinearity ----------------------------------------------

cors_env <- algatr::check_env(env) # 31 pairs of vars had correlation coefficients > 0.7.
dist <- as.data.frame(cors_env$cor_matrix)
write_csv(dist, file = here("data", "forreri_cors_env.csv"), col_names = TRUE)

# Take a look at the collinearity among layers
dist %>%
  tibble::rownames_to_column("sample") %>%
  tidyr::gather("sample_comp", "dist", -"sample") %>%
  ggplot2::ggplot(ggplot2::aes(x = sample, y = sample_comp, fill = dist)) +
  ggplot2::geom_tile() +
  ggplot2::coord_equal() +
  viridis::scale_fill_viridis(option = "magma") +
  ggplot2::xlab("Sample") +
  ggplot2::ylab("Sample") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))


# (4) Raster PCA ----------------------------------------------------------

env_pcs <- RStoolbox::rasterPCA(env, spca = TRUE)

# Let's take a look at the results for the top three PCs
plots <- lapply(1:3, function(x) RStoolbox::ggR(env_pcs$map, x, geom_raster = TRUE))
plots[[1]] 
plots[[2]] 
plots[[3]]

# Plot the PCA itself
pcs <- terra::as.data.frame(env_pcs$map, xy = TRUE)

# How much variance is explained by each env PC?
summary(env_pcs$model) # <- made this into a txt file externally called "forreri_envpcs_results.txt"

dat <- readr::read_tsv(here("data", "forreri_envpcs_results.txt"), col_names = TRUE) %>% 
  tidyr::pivot_longer(cols = Comp.1:Comp.19,
                      names_to = "comp",
                      values_to = "value")

dat$comp <- as.character(dat$comp)
dat$comp <- factor(dat$comp, levels = c("Comp.1", "Comp.2", "Comp.3",
                                        "Comp.4", "Comp.5", "Comp.6",
                                        "Comp.7", "Comp.8", "Comp.9",
                                        "Comp.10", "Comp.11", "Comp.12",
                                        "Comp.13", "Comp.14", "Comp.15",
                                        "Comp.16", "Comp.17", "Comp.18",
                                        "Comp.19"))
# Visualize results
dat %>% 
  dplyr::filter(Statistic == "Prop_variance") %>% 
  ggplot(aes(x = comp, y = value)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90)) +
  xlab("Env PC comparison") +
  ylab("Proportion of variance")

# Get statistics on proportion of variance explained
dat %>% 
  filter(Statistic=="Prop_variance") %>% 
  filter(comp=="Comp.1" | comp=="Comp.2" | comp=="Comp.3") %>% 
  summarize(top3 = sum(value)) # 85% explained by top 3 PCs


# (5) Save data -----------------------------------------------------------

# Save env layers
terra::writeRaster(env, file = here("data", "forreri_envlayers.tif"),
                   overwrite = TRUE)

# Save top 3 env PCs
lapply(1:3, function(x) terra::writeRaster(env_pcs$map[[x]],
                                           file = paste0(here("data", "PC_layers"), "/forreri_PC", x, ".tif"),
                                           overwrite = TRUE))

# Stack the top 3 env PCs and save as a single raster
forr_env <- raster::stack(list.files(here("data", "PC_layers"), full.names = TRUE))
raster::writeRaster(forr_env, here("data", "forreri_PCenv.tif"), overwrite = TRUE)


# (6) Plot raster PCA results ---------------------------------------------

# If you haven't run the above, can just start with:
# forr_env <- raster::stack(list.files(here("data", "PC_layers"), full.names = TRUE))

# Take a look at PCs 1 and 2 plotted
# pcs <- as.data.frame(forr_env, xy = TRUE)
# p_hex12 <-
#   pcs %>% 
#   ggplot(aes(x = PC1, y = PC2)) +
#   geom_hex() +
#   coord_equal() +
#   scale_fill_viridis(option = "A") # export 8x6

# Read in borders from shapefile
border <- sf::st_read(here("data"), layer = "mx_centam_merged")
sf::st_crs(border) = 4326 # WGS 84

pal <- MVZ_palette("LifeHistories",  type = "discrete")
palette = grDevices::colorRampPalette(pal)(100)

forr_env_r <- terra::rast(forr_env)
xlim <- c(xmin(forr_env_r), xmax(forr_env_r))
ylim <- c(ymin(forr_env_r), ymax(forr_env_r))

data("mxstate.map")

ggplot() +
  geom_spatraster(data = forr_env_r[[3]], aes(fill = PC3)) +
  scale_fill_gradientn(colors = palette, na.value = NA, name = "PC3") +
  theme_map() +
  geom_sf(data = border, fill = NA, color = "white", size = 0.5) +
  coord_sf(xlim = xlim, ylim = ylim) +
  geom_polygon(data = mxstate.map, aes(x = long, y = lat, group = group), color = "white", fill = NA, linewidth = 0.25) +
  geom_point(data = samps, aes(x = x, y = y))
ggsave(here("plots", "forr_PC3_map.pdf"), width = 10, height = 8, units = "in")


# (7) Bioclimatic vars in each PC -----------------------------------------

l <- loadings(env_pcs$model)
loaddf <- data.frame(matrix(as.numeric(l), attributes(l)$dim, dimnames = attributes(l)$dimnames))
bioclim <- read_tsv(here("data", "bioclim_vars.txt"), col_names = c("BIO", "Description"))

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
ggsave(here("plots", "forr_PC3_loadings.pdf"), width = 2.5, units = "in")
