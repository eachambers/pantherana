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

##    FILES GENERATED:
##          forreri_envlayers.tif
##          forreri_envpcs_results.txt


# (1) Retrieve metadata ---------------------------------------------------

samps <- read_tsv(here("data", "2_Data_processing", "data_files_input_into_scripts", "forreri_metadata.txt"),
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
write_csv(dist, file = here("data", "3_Analyses", "4_landgen", "forreri_cors_env.csv"), col_names = TRUE)

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

dat <- readr::read_tsv(here("data", "3_Analyses", "4_landgen", "forreri_envpcs_results.txt"), col_names = TRUE) %>% 
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
terra::writeRaster(env, file = here("data", "3_Analyses", "4_landgen", "forreri_envlayers.tif"),
                   overwrite = TRUE)

# Save top 3 env PCs
# lapply(1:3, function(x) terra::writeRaster(env_pcs$map[[x]],
#                                            file = paste0(here("data", "PC_layers"), "/forreri_PC", x, ".tif"),
#                                            overwrite = TRUE))

# Save env PCs as a single raster
raster::writeRaster(env_pcs$map, here("data", "3_Analyses", "4_landgen", "forreri_PCenv.tif"), overwrite = TRUE)
