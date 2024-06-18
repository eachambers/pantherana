library(raster)
library(rgeos)
library(dismo)
library(here)
library(rJava)
library(tidyverse)
library(terra)
library(tidyterra)
library(cowplot)
library(algatr)
library(RStoolbox)
theme_set(theme_cowplot())

utils::download.file(url = "https://github.com/mrmaxent/Maxent/blob/master/ArchivedReleases/3.4.4/maxent.jar",
                     destfile = paste0(system.file("java", package = "dismo"), 
                                       "/maxent.jar"), mode = "wb")  ## wb for binary file, otherwise maxent.jar can not execute

# Build ENMs with sampling locations only ---------------------------------

# metadata <- read_tsv(here("data", "forreri_metadata.txt"))
# Import admixture results for K=5; X1 is forreri and X5 is floresi
admix <- read_tsv(here("data", "forreri_K5.txt")) %>% 
  filter(max_K %in% c("X1", "X5")) %>% 
  mutate(species = case_when(max_K == "X1" ~ "forreri",
                             max_K == "X5" ~ "floresi")) %>% 
  dplyr::select(lat, long, species) %>% 
  distinct() %>% 
  filter(!is.na(lat))

# Add in type localities
typelocs <- read_tsv(here("data", "forreri_typelocalities.txt")) %>% 
  rename(lat = Latitude, long = Longitude) %>% 
  filter(Species %in% c("adleri", "cora", "forreri", "floresi")) %>% 
  rename(species = Species)
# Join types with admixture, retaining only coords, removing dupes
coords <- bind_rows(admix, typelocs) # 19 observations

##### Re-run raster PCA based on all localities (increasing northward)
# Retrieve worldclim tiles for coords
env <- algatr::get_worldclim(coords = coords %>% 
                               rename(x = long, y = lat) %>% 
                               dplyr::select(x, y),
                             res = 0.5,
                             buff = 0.01,
                             save_output = TRUE)

# Take a look at WorldClim tiles that were retrieved
plot(env[[1]], col = magma(100), axes = FALSE)
points(coords %>% 
         rename(x = long, y = lat) %>% 
         dplyr::select(x, y), 
       pch = 19)

env_pcs <- RStoolbox::rasterPCA(env, spca = TRUE)

# Save top 3 env PCs
lapply(1:3, function(x) terra::writeRaster(env_pcs$map[[x]],
                                           file = paste0("ENM_PC_layers/ENM_forreri_PC", x, ".tif"),
                                           overwrite = TRUE))

# Read in resulting rasterPCA as a stack
# forr_env <- raster::stack(list.files(here("data", "ENM_PC_layers/"), full.names = TRUE))

rast <- rast(forr_env)

ggplot() +
  geom_spatraster(data = rast[[1]]) +
  geom_point(data = admix, aes(x = long, y = lat, color = species)) +
  scale_fill_viridis_c(na.value = NA)

# Make into spdf object
coordinates(admix) <- ~long + lat



# Build ENMs with occurrence data from GBIF -------------------------------

# GBIF.org (17 June 2024) GBIF Occurrence Download https://doi.org/10.15468/dl.zqykdz
# TODO We need to categorize samples as being floresi or forreri, etc.

occ <- read_tsv(here("data", "forreri_gbif.txt"))

## Clean data
occ_clean <- occ %>% 
  dplyr::select(decimalLatitude, decimalLongitude) %>% 
  # Remove samples without sampling coords
  filter(!is.na(decimalLatitude)) %>% 
  # Remove duplicated data
  distinct()

# Plot remaining observations
plot(forr_env[[1]])
plot(occ_clean, col = "red")  # the 'add=T' tells R to put the incoming data on the existing layer

rast <- rast(forr_env)

ggplot() +
  geom_spatraster(data = rast[[1]]) +
  geom_point(data = occ_clean, aes(x = decimalLongitude, y = decimalLatitude)) +
  scale_fill_viridis_c(na.value = NA)

# Remove sample in Baja?
