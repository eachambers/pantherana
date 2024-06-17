library(terra)
library(raster)
library(tidyverse)
library(here)
library(cowplot)
library(elevatr)
library(mxmaps)
library(tidyterra)
library(rgdal)
theme_set(theme_cowplot())
setwd("~/Box Sync/Rana project/ddRADseq/ALL_RANA/Landgen")

## The following file generates population assignment files (Imap files) for
## BPP and HHSD analyses using the results from ADMIXTURE.
##    (1) Gathers loci numbers from the lmiss Plink report
##    (2) Generate Imap file from ADMIXTURE results

##    FILES REQUIRED:
##            *_metadata.txt # metadata for each of the four separate assemblies
##            PC_layers/ # directory with PC layers (see `Env_data.R` script to generate)


# Import coordinates ------------------------------------------------------

files <- list.files(here("data"), pattern = "metadata.txt", full.names = TRUE)
fs <- 1:length(files)
coords <- 
  fs %>% 
  lapply(function(x) {
    dat <- readr::read_tsv(files[x], col_names = TRUE) %>% 
      dplyr::rename(x = long, y = lat) %>% 
      drop_na() %>% 
      as.data.frame()
    return(dat)
  }) %>% 
  dplyr::bind_rows() %>% 
  dplyr::select(-dataset) %>% 
  distinct() # should be 507 samples

# Check that there are no duplicates
coords %>% 
  filter(duplicated(.[["Bioinformatics_ID"]]))


# Import map data ---------------------------------------------------------

# 33.486042, -119.324656
x = c(-119.32, max(coords$x), -119.32, max(coords$x))
y = c(33.486, 33.486, min(coords$y), min(coords$y))
bounds <- data.frame(x, y)

mxx = c(-119.32, max(coords$x), -119.32, -86.261)
mxy = c(33.486, 33.486, 14.1596, 14.1596)
mxbounds <- data.frame(x = mxx, y = mxy)

# Import DEM using elevatr
dem <- elevatr::get_elev_raster(bounds, prj = "WGS84", z = 6)
demmx <- elevatr::get_elev_raster(mxbounds, prj = "WGS84", z = 6)



# Define min value (otherwise it'll go negative) 
min <- 0
dem2 <- dem
dem2[dem <= min] <- NA
# Convert to SpatRaster for plotting
dem <- terra::rast(dem2)

demmx2 <- demmx
demmx2[demmx <= min] <- NA
# Convert to SpatRaster for plotting
demmx <- terra::rast(demmx2)

# Create color palette for DEM
pal = c("#844836", "#b07156", "#b99364", "#dac3a9", "#8fa997") # , "#6baeb0", "#376988"
palette = grDevices::colorRampPalette(pal)(100)
palette <- rev(palette)

## Build background map
world <- map_data("world")
centamer <- filter(world, region == "Mexico" | region == "Belize" | region == "Guatemala" | region == "Honduras" | region == "Nicaragua" | region == "Costa Rica" | region == "Panama")

# Load Mexico state map using mxmaps package
data("mxstate.map")
mxstate.map$group <- as.numeric(mxstate.map$group)

map_data <-
  bind_rows(centamer, mxstate.map)

# Load biogeographic provinces of Mexico which were downloaded here: http://geoportal.conabio.gob.mx/metadatos/doc/html/rbiog4mgw.html
provinces <- readOGR(dsn = "~/Box Sync/Rana project/Range maps/biogeo_provinces_mx/rbiog4mgw.shp", stringsAsFactors = F)


# Build plots -------------------------------------------------------------

ggplot() +
  geom_spatraster(data = dem) +
  scale_fill_gradientn(colors = palette, na.value = NA) +
  theme_map() +
  # geom_polygon(data = mxstate.map, aes(x = long, y = lat, group = group), color = "white", fill = NA, linewidth = 0.25) +
  geom_polygon(data = centamer, aes(x = long, y = lat, group = group), color = "black", fill = NA, linewidth = 0.25) +
  # geom_point(data = coords, aes(x = x, y = y), color = "white") +
  geom_polygon(data = provinces, aes(x = long, y = lat, group = group), color = "black", fill = NA, linewidth = 0.25, linetype = "dotted")

# Plot Mexican biogeographic provinces (export 12x8)
ggplot() +
  geom_spatraster(data = demmx) +
  scale_fill_gradientn(colors = palette, na.value = NA) +
  theme_map() +
  # geom_polygon(data = mxstate.map, aes(x = long, y = lat, group = group), color = "white", fill = NA, linewidth = 0.25) +
  # geom_point(data = coords, aes(x = x, y = y), color = "white") +
  geom_polygon(data = provinces, aes(x = long, y = lat, group = group), color = "white", fill = NA, linewidth = 0.25) + # linetype = "dashed"
  geom_polygon(data = centamer %>% dplyr::filter(region == "Mexico"), aes(x = long, y = lat, group = group), color = "black", fill = NA, linewidth = 0.5)


# Map environmental PCs ---------------------------------------------------

library(viridis)

# Load env data
forr_env <- raster::stack(list.files("PC_layers/", full.names = TRUE))
forr_env <- raster::readAll(forr_env)

# Simple plotting:
raster::plot(forr_env[[3]])

# ggplot plotting:
pc3 <- terra::as.data.frame(forr_env[[3]], xy = TRUE)

pc3 %>% 
  # drop_na(forreri_PC3) %>% 
  ggplot() +
  geom_raster(data = pc3, aes(x = x, y = y, fill = forreri_PC3)) +
  theme_map() +
  # scale_fill_viridis(option = "magma", na.value = "white") +
  scale_fill_viridis_c(na.value = "white") +
  coord_fixed()

library(stars)
tif = read_stars("MexicanBiogeographicProvinces.tif", package = "stars")
sf = st_as_sf(tif)


# forreri type localities -------------------------------------------------

# Read in type localities
types <- read_tsv(here("data", "forreri_type_localities.txt"), col_names = TRUE)

## Build background map
world <- map_data("world")
centamer <- filter(world, region == "Belize" | region == "Guatemala" | region == "Honduras" | region == "Nicaragua" | region == "Costa Rica")
# Load Mexico state map using mxmaps package
data("mxstate.map")
mxstate.map$group <- as.numeric(mxstate.map$group)
# centamer$group <- as.factor(centamer$group)

map_data <-
  bind_rows(centamer, mxstate.map)

base_map <-
  map_data %>% 
  ggplot() +
  geom_polygon(aes(x = long, y = lat, group = group), color = "white", fill = "#e6e6e6", linewidth = 0.25) +
  theme_map() +
  geom_point(data = final, aes(long, lat), size = 2) +
  geom_point(data = types, aes(x = longitude, y = latitude), size = 2, color = "red") +
  coord_fixed()

base_map +
  geom_point(data = types, aes(x = longitude, y = latitude), size = 2)
