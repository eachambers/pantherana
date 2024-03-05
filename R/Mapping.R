library(terra)
library(raster)
library(tidyverse)
library(here)
library(cowplot)
library(elevatr)
theme_set(theme_cowplot())
setwd("~/Box Sync/Rana project/ddRADseq/ALL_RANA/Landgen")

## The following file generates population assignment files (Imap files) for
## BPP and HHSD analyses using the results from ADMIXTURE.
##    (1) Gathers loci numbers from the lmiss Plink report
##    (2) Generate Imap file from ADMIXTURE results

##    FILES REQUIRED:
##            forreri.ustr # n=104

coords <- read_tsv(here("data", "forreri_coords.txt"), col_names = TRUE) %>% 
  rename(x = long, y = lat) %>% 
  dplyr::select(x,y) %>% 
  as.data.frame()
dem <- get_elev_raster(coords, prj = "WGS84", z = 4)

# Define min value (otherwise it'll go negative) 
min <- 0
dem2 <- dem
dem2[dem <= min] <- NA

pal = c("#844836", "#b07156", "#b99364", "#dac3a9", "#8fa997") # , "#6baeb0", "#376988"
palette = grDevices::colorRampPalette(pal)(100)
palette <- rev(palette)
raster::plot(dem2, col = palette, axes = FALSE) # export 12x8
points(types)

ggplot() +
  geom_raster(data = dem2, aes(x = x, y = y, fill = file1246818cc0fe)) +
  # scale_fill_manual(values = palette, na.value = "white") +
  theme_map()


# Import tif of provinces
provinces <- raster("MexicanBiogeographicProvinces.tif")
raster::plot(provinces)

rast <- terra::rast("MexicanBiogeographicProvinces.tif")


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
