library(terra)
library(raster)
library(tidyverse)
library(here)
library(cowplot)
library(elevatr)
library(mxmaps)
library(tidyterra)
# library(rgdal)
library(viridis)
theme_set(theme_cowplot())

source(here("data_viz", "mapping_functions.R"))
source(here("data_viz", "pop_gen_functions.R"))
source(here("analysis", "4_hhsd", "hhsd_functions.R"))
source(here("data_viz", "rana_colors.R"))

## The following file:
##    (1) Imports coordinates
##    (2) Imports map data
##    (3) Range maps: Figs. 1C, 2A, and 3A
##    (4) Fig. S3: Builds map of elevation with MX biogeographic regions and state boundaries
##    (4) Map environmental PCs
##    (5) forreri type localities
##    (7) Plots the Plateau of Guatemala (for R. macroglossa type locality)

##    FILES REQUIRED:
##            *_metadata.txt # metadata for each of the four separate assemblies
##            PC_layers/ # Environmental PC layers (generated using `Env_data.R`)


# (1) Import coordinates --------------------------------------------------

files <- list.files(here("data", "2_Data_processing", "data_files_input_into_scripts"), pattern = "_metadata.txt", full.names = TRUE)
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
  filter(duplicated(.[["Bioinformatics_ID"]])) # should be empty


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
centamer <- filter(world, region == "Mexico" | region == "Belize" | region == "Guatemala" | region == "Honduras" | region == "Nicaragua" | region == "Costa Rica" | region == "Panama" | region == "El Salvador")

# Load Mexico state map using mxmaps package
data("mxstate.map")
mxstate.map$group <- as.numeric(mxstate.map$group)

map_data <-
  bind_rows(centamer, mxstate.map)

usa <- map_data("state")


# Range maps --------------------------------------------------------------

# Read in all range map shapefiles
ranges <- import_range_maps(path = here("data", "4_Data_visualization", "data_files_input_into_scripts", "range_maps"))
# Read in type localities
tls <- read_tsv(here("data", "4_Data_visualization", "data_files_input_into_scripts", "type_localities.txt"))
# Get colors
spp_colors <- mapping_colors(here("data", "4_Data_visualization", "data_files_input_into_scripts", "type_localities_rec.txt"))

# Make range map
p_ranges <-
  ggplot() +
  geom_polygon(data = usa, aes(x = long, y = lat, group = group), color = NA, fill = "#e6e6e6") + # color = "#969696" if CENTAM
  geom_polygon(data = map_data, aes(x = long, y = lat, group = group), color = NA, fill = "#e6e6e6") +
  theme_map() +
  # forreri
  geom_sf(data = ranges$forr, fill = spp_colors$forr, color = NA, alpha = 0.75) +
  # foothills species
  geom_sf(data = ranges$yava, color = NA, fill = spp_colors$yava, alpha = 0.75) +
  geom_sf(data = ranges$magn, color = NA, fill = spp_colors$magn, alpha = 0.75) +
  geom_sf(data = ranges$omil, color = NA, fill = spp_colors$omil, alpha = 0.75) +
  geom_point(data = ranges$aten_long, aes(x = Longitude, y = Latitude), fill = "#428b9b", color = "black", pch = 21) +
  geom_point(data = ranges$aten_short, aes(x = Longitude, y = Latitude), fill = "#c0a06f", color = "black", pch = 21) +
  # ATL_MXPL species
  geom_sf(data = ranges$berl, color = NA, fill = spp_colors$berl, alpha = 0.75) +
  geom_sf(data = ranges$macro1, color = NA, fill = spp_colors$macro, alpha = 0.75) +
  geom_sf(data = ranges$macro2, color = NA, fill = spp_colors$macro, alpha = 0.75) +
  geom_sf(data = ranges$spec, color = NA, fill = spp_colors$spec, alpha = 0.75) +
  # CENTAM species
  geom_sf(data = ranges$spnov, color = NA, fill = spp_colors$spnov, alpha = 0.75) +
  geom_sf(data = ranges$lenca, color = NA, fill = spp_colors$lenca, alpha = 0.75) +
  coord_sf(ylim = c(min(centamer$lat), 36.97)) + # maxx = -82.13
  geom_polygon(data = usa, aes(x = long, y = lat, group = group), color = "white", fill = NA, linewidth = 0.25) + # color = "#969696" if CENTAM
  geom_polygon(data = map_data, aes(x = long, y = lat, group = group), color = "white", fill = NA, linewidth = 0.25)

# Add type localities of recognized and unrecognized species on to map
p_ranges +
  geom_point(data = tls %>% filter(Status == "rec"), aes(Longitude, Latitude, color = Species), size = 3) +
  scale_color_manual(values = c("Arcelia" = "#d2a3a6", "Atenquique_long" = "#428b9b", "Atenquique_short" = "#c0a06f", 
                                "berlandieri" = "#984625", "forreri" = "#73806d",
                                "lenca" = "#924c62", "macroglossa" = "#e0895a", "magnaocularis" = "#f7cd5e",
                                "omiltemana" = "#ba94a6", "spectabilis" = "#80a4bc", "yavapaiensis" = "#054051")) +
  geom_point(data = tls %>% filter(Status == "rec"), aes(Longitude, Latitude), shape = 1, color = "black", size = 3) +
  geom_point(data = tls %>% filter(Status == "notrec"), aes(Longitude, Latitude), color = "white", size = 2) +
  geom_point(data = tls %>% filter(Status == "notrec"), aes(Longitude, Latitude), shape = 1, color = "black", size = 2) +
  theme(legend.position = "none") # export 10x8

# Instead of type localities, add sampling coordinates, colorized by species
# allsamps <- read_tsv(here("data", "all_samps.txt")) %>%
#   na.omit()
# p_ranges +
#   geom_point(data = allsamps, aes(x = Longitude, y = Latitude, color = Species), size = 1) +
#   scale_color_manual(values = c("aten_long" = "#428b9b", "aten_short" = "#c0a06f",
#                               "berlandieri" = "#984625", "forreri" = "#73806d",
#                               "lenca" = "#924c62", "macroglossa" = "#e0895a", "magnaocularis" = "#f7cd5e",
#                               "omiltemana" = "#ba94a6", "spectabilis" = "#80a4bc", "yavapaiensis" = "#054051")) +
#   geom_point(data = allsamps, aes(x = Longitude, y = Latitude), shape = 1, color = "black", size = 1)
  

# Fig. 4: HHSD sampling ---------------------------------------------------

dataset_name = "forreri" # forreri / foothills / mxpl

if (dataset_name == "foothills") metadata <- read_tsv(here("data", "2_Data_processing", "data_files_input_into_scripts", "PACMX_metadata.txt"), col_names = TRUE)
if (dataset_name == "forreri") metadata <- read_tsv(here("data", "2_Data_processing", "data_files_input_into_scripts", "forreri_metadata.txt"), col_names = TRUE)
if (dataset_name == "mxpl") metadata <- read_tsv(here("data", "2_Data_processing", "data_files_input_into_scripts", "ATL_MXPL_metadata.txt"), col_names = TRUE)
if (dataset_name == "forreri") colvals = c("hilli" = "#428b9b", "miad" = "#73806d", "arce" = "#d2a3a6", "forr" = "gray74", "flor" = "gray30")
if (dataset_name == "foothills") colvals = c("omig" = "#a8a2ca", "magn" = "#f7cd5e", "omio" = "#ba94a6", 
                                             "yava" = "#054051", "ates" = "#c0a06f", "atel" = "#428b9b")
if (dataset_name == "mxpl") colvals = c("spec" = "#80a4bc", "chic" = "#6f82b7", "berl" = "#984625", "neov" = "#e97490")

dat <- retrieve_hhsd_coding(dataset_name, save_imap = TRUE)
final <- left_join(dat, metadata, by = "Bioinformatics_ID") %>% na.omit()
states <- map_data("state")

p_base <-
  ggplot() +
  geom_polygon(data = map_data, aes(x = long, y = lat, group = group), color = "white", fill = "#e6e6e6", linewidth = 0.5) + 
  geom_polygon(data = states, aes(x = long, y = lat, group = group), color = "white", fill = "#e6e6e6", linewidth = 0.5) +
  geom_point(data = final, aes(long, lat, color = pop), size = 4) + # if mxpl, size = 3
  geom_point(data = final, aes(long, lat), shape = 1, color = "black", size = 4) + # if mxpl, size = 3
  scale_color_manual(values = colvals) +
  theme_map() +
  theme(legend.position = "none") # export 6x4

# Fix the limits depending on the dataset
ymin <- min(final$lat)-2
ymax <- max(final$lat)+2
xmin <- min(final$long)-2
xmax <- max(final$long)+2

p_base +
  coord_fixed(xlim = c(xmin, xmax), ylim = c(ymin, ymax))


# Fig. S3: DEM with biogeo provinces --------------------------------------

# Import biogeographic provinces of Mexico which were downloaded here: http://geoportal.conabio.gob.mx/metadatos/doc/html/rbiog4mgw.html
provinces <- sf::st_read(here("data", "4_Data_visualization", "data_files_input_into_scripts", "biogeo_provinces_mx", "rbiog4mgw.shp"))

# Plot Mexican biogeographic provinces (export 12x8)
ggplot() +
  geom_spatraster(data = demmx) +
  scale_fill_gradientn(colors = palette, na.value = NA) +
  theme_map() +
  geom_polygon(data = mxstate.map, aes(x = long, y = lat, group = group), color = "grey40", fill = NA, linewidth = 0.25, linetype = "dashed") +
  # geom_point(data = coords, aes(x = x, y = y), color = "pink") +
  geom_sf(data = provinces, color = "white", fill = NA, linewidth = 0.25) + # linetype = "dashed"
  geom_polygon(data = centamer %>% dplyr::filter(region == "Mexico"), aes(x = long, y = lat, group = group), color = "black", fill = NA, linewidth = 0.5)


# forreri type localities -------------------------------------------------

# Read in type localities
types <- read_tsv(here("data", "4_Data_visualization", "data_files_input_into_scripts", "forreri_typelocalities.txt"), col_names = TRUE)

## Build background map
world <- map_data("world")
usa <- map_data("state")
centamer <- filter(world, region == "Belize" | region == "Guatemala" | 
                     region == "Honduras" | region == "Nicaragua" | 
                     region == "Costa Rica" | region == "Panama")
# Load Mexico state map using mxmaps package
data("mxstate.map")
mxstate.map$group <- as.numeric(mxstate.map$group)

map_data <-
  bind_rows(centamer, mxstate.map) # keep usa separate because order is screwing up lines

base_map <-
  map_data %>% 
  ggplot() +
  geom_polygon(aes(x = long, y = lat, group = group), color = "white", fill = "#e6e6e6", linewidth = 0.25) +
  theme_map() +
  geom_point(data = final, aes(long, lat), size = 2) +
  geom_point(data = types, aes(x = Longitude, y = Latitude), size = 2, color = "red") +
  coord_fixed()


# Plateau of Guatemala ----------------------------------------------------

guat <- filter(world, region == "Guatemala")
xguat = c(min(guat$long)-1, max(guat$long)+1, min(guat$long)-1, max(guat$long)+1)
yguat = c(max(guat$lat)+1, max(guat$lat)+1, min(guat$lat)-1, min(guat$lat-1))
guat_bounds <- data.frame(x = xguat, y = yguat)
demguat <- elevatr::get_elev_raster(guat_bounds, prj = "WGS84", z = 8)
min <- 0
demguat2 <- demguat
demguat2[demguat <= min] <- NA
# Convert to SpatRaster for plotting
demguat <- terra::rast(demguat2)
  
plat_guat <- sf::st_read(here("data", "4_Data_visualization", "data_files_input_into_scripts", "range_maps", "plateau_de_guatemala.shp"))
guat_cities <- read_tsv(here("data", "4_Data_visualization", "data_files_input_into_scripts", "guatemala_bocourt.txt"))

ggplot() +
  geom_spatraster(data = demguat) +
  scale_fill_gradientn(colors = palette, na.value = NA) +
  theme_map() +
  geom_polygon(data = guat, aes(x = long, y = lat, group = group), color = "black", fill = NA, linewidth = 0.5) +
  # geom_polygon(data = centamer %>% dplyr::filter(region == "Mexico" | region == "Guatemala" | region == "Honduras" | region == "El Salvador"), aes(x = long, y = lat, group = group), color = "black", fill = NA, linewidth = 0.25) +
  geom_sf(data = plat_guat, fill = "darkorange", color = "white", alpha = 0.4, linetype = "dashed", linewidth = 1) +
  # geom_point(data = guat_cities, aes(x = long, y = lat)) +
  coord_sf(ylim = c(min(guat$lat)-0.5, max(guat$lat)+0.5),
           xlim = c(min(guat$long)-0.5, max(guat$long)+0.5))  
ggsave(here("plots", "guatemala.pdf"), height = 8, width = 8)
