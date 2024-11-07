library(raster)
# library(rgeos)
library(dismo) # where maxent function lives
library(here)
library(rJava)
library(tidyverse)
library(terra)
library(tidyterra)
library(cowplot)
library(algatr)
library(RStoolbox)
library(viridis)
library(countrycode)
library(CoordinateCleaner) # install_github("ropensci/CoordinateCleaner")
library(rgbif)
library(sf)
library(rnaturalearth)
theme_set(theme_cowplot())

utils::download.file(url = "https://github.com/mrmaxent/Maxent/blob/master/ArchivedReleases/3.4.3/maxent.jar",
                     destfile = paste0(system.file("java", package = "dismo"), 
                                       "/maxent.jar"), mode = "wb")  ## wb for binary file, otherwise maxent.jar cannot execute

# GBIF.org (17 June 2024) GBIF Occurrence Download https://doi.org/10.15468/dl.zqykdz

# Gather and clean GBIF data ----------------------------------------------

# Observations for R. forreri (n=1848 obs)
dat <- occ_search(scientificName = "Lithobates forreri", 
                  limit = 5000, 
                  hasCoordinate = TRUE)
dat <- dat$data

# Subset out relevant cols
dat <- dat %>%
  dplyr::select(species, decimalLongitude, 
                decimalLatitude, countryCode, individualCount,
                gbifID, family, taxonRank, coordinateUncertaintyInMeters,
                year, basisOfRecord, institutionCode, datasetName)

# Convert country code from ISO2c to ISO3c
dat$countryCode <- countrycode(dat$countryCode, 
                               origin =  'iso2c',
                               destination = 'iso3c')

# Clean coordinates; 31 records were flagged (1817 remain)
dat <- data.frame(dat)
flags <- clean_coordinates(x = dat, 
                           lon = "decimalLongitude", 
                           lat = "decimalLatitude",
                           countries = "countryCode",
                           species = "species",
                           tests = c("capitals", "centroids",
                                     "equal", "zeros", "countries"))

# Exclude flagged records
dat_cl <- dat[flags$.summary,]

# Further data cleaning
dat_cl <- dat_cl %>%
  # Remove records with low coordinate precision (uncertainty > 100 km)
  filter(coordinateUncertaintyInMeters / 1000 <= 100 | 
           is.na(coordinateUncertaintyInMeters)) %>% 
  # Only use museum collection data
  filter(basisOfRecord == "PRESERVED_SPECIMEN") %>% 
  # Remove localities with > 99 samples
  filter(individualCount < 99 | is.na(individualCount)) %>% 
  # Remove observations from before 1945
  filter(year > 1945)

# Finally, prune out observations based on species range
filepath <- here("data", "forreri_range", "data_0.shp")
wgs84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
nat_range <- terra::vect(filepath, crs = wgs84)
nat_range$species <- "Lithobates forreri" # species name must be consistent with GBIF data

range_flags <- cc_iucn(x = dat_cl,
                       range = nat_range,
                       lon = "decimalLongitude",
                       lat = "decimalLatitude",
                       value = "flagged",
                       buffer = 30000) # 30 km buffer around range
dat_fin <- dat_cl[range_flags, ]


# Plot results ------------------------------------------------------------

## Build background map
world <- map_data("world")
centamer <- filter(world, region == "Mexico" | region == "Belize" | region == "Guatemala" | region == "Honduras" | region == "Nicaragua" | region == "Costa Rica" | region == "Panama" | region == "El Salvador")

# Load Mexico state map using mxmaps package
data("mxstate.map")
mxstate.map$group <- as.numeric(mxstate.map$group)

map_data <-
  bind_rows(centamer, mxstate.map)

plo <- sf::st_as_sf(nat_range)

ggplot() +
  theme_map() +
  geom_polygon(data = map_data, aes(x = long, y = lat, group = group), 
               color = "black", fill = NA, linewidth = 0.25) +
  coord_equal() +
  geom_sf(data = plo, aes(fill = species), alpha = 0.5) +
  scale_fill_manual(values = c("#73806d")) +
  theme(legend.position = "none",
        axis.title = element_blank()) +
  geom_point(data = dat_fin,
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "darkred",
             size = 0.5)


# Divide into different species -------------------------------------------

# Boundary of forreri is between the states of Jalisco and Nayarit
dat_forr <-
  dat_fin %>% 
  filter(decimalLatitude > 20.549714)

# Southern boundary of floresi is between the states of Michoacan and Guerrero
dat_flor <-
  dat_fin %>% 
  filter(decimalLatitude < 20.549714) %>% 
  filter(decimalLatitude > 17.989559 & decimalLongitude < -101.862222)

# Plot them out together
ggplot() +
  theme_map() +
  geom_polygon(data = map_data, aes(x = long, y = lat, group = group), 
               color = "black", fill = NA, linewidth = 0.25) +
  coord_equal() +
  geom_sf(data = plo, aes(fill = species), alpha = 0.5) +
  scale_fill_manual(values = c("#73806d")) +
  theme(legend.position = "none",
        axis.title = element_blank()) +
  geom_point(data = dat_forr,
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "gray50",
             size = 0.5) +
  geom_point(data = dat_flor,
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "gray20",
             size = 0.5)


# My own data -------------------------------------------------------------

# Add in my own samples and add species names based on admixture results
forreri_eac <- read_tsv(here("data", "forreri_K5.txt")) %>% 
  mutate(species = case_when(max_K == "X1" ~ "forreri",
                             max_K == "X2" ~ "miadis",
                             max_K == "X3" ~ "hillisi",
                             max_K == "X4" ~ "arcelia",
                             max_K == "X5" ~ "floresi")) %>% 
  dplyr::select(Bioinformatics_ID, species, lat, long) %>% 
  filter(!is.na(lat))

# Combine GBIF with my data
enm_dat_forr <- bind_rows(forreri_eac %>% 
                            filter(species == "forreri") %>% 
                            rename(ID = Bioinformatics_ID) %>% 
                            dplyr::select(-species),
                          dat_forr %>% 
                            rename(ID = gbifID,
                                   lat = decimalLatitude,
                                   long = decimalLongitude) %>% 
                            dplyr::select(ID, lat, long))

enm_dat_flor <- bind_rows(forreri_eac %>% 
                            filter(species == "floresi") %>% 
                            rename(ID = Bioinformatics_ID) %>% 
                            dplyr::select(-species),
                          dat_flor %>% 
                            rename(ID = gbifID,
                                   lat = decimalLatitude,
                                   long = decimalLongitude) %>% 
                            dplyr::select(ID, lat, long))


# Enviromental data -------------------------------------------------------

# Retrieve worldclim tiles for coords
env_forr <- algatr::get_worldclim(coords = enm_dat_forr %>% 
                                    rename(x = long, y = lat) %>% 
                                    dplyr::select(x, y),
                                  res = 0.5,
                                  buff = 0.01,
                                  save_output = TRUE)

env_flor <- algatr::get_worldclim(coords = enm_dat_flor %>% 
                                    rename(x = long, y = lat) %>% 
                                    dplyr::select(x, y),
                                  res = 0.5,
                                  buff = 0.01,
                                  save_output = TRUE)

# Take a look at WorldClim tiles that were retrieved
plot(env_forr[[1]], col = magma(100), axes = FALSE)
plot(env_flor[[1]], col = magma(100), axes = FALSE)
points(enm_dat_flor %>% 
         rename(x = long, y = lat) %>% 
         dplyr::select(x, y), 
       pch = 19)

env_pcs_forr <- RStoolbox::rasterPCA(env_forr, spca = TRUE)
env_pcs_flor <- RStoolbox::rasterPCA(env_flor, spca = TRUE)

# Save top 3 env PCs
purrr::map(1:3, ~terra::writeRaster(env_pcs_forr$map[[.x]],
                                    file = paste0(here("data", "ENM_PC_layers"), "/ENM_forreri_PC", .x, ".tif"),
                                    overwrite = TRUE))
purrr::map(1:3, ~terra::writeRaster(env_pcs_flor$map[[.x]],
                                    file = paste0(here("data", "ENM_PC_layers"), "/ENM_floresi_PC", .x, ".tif"),
                                    overwrite = TRUE))

# Read in resulting rasterPCAs stacks
env_pcs_forr <- raster::stack(list.files(here("data", "ENM_PC_layers"), pattern = "forreri", full.names = TRUE))
env_pcs_flor <- raster::stack(list.files(here("data", "ENM_PC_layers"), pattern = "floresi", full.names = TRUE))

# ggplot() +
#   geom_spatraster(data = env_pcs_flor[[1]]) +
#   geom_point(data = enm_dat_flor, aes(x = long, y = lat)) +
#   scale_fill_viridis_c(na.value = NA)


# Clip layers -------------------------------------------------------------

coordinates(enm_dat_forr) <- ~long + lat
# Create a 50m buffer around the occurrence data
enm_dat_forr_buff <- terra::buffer(enm_dat_forr, width = 50000)
# Crop study area (rectangle)
study_area <- terra::crop(rast(env_pcs_forr), extent(enm_dat_forr_buff))  
# the 'study area' created by extracting the buffer area from the raster stack
vec <- terra::vect(enm_dat_forr_buff)
study_area <- terra::mask(study_area, vec)
# Check it
plot(study_area)

terra::writeRaster(study_area,
                   filename = paste0(here("data", "study_area"), "/", names(study_area), "_forr.asc"))


# Select background points ------------------------------------------------

# select background points from this buffered area; when the number provided 
# to set.seed() function, the same random sample will be selected in the next line			
# use this code before the sampleRandom function every time, if you want to get
# the same "random samples"
set.seed(808)
bg <- raster::sampleRandom(x = raster(study_area),
                   size = 20000,
                   na.rm = T,
                   sp = T) # return SpatialPoints 

plot(study_area[[1]])
# add the background points to the plotted raster
plot(bg, add = T) 
# add the occurrence data to the plotted raster
coordinates(dat_forr) <- ~decimalLongitude + decimalLatitude
plot(dat_forr, add = T, col = "red")


# Split into training and testing -----------------------------------------

set.seed(808)
# Randomly select 30% for training
selected <- sample(1:nrow(dat_forr), nrow(dat_forr) * 0.3)

train <- dat_forr[selected, ]
test <- dat_forr[-selected, ]


# Format data for Maxent --------------------------------------------------

# Extract environmental values using the training dataset; PRESENCE
p <- extract(env_pcs_forr, train)
# Extract environmental values using the testing dataset
p_test <- extract(env_pcs_forr, test)
# Extract environmental values using the background (n=20000); ABSENCE
a <- extract(env_pcs_forr, bg)

# Maxent requires presence/absence data, encoded as 0s and 1s
pa <- c(rep(1, nrow(p)), rep(0, nrow(a)))
pder <- as.data.frame(rbind(p, a))


# Run Maxent --------------------------------------------------------------

forr_pa <- read.table(here("data", "maxent", "forr_pa.txt"), header = TRUE)
forr_pder <- read.table(here("data", "maxent", "forr_pder.txt"), header = TRUE)

cat(class(forr_pder),"  ",  class(forr_pa))

# Train Maxent with tabular data
mod <- dismo::maxent(x = forr_pder, # predictors: env conditions
                     p = forr_pa,
                     args = c('hinge = false', 'threshold = false'),
                     silent = F)

mod

