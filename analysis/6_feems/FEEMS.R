library(tidyverse)
library(here)
library(vcfR)
library(dggridR)
library(sf)
library(raster)
library(terra)

## The following file:
##    (1) Generates input files for FEEMS
##    (2) Creates a discrete global grid for FEEMS

##    FILES REQUIRED:
##            "forreri_metadata.txt" and "forreri_0.25miss_ldp_n103.recode.vcf"
##            Shapefile with Mexico + CentAm: "mx_centam_merged.shp" for DGG


# (1) Generate input files for FEEMS --------------------------------------

forreri_metadata <- read_tsv(here("data", "forreri_metadata.txt"))
coords <- forreri_metadata %>% na.omit()

# Verify ordering with vcf
vcf <- read.vcfR(here("data", "forreri_0.25miss_ldp_n103.recode.vcf"))
inds <- data.frame(Bioinformatics_ID = colnames(vcf@gt[, -1]))
coords <- left_join(inds, coords) %>% 
  dplyr::rename(x = long, y = lat) %>% 
  dplyr::relocate(y, .after = x)

all(coords$Bioinformatics_ID == colnames(vcf@gt[, -1]))

# write_tsv(coords, here("data", "forreri_coords.txt"), col_names = FALSE)


# (2) Create DGG for FEEMS ------------------------------------------------

coords <- read_tsv(here("data", "forreri_coords.txt"), col_names = FALSE) # generated above
dggs <- dgconstruct(res = 8, projection = "ISEA", aperture = 4, topology = "TRIANGLE")

# Read in borders from shapefile
border <- sf::st_read(here("data"), layer = "mx_centam_merged")
st_crs(border) = 4326 # WGS 84

# Get a grid covering Mexico
forr_grid <- dgshptogrid(dggs, here("data", "mx_centam_merged.shp")) # path to shapefile?

# Plot to check
ggplot() +
  geom_sf(data = border, fill = NA, color = "black")   +
  geom_sf(data = forr_grid, fill = alpha("blue", 0.4), color = alpha("white", 0.4)) +
  geom_point(data = coords, aes(x = x, y = y), color = "red")

# Save grid
# st_crs(forr_grid) <- NA
sf::st_write(forr_grid, here("data", "forr_grid.shp"), append = FALSE)
