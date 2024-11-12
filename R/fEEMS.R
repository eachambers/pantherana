library(tidyverse)
library(here)
library(vcfR)
library(dggridR)
library(sf)
library(raster)
library(terra)
library(mxmaps)
library(cowplot)

## The following file:
##    (1) Generates input files for FEEMS
##    (2) Creates a discrete global grid for FEEMS
##    (3) Imports FEEMS output files
##    (4) Makes FEEMS map (Fig. 5C)
##    (5) Plots FEEMS cross-validation analysis results (Fig. S8)

##    FILES REQUIRED:
##            "forreri_metadata.txt" and "forreri_0.25miss_ldp_n103.recode.vcf"
##            Shapefile with Mexico + CentAm: "mx_centam_merged.shp" for DGG
##            FEEMS output files


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


# (3) Import FEEMS output files -------------------------------------------
# Adapted from https://github.com/karolisr/pitcairnia-dr-nrv/blob/93534b6906d2a59a821c19deaef649eb1a3fb283/25-geo-dist.R

feems_nodes <- read_csv(here("data", "feems_nodes.csv"), col_names = "node_id", col_types = "i")
feems_edges <- read_csv(here("data", "feems_edges.csv"), col_names = c("n1", "n2"), col_types = "ii")
feems_w <- read_csv(here("data", "feems_w.csv"), col_names = "w", col_types = "d")
feems_node_pos <- read_csv(here("data", "feems_node_pos_T.csv"), 
                           col_names = c("lon", "lat", "nsamp"), 
                           col_types = "dd")

# Combine node and edge data
nodes <- add_column(feems_nodes, feems_node_pos)
edges <- add_column(feems_edges, feems_w)

# Add coords for nodes
edges_n1 <- left_join(edges, nodes, by = c("n1" = "node_id"))
edges_n2 <- left_join(edges_n1, nodes, by = c("n2" = "node_id"))
feems_edges <- dplyr::select(edges_n2, lon1 = lon.x, lat1 = lat.x, lon2 = lon.y, lat2 = lat.y, w)

# From original FEEMS paper: the effective migration rate on the log10 scale, relative to 
# the overall migration rate across the habitat
# Add log-transformed weights
feems_edges <- feems_edges %>% 
  mutate(logw = log10(w),
         logwmean = log10(w) - mean(log10(w)))

projcrs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
df1 <- st_as_sf(x = feems_edges,                         
               coords = c("lon1", "lat1"),
               crs = projcrs) %>% 
  dplyr::select(-lat2, -lon2)
df2 <- st_as_sf(x = feems_edges,                         
                coords = c("lon2", "lat2"),
                crs = projcrs) %>% 
  dplyr::select(-lat1, -lon1)

combined_coords = cbind(st_coordinates(df1$geometry), st_coordinates(df2$geometry))
# Construct linestrings by row:
geometry = st_sfc(
  lapply(1:nrow(combined_coords),
         function(i) {
           st_linestring(matrix(combined_coords[i,], ncol = 2, byrow = TRUE))
         }))

feems_new_edges <- cbind(feems_edges, geometry) %>% 
  dplyr::select(-c(lat1, lat2, lon1, lon2)) %>% 
  st_as_sf() %>%
  st_set_crs(4326)

# Crop FEEMS edges to borders of MX and CentAm
feems_cropped <- st_intersection(feems_new_edges, border)


# (4) Make FEEMS map ------------------------------------------------------

# feems.colors <- c('#994000', '#CC5800', '#FF8F33', '#FFAD66',
#                   '#FFCA99', '#FFE6CC', '#FBFBFB', '#CCFDFF',
#                   '#99F8FF', '#66F0FF', '#33E4FF', '#00AACC',
#                   '#007A99')

outer <- read_delim(here("data", "forreri_outer.txt"), col_names = c("long", "lat"))

data("mxstate.map")
mxstate.map$group <- as.numeric(mxstate.map$group)

ggplot() +
  geom_sf(data = border, fill = "grey90", color = NA) + 
  geom_sf(data = feems_cropped, aes(color = logwmean), alpha = 0.75) +
  # scale_color_gradientn(colors = feems.colors) +
  scale_color_gradient2(low = "#FF8F33",
                        mid = "white",
                        high = "#33E4FF") +
  geom_sf(data = border, fill = NA, color = "black") +
  coord_sf(xlim = unique(outer$long), ylim = unique(outer$lat)) +
  theme_map() +
  geom_polygon(data = mxstate.map, aes(x = long, y = lat, group = group), color = "gray48", fill = NA, linewidth = 0.25) +
  geom_point(data = nodes %>% filter(nsamp > 0), aes(x = lon, y = lat, size = nsamp), color = "#a8a2ca") +
  geom_point(data = coords, aes(x = x, y = y), color = "black")

ggsave(paste0(here("plots"), "/forreri_feems_logwmean.pdf"), width = 10, height = 8, units = "in")


# (5) FEEMS cross-validation ----------------------------------------------

cv <- read_csv(here("data", "mean_cv_err.csv"), col_names = "CV_error")
lamb <- read_csv(here("data", "lamb_grid.csv"), col_names = "lambda")
cv <- bind_cols(cv, lamb)

theme_set(theme_cowplot())
cv %>% 
  ggplot(aes(x = log10(lambda), y = CV_error)) +
  geom_vline(xintercept = log10(20), color = "red", linetype = "dashed") +
  geom_point() +
  ylab("Cross-validation error")

ggsave(here("plots", "FEEMS_cv_plot.pdf"), width = 5, height = 4)
