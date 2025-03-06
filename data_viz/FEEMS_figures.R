library(tidyverse)
library(here)
library(sf)
library(raster)
library(terra)
library(mxmaps)
library(cowplot)

## The following file:
##    (1) Imports FEEMS output files
##    (2) Makes FEEMS map (Fig. 5C)
##    (3) Plots FEEMS cross-validation analysis results (Fig. S9)

##    FILES REQUIRED:
##            FEEMS output files
##            forreri sampling coordinates


# (1) Import FEEMS output files -------------------------------------------

coords <- read_tsv(here("data", "forreri_metadata.txt")) %>% na.omit()

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

# Read in borders from shapefile
border <- sf::st_read(here("data"), layer = "mx_centam_merged")
st_crs(border) = 4326 # WGS 84

# Crop FEEMS edges to borders of MX and CentAm
feems_cropped <- st_intersection(feems_new_edges, border)


# (2) Make FEEMS map ------------------------------------------------------

outer <- read_delim(here("data", "forreri_outer.txt"), col_names = c("long", "lat"))

data("mxstate.map")
mxstate.map$group <- as.numeric(mxstate.map$group)

ggplot() +
  geom_sf(data = border, fill = "grey90", color = NA) +
  geom_sf(data = feems_cropped, aes(color = logwmean, linewidth = logwmean), alpha = 0.75) +
  scale_linewidth(range = c(1, 0.25)) +
  scale_size_area(max_size = 3) + # custom size scale
  scale_color_gradient2(low = "#FF8F33",
                        mid = "white",
                        high = "#33E4FF",
                        name = "Effective migration\nrate, log(w)") +
  geom_sf(data = border, fill = NA, color = "black") +
  coord_sf(xlim = unique(outer$long), ylim = unique(outer$lat)) +
  theme_map() +
  # Add MX state borders
  geom_polygon(data = mxstate.map, aes(x = long, y = lat, group = group), color = "gray48", fill = NA, linewidth = 0.25) +
  # Add sampling coordinates and nodes (i.e., sampling coords snapped onto lattice)
  geom_point(data = nodes %>% filter(nsamp > 0), aes(x = lon, y = lat, size = nsamp), color = "#a8a2ca") +
  scale_size(range = c(1, 6)) +
  guides(linewidth = "none") + # get rid of linewidth legend
  geom_point(data = coords, aes(x = long, y = lat), color = "black")

ggsave(paste0(here("plots"), "/forreri_feems_logwmean.pdf"), width = 10, height = 8, units = "in")


# (3) FEEMS cross-validation ----------------------------------------------

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