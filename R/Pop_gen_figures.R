library(tidyverse)
library(here)
library(ggplot2)
library(cowplot)
library(maps)
library(mxmaps) # devtools::install_github("diegovalle/mxmaps")
library(ggrepel)
library(repel) # devtools::install_github("hms-dbmi/repel")
library(usmap)
library(sf)
library(sfheaders)
library(maptools)

theme_set(theme_cowplot())
setwd("~/Box Sync/Rana project/ddRADseq/ALL_RANA/Separate_assemblies/ADMIXTURE/")

## The following file:
##    (1) Plots CV error from a set of K-values from ADMIXTURE results using the `cv_error()` function
##    (2) Imports results from ADMIXTURE (using `*.Q` and `*.fam` files) using the `import_admix_data()` function
##    (3) Builds structure plots from ADMIXTURE results using the `build_str_plot()` and `retrieve_kcols()` functions
##    (4) Plot pie charts on a map

##    FILES REQUIRED:
##            CV_error.txt

source(here("R", "Pop_gen.R"))
source(here("R", "rana_colors.R"))


# (1) Examine CV error from ADMIXTURE -------------------------------------

# Specify the dataset name
dataset_name = "forreri" # "ATL_MXPL" "forreri" "PACMX" "CENTAM"

# Build all together
cv_error(dat = here("data", "CV_error.txt"), faceted = TRUE)

# Alternatively, plot for one of the assemblies:
cv <- read_tsv(here("data", "CV_error.txt"), col_names = TRUE) %>% 
  filter(dataset == dataset_name)
maxy = max(cv$CV_error)

if (dataset_name == "ATL_MXPL" | dataset_name == "CENTAM") {
  range = c(3:15)
  hilight = c(3:7)
  if (dataset_name == "ATL_MXPL") xint = 6
  if (dataset_name == "CENTAM") xint = 4
}
if (dataset_name == "forreri") {
  range = c(1:10)
  hilight = c(2:6)
  xint = 5
}
if (dataset_name == "PACMX") {
  range = c(1:15)
  hilight = c(4:10)
  xint = 8
}

cv %>% 
  ggplot(aes(x = K_value, y = CV_error)) +
  scale_x_continuous(breaks = range) +
  annotate("rect", xmin = min(hilight), xmax = max(hilight), ymin = 0, ymax = maxy, # 3-7 CENTAM & ATL_MXPL; 4-10 PACMX; 2-6 forreri
           alpha = 0.5, fill = "#88a1c6") +
  geom_vline(xintercept = xint, color = "red") + # 4 CENTAM; 8 PACMX; 6 ATL_MXPL; 5 forreri
  geom_point() +
  geom_line() +
  xlab("K value") +
  ylab("Cross-validation error") # export 4x3


# (2) Import ADMIXTURE results --------------------------------------------

# The following will produce a list with each element being a K-value specified 
# using the `K_values` argument. This function will import Q-matrices using the 
# `*.Q` files, it will append sample IDs based on the `*.fam` files, and it will 
# determine which K-values are the highest (`max_K` column) and what those values 
# are (`max_K_val` column) for each sample. The `order` column is for subsequent plotting.

atlmx <- import_admix_data(path = ".", file = "ATL_MXPL/ATL_MXPL_relaxed_0.25miss_ldp", K_values = c(3, 4, 5, 6, 7))
atlmx <- list(K3 = atlmx[[1]], K4 = atlmx[[2]], K5 = atlmx[[3]], K6 = atlmx[[4]], K7 = atlmx[[5]])

centam <- import_admix_data(path = ".", file = "CENTAM/new_CENTAM_relaxed_0.25miss_ldp", K_values = c(3, 4, 5, 6, 7))
centam <- list(K3 = centam[[1]], K4 = centam[[2]], K5 = centam[[3]], K6 = centam[[4]], K7 = centam[[5]])

pacmx <- import_admix_data(path = ".", file = "PACMX/new_PACMX_relaxed_0.25miss_ldp", K_values = c(4, 5, 6, 7, 8, 9, 10))
pacmx <- list(K4 = pacmx[[1]], K5 = pacmx[[2]], K6 = pacmx[[3]], K7 = pacmx[[4]],
              K8 = pacmx[[5]], K9 = pacmx[[6]], K10 = pacmx[[7]])

forreri <- import_admix_data(path = ".", file = "forreri/forreri_0.25miss_ldp", K_values = c(2, 3, 4, 5, 6))
forreri <- list(K2 = forreri[[1]], K3 = forreri[[2]], K4 = forreri[[3]], K5 = forreri[[4]], K6 = forreri[[5]])


# (3) Structure plots of ADMIXTURE results --------------------------------

# Set which element of the list produced above (K-value) you want plotted
# for a given assembly:
# dat = atlmx$K7
# dat = centam$K4
dat = pacmx$K6
# dat = forreri$K5

if (dataset_name == "ATL_MXPL") dat <- dat %>% dplyr::filter(Bioinformatics_ID != "IRL57_LCA")

K_value = ncol(dat)-5
kcols <- retrieve_kcols(K_value, dataset = "admixture", dataset_name)
indorder <- read_tsv(paste0(here("data"), "/", dataset_name, "_order.txt"), col_names = TRUE)

p <- build_str_plot(dat = dat, 
                    K_value, 
                    order = "manual",
                    man_order = indorder$original_order,
                    kcols = kcols, 
                    xaxis_labels = TRUE,
                    write_output = TRUE, 
                    export_plot = TRUE,
                    metadata_path = paste0(here("data"), "/", dataset_name, "_metadata.txt"),
                    output_name = dataset_name) # export 12x7


# (4) Plot pie charts on map ----------------------------------------------

# Import data -------------------------------------------------------------

# dat = atlmx$K6
# dat = centam$K4
# dat = pacmx$K6
dat = forreri$K5

# Also provide the dataset name
dataset_name = "forreri" # "forreri" "PACMX" "CENTAM" "ATL_MXPL"

# Get color values
kcols <- retrieve_kcols(K_value = ncol(dat)-5, dataset = "admixture", dataset_name)
# Import metadata and join with ADMIXTURE results
metadata <- read_tsv(paste0(here("data"), "/", dataset_name, "_metadata.txt"), col_names = TRUE)
final <- left_join(dat, metadata, by = "Bioinformatics_ID")


# Retrieve data for background map ----------------------------------------

world <- map_data("world")
centamer <- filter(world, region == "Belize" | region == "Guatemala" | region == "Honduras" | 
                     region == "Nicaragua" | region == "Costa Rica" | region == "El Salvador" | 
                     region == "Panama")

if (dataset_name == "ATL_MXPL" | dataset_name == "forreri" | dataset_name == "PACMX") {
  data("mxstate.map")
  mxstate.map$group <- as.numeric(mxstate.map$group)
  map_data <-
    bind_rows(centamer, mxstate.map)
  
  if (dataset_name == "PACMX") {
    # Include some U.S. states
    include <- c("california", "arizona", "new mexico", "nevada", "texas")
    extended <- c("california", "arizona", "new mexico", "nevada", "texas",
                  "louisiana", "arkansas", "oklahoma", "mississippi", "tennessee",
                  "alabama", "georgia", "florida", "missouri", "south carolina", "north carolina")
    # Get polygons for mapping
    states <- map_data("state") %>% 
      filter(region %in% include)
    map_data <-
      bind_rows(map_data, states)
  }
}
if (dataset_name == "CENTAM") map_data <- centamer


# Add jittering to coords -------------------------------------------------

# The `repel_coords()` function uses the `repel_text()` function in the repel 
# package and will append two new columns to the data which correspond to 
# jittered/repelled lat and long. These will be used for the final locations 
# of the pie charts on the map (and will be the end points of lines drawn from 
# the sampling location to the pies). Within this function, you can also choose 
# to not repel certain points (i.e., the "new_lat" and "new_long" values will 
# stay the same) using the `ignore_repel` argument and providing a list of samples 
# to ignore with the `ignore_repel_list` argument. This function requires that your 
# coordinate columns are named "lat" and "long"; it will also remove NAs.
final <- repel_coords(dat = final, ignore_repel = FALSE)


# Make pies separately ----------------------------------------------------

# Draw the individual pies using the `make_pies()` and `draw_pies()` functions
kvalmaxcol <- paste0("X", ncol(dat)-5)

pies_to_add <- final %>% 
  # Nest a tibble in the data column with the admixture values
  nest(data = c(X1:kvalmaxcol)) %>% 
  # Make a new column with a pie chart for the data in each row
  mutate(pies = purrr::map(data, make_pies, kvalmaxcol = kvalmaxcol, kcols, alpha = 1)) %>% 
  # draw_pies uses cowplot::draw_plot which builds inset plots; the following adds another col with drawn pies
  mutate(drawn_pies = pmap(tibble(plot = pies, lat = new_lat, long = new_long), draw_pies, height = 0.75, width = 0.75))


# Make base map with points and lines -------------------------------------

# The base map will include the background map as well as the sampling localities 
# and lines that will eventually connect to the pies.
base_map <-
  map_data %>% 
  ggplot() +
  geom_polygon(aes(x = long, y = lat, group = group), color = "white", fill = "#e6e6e6", linewidth = 0.25) + # color = "#969696" if CENTAM
  theme_map() +
  geom_segment(data = final, aes(x = long, y = lat, xend = new_long, yend = new_lat), linewidth = 0.75) +
  geom_point(data = final, aes(long, lat), size = 2) +
  coord_fixed()


# Plot the pies -----------------------------------------------------------

# Finally, plot the pies you've drawn on the map you've drawn
list(base_map, pies_to_add$drawn_pies) %>% 
  reduce(.f = `+`) # export 12x8



# -------------------------------------------------------------------------
# (5) Plot pie charts with average site values ----------------------------

# To plot a locality-based average of the ancestry coefficients, summarize data:
sites <-
  final %>% 
  pivot_longer(cols = X1:kvalmaxcol, names_to = "K_value", values_to = "proportion") %>% 
  group_by(lat, long, K_value) %>% 
  summarize(site_prop = mean(proportion)) %>% 
  pivot_wider(names_from = K_value, values_from = site_prop)

# Add repelled coordinates:
sites <- repel_coords(dat = sites, ignore_repel = FALSE)

# Now get pies:
pies_to_add <- sites %>% 
  # Nest a tibble in the data column with the admixture values
  nest(data = c(X1:kvalmaxcol)) %>% 
  # Make a new column with a pie chart for the data in each row
  mutate(pies = purrr::map(data, make_pies, kvalmaxcol = kvalmaxcol, kcols, alpha = 1)) %>% 
  # draw_pies uses cowplot::draw_plot which builds inset plots; the following adds another col with drawn pies
  mutate(drawn_pies = pmap(tibble(plot = pies, lat = new_lat, long = new_long), draw_pies, height = 0.75, width = 0.75))

# Build base map:
base_map_sites <-
  map_data %>%
  ggplot() +
  geom_polygon(aes(x = long, y = lat, group = group), color = "white", fill = "#e6e6e6", linewidth = 0.25) + # color = "#969696" if CENTAM
  theme_map() +
  geom_segment(data = sites, aes(x = long, y = lat, xend = new_long, yend = new_lat), linewidth = 0.75) +
  geom_point(data = sites, aes(long, lat), size = 2) +
  coord_fixed()

# Plot pies on base map:
list(base_map_sites, pies_to_add$drawn_pies) %>% 
  reduce(.f = `+`) # export 12x8


# Add observation numbers on pies -----------------------------------------

obs <-
  final %>% 
  dplyr::select(lat, long) %>% 
  group_by(lat, long) %>% 
  summarize(no_inds = n())

map_data %>%
  ggplot() +
  geom_polygon(aes(x = long, y = lat, group = group), color = "white", fill = "#e6e6e6", linewidth = 0.25) + # color = "#969696" if CENTAM
  theme_map() +
  geom_point(data = sites, aes(long, lat), size = 2) +
  coord_fixed() +
  geom_text(data = obs, aes(x = long, y = lat, label = no_inds), color = "white", size = 1.5)


# -------------------------------------------------------------------------
# Plot forreri with range + type localities -------------------------------

shp <- maptools::readShapePoly("forreri_58599/species_58599.shp")
forrrange <- fortify(shp)
sp6 <- maptools::readShapePoly("sp6/species d (sp.6).shp")
sp6range <- fortify(sp6)
tls <- read_tsv(here("data", "forreri_typelocalities.txt"), col_names = TRUE)

# Build base map:
base_map_sites <-
  map_data %>%
  ggplot() +
  geom_polygon(aes(x = long, y = lat, group = group), color = "white", fill = "#e6e6e6", linewidth = 0.25) + # color = "#969696" if CENTAM
  theme_map() +
  geom_polygon(data = forrrange, aes(long, lat, group = group), fill = "#73806d", alpha = 0.5) +
  geom_polygon(data = sp6range, aes(long, lat, group = group), fill = "#73806d", alpha = 0.5) +
  geom_segment(data = sites, aes(x = long, y = lat, xend = new_long, yend = new_lat), linewidth = 0.75) +
  geom_point(data = sites, aes(long, lat), size = 2) +
  geom_point(data = tls, aes(Longitude, Latitude), size = 2, color = "red")
  coord_fixed()

# Plot pies on base map:
list(base_map_sites, pies_to_add$drawn_pies) %>% 
  reduce(.f = `+`) # export 12x8


# -------------------------------------------------------------------------
# Plot PACMX with ranges --------------------------------------------------

yava <- maptools::readShapePoly("pacmx_ranges/species_19181.shp")
yava <- fortify(yava)
magn <- maptools::readShapePoly("pacmx_ranges/species_58656.shp")
magn <- fortify(magn)
forr <- maptools::readShapePoly("pacmx_ranges/species_58599.shp")
forr <- fortify(forr)
sp6 <- maptools::readShapePoly("pacmx_ranges/species d (sp.6).shp")
sp6range <- fortify(sp6)
omil <- maptools::readShapePoly("pacmx_ranges/species_58687.shp")
omil <- fortify(omil)
spec <- maptools::readShapePoly("pacmx_ranges/species_58722.shp")
spec <- fortify(spec)
tayl <- maptools::readShapePoly("pacmx_ranges/species_58732.shp")
tayl <- fortify(tayl)
macro <- maptools::readShapePoly("pacmx_ranges/data_0.shp")
macro <- fortify(macro)
brown <- maptools::readShapePoly("pacmx_ranges/browndata_0.shp")
brown <- fortify(brown)

tls <- read_tsv(here("data", "PACMX_typelocalities.txt"), col_names = TRUE)

ggplot() +
  geom_polygon(data = map_data, aes(x = long, y = lat, group = group), color = "white", fill = "#e6e6e6", linewidth = 0.25) + # color = "#969696" if CENTAM
  theme_map() +
  geom_polygon(data = forr, aes(long, lat, group = group), fill = "#73806d", alpha = 0.75) +
  geom_polygon(data = sp6range, aes(long, lat, group = group), fill = "#73806d", alpha = 0.75) +
  geom_polygon(data = yava, aes(long, lat, group = group), fill = "#054051", alpha = 0.75) +
  geom_polygon(data = magn, aes(long, lat, group = group), fill = "#f7cd5e", alpha = 0.75) +
  geom_polygon(data = omil, aes(long, lat, group = group), fill = "#ba94a6", alpha = 0.75) +
  geom_polygon(data = spec, aes(long, lat, group = group), fill = "#80a4bc", alpha = 0.75) +
  geom_polygon(data = macro, aes(long, lat, group = group), fill = "#e0895a", alpha = 0.75) +  
  geom_polygon(data = brown, aes(long, lat, group = group), fill = "#e0895a", alpha = 0.75) +  
  geom_point(data = tls, aes(Longitude, Latitude), size = 2, color = "red") +
  coord_fixed() +
  geom_polygon(data = map_data, aes(x = long, y = lat, group = group), color = "white", fill = NA, linewidth = 0.25)
