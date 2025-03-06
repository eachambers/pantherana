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
library(algatr)
theme_set(theme_cowplot())

## The following file:
##    (1) Imports results from ADMIXTURE (using "*.Q", "*.fam", and "*.out" files) using the `import_admix_data()` function
##    (2) Fig. S5 (top): Plots CV error from a set of K-values from ADMIXTURE results using the `cv_error()` function
##    (3) Figs. 2B, 2C, 3B, 3C, and S5 (bottom): Builds structure plots from ADMIXTURE results using the `build_str_plot()` and `retrieve_kcols()` functions
##    (4) Plot pie charts on a map
##    (5) Figs. 2B, 2C, 3B, and 3C: Plot pie charts with average ancestry values
##    (6) Identify sympatric localities and calculate distances between them
##    (7) Figs. 2C and 3C: Build inset sympatric locality maps

##    FILES REQUIRED:
##            "*.Q", "*.fam", and "*.out" files from admixture
##            *_metadata files for assemblies
##            *_order.txt for assemblies which specifies how inds are ordered in barplots

source(here("data_viz", "Pop_gen_functions.R"))
source(here("data_viz", "rana_colors.R"))


# (1) Import ADMIXTURE results --------------------------------------------

# The following will produce a list with each element being a K-value specified 
# using the `K_values` argument. This function will import Q-matrices using the 
# `*.Q` files, it will append sample IDs based on the `*.fam` files, and it will 
# determine which K-values are the highest (`max_K` column) and what those values 
# are (`max_K_val` column) for each sample.

atlmx <- import_admix_data(path = here("data", "admixture"), 
                           prefix = "ATL_MXPL_relaxed_0.25miss_ldp",
                           K_values = c(3, 4, 5, 6, 7))

centam <- import_admix_data(path = here("data", "admixture"), 
                            prefix = "new_CENTAM_relaxed_0.25miss_ldp", 
                            K_values = c(3, 4, 5, 6, 7))

pacmx <- import_admix_data(path = here("data", "admixture"), 
                           prefix = "new_PACMX_relaxed_0.25miss_ldp", 
                           K_values = c(4, 5, 6, 7, 8, 9, 10))

forreri <- import_admix_data(path = here("data", "admixture"), 
                             prefix = "forreri_0.25miss_ldp", 
                             K_values = c(2, 3, 4, 5, 6))


# (2) Examine CV error from ADMIXTURE -------------------------------------

plot_cv_error(cv_scores = atlmx$cv_scores,
              hilite = c(3:7),
              bestk = 6)
ggsave(here("plots", "atlmx_cv.pdf"), width = 4, height = 3, units = "in")

plot_cv_error(cv_scores = centam$cv_scores,
              hilite = c(3:7),
              bestk = 4)
ggsave(here("plots", "centam_cv.pdf"), width = 4, height = 3, units = "in")

plot_cv_error(cv_scores = pacmx$cv_scores,
              hilite = c(4:10),
              bestk = 9)
ggsave(here("plots", "pacmx_cv.pdf"), width = 4, height = 3, units = "in")

plot_cv_error(cv_scores = forreri$cv_scores,
              hilite = c(2:6),
              bestk = 5)
ggsave(here("plots", "forreri_cv.pdf"), width = 4, height = 3, units = "in")


# (3) Structure plots of ADMIXTURE results --------------------------------

dataset_name = "PACMX" # "forreri" / "ATL_MXPL" / "CENTAM" / "PACMX"

# Set which element of the list produced above (K-value) you want plotted
# for a given assembly:
# dat = atlmx$dat$K7
# dat = centam$dat$K4
dat = pacmx$dat$K9
# dat = forreri$dat$K5
# The `order` column is for subsequent plotting.

if (dataset_name == "ATL_MXPL") dat <- dat %>% dplyr::filter(Bioinformatics_ID != "IRL57_LCA")

K_value = ncol(dat)-5
kcols <- retrieve_kcols(K_value, dataset = "admixture", dataset_name)
indorder <- read_tsv(paste0(here("data"), "/", dataset_name, "_order.txt"), col_names = TRUE)

if (dataset_name == "ATL_MXPL") indorder <- indorder %>% dplyr::filter(Bioinformatics_ID != "IRL57_LCA")

build_str_plot(dat = dat, 
               K_value, 
               order = "manual",
               man_order = indorder$original_order,
               kcols = kcols, 
               xaxis_labels = TRUE,
               write_output = FALSE, 
               export_plot = FALSE,
               metadata_path = paste0(here("data"), "/", dataset_name, "_metadata.txt"),
               output_path = here("plots")) # export 12x7
ggsave(paste0(here("plots"), "/", dataset_name, "_strplot.pdf"), width = 12, height = 7, units = "in")

# (4) Plot pie charts on map ----------------------------------------------

# Import data -------------------------------------------------------------

# dat = atlmx$dat$K6
# dat = centam$dat$K4
dat = pacmx$dat$K9
# dat = forreri$dat$K5

# Also provide the dataset name
dataset_name = "PACMX" # "forreri" "PACMX" "CENTAM" "ATL_MXPL"

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


# Plot PACMX with ranges --------------------------------------------------

# TODO switch below to using function
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


# (6) Identify sympatric localities ---------------------------------------

dat_centam = centam$dat$K4
dat_pacmx = pacmx$dat$K9

# Pull in and bind with locality info for each
dat <- dat_pacmx
dataset_name = "PACMX"

kcols <- retrieve_kcols(ncol(dat)-5, dataset = "admixture", dataset_name)
metadata <- read_tsv(paste0(here("data"), "/", dataset_name, "_metadata.txt"), col_names = TRUE)
final <- left_join(dat, metadata, by = "Bioinformatics_ID")

yava_forr_naya_inds <- c("V2790_PAC", "V2791_PAC", "V2792_PAC", "V2793_PAC", "V2794_PAC", "V2798_PAC", "V2800_PAC", "V2799_PAC")
yava_forr_jali_inds <- c("T2020_Aten_PAC", "T2022_Aten_PAC", "T3437_PAC")
omil_forr_inds <- c("V2593_pap_PAC", "V2655_pap_PAC", "V2606_pap_PAC", "V2602_pap_PAC", "T14192_papa_PAC", "V2598_pap_PAC")
spec_omil_inds <- c("TF8695_Rspe_CMX", "TF8696_Rspe_CMX", "TF8697_Rspe_CMX", "TF8698_Rspe_CMX", "TF8699_Rspe_CMX",
                    "TF7120_Rfor_PAC", "TF7121_Rfor_PAC", "TF7122_Rfor_PAC", "TF7123_Rfor_PAC",
                    "TF7124_Rfor_PAC", "TF7125_Rfor_PAC", "TF7126_Rfor_PAC", "TF7127_Rfor_PAC",
                    "TF7129_Rfor_PAC", "TF7128_Rfor_PAC")
forr_macro_inds <- c("TF8646_Rspe_CMX", "TF8647_Rspe_CMX", "TF8648_Rspe_CMX", "TF8649_Rspe_CMX",
                     "TF8650_Rspe_CMX", "TF8651_Rspe_CMX", "TF8652_Rspe_CMX", "TF8653_Rspe_CMX",
                     "TF8635_Rspe_CMX", "TF8636_Rspe_CMX", "TF8637_Rspe_CMX", "TF8638_Rspe_CMX",
                     "TF8639_Rspe_CMX", "TF8640_Rspe_CMX", "TF8641_Rspe_CMX", "TF8642_Rspe_CMX", 
                     "TF8643_Rspe_CMX", "TF8644_Rspe_CMX", "TF8645_Rspe_CMX")

forr_macro <- final %>% 
  filter(Bioinformatics_ID %in% forr_macro_inds)
geodists <- geo_dist(forr_macro %>% dplyr::select(lat, long) %>% rename(x = long, y = lat))
colnames(geodists) <- forr_macro$Bioinformatics_ID
geodists <- geodists %>% as.data.frame() %>% mutate(Bioinformatics_ID = forr_macro$Bioinformatics_ID,
                                                    max_K = forr_macro$max_K)
min(forr_macro$max_K_val)

# Counted as sympatric:
# omil_forr: 1.3337578, 0.8037569, 6.0612121, 9.4624518, 6.061212, 6.583978, 5.426628 at min of 0.998966 ancestry
# forr_macro: 0, 6.887742km @ min of 0.998093
# yava_forr_jali: 8.765118km between X6 and X7 @ 0.99 ancestry

## Not counted as sympatric:
# yava_forr_naya: 15.83702 & 16.35864 km dists between X1 and X2 @ 0.9831 ancestry
# omil_forr: 15.8070380, 13.337578, 15.80704, 15.96336, 22.74879
# at min of 0.998966 ancestry
# spec_omil: 19.07065km @ 0.999906 min


# (7) Build inset sympatric locality structure plots ----------------------

# TODO check below stuff is consistent with above

dataset_name = "PACMX" # "PACMX" / "CENTAM"

# if (dataset_name == "ATL_MXPL") dat <- dat %>% dplyr::filter(Bioinformatics_ID != "IRL57_LCA")
K_value = ncol(dat)-5
kcols <- retrieve_kcols(K_value, dataset = "admixture", dataset_name)
metadata <- read_tsv(paste0(here("data"), "/", dataset_name, "_metadata.txt"), col_names = TRUE)
sites <- read_tsv(paste0(here("data"), "/", dataset_name, "_sites.txt"), col_names = TRUE)

final <- left_join(dat, metadata, by = "Bioinformatics_ID")
final <- left_join(final, sites)

# Sympatric sites are those that have 2 max_Ks
check <-
  final %>%
  group_by(site_name, max_K) %>%
  count()

# Extract sympatric localities only; N.B.: RioelAmate and RioPapagayo are identical
if (dataset_name == "PACMX") symps <- final %>% dplyr::filter(site_name == "PuenteTehuantepec" | site_name == "RiolasTejas" | 
                                                                site_name == "RioelAmate" | site_name == "RioPapagayo") %>% dplyr::arrange(site_name, max_K) %>% rownames_to_column(var = "new_order")
if (dataset_name == "CENTAM") symps <- final %>% dplyr::filter(site_name == "CerroSPLaLoma") %>% dplyr::arrange(site_name, max_K) %>% rownames_to_column(var = "new_order")

symps$new_order <- factor(symps$new_order, levels = symps$new_order)

symps <-
  symps %>% 
  tidyr::pivot_longer(names_to = "cluster", values_to = "proportion", 
                      -c(Bioinformatics_ID, K_val, order, new_order, max_K, max_K_val, lat, long, site_name, dataset, state))
p <- symps %>% 
  ggplot2::ggplot(aes(x = new_order, y = proportion, fill = cluster)) +
  ggplot2::geom_bar(stat = "identity") + 
  panel_border() + 
  scale_y_continuous(expand = c(0,0)) +
  # scale_x_discrete(expand = c(0,0)) +
  scale_fill_manual(values = kcols) +
  theme_cowplot() %+replace% theme(axis.line = element_line(colour = "black"),
                                   axis.text.x = element_blank(),
                                   axis.text.y = element_blank(),
                                   axis.ticks = element_blank(),
                                   axis.title.x = element_blank(),
                                   axis.title.y = element_blank(),
                                   legend.position = "none",
                                   panel.border = element_rect(fill = NA, colour = "black", linetype = "solid", linewidth = 1.5),
                                   strip.text.y = element_text(size = 30, face = "bold"),
                                   strip.background = element_rect(colour = "white", fill = "white"),
                                   panel.spacing = unit(-0.1, "lines"))

ids <- as.list(unique(symps$Bioinformatics_ID))
p <- p + 
  scale_x_discrete(expand = c(0,0), labels = ids) +
  theme(axis.text.x = element_text(angle = 90, size = 5)) # 10x4 PACMX / 2.3x3 CENTAM


# Guerrero samples for PACMX ----------------------------------------------

data("mxstate.map")
mxstate.map$group <- as.numeric(mxstate.map$group)
gueoax <- mxstate.map %>% dplyr::filter(id == 12 | id == 20)

dataset_name = "PACMX"
K_value = ncol(dat)-5
kcols <- retrieve_kcols(K_value, dataset = "admixture", dataset_name)
metadata <- read_tsv(paste0(here("data"), "/", dataset_name, "_metadata.txt"), col_names = TRUE)
sites <- read_tsv(paste0(here("data"), "/", dataset_name, "_sites.txt"), col_names = TRUE)

final <- left_join(dat, metadata, by = "Bioinformatics_ID")
gro <- left_join(final, sites) %>% 
  dplyr::filter(state == "Guerrero" | state == "Oaxaca") %>% 
  dplyr::arrange(site_name, max_K) %>% 
  rownames_to_column(var = "new_order")
gro$new_order <- factor(gro$new_order, levels = gro$new_order)

# Get observation numbers (optional)
obs <-
  gro %>% 
  # filter(state == "Oaxaca") %>% 
  dplyr::select(site_name) %>% 
  group_by(site_name) %>% 
  summarize(no_inds = n())

gro <-
  gro %>% 
  tidyr::pivot_longer(names_to = "cluster", values_to = "proportion", 
                      -c(Bioinformatics_ID, K_val, order, new_order, max_K, max_K_val, lat, long, site_name, dataset, state))

# Plot structure plot, ordered by site
gro %>% 
  filter(state == "Oaxaca") %>%
  ggplot2::ggplot(aes(x = new_order, y = proportion, fill = cluster)) +
  ggplot2::geom_bar(stat = "identity") + 
  panel_border() + 
  scale_y_continuous(expand = c(0,0)) +
  # scale_x_discrete(expand = c(0,0)) +
  scale_fill_manual(values = kcols) +
  theme_cowplot() %+replace% theme(axis.line = element_line(colour = "black"),
                                   axis.text.x = element_blank(),
                                   axis.text.y = element_blank(),
                                   axis.ticks = element_blank(),
                                   axis.title.x = element_blank(),
                                   axis.title.y = element_blank(),
                                   legend.position = "none",
                                   panel.border = element_rect(fill = NA, colour = "black", linetype = "solid", linewidth = 1.5),
                                   strip.text.y = element_text(size = 30, face = "bold"),
                                   strip.background = element_rect(colour = "white", fill = "white"),
                                   panel.spacing = unit(-0.1, "lines")) # export 12x8

# Make the map with site localities
ggplot() +
  geom_polygon(data = gueoax, aes(x = long, y = lat, group = group), color = "white", fill = "#e6e6e6", linewidth = 0.25) + # color = "#969696" if CENTAM
  theme_map() +
  coord_fixed() +
  geom_point(data = gro, aes(long, lat), size = 3) +
  geom_text(data = gro, aes(x = long, y = lat, label = site_name), color = "black", size = 3) # export 10x8


