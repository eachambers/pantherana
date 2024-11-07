library(algatr)
library(terra)
library(raster)
library(tidyverse)
library(geodata)
library(corrplot)
library(vegan)
library(gdistance)
library(topoDistance)
library(rmapshaper)
library(viridis)
library(vcfR)
library(here)
library(cowplot)
library(SNPRelate)
library(tidyterra)
theme_set(theme_cowplot())

source(here("R", "Land_gen_functions.R"))

## This code performs landscape genomic analyses related to testing for IBD in the R. forreri complex:
##     (1) Calculate site-based genetic distances
##     (2) Calculate site-based geographic distances
##     (3) Run a Mantel test
##     (4) Build Fig. S6: Mantel test results
##     (5) Run MMRR
##     (6) Build Fig. 5A: MMRR results
##     (7) Run GDM
##     (8) Build Fig. 5B: GDM results

##    FILES REQUIRED:
##          LD-pruned vcf (forreri_0.25miss_ldp.recode.vcf)
##          Sample IDs and coords (forreri_metadata.txt)
##          Unique sampling sites (forreri_sites.txt)
##          files in PC_layers/ - Environmental PCs (from raster PCA in Env_data.R)


# (1) Calculate site-based genetic distances ------------------------------

# Import site data
sites <- read_tsv(here("data", "forreri_sites.txt"), col_names = TRUE)

# Import locality information for each individual
samps <- read_tsv(here("data", "forreri_metadata.txt"),
                         col_names = TRUE) %>% 
  dplyr::rename(x = long,
                y = lat) %>% 
  na.omit() # V2797 does not have coords

# Import genetic data; 65000 variants
vcf <- read.vcfR(file = here("data", "forreri_0.25miss_ldp.recode.vcf"))

# Subset vcf to match coords
index <- colnames(vcf@gt) %in% samps$Bioinformatics_ID
# First col is format col
index[1] <- TRUE
# Subset vcf to match coords
vcf <- vcf[, index]
colnames(vcf@gt) # should be 104 (because FORMAT is a col so 103 inds)

# For several localities, we have genetic data for >1 individual.
# To not bias our landscape genomic analyses, we need to calculate
# mean allele frequencies for each locality. To do so, we'll first 
# convert the vcf to a dosage matrix.
dos <- algatr::vcf_to_dosage(vcf) # 951 loci omitted from genlight object because >2 alleles; dosage has coding 0, 1, 2
freqs <- dos/2 # coding is now 0, 0.5, 1 as allele frequencies

# Check ordering consistency between vcf and metadata
all(colnames(vcf@gt[,-1]) == samps$Bioinformatics_ID) # should be all TRUE, if not run the following lines
vcford <- colnames(vcf@gt[,-1])
# Combine lat/longs with site names, arranging so it's ordered the same as vcf - REQUIRED!!!
sites <- left_join(samps, sites) %>% 
  arrange(Bioinformatics_ID, vcford)
all(colnames(vcf@gt[,-1]) == sites$Bioinformatics_ID) # should be TRUE

# Now, calculate mean site-based allele frequencies for each locus
# !BE AWARE: the following changes the ordering of the sites; ensure these are consistent
# between gendist and geodist matrices
site_freqs <- aggregate(freqs, list(sites$site_name), FUN = mean, na.rm = TRUE) # 38 obs of 64050 vars
write_tsv(site_freqs, here("data", "site_freqs.txt"))

# Finally, calculate genetic distances and export:
gendist <- as.matrix(dist(site_freqs, method = "euclidean"))
unique_sites <- site_freqs[1]
colnames(gendist) <- unique_sites$Group.1
write_tsv(as.data.frame(gendist), file = here("data", "forreri_gendist_ldp.txt"), col_names = TRUE)


# (2) Calculate site-based geographic distances ---------------------------

# There are two sets of inds from Ticuizitan; they're all from the same site so 
# should be considered the same. For ease for geodist calcs we need to remove them:
tic_inds <- c("V2717_PAC", "V2718_PAC", "V2719_PAC", "V2721_PAC", 
              "V2722_PAC", "V2724_PAC", "V2728_PAC") # y = 19.183056

subsites <- sites %>% 
  filter(!Bioinformatics_ID %in% tic_inds) %>% 
  dplyr::select(-Bioinformatics_ID) %>% 
  filter(x != "-82.97517" & x != "-82.97519") %>% 
  distinct() %>% 
  arrange(site_name)

unique_sites <- as.data.frame(site_freqs[1]) %>% 
  dplyr::rename(site_name = "Group.1")
subsites <- left_join(unique_sites, subsites)

# Check that ordering is correct:
all(unique_sites$site_name == subsites$site_name)

# Should be 38 sites
# Calculate geographic distances and export
geodist <- geo_dist(subsites %>% dplyr::select(x, y)) # you will get an NA message if you have more than 2 cols
colnames(geodist) <- subsites$site_name
write_tsv(as.data.frame(geodist), file = "forreri_geodist_ldp.txt", col_names = TRUE)


# (3) Run Mantel test -----------------------------------------------------

geodist[upper.tri(geodist, diag = FALSE)] <- NA
gendist[upper.tri(gendist, diag = FALSE)] <- NA

melt_geo <-
  geodist %>%
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "comparison") %>%
  tidyr::pivot_longer(cols = -(comparison)) %>%
  na.omit() %>%
  dplyr::filter(comparison != name) %>%
  dplyr::rename(geodist := value)

melt_gen <- 
  gendist %>%
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "comparison") %>%
  tidyr::pivot_longer(cols = -(comparison)) %>%
  na.omit() %>%
  dplyr::filter(comparison != name) %>%
  dplyr::rename(gendist := value)

dat <- full_join(melt_gen, melt_geo) %>% 
  mutate(geodist_km = geodist*100) # joining by comparison, name; should be 741 obs (38*37/2 / 2 + 38)

##### Run Mantel test
# ibd <- vegan::mantel(gendist, geodist, method = "pearson")
# ibd # Mantel stat r = 0.6686, sig=0.001
# spear <- vegan::mantel(gendist, geodist, method = "spearman")
# spear # Mantel stat r = 0.7411, sig=0.001
IBD <- vegan::mantel(gendist, log(geodist), method = "pearson")
IBD # Mantel stat r = 0.7699, sig=0.001


# (4) Visualize Mantel test -----------------------------------------------

# First, let's just see what this looks like plotted
dat %>% 
  ggplot(aes(x = geodist_km, y = gendist)) +
  geom_point() +
  xlab("Geographic distance (km)") +
  ylab("Genetic distance")

# Generalized additive model, density fill
dat %>% 
  filter(gendist > 0) %>% 
  ggplot(aes(x = geodist_km, y = gendist)) +
  stat_density_2d(aes(fill = after_stat(density), alpha = ..density..), geom = "raster", contour = FALSE, na.rm = TRUE) +
  scale_fill_distiller(palette = "Spectral", direction = -1) +
  scale_alpha_continuous(range = c(0, 10), guide = guide_none()) +
  theme(legend.position = "none") +
  geom_point(alpha = 0.5) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("Geographic distance (km)") +
  ylab("Genetic distance") +
  ggtitle("Site-based comparisons (generalized additive model)") +
  # geom_smooth(se = FALSE, method = "glm")
  geom_smooth(se = FALSE, method = "gam", formula = y ~ s(log(x)), color = "red", linewidth = 1) # export 8x6

dat %>% 
  filter(gendist > 0) %>% 
  filter(geodist > 0) %>% 
  ggplot(aes(x = log(geodist_km), y = gendist)) +
  stat_density_2d(aes(fill = ..density.., alpha = ..density..), geom = "raster", contour = FALSE, na.rm = TRUE) +
  scale_fill_distiller(palette = "Spectral", direction = -1) +
  scale_alpha_continuous(range = c(0, 10), guide = guide_none()) +
  theme(legend.position = "none") +
  geom_point(alpha = 0.5) +
  geom_smooth(method = lm, formula = y~x-1, color = "darkgrey") +
  xlab("log(Geographic distance (km))") +
  ylab("Genetic distance") +
  ggtitle("Site-based comparisons")

ggsave(here("plots", "Mantel_test.pdf"), width = 8, height = 6, units = "in")


# (5) Run MMRR ------------------------------------------------------------

# Load env data
forr_env <- raster::stack(list.files(here("data", "PC_layers/"), full.names = TRUE))
forr_env <- raster::readAll(forr_env)

# Extract enviro vars (3 enviro PCs for each site coordinate)
env <- raster::extract(forr_env, subsites %>% 
                         dplyr::select(x, y))

# Calculate environmental distances
X <- env_dist(env)

# Add geographic distance to X from the Mantel test
X[["geodist"]] <- geodist

Y <- as.matrix(gendist)

results_best <- mmrr_run(Y, X, nperm = 99, stdz = TRUE, model = "best")

# Look at plots of results
mmrr_plot_fitted(mod = results_best$mod, Y, X, stdz = TRUE)
mmrr_plot_cov(X, stdz = TRUE)
mmrr_plot_vars(Y, X, stdz = TRUE)

# Get summary stats
mmrr_table(results_best, digits = 2, summary_stats = TRUE)


# (6) Visualize MMRR results ----------------------------------------------

# Build MMRR plot: Fig. 5A
stdz = TRUE
# Unfold X and Y
y <- unfold(Y, scale = stdz)
dfX <- purrr::map_dfc(X, unfold, scale = stdz) %>% purrr::map_dfc(as.numeric)

# Make single variable dataframe
df <- dfX %>%
  dplyr::mutate(Y = y) %>%
  tidyr::gather("var", "X", -Y)

# Build plot for significant vars only
df %>% 
  filter(var == "PC3" | var == "geodist") %>% 
  ggplot(aes(x = X, y = Y)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", formula = y ~ x, color = "#428b9b") +
  facet_wrap(~var, scales = "free", nrow = 1) +
  xlab("Variable distance") +
  ylab("Genetic distance") +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold", hjust = 0),
        axis.title = element_text(size = 14, face = "bold"))
ggsave(here("plots", "MMRR_results.pdf"), width = 9, height = 3.5, units = "in")


# (7) Run GDM -------------------------------------------------------------

# Extract vars, same as MMRR
env <- raster::extract(forr_env, subsites %>% 
                         dplyr::select(x, y))

gdm_best <- gdm_run(gendist, 
                    coords = subsites %>% 
                      dplyr::select(x, y), 
                    env = env, 
                    model = "best", 
                    scale_gendist = TRUE,
                    nperm = 500, 
                    sig = 0.05)

# Look at p-values
gdm_best$pvalues
gdm_best$varimp

# Plot the results
gdm_plot_diss(gdm_best$model)
gdm_plot_isplines(gdm_best$model)
gdm_table(gdm_best)


# (8) Visualize GDM results -----------------------------------------------

# Plot the GDM map
gdm_map(gdm_best$model, 
        forr_env, 
        subsites %>% 
          dplyr::select(x, y)) # export 10x8

# Extract the GDM map from the GDM model object
map <- build_gdm_map(gdm_best$model, 
                     forr_env, 
                     subsites %>% 
                       dplyr::select(x, y), 
                     plot_vars = FALSE)

# Now, use `extrap_mask()` to do masking based on SD, a buffer, or a range
map_mask_sd <- extrap_mask(subsites %>% 
                             dplyr::select(x, y), 
                           map$pcaRastRGB, 
                           method = "sd",
                           nsd = 2)

# Plot the resulting masked map
ggplot() +
  geom_spatraster_rgb(data = map$pcaRastRGB) +
  geom_spatraster(data = map_mask_sd, aes(fill = sum), alpha = 0.6) +
  scale_fill_gradientn(colors = "white", na.value = NA) +
  geom_point(data = subsites, aes(x = x, y = y), size = 4, colour = "black", pch = 21) +
  theme_map() +
  theme(legend.position = "none")
ggsave(here("plots", "forreri_GDM.pdf"), width = 10, height = 8, units = "in")

# Plot vector loadings for sig vars
gdm_plot_vars(map$pcaSamp, 
              map$pcaRast, 
              map$pcaRastRGB,
              coords = subsites %>% 
                dplyr::select(x, y), 
              display_axes = TRUE, 
              x = "PC1", 
              y = "PC2")
ggsave(here("plots", "forreri_GDM_loadings.pdf"), width = 4, height = 4, units = "in")

