library(devtools)
library(conStruct) # install_github("gbradburd/conStruct", ref="covariance_fix")
library(tidyverse)
library(fields)
library(here)
library(cowplot)
library(algatr)
theme_set(theme_cowplot())
setwd("~/Box Sync/Rana project/ddRADseq/ALL_RANA/Landgen")

## The following file:
##    (1) Processes input data by removing loci and inds with high missing data
##    (2) Converts a Structure file (.ustr) to conStruct file (allele freq data)
##    (3) Calculates geographic distances between samples
##    (4) Runs conStruct
##    (5) Runs cross-validation analysis comparing spatial and non-spatial models

##    FILES REQUIRED:
##            forreri.ustr # n=104
##            forreri_noNAs.ustr # n=XX, created below
##            forreri_md90.ustr
##            forreri_md75.ustr
##            forreri_md50.ustr
##      Geographic sampling coordinates (long, lat) in tsv format:
##            forreri_metadata.txt # n=XX


# (1) Process input data --------------------------------------------------

# Read in ustr file
str <- read_tsv("../Separate_assemblies/iPyrad/forreri/forreri.ustr", col_names = FALSE) # 66,037 SNPs

str <-
  str %>% 
  dplyr::select(-c(X2:X5)) %>% 
  rename(INDV = X1) %>% 
  dplyr::filter(INDV != "V2797_PAC") # remove sample that has no coordinate data, leaving 103 individuals (206 obs)

# conStruct has trouble running with samples that have lots of missing data, so
# let's remove some of these individuals.
md90 <- c('MVZ149013_LCA', 'V2832_PAC', 'T3042_LCA', 'N918_LCA', 'T3762_LCA') # 5 inds; 98 inds remain

md75 <- c('MVZ149013_LCA', 'V2832_PAC', 'T3042_LCA', 'N918_LCA', 'T3762_LCA',
          'V1_PAC', 'T14192_papa_PAC', 'V5_PAC', 'V12_PAC', 'V2695_Rfor_PAC', 'MVZ269990_LCA', 
          'JA24801_PAC', 'V15_PAC', 'V2906_PAC', 'UTA60803_Rfor_PAC', 'MVZ269989_LCA', 'V1311_PAC', 
          'KEN592_LCA', 'V2831_PAC', 'KEN590_LCA', 'T14199_PAC', 'KEN591_LCA', 'V13_PAC', 'V14_PAC') # 24 inds; 79 inds remain

md50 <- c('MVZ149013_LCA', 'V2832_PAC', 'T3042_LCA', 'N918_LCA', 'T3762_LCA',
               'V1_PAC', 'T14192_papa_PAC', 'V5_PAC', 'V12_PAC', 'V2695_Rfor_PAC', 'MVZ269990_LCA', 
               'JA24801_PAC', 'V15_PAC', 'V2906_PAC', 'UTA60803_Rfor_PAC', 'MVZ269989_LCA', 'V1311_PAC', 
               'KEN592_LCA', 'V2831_PAC', 'KEN590_LCA', 'T14199_PAC', 'KEN591_LCA', 'V13_PAC', 'V14_PAC',
               'V2725_PAC', 'V16_PAC', 'V2797_PAC', 'V8_PAC', 'V1312_PAC', 'V2723_PAC', 
                'V2707_Rfor_PAC', 'MVZ269991_LCA', 'V2794_PAC', 'V2800_PAC', 'V2728_PAC', 'TF8644_Rspe_CMX', 
                'V2706_Rfor_PAC', 'TF8647_Rspe_CMX', 'T3703_LCA', 'V2702_Rfor_PAC', 'V2692_Rfor_PAC', 'TF8649_Rspe_CMX', 
               'V7_PAC', 'MVZ207327_LCA', 'TF8648_Rspe_CMX', 'V2697_Rfor_PAC', 'MVZ175967_LCA', 'T3715_LCA', 
               'TF8650_Rspe_CMX', 'T4029_bals_PAC', 'T3588_LCA', 'V2726_PAC', 'V2686_Rfor_PAC', 'V2685_Rfor_PAC', 
               'V2701_Rfor_PAC', 'V2687_Rfor_PAC') # 56 inds; 47 inds remain

md90_dat <- str %>% 
  dplyr::filter(!INDV %in% md90) # 196 rows (98 inds)

md75_dat <- str %>% 
  dplyr::filter(!INDV %in% md75) # 158 rows (79 inds)

md50_dat <- str %>% 
  dplyr::filter(!INDV %in% md50) # 96 rows (48 inds)

# Check for loci that have all missing values (-9s)
full <-
  str %>% 
  select(-(which((str %>% select(-INDV) %>% colSums()) == 
                   (nrow(str))*-9) + 1)) # 66,033 SNPs
md90_str <-
  md90_dat %>% 
  select(-(which((md90_dat %>% select(-INDV) %>% colSums()) == 
                   (nrow(md90_dat))*-9) + 1)) # 66,033 SNPs
md75_str <-
  md75_dat %>% 
  select(-(which((md75_dat %>% select(-INDV) %>% colSums()) == 
                   (nrow(md75_dat))*-9) + 1)) # 66,032 SNPs
md50_str <-
  md50_dat %>% 
  select(-(which((md50_dat %>% select(-INDV) %>% colSums()) == 
                   (nrow(md50_dat))*-9) + 1)) # 66,021 SNPs

# Export files
write_tsv(full, "forreri_noNAs.ustr", col_names = FALSE)
write_tsv(md90_str, "forreri_md90.ustr", col_names = FALSE)
write_tsv(md75_str, "forreri_md75.ustr", col_names = FALSE)
write_tsv(md50_str, "forreri_md50.ustr", col_names = FALSE)


# Check if distinct -------------------------------------------------------

# We need to see if there are any samples that are identical
md50_str %>% 
  dplyr::select(-INDV) %>% 
  duplicated()


# (2) Convert ustr to conStruct file ------------------------------------------

# Saves an RData file in working directory
data <- conStruct::structure2conStruct(infile = here("forreri_noNAs.ustr"),
                                       start.loci = 2, # first col that contains data
                                       start.samples = 1, # first row that contains samples
                                       onerowperind = FALSE, # two rows per individual
                                       missing.datum = -9, # missing data
                                       outfile = "forreri_noNAs")

data_md90 <- conStruct::structure2conStruct(infile = here("forreri_md90.ustr"),
                                            start.loci = 2, # first col that contains data
                                            start.samples = 1, # first row that contains samples
                                            onerowperind = FALSE, # two rows per individual
                                            missing.datum = -9, # missing data
                                            outfile = "forreri_md90")

data_md75 <- conStruct::structure2conStruct(infile = here("forreri_md75.ustr"),
                                            start.loci = 2, # first col that contains data
                                            start.samples = 1, # first row that contains samples
                                            onerowperind = FALSE, # two rows per individual
                                            missing.datum = -9, # missing data
                                            outfile = "forreri_md75")

data_md50 <- conStruct::structure2conStruct(infile = here("forreri_md50.ustr"),
                                            start.loci = 2, # first col that contains data
                                            start.samples = 1, # first row that contains samples
                                            onerowperind = FALSE, # two rows per individual
                                            missing.datum = -9, # missing data
                                            outfile = "forreri_md50")

load("forreri_md50.RData")


# (3) Calculate geographic distances --------------------------------------

dat <- md50_str

samps <- dat %>% dplyr::select(INDV) %>% distinct()

# Extract coords from metadata, making sure it's a) long then lat and b) in the same order as the ustr file
coords <- read_tsv(here("data", "forreri_metadata.txt"), col_names = TRUE) %>% 
  rename(x = long, y = lat, INDV = Bioinformatics_ID)
coords <- left_join(samps, coords) # This will deal with ordering as well as any mismatched samples that are present in coords and not in dat

geodist <- algatr::geo_dist(coords %>% 
                      dplyr::select(x, y))
colnames(geodist) <- coords$INDV
rownames(geodist) <- coords$INDV

# Check the order of individuals in coords is the same as in ustr file; CHANGE DEPENDING ON DATASET
rownames(data_md50) == rownames(geodist)

# Remove all other columns from the coordinate data
coords <- coords %>% dplyr::select(x, y) %>% as.matrix()


# (4) Run conStruct -------------------------------------------------------

forr_k4 <- conStruct(spatial = TRUE,
                     K = 4,
                     freqs = data_md75,
                     geoDist = geodist,
                     coords = coords,
                     prefix = "forreri_k4",
                     control = setNames(list(0.9), "adapt_delta"),
                     n.chains = 2,
                     n.iter = 15000,
                     make.figs = TRUE,
                     save.files = TRUE)

# Save results
save(forr_k4, file = "md50_k4_10k_conStruct.Robj")

### Export results for plotting
# load("../conStruct/forreri_k4_conStruct.results.Robj")
results <- conStruct.results$chain_2$MAP$admix.proportions %>% 
  as.data.frame() %>% 
  rename(k1 = V1,
         k2 = V2,
         k3 = V3,
         k4 = V4) %>% 
  cbind(coords)

write_tsv(results, "conStruct_md50_k4_results.txt")
