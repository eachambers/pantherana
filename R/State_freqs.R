library(tidyverse)
library(here)

## This code does the following:
##     1. Calculates state frequencies for pooled assemblies to run RAxML-ng with an ascertainment bias correction
##        on each of the two pooled assemblies (min_500, max_80p)

##    FILES REQUIRED:
##          pooled_missing.txt # info on missing data and assignment of samples into various assemblies (from PAUP*; n=632)
##          pooled_statefreqs.txt # state frequencies obtained with PAUP* for pooled assembly (n=632)


# Import and join data ----------------------------------------------------

miss <- read_tsv(here("data", "pooled_missing.txt"), col_names = TRUE)
freqs <- read_tsv(here("data", "pooled_statefreqs.txt"), col_names = TRUE)

freqs <- read_tsv(here("data", "ATL_MXPL_relaxed_statefreqs.txt"), col_names = TRUE)

# Join datasets by sample ID and remove samples that have way too much missing data
# Some of the basefreqs have NA values which is why we'll do a left_join
dat <- full_join(miss, freqs, by = "sample_ID") %>%
  dplyr::filter(pooled_dataset != "remove") # 595 obs, 11 vars

# Calculate state frequencies for each dataset ----------------------------

# min_500 dataset includes all sites in dat
# atl_mxpl has 6583661 sites total (invar + var)
dat %>% 
  summarize(mean_A = mean(baseA),
            mean_C = mean(baseC),
            mean_G = mean(baseG),
            mean_T = mean(baseT),
            sites_A = (mean(baseA))*6583661,
            sites_C = mean(baseC)*6583661,
            sites_G = mean(baseG)*6583661,
            sites_T = mean(baseT)*6583661)

dat %>% 
  summarize()

# max_80p includes only max_80p dataset
dat %>% 
  dplyr::filter(pooled_dataset == "max_80p") %>%
  summarize(mean_A = mean(baseA),
            mean_C = mean(baseC),
            mean_G = mean(baseG),
            mean_T = mean(baseT),
            no_inds = n())
