library(tidyverse)
library(here)

## This code does the following:
##     1. Calculates state frequencies for pooled assemblies to run RAxML-ng with an ascertainment bias correction
##        on the pooled assembly

##    FILES REQUIRED:
##          pooled_missing.txt # info on missing data and assignment of samples into various assemblies (from PAUP*; n=632)
##          pooled_statefreqs.txt # state frequencies obtained with PAUP* for pooled assembly (n=632)


# Import and join data ----------------------------------------------------

miss <- read_tsv(here("data", "2_Data_processing", "data_files_input_into_scripts", "pooled_missing.txt"), col_names = TRUE)
freqs <- read_tsv(here("data", "2_Data_processing", "data_files_input_into_scripts", "pooled_statefreqs.txt"), col_names = TRUE)

# Join datasets by sample ID and remove samples that have way too much missing data
# Some of the basefreqs have NA values which is why we'll do a left_join
dat <- full_join(miss, freqs, by = "sample_ID") %>%
  dplyr::filter(pooled_dataset != "remove") # 595 obs, 11 vars


# Calculate state frequencies ---------------------------------------------

# min_10K refers to dataset that retains samples with at least 10K SNPs
dat %>% 
  dplyr::filter(pooled_dataset == "min_10K") %>%
  summarize(mean_A = mean(baseA),
            mean_C = mean(baseC),
            mean_G = mean(baseG),
            mean_T = mean(baseT), # sites_T = mean(baseT)*length
            no_inds = n())
