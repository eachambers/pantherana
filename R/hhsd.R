library(tidyverse)
library(here)
library(cowplot)
library(maps)
library(mxmaps) # devtools::install_github("diegovalle/mxmaps")
library(usmap)
theme_set(theme_cowplot())

source(here("R", "hhsd_functions.R"))

## The following file generates population assignment files (Imap files) for
## HHSD analysis using the results from ADMIXTURE.
##    (1) Gathers loci numbers from the lmiss Plink report
##    (2) Generates Imap file from ADMIXTURE results

# (1) Retrieve loci numbers -----------------------------------------------

# Specify dataset
lmiss_prefix = "spp_delim/forreri_0.6miss_ldp" # forreri_0.75miss_ldp / foothills / forreri_0.6miss_ldp / mxpl

# Get the missing data report generated
lmiss <- read.table(paste0(lmiss_prefix, ".lmiss"), header = TRUE) %>% 
  dplyr::select(SNP) %>% 
  separate(SNP, into = c("locus", "position"), sep = "_") %>% 
  dplyr::select(locus)
loci <- as.data.frame(str_remove(lmiss$locus, "loc"))

# Generate file with locus numbers
write_tsv(loci, paste0(lmiss_prefix, "_locusnames.txt"), col_names = FALSE)

# The above file will then be used in the process_loci_files.py script.


# (2) Generate Imap file --------------------------------------------------

dataset_name = "mxpl" # forreri / foothills / mxpl

if (dataset_name == "foothills") metadata <- read_tsv(paste0(here("data"), "/PACMX_metadata.txt"), col_names = TRUE)
if (dataset_name == "forreri") metadata <- read_tsv(paste0(here("data"), "/forreri_metadata.txt"), col_names = TRUE)
if (dataset_name == "mxpl") metadata <- read_tsv(paste0(here("data"), "/ATL_MXPL_metadata.txt"), col_names = TRUE)

dat <- retrieve_hhsd_coding(dataset_name, save_imap = TRUE)
final <- left_join(dat, metadata, by = "Bioinformatics_ID")

# To export lists of individuals to keep and to remove:
write_tsv(dat %>% dplyr::select(Bioinformatics_ID), paste0(here("data", "hhsd"), "/", dataset_name, "_indv.txt"), col_names = FALSE) # list to keep
# List to remove
toremove <- setdiff(metadata %>% dplyr::select(Bioinformatics_ID),
                    dat %>% dplyr::select(Bioinformatics_ID)) # lapply(function(x) paste("^", x, sep = ""))

write_tsv(toremove, paste0(here("data", "hhsd"), "/", dataset_name, "_remove.txt"), col_names = FALSE)
