library(tidyverse)
library(here)

source(here("analysis", "4_hhsd", "hhsd_functions.R"))
source(here("data_viz", "pop_gen_functions.R"))

# example call: `Rscript hhsd.R "forreri_0.6miss_ldp" "forreri"`

#!/usr/bin/env Rscript # leave line commented
args = commandArgs(trailingOnly = TRUE)
lmiss_prefix = args[1]
dataset_name = args[2]

# Alternatively, if you want to run the script from within R, do:
lmiss_prefix = here("data", "3_Analyses", "3_hhsd", "foothills", "input_files", "foothills_06") # "forreri_0.6miss_ldp" # foothills_06 / mxpl_06
dataset_name = "foothills" # forreri / foothills / mxpl

## The following file generates population assignment files (Imap files) for
## HHSD analysis using the results from ADMIXTURE.
##    (1) Gathers loci numbers from the lmiss Plink report
##    (2) Generates Imap file from ADMIXTURE results


# (1) Retrieve loci numbers -----------------------------------------------

miss <- read.table(paste0(lmiss_prefix, ".lmiss"), header = TRUE)

if (max(miss$N_MISS) >= unique(miss$N_GENO)-2) print("There are some SNPs with no (or only 1) ind with data...")

# Get the missing data report generated
lmiss <- miss %>% 
  dplyr::select(SNP) %>% 
  separate(SNP, into = c("locus", "position"), sep = "_") %>% 
  dplyr::select(locus)
loci <- as.data.frame(str_remove(lmiss$locus, "loc"))

# Generate file with locus numbers
write_tsv(loci, paste0(lmiss_prefix, "_locusnames.txt"), col_names = FALSE)

# The above file will then be used in the process_loci_files.py script.


# (2) Generate Imap file --------------------------------------------------

if (dataset_name == "foothills") metadata <- read_tsv(paste0(here("data", "2_Data_processing", "data_files_input_into_scripts"), "/PACMX_metadata.txt"), col_names = TRUE)
if (dataset_name == "forreri") metadata <- read_tsv(paste0(here("data", "2_Data_processing", "data_files_input_into_scripts"), "/forreri_metadata.txt"), col_names = TRUE)
if (dataset_name == "mxpl") metadata <- read_tsv(paste0(here("data", "2_Data_processing", "data_files_input_into_scripts"), "/ATL_MXPL_metadata.txt"), col_names = TRUE)

dat <- retrieve_hhsd_coding(dataset_name, save_imap = TRUE)
final <- left_join(dat, metadata, by = "Bioinformatics_ID")

# To export lists of individuals to keep and to remove:
write_tsv(dat %>% dplyr::select(Bioinformatics_ID), paste0(here("data", "3_Analyses", "3_hhsd"), "/", dataset_name, "/input_files/", dataset_name, "_indv.txt"), col_names = FALSE) # list to keep
# List to remove
toremove <- setdiff(metadata %>% dplyr::select(Bioinformatics_ID),
                    dat %>% dplyr::select(Bioinformatics_ID)) # lapply(function(x) paste("^", x, sep = ""))

write_tsv(toremove, paste0(here("data", "3_Analyses", "3_hhsd"), "/", dataset_name, "/input_files/", dataset_name, "_remove.txt"), col_names = FALSE)
