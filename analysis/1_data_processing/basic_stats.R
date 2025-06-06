library(tidyverse)
library(here)

## This code does the following:
##    (1) Calculate average read depth from iPyrad stats files

##    FILES REQUIRED:
##          pooled_dataset_assignments.txt
##          All s3 stats files from iPyrad output

##  Although there is only a single pooled assembly, some individuals needed to be removed
##  from the pooled assembly based on how much missing data they contained.
##      min_10K dataset contains samples that have at least 10K SNPs (n=555)


# Calculate average depth -------------------------------------------------

#' Helper function to read tables and make sample ID a column
#'
#' @param x list of files to read in
#'
#' @return
#' @export
#'
#' @examples
read_in_table <- function(x){
  read.table(file = x, header = TRUE) %>% 
    rownames_to_column(var = "Bioinformatics_ID")
}

# Run above function
filenames = list.files(path = here("data", "1_Bioinformatics", "iPyrad_output_files", "stats_files"), 
                       pattern = "s3_cluster_stats.txt", full.names = FALSE)
datalist <- map(filenames, read_in_table)
dat <- bind_rows(datalist) # should be 633 obs of 10 vars

### Calculate average stat read depth

dat %>% summarize(mean(avg_depth_stat)) # 17.56313

### Calculate average stat read depth for each subset dataset

# Read in file with individual assignments into datasets
datasets <- read_tsv(here("data", "2_Data_processing", "data_files_input_into_scripts", "pooled_dataset_assignments.txt"), 
                     col_names = TRUE)

# If you read through our bioinformatics walkthrough, you'll know that after step 3,
# no clusters were found in a single sample (TF8608_Rber_CMX) so it was removed from 
# the assembly. Let's verify that now:
setdiff(dat$Bioinformatics_ID, datasets$Bioinformatics_ID) # TF8608_Rber_CMX

# Based on the above, let's remove that sample so it doesn't mess anything up.
dat <- dat %>% dplyr::filter(Bioinformatics_ID != "TF8608_Rber_CMX")

# Join together the two dfs
final <- left_join(dat, datasets, by = "Bioinformatics_ID")

# Now, let's calculate read depth for the min500 dataset
# This dataset includes all samples other than `remove` 
final %>% 
  dplyr::filter(pooled_dataset != "remove") %>% 
  summarize(mean(avg_depth_stat),
            no_inds = n())

# And for the dataset containing only samples that have min 10K SNPs
# This dataset ONLY includes rows that say `min_10K`
final %>% 
  dplyr::filter(pooled_dataset == "min_10K") %>% 
  summarize(mean(avg_depth_stat),
            no_inds = n())
