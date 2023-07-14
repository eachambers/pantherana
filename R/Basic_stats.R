library(tidyverse)
library(here)

theme_set(theme_cowplot())
setwd("~/Box Sync/Rana project/ddRADseq/ALL_RANA/Pooled_assembly/iPyrad/outfiles/stats_files")

## This code does the following:
##     1. Calculate proportions of missing data from PAUP* output file
##     2. Calculate average read depth from iPyrad stats files

##    FILES REQUIRED:
##          *_missing.txt # for pooled, ATL_MXPL, PAC_FOOT, and PAC_XX
##          pooled_dataset_assignments.txt
##          All s3 stats files from iPyrad output (located: Pooled_assembly/iPyrad/outfiles/stats_files/; not on Github)

##  Although there is only a single pooled assembly, some individuals needed to be removed
##  from the pooled assembly based on how much missing data they contained. Two subsets of
##  samples were created from the pooled assembly:
##      min_500 dataset contains samples that have at least 500 SNPs (n=595)
##      80p dataset contains samples that contain at most 80% missing data (n=414)


# Calculate average missing data ------------------------------------------

### FOR POOLED DATASETS:

miss <- read_tsv(here("data", "pooled_missing.txt"), col_names = TRUE)

# Average missing data across all 632 samples
miss %>% summarize(mean(percent_missing))

# For min_500 dataset:
miss %>% 
  dplyr::filter(pooled_dataset != "remove") %>% 
  summarize(mean_snps = mean(no_sites_present),
            min_snps = min(no_sites_present),
            max_snps = max(no_sites_present),
            mean(percent_missing),
            min_miss = min(percent_missing),
            max_miss = max(percent_missing),
            no_inds = n())

# For max_80p dataset:
miss %>% 
  dplyr::filter(pooled_dataset == "max_80p") %>% 
  summarize(mean_snps = mean(no_sites_present),
            min_snps = min(no_sites_present),
            max_snps = max(no_sites_present),
            mean(percent_missing),
            min_miss = min(percent_missing),
            max_miss = max(percent_missing),
            no_inds = n())

### FOR SEPARATE ASSEMBLIES:

miss_atl <- read_tsv(here("data", "ATL_MXPL_missing.txt"), col_names = TRUE)
miss_atl %>% 
  summarize(avg_perc_miss = mean(percent_missing),
            avg_snps = mean(no_sites_present),
            min_snps = min(no_sites_present),
            max_snps = max(no_sites_present))



# 2. Calculate average depth -------------------------------------------------

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
filenames = list.files(path = ".", pattern = "s3_cluster_stats.txt", full.names = FALSE)
datalist <- map(filenames, read_in_table)
dat <- bind_rows(datalist) # should be 633 obs of 10 vars

### Calculate average stat read depth

dat %>% summarize(mean(avg_depth_stat)) # 17.56313

### Calculate average stat read depth for each subset dataset

# Read in file with individual assignments into datasets
datasets <- read_tsv(here("data", "pooled_dataset_assignments.txt"), col_names = TRUE)

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

# And for the max80p dataset
# This dataset ONLY includes rows that say `max_80p`
final %>% 
  dplyr::filter(pooled_dataset == "max_80p") %>% 
  summarize(mean(avg_depth_stat),
            no_inds = n())




# -------------------------------------------------------------------------
# GRAVEYARD ---------------------------------------------------------------

# Process input data ------------------------------------------------------

# Retrieve metadata; also contains "order" which is the order of inds in the vcf (relevant later)
strata <- read_tsv("PAC_min10K_metadata.txt", col_names = TRUE) %>% 
  rename(order = "order_vcf")

# Remove individuals with lots of missing data using vcftools

# Run NJ tree -------------------------------------------------------------

# Process plink genetic distances
gendist <- as.data.frame(readr::read_tsv("PAC_project_outfiles/PAC_min10K_plink.dist", col_names = FALSE))
plink_names <- readr::read_tsv("PAC_project_outfiles/PAC_min10K_plink.dist.id", col_names = FALSE) %>%
  dplyr::select(-`X1`) %>%
  as.matrix()

# Assign row and col names according to sampleID
rownames(gendist) <- plink_names
colnames(gendist) <- plink_names

# Use ape to generate NJ tree


# Run SNMF ----------------------------------------------------------------

# Convert vcf to lfmm object
# Number of detected inds: 325
# Number of detected loci: 566177
# 37118 lines were removed because these are not SNPs.
lfmm <- 
  LEA::vcf2lfmm(input.file = "PAC_project_outfiles/PAC_min10K.vcf.recode.vcf", 
                output.file = "PAC_project_outfiles/PAC_min10K.lfmm", 
                force = TRUE)

# Run sNMF
PAC_min10K_2runs = snmf(lfmm,
                        iterations=1000, K = 1:25, rep = 2,
                        entropy = T, CPU = 8, ploidy = 2) # rep is no. of runs for each value of K
# above should be 10,000 iterations with 5 reps but needed it to go faster

# PAC_min10K_2runs <- load.snmfProject("PAC_min10K_2runs.snmfProject")

# figure out which of the runs for each K is the best (minimized)
which.min(cross.entropy(PAC_min10K_2runs, K = 25)) # do for K 1 through 25

## Take a look at the entropy values for each K
p_min10K <- plot(PAC_min10K_2runs, col = "dodgerblue4", pch = 19, cex = 1.2)

######### Do same but for 80p dataset now and fewer K range

# lfmm_80p <- 
#   LEA::vcf2lfmm(input.file = "PAC_project_outfiles/PAC_80p.vcf.recode.vcf", 
#                 output.file = "PAC_project_outfiles/PAC_80p.lfmm", 
#                 force = TRUE)
# 
# PAC_min10K_2runs = snmf(lfmm,
#                         iterations=1000, K = 12:25, rep = 1,
#                         entropy = T, CPU = 8, ploidy = 2)


# Visualize snmf results --------------------------------------------------

library(RColorBrewer)
colourCount = 20
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

# Let's do increments of 3 starting at 12
barchart(PAC_min10K_2runs, K = 25, run = 2, 
         border = NA, space = 0, col = getPalette(colourCount),
         xlab = "Individuals", ylab = "Ancestry proportions", 
         main = "Ancestry matrix", sort.by.Q = TRUE) -> bp

test <- as.data.frame(bp$order) %>% 
  rename(order = "bp$order") %>% 
  full_join(strata)

axis(1, at = 1:length(bp$order), 
     labels = test$sampleID, las = 3, 
     cex.axis = .2)


# Process PCA data --------------------------------------------------------

# Read in vcf and convert to genind object
# vcf <- vcfR::read.vcfR("PAC_project_outfiles/PAC_min10K.vcf.recode.vcf") # 603295 variants, 334 cols
# gen <- vcfR2genind(vcf) # takes a little while to convert
# save(gen, file = "genind_min10K.RDA")
load(file = "genind_min10K.RDA")

# Get all indNames; should be 325 inds
indNames(gen)

# Convert genind to genlight for dartR
gl <- gi2gl(gen)
# save(gl, file = "genlight_min10K.RDA")
load(file = "genlight_min10K.RDA")

# Coastal Mexico samples only:
forreri <- strata %>% 
  filter(tree_set == "forreri") %>% 
  dplyr::select(sampleID) %>% 
  as.list()

omiltemana <- strata %>% 
  filter(tree_set == "omiltemana") %>% 
  dplyr::select(sampleID) %>% 
  as.list()

macroglossa <- strata %>% 
  filter(tree_set == "macroglossa") %>% 
  dplyr::select(sampleID) %>% 
  as.list()

spectabilis <- strata %>% 
  filter(tree_set == "spectabilis") %>% 
  dplyr::select(sampleID) %>% 
  as.list()

foothills <- strata %>% 
  filter(tree_set == "foothills") %>% 
  dplyr::select(sampleID) %>% 
  as.list()

removeInds <- foothills$sampleID
gl_foot <- gl.drop.ind(gl, removeInds, recalc = FALSE, mono.rm = FALSE, verbose = NULL)
gl_foot$ind.names

# Assign strata to strata within the genind object
# I've assigned six strata: locality, fieldID, species, lat, long, state
strata(gen) <- strata

# Can take a look at separate parts of the strata by doing:
head(strata(gen, ~tree_set))

# Now, let's take one of these strata columns (species) and assign it as the pop
setPop(gen) <- ~tree_set

# Check that it worked
popNames(gen) # 25 diff species

# gen_foot <- gen[!row.names(gen@tab) %in% foothills$sampleID]
# gen_foot <- gen[indNames(gen) != foothills$sampleID] # doesn't work
# gen_foot <- gen[!indNames(gen) %in% foothills$sampleID]


# Run PCA -----------------------------------------------------------------

sum(is.na(gen$tab)) # 244245829

# Scale values to mean; 403690950 elements
X_full <- scaleGen(gen, NA.method = "mean") # takes a while

# Verify that it worked:
X_full[1:20, 1:20]

pca1_full <- dudi.pca(X_full, cent = FALSE, scale = FALSE, scannf = FALSE, nf = 3)

save(pca1_full, file = "PCA_min10K.RDA")

# Let's pull out samples that occur along coastal states of MX:
dat <- bind_cols(pca1_full$li, strata)


# Build PCA plots ---------------------------------------------------------

### Build PCA plots
palette <- c("KS"="#b67431", "MO"="#a8a2ca", "alterna"="gray74", "elapsoides"="gray30")

"#80a4bc"
"#eed2a7"
"#cd7a3f"

palette <- c("spnov_Panama", "cf_forreri" = "#73806d", "papagayo" = "#ba94a6", "macroglossa?" = "#9f9f9f", "magnaocularis", "sp6", 
             "macro_forreri" = "#9f9f9f", "forreri" = "#73806d", "miadis", "spectabilis" = "#6f82b7", "yavapaiensis", "spectabilis?" = "#6f82b7",
             "sp5", "lenca", "spnov_CR", "omiltemana", "atenquique_sp7", "sp8", "sp4", "brown_macro" = "#984625", "macro_papa" = "#9f9f9f",
             "macroglossa" = "#e0895a", "sp7", "balsas_aten", "macroglossa?forreri?" = "#9f9f9f")

coast %>% 
    ggplot(aes(x = Axis1, y = Axis2, group = state, color = state)) +
    geom_vline(xintercept = 0, linewidth = 1, color = "lightgrey") +
    geom_hline(yintercept = 0, linewidth = 1, color = "lightgrey") +
    geom_point(alpha = 0.5, size = 5) +
    theme(legend.position = "none",
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20)) +
    xlab("PC 1") +
    ylab("PC 2")


ggsave("FigS4EF_KSMO_PCA.pdf", width = 11, height = 6)



