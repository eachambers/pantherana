library(tidyverse)
library(here)

## This code does the following:
##     1. Selects one SNP per locus for LD-pruning prior to running population genetic
##        structure analyses. The code does so with the same rules as how iPyrad generates
##        a usnps file: for each RAD tag (locus), it selects a SNP with the least amount of variant-
##        based missing data; if there is more than one SNP with the same amount of missing
##        data within the same RAD tag (locus), it will randomly select a SNP.


##    FILES REQUIRED (all generated using Plink 1.9 --missing):
##          The .lmiss files for relevant vcfs, which has per-variant missing data proportions.


# Define function ---------------------------------------------------------

#' Reads in Plink lmiss file and provides list of SNPs to retain
#'
#' @param path path to .lmiss file
#' @param name name of .lmiss file (without lmiss suffix)
#' @param method how lmiss was generated (`"plink"` is default; alternatively `"vcftools"`)
#' @param keepsmall whether to keep the max number of SNPs at 65K (for Plink)
#'
#' @return list of CHR and SNP coding to be retained
#' @export
#'
#' @examples
snps_to_ldprune <- function(path, name, method = "plink", keepsmall = FALSE){
  dat <- read.table(paste0(path, "/", name, ".lmiss"), header = TRUE)
  
  if (method == "vcftools") {
    dat <-
      dat %>% 
      dplyr::select(CHR, POS, N_MISS, F_MISS) %>% 
      mutate(SNP = paste0(CHR, "_pos", POS)) %>% 
      mutate(across("SNP", str_replace, "RAD_", "loc")) %>% 
      dplyr::select(-POS)
    
    ldpruned <-
      dat %>% 
      # Only consider one locus (CHR) at a time
      group_by(CHR) %>% 
      # Filter on minimum missing data prop
      dplyr::filter(F_MISS == min(F_MISS)) %>% 
      # Randomly select one entry per locus (CHR); should be (roughly) the same number of RAD tags
      slice_sample(n = 1)
    
    ldpruned <-
      ldpruned %>% 
      dplyr::select(CHR, SNP) %>% 
      # Split on _pos
      separate(SNP, into = c("locus", "POS"), sep = "_pos") %>% 
      # Rename columns so they're standardized with vcf
      rename(CHROM = CHR) %>% 
      # Remove unwanted columns
      dplyr::select(CHROM, POS)
  }
  
  if (method == "plink") {
    ldpruned <-
      dat %>% 
      # Only consider one locus (CHR) at a time
      group_by(CHR) %>% 
      # Filter on minimum missing data prop
      dplyr::filter(F_MISS == min(F_MISS)) %>% 
      # Randomly select one entry per locus (CHR); should be (roughly) the same number of RAD tags
      slice_sample(n = 1)
    
    # The vcftools command to subsample SNPs (--positions) needs two cols:
    # CHROM and POS (no header). Importantly, POS is not the same as SNP
    # in the current ldpruned object.
    ldpruned <-
      ldpruned %>% 
      dplyr::select(CHR, SNP) %>% 
      # Split on _pos
      separate(SNP, into = c("locus", "oldPOS"), sep = "_pos") %>% 
      # Add new col with old POS + 1
      mutate(POS = as.numeric(oldPOS) + 1) %>%
      # Rename columns so they're standardized with vcf
      rename(CHROM = CHR) %>% 
      # Remove unwanted columns
      dplyr::select(CHROM, POS)
    
    if (nrow(ldpruned) > 65000 & keepsmall) {
      ldpruned <-
        ldpruned %>%
        ungroup() %>% 
        sample_n(65000)
    }
  }
    
  write_delim(ldpruned, paste0(path, "/", name, "_ldpruned.txt"), col_names = FALSE)
}


# Run above function ------------------------------------------------------

# Should end up with the same number of SNPs as loci
snps_to_ldprune(path = here("data"), name = "ATL_MXPL_relaxed_0.25miss")
snps_to_ldprune(path = here("data"), name = "new_CENTAM_relaxed_0.25miss")
snps_to_ldprune(path = here("data"), name = "new_PACMX_relaxed_0.25miss")
snps_to_ldprune(path = here("data"), name = "forreri_0.25miss", method = "vcftools", keepsmall = TRUE)

