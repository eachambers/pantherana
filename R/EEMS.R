library(devtools)
library(tidyverse)
library(here)
library(cowplot)
library(algatr)
library(maps)
library(Rcpp)
library(RcppEigen)
library(raster)
# library(rgeos)
library(sp)
library(reemsplots2) # devtools::install_github("dipetkov/reemsplots2")
library(rworldxtra)
theme_set(theme_cowplot())
setwd("~/Box Sync/Rana project/ddRADseq/ALL_RANA/Landgen/EEMS/")

# We first take our LD-pruned vcf file which only contains samples for which we have coordinates "forreri_n92_ldp.recode.vcf"
# Then we need to create a diffs file which contains a matrix of average pairwise genetic differences
# To do so, we can run Plink on our vcf as follows:
# `plink --vcf forreri_n92_ldp.recode.vcf --distance square 1-ibs --const-fid --allow-extra-chr --out forreri_n92_ldp`
# The above will generate two files: .dist which contains the distance matrix itself and .dist.id which contains the corresponding sample IDs

coords <- read_tsv("forreri_coords.txt", col_names = FALSE)
diffs <- read_tsv("forreri_n92_ldp.mdist", col_names = FALSE)
outer <- read_csv("forreri_broad_outer.csv", col_names = FALSE)
