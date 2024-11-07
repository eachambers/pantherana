# Distinguishing species boundaries from geographic variation
Phylogenetics of leopard frogs (subgenus *Pantherana*)

# Distinguishing species boundaries from geographic variation

This repository contains code for analyses and data visualization from:

**Chambers EA, Lara-Tufiño JD, Campillo-García G, Cisneros-Bernal AY, Dudek Jr. DJ, León-Règagnon V, Townsend JH, Flores-Villela O, & Hillis DM. (XXX) Distinguishing species boundaries from geographic variation. *XXX*, [DOI](LINKTODOI).**

This paper provided a suggested framework for testing whether genetic variation arises through species boundaries or intraspecific geographic variation and focused on species in the leopard frog species complex (*Rana pipiens*) in Mexico and Central America. All raw data and input files can be found on Dryad [here](DRYADLINK), and raw, demultiplexed fastqs can be found on SRA (BioProject ACCESSIONNUMBER).

**Assembly coding:**

Analyses were performed on different bioinformatics assemblies of our data. They are as follows:
* **Pooled assembly:** 93 samples (KS-MO transect samples (n=85), *L. alterna* (n=5), *L. elapsoides* (n=3))
* **KS-MO + alterna dataset:** 90 samples (KS-MO transect samples (n=85), *L. alterna* (n=5))
* **KS-MO transect dataset:** 85 samples (KS-MO transect samples)

## Data processing
* [Missing data calculations script](https://github.com/eachambers/ksmo_lampro/blob/main/R/data_analysis/Missing_data_calcs.R)
    1. Takes PAUP* output and calculates missing data per sample and on average
* [Missing data and loci in assembly visualization](https://github.com/eachambers/ksmo_lampro/blob/main/R/data_visualization/FigS1_Missing_data.R)
    1. Plots missing data and loci in assembly per sample (**Fig. S1**)

## Phylogenetic tree
* [Data visualization script](https://github.com/eachambers/ksmo_lampro/blob/main/R/data_visualization/FigS2_phylotree.R)
    1. Builds tree figure with node support values colorized according to bootstrap support
    2. Colorizes terminal edges according to group assignment (**Fig. S2**)

## Population Genetics Analyses

### PCA and correlation test
* [Analysis & visualization script](https://github.com/eachambers/ksmo_lampro/blob/main/R/data_visualization/FigS4_PCA.R)
    1. Runs a PCA on three subsets of the data: all samples, KS-MO + alterna, KS-MO transect only
    2. Visualizes PCA for PC1, PC2, and PC3 for each subset **Fig. S4**

### sNMF analysis
* [Analysis script](https://github.com/eachambers/ksmo_lampro/blob/main/R/data_analysis/sNMF_analyses.R)
    1. Runs sNMF on all three data subsets

### conStruct analysis
* [Analysis script](https://github.com/eachambers/ksmo_lampro/blob/main/R/data_analysis/conStruct_analysis.R)
    1. Runs conStruct on KS-MO transect samples

### Structure-based data visualization scripts
* [Piecharts, repelled (**Figs. 2a, S7**)](https://github.com/eachambers/ksmo_lampro/blob/main/R/data_visualization/Fig2A_S7_piecharts_repelled.R)
* [STRUCTURE-style plots (**Figs. 2b, S6**)](https://github.com/eachambers/ksmo_lampro/blob/main/R/data_visualization/Fig2B_S6_popgen.R)

## Other Analyses

### Admixture index analysis
* [Analysis script](https://github.com/eachambers/ksmo_lampro/blob/main/R/data_analysis/Admixture_index_analysis.R)
    1. Calculates diagnostic differences (given a user-specified threshold) between reference groups (east and west at KS-MO transect)
    2. Based on diagnostic diffs between reference groups, calculates admixture index values for remaining groups across transect
* [Data visualization script](https://github.com/eachambers/ksmo_lampro/blob/main/R/data_visualization/Fig3_admixture_index.R)
    1. Build histograms of admixture index values per group (**Fig. 3b**)
    2. Plots per-locus admixture index values (**Fig. 3c**)
    3. Plots average (across all loci) admixture index values over longitude (**Fig. 3d**)

### Fixed difference analysis
* [Analysis script](https://github.com/eachambers/ksmo_lampro/blob/main/R/data_analysis/Fixed_diff_analysis.R)
* [Data visualization script](https://github.com/eachambers/ksmo_lampro/blob/main/R/data_visualization/Fig4_S9_S10_Fixed_diff_analysis.R)
    1. Builds PCoA plots for pre- and post-collapsed groups (**Fig. 4**)
    2. Builds heat map with raw fixed difference values among individuals (**Fig. S9**)
    3. Builds heat map with raw fixed difference values among groups (**Fig. S10**)

### HZAR cline fitting
* [Generate HZAR input file from admixture index values](https://github.com/eachambers/ksmo_lampro/blob/main/R/data_analysis/Generate_HZAR_input_file.R)
* [Run HZAR](https://github.com/eachambers/ksmo_lampro/blob/main/R/data_analysis/HZAR_analysis.R)
* [Data visualization script](https://github.com/eachambers/ksmo_lampro/blob/main/R/data_visualization/FigS8_HZAR.R)
    1. Builds HZAR cline plot (**Fig. S8**)

### Morphological re-assessment of Armstrong et al. (2001)
* [Analysis and data visualization script](https://github.com/eachambers/ksmo_lampro/blob/main/R/data_visualization/Fig5_Armstrong_analysis.R)
    1. Takes in canonical axis 1 data from Armstrong et al. (2001)
    2. Builds histograms for reference and transect samples (**Fig. 5**)

### Map figure
* [Samples into group assignments on map (**Fig. 3a**)](https://github.com/eachambers/ksmo_lampro/blob/main/R/data_visualization/Fig3A_map_groups.R)
