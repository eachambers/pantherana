# Distinguishing species boundaries from geographic variation

***This manuscript is currently in review; therefore, we discourage the use of the following code at this time.***

This repository contains code for analyses and data visualization from:

**Chambers EA, Lara-Tufiño JD, Campillo-García G, Cisneros-Bernal AY, Dudek Jr. DJ, León-Règagnon V, Townsend JH, Flores-Villela O, & Hillis DM. (XXX) Distinguishing species boundaries from geographic variation. *XXX*, [DOI](LINKTODOI).**

This paper provided a suggested framework for testing whether genetic variation arises through species boundaries or intraspecific geographic variation and focused on species in the leopard frog species complex (*Rana pipiens*) in Mexico and Central America. All raw data and input files can be found on Dryad [here](https://doi.org/10.5061/dryad.cjsxksnhc), and raw, demultiplexed fastqs can be found on SRA (BioProject PRJNA1233814).

To run code in this repository, download the supplementary data files available on Dryad and dump them all into the [data](https://github.com/eachambers/pantherana/tree/main/data) directory. Any resulting plots from these scripts will be dumped into the `plots` directory.

**Assembly nomenclature:**

Analyses were performed on different bioinformatics assemblies of our data: a pooled assembly, which included all individuals, and four separate assemblies, which were various subsets of individuals from the pooled assembly. They are as follows:

| Assembly | No. inds | Analyses
| -------- | ------------- | ------------ |
| Pooled | 479 | - Phylogenetic inference |
| ATL_MXPL | 189 | - Admixture <br>- HHSD (*on mxpl, which excluded R. macroglossa*) |
| PACMX | 245 | - Admixture <br>- HHSD (*on foothills, which excluded R. forreri*) |
| CENTAM | 140 | - Admixture |
| forreri | 104 | - Admixture <br>- HHSD <br>- Landscape genomics <br>- FEEMS |

## Bioinformatics pipeline and data processing
Files are available within the `analysis/1_data_processing` directory [here](https://github.com/eachambers/pantherana/tree/main/analysis/1_data_processing).
* [Bioinformatics pipeline script](https://github.com/eachambers/pantherana/blob/main/R/bioinformatics_processing.sh)
    - Runs iPyrad to obtain pooled assembly and four separate assemblies (see table above for details)
* [Post-processing script](https://github.com/eachambers/pantherana/blob/main/R/basic_data_characteristics.sh)
    - Gets basic data characteristics from pooled assembly
* [Basic statistics script](https://github.com/eachambers/pantherana/blob/main/R/Basic_stats.R)
    - Calculates average read depth from iPyrad stats files

## Phylogenetic tree inference
Files are available within the `analysis/2_phylo` directory [here](https://github.com/eachambers/pantherana/tree/main/analysis/2_phylo).
* [RAxML-ng analysis](https://github.com/eachambers/pantherana/blob/main/analysis/2_phylo/RAXML-ng.sh)
    - Uses state frequency calculations from [`state_freqs.R`](https://github.com/eachambers/pantherana/blob/main/analysis/2_phylo/state_freqs.R) for running RAXML-ng with an ascertainment bias correction
    - Removes ambiguous sites that may resolve to being invariant using [`remove_invariant_sites.R`](https://github.com/eachambers/pantherana/blob/main/analysis/2_phylo/remove_invariant_sites.R) script
    - Runs RAxML-ng on pooled assembly

## Population structure analyses
Files are available within the `analysis/3_popgen` directory [here](https://github.com/eachambers/pantherana/tree/main/analysis/3_popgen).
* [Population structure processing script](https://github.com/eachambers/pantherana/blob/main/analysis/3_popgen/population_structure.sh)
    - Further filters separate assemblies based on missing data
    - Generates missing data reports for separate assemblies
    - Removes SNPs based on output from [`LD-pruning.R`](https://github.com/eachambers/pantherana/blob/main/analysis/3_popgen/LD-pruning.R) script, which performs LD-pruning by selecting one random SNP (the least amount of missing data) for each RAD tag
    - Runs admixture

## Model-based species delimitation
Files are available within the `analysis/4_hhsd` directory [here](https://github.com/eachambers/pantherana/tree/main/analysis/4_hhsd).
* [Species delimitation script](https://github.com/eachambers/pantherana/blob/main/analysis/4_hhsd/HHSD.sh)
    - Install HHSD
    - Further prune HHSD datasets on the basis of missing data
    - Subset vcfs containing only individuals relevant to HHSD analyses
    - Retrieve full loci (as these are required for HHSD) based on output from `hhsd.R` and `process_loci_file.py` scripts (see below)
    - Generate a Phylip file with SNPs and invariant sites
    - Run HHSD
* [HHSD script](https://github.com/eachambers/pantherana/blob/main/analysis/4_hhsd/hhsd.R); uses functions imported from [`hhsd_functions.R`](https://github.com/eachambers/pantherana/blob/main/analysis/4_hhsd/hhsd_functions.R)
    - Uses missing data report to generate locus names
    - Generates Imap files required for running HHSD
* [Process iPyrad loci files](https://github.com/eachambers/pantherana/blob/main/analysis/4_hhsd/process_loci_file.py)
    - Generates separate txt files for each relevant locus from the iPyrad .loci files for each separate assembly
 
## Landscape genomic analyses
Files are available within the `analysis/5_landgen` directory [here](https://github.com/eachambers/pantherana/tree/main/analysis/5_landgen).
* [Environmental data processing](https://github.com/eachambers/pantherana/blob/main/analysis/5_landgen/env_data.R)
    - Retrieve environmental layers using WorldClim
    - Check for collinearity among environmental variables
    - Run a raster PCA for dimensional reduction of environmental variables (to deal with collinearity)
* [Landscape genomics script](https://github.com/eachambers/pantherana/blob/main/analysis/5_landgen/land_gen.R); uses functions imported from [`land_gen_functions.R`](https://github.com/eachambers/pantherana/blob/main/analysis/5_landgen/land_gen_functions.R)
    - Run a Mantel test, MMRR, GDM, and a PCA

## FEEMS
Files are available within the `analysis/6_feems` directory [here](https://github.com/eachambers/pantherana/tree/main/analysis/6_feems).
* [FEEMS pre-processing script](https://github.com/eachambers/pantherana/blob/main/analysis/6_feems/FEEMS_preprocessing.sh)
    - Further prune based on missing data and remove single individual lacking coordinates
    - Set up FEEMS environment
* [FEEMS input files](https://github.com/eachambers/pantherana/blob/main/analysis/6_feems/FEEMS.R)
    - Generate input data files for FEEMS
* [Run FEEMS](https://github.com/eachambers/pantherana/blob/main/analysis/6_feems/run_feems.py)
    - Script to generate the spatial graph and run the cross-validation analysis of FEEMS

## Data visualization
* [Making maps](https://github.com/eachambers/pantherana/blob/main/R/Mapping.R); uses functions imported from [`Mapping_functions.R`](https://github.com/eachambers/pantherana/blob/main/R/mapping_functions.R), [`rana_colors.R`](https://github.com/eachambers/pantherana/blob/main/R/rana_colors.R), and [`hhsd_functions.R`](https://github.com/eachambers/pantherana/blob/main/R/hhsd_functions.R)
    - Make range maps, with and without type localities indicated. Used for **Figs. 1C, 2A,** and **3A**
    - Make map of localities for HHSD analysis, color-coded to Imap population assignment: **Fig. 4**
    - Make DEM map of Mexico with biogeographic provinces plotted: **Fig. S3**
* [Population structure results](https://github.com/eachambers/pantherana/blob/main/R/Pop_gen_figures.R); uses functions imported from [`Pop_gen_functions.R`](https://github.com/eachambers/pantherana/blob/main/R/Pop_gen_functions.R) and [`rana_colors.R`](https://github.com/eachambers/pantherana/blob/main/R/rana_colors.R)
    - Plot CV error: **Fig. S5 (top)**
    - Build structure-style plots for all K values: **Figs. 2B, 2C, 3B, 3C,** and **S5 (bottom)**
    - Plot pie charts of admixture results on maps: **Figs. 2B, 2C, 3B,** and **3C**
    - Identify sympatric localities and calculate distances between them
    - Build structure-style plots for sympatric localities: **Figs. 2C** and **3C**
* [Landscape genomics results](https://github.com/eachambers/pantherana/blob/main/R/Land_gen.R); uses functions imported from [`Land_gen_functions.R`](https://github.com/eachambers/pantherana/blob/main/R/Land_gen_functions.R)
    - Plot Mantel test results: **Fig. S6**
    - Plot MMRR results: **Fig. 5A**
    - Plot GDM results: **Fig. 5B**
* [Environmental data results](https://github.com/eachambers/pantherana/blob/main/R/Env_data.R)
    - Plot environmental PC3: **Fig. S7A**
    - Plot loadings for environmental PC3: **Fig. S7B**
* [Gene flow results](https://github.com/eachambers/pantherana/blob/main/R/fEEMS.R)
    - Plot FEEMS results: **Fig. 5C**
    - Plot FEEMS cross-validation analysis results: **Fig. S8**
* [Rana taxonomic history](https://github.com/eachambers/pantherana/blob/main/R/rana_taxonomy.R)
    - Plot taxonomic history of the leopard frog complex: **Fig. S2**
