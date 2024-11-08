# Distinguishing species boundaries from geographic variation

This repository contains code for analyses and data visualization from:

**Chambers EA, Lara-Tufiño JD, Campillo-García G, Cisneros-Bernal AY, Dudek Jr. DJ, León-Règagnon V, Townsend JH, Flores-Villela O, & Hillis DM. (XXX) Distinguishing species boundaries from geographic variation. *XXX*, [DOI](LINKTODOI).**

This paper provided a suggested framework for testing whether genetic variation arises through species boundaries or intraspecific geographic variation and focused on species in the leopard frog species complex (*Rana pipiens*) in Mexico and Central America. All raw data and input files can be found on Dryad [here](DRYADLINK), and raw, demultiplexed fastqs can be found on SRA (BioProject ACCESSIONNUMBER).

**Assembly coding:**

Analyses were performed on different bioinformatics assemblies of our data. They are as follows:

| Assembly | No. inds | Analyses
| -------- | ------------- | ------------ |
| Pooled | 479 | - Phylogenetic inference |
| ATL_MXPL | 189 | - Admixture <br>- HHSD (*on mxpl, which excluded R. macroglossa*) |
| PACMX | 245 | - Admixture <br>- HHSD (*on foothills, which excluded R. forreri*) |
| CENTAM | 140 | - Admixture |
| forreri | 104 | - Admixture <br>- HHSD <br>- Landscape genomics <br>- FEEMS |

## Bioinformatics pipeline and data processing
* [Bioinformatics pipeline script](https://github.com/eachambers/pantherana/blob/main/R/bioinformatics_processing.sh)
    1. Runs iPyrad to obtain pooled assembly and four separate assemblies (see table above for details)
* [Post-processing script](https://github.com/eachambers/pantherana/blob/main/R/basic_data_characteristics.sh)
    1. Gets basic data characteristics from pooled assembly
* [Basic statistics script](https://github.com/eachambers/pantherana/blob/main/R/Basic_stats.R)
    1. Calculates average read depth from iPyrad stats files

## Phylogenetic tree inference
* [RAxML-ng analysis](https://github.com/eachambers/pantherana/blob/main/R/RAXML-ng.sh)
    1. Uses state frequency calculations from [`State_freqs.R`](https://github.com/eachambers/pantherana/blob/main/R/State_freqs.R) for running RAXML-ng with an ascertainment bias correction
    2. Removes ambiguous sites that may resolve to being invariant using [`Remove_invariant_sites.R`](https://github.com/eachambers/pantherana/blob/main/R/Remove_invariant_sites.R) script
    3. Runs RAxML-ng on pooled assembly

## Population structure analyses
* [Population structure processing script](https://github.com/eachambers/pantherana/blob/main/R/population_structure.sh)
    1. Further filters separate assemblies based on missing data
    2. Generates missing data reports for separate assemblies
    3. Removes SNPs based on output from [`LD-pruning.R`](https://github.com/eachambers/pantherana/blob/main/R/LD-pruning.R) script, which performs LD-pruning by selecting one random SNP (the least amount of missing data) for each RAD tag
    4. Runs admixture

## Landscape genomic analyses
### Environmental data
* [Environmental data processing](https://github.com/eachambers/pantherana/blob/main/R/Env_data.R)
    1. Retrieves environmental layers using WorldClim
    2. Checks for collinearity among environmental variables
    3. Runs a raster PCA for dimensional reduction of environmental variables (to deal with collinearity)
    4. Misc. visualizations of resulting raster PCA

### FEEMS
* [FEEMS pre-processing script](https://github.com/eachambers/pantherana/blob/main/R/FEEMS_preprocessing.sh)
    1. Further prune based on missing data and remove single individual lacking coordinates
    2. Set up FEEMS environment
* [FEEMS input files](https://github.com/eachambers/pantherana/blob/main/R/fEEMS.R)
    1. Generates input data files for FEEMS
* [Run FEEMS](https://github.com/eachambers/pantherana/blob/main/R/run_feems.py)
    1. Script to generate the spatial graph and run the cross-validation analysis of FEEMS

### Remaining LG analyses
* [Landscape genomics script](https://github.com/eachambers/pantherana/blob/main/R/Land_gen.R); uses functions imported from [`Land_gen_functions.R`](https://github.com/eachambers/pantherana/blob/main/R/Land_gen_functions.R)
    1. Runs a Mantel test, MMRR, GDM, and a PCA

## Model-based species delimitation
* [Species delimitation script](https://github.com/eachambers/pantherana/blob/main/R/HHSD.sh)
    1. Install HHSD
    2. Further prune HHSD datasets on the basis of missing data
    3. Subset vcfs containing only individuals relevant to HHSD analyses
    4. Retrieve full loci (as these are required for HHSD) based on output from `hhsd.R` and `process_loci_file.py` scripts
    5. Generate a Phylip file with SNPs and invariant sites
    6. Run HHSD using the Imap file generated using the `BPP.R` script
* [HHSD script](https://github.com/eachambers/pantherana/blob/main/R/hhsd.R); uses functions imported from [`hhsd_functions.R`](https://github.com/eachambers/pantherana/blob/main/R/hhsd_functions.R)
    1. Uses missing data report to generate locus names
    2. Generates Imap files required for running HHSD
* [Process iPyrad loci files](XXX)
    1. Generates separate txt files for each relevant locus from the iPyrad .loci files for each separate assembly

## Data visualization
* [Making maps](https://github.com/eachambers/pantherana/blob/main/R/Mapping.R); uses functions imported from [`Mapping_functions.R`](https://github.com/eachambers/pantherana/blob/main/R/mapping_functions.R), [`rana_colors.R`](https://github.com/eachambers/pantherana/blob/main/R/rana_colors.R), and [`hhsd_functions.R`](https://github.com/eachambers/pantherana/blob/main/R/hhsd_functions.R)
    1. Make range maps, with and without type localities indicated. Used for **Figs. 1C, 2A,** and **3A**
    2. Make map of localities for HHSD analysis, color-coded to Imap population assignment: **Fig. 4**
    3. Make DEM map of Mexico with biogeographic provinces plotted: **Fig. S3**
* [Population structure results](https://github.com/eachambers/pantherana/blob/main/R/Pop_gen_figures.R); uses functions imported from [`Pop_gen_functions.R`](https://github.com/eachambers/pantherana/blob/main/R/Pop_gen_functions.R) and [`rana_colors.R`](https://github.com/eachambers/pantherana/blob/main/R/rana_colors.R)
    1. Plots CV error: **Fig. S5 (top)**
    2. Builds structure-style plots for all K values: **Figs. 2B, 2C, 3B, 3C,** and **S5 (bottom)**
    3. Plots pie charts of admixture results on maps: **Figs. 2B, 2C, 3B,** and **3C**
    4. Identify sympatric localities and calculate distances between them
    5. Build structure-style plots for sympatric localities: **Figs. 2C** and **3C**
* [Landscape genomics results](https://github.com/eachambers/pantherana/blob/main/R/Land_gen.R); uses functions imported from [`Land_gen_functions.R`](https://github.com/eachambers/pantherana/blob/main/R/Land_gen_functions.R)
    1. Plots Mantel test results: **Fig. S6**
    2. Plots MMRR results: **Fig. 5A**
    3. Plots GDM results: **Fig. 5B**
* [Environmental data results](https://github.com/eachambers/pantherana/blob/main/R/Env_data.R)
    1. Plots environmental PC3: **Fig. S7A**
    2. Plots loadings for environmental PC3: **Fig. S7B**
* [Gene flow results](https://github.com/eachambers/pantherana/blob/main/R/fEEMS.R)
    1. Plots FEEMS results: **Fig. 5C**
* [Rana taxonomic history](https://github.com/eachambers/pantherana/blob/main/R/rana_taxonomy.R)
    1. Plots taxonomic history of the leopard frog complex: **Fig. S2**
