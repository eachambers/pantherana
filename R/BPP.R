library(tidyverse)
library(here)
library(cowplot)
theme_set(theme_cowplot())
setwd("~/Box Sync/Rana project/ddRADseq/ALL_RANA/Landgen")

## The following file generates population assignment files (Imap files) for
## BPP and HHSD analyses using the results from ADMIXTURE.
##    (1) Gathers loci numbers from the lmiss Plink report
##    (2) Generate Imap file from ADMIXTURE results

##    FILES REQUIRED:
##            forreri.ustr # n=104
##            forreri_noNAs.ustr # n=XX, created below
##            forreri_md90.ustr
##            forreri_md75.ustr
##            forreri_md50.ustr
##      Geographic sampling coordinates (long, lat) in tsv format:
##            forreri_metadata.txt # n=XX


# (1) Retrieve loci numbers -----------------------------------------------

# Specify dataset
lmiss_prefix = "spp_delim/spectabilis" # forreri_0.75miss_ldp / foothills / spectabilis / forreri_0.25miss_ldp / forreri_0.5miss_ldp

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

source(here("R", "Pop_gen.R"))
source(here("R", "rana_colors.R"))
setwd("~/Box Sync/Rana project/ddRADseq/ALL_RANA/Separate_assemblies/ADMIXTURE/")

dataset_name = "spectabilis" # forreri / foothills / spectabilis / berlandieri

if (dataset_name == "foothills") metadata <- read_tsv(paste0(here("data"), "/PACMX_metadata.txt"), col_names = TRUE)
if (dataset_name == "forreri") metadata <- read_tsv(paste0(here("data"), "/forreri_metadata.txt"), col_names = TRUE)
if (dataset_name == "spectabilis" | dataset_name == "berlandieri") metadata <- read_tsv(paste0(here("data"), "/ATL_MXPL_metadata.txt"), col_names = TRUE)
if (dataset_name == "forreri") colvals = c("hilli" = "#428b9b", "miad" = "#73806d", "arce" = "#d2a3a6", "forr" = "gray74", "flor" = "gray30")
if (dataset_name == "foothills") colvals = c("omig" = "#a8a2ca", "magn" = "#f7cd5e", "omio" = "#ba94a6", 
                                             "yava" = "#054051", "ates" = "#c0a06f", "atel" = "#428b9b")
if (dataset_name == "spectabilis") colvals = c("spec" = "#80a4bc", "chic" = "#6f82b7")
if (dataset_name == "berlandieri") colvals = c("berl" = "#984625", "neov" = "#e97490")

dat <- retrieve_kcols(dataset = "spp_delim", dataset_name, save_imap = TRUE)
final <- left_join(dat, metadata, by = "Bioinformatics_ID")

# To export lists of individuals to keep and to remove:
write_tsv(dat %>% dplyr::select(Bioinformatics_ID), paste0("../../Landgen/spp_delim/", dataset_name, "_indv.txt"), col_names = FALSE) # list to keep
# List to remove
toremove <- setdiff(metadata %>% dplyr::select(Bioinformatics_ID),
                    dat %>% dplyr::select(Bioinformatics_ID)) # lapply(function(x) paste("^", x, sep = ""))

write_tsv(toremove, paste0("../../Landgen/spp_delim/", dataset_name, "_remove.txt"), col_names = FALSE)


# -------------------------------------------------------------------------
# Generate maps -----------------------------------------------------------

# Retrieve data for background map ----------------------------------------

world <- map_data("world")
centamer <- filter(world, region == "Belize" | region == "Guatemala" | region == "Honduras" | 
                     region == "Nicaragua" | region == "Costa Rica" | region == "El Salvador" | 
                     region == "Panama")

data("mxstate.map")
mxstate.map$group <- as.numeric(mxstate.map$group)
map_data <- bind_rows(centamer, mxstate.map)
  
if (dataset_name == "foothills") {
  # Include some U.S. states
  include <- c("california", "arizona", "new mexico", "nevada", "texas")
  # Get polygons for mapping
  states <- map_data("state") %>% 
    filter(region %in% include)
  map_data <-
    bind_rows(map_data, states)
}

map_data %>% 
  ggplot() +
  geom_polygon(aes(x = long, y = lat, group = group), color = "white", fill = "#e6e6e6", linewidth = 0.5) + 
  theme_map() +
  geom_point(data = final, aes(long, lat, color = pop), size = 4) +
  scale_color_manual(values = colvals) +
  coord_fixed() +
  theme(legend.position = "none") # export 6x4


# Build chord diagrams ----------------------------------------------------

library(circlize)

estM <- read_csv(here("data", "forreri_K5_merge_estM.txt")) %>% 
  dplyr::select(-`...4`)

# grid.col is for the actual ring colour; col is for the link color; order is obvious
gridcol = c(miad = "#73806d", flor = "gray30", hilli = "#428b9b", forr = "gray74", arce = "#d2a3a6")
order = c("forr", "flor", "hilli", "arce", "miad")

# Sort data based on whichever factor you're colorizing, and then adjust the 
# breaks depending on the balance of color in resulting arrows
colors <- colorRampPalette(brewer.pal(11, "Spectral"))(10)

chordDiagramFromDataFrame(estM, directional = 1, 
                          grid.col = gridcol, order = order,
                          col = colors, link.border = "black",
                          link.arr.length = 0.07,
                          link.arr.type = "big.arrow", link.sort = TRUE, 
                          direction.type = c("arrows", "diffHeight"), diffHeight = -0.02,
                          transparency = 0.25)


