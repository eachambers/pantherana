# TRees

library(ape)
setwd("~/Box Sync/Rana project/ddRADseq/ALL_RANA/Pooled_assembly/RAxML-ng/")

tree <- read.tree("min10K.raxml.support.tree")

remove <- read_tsv(here("data", "to_prune.txt"), col_names = FALSE)
remove <- remove$X1

pruned <- drop.tip(tree, tree$tip.label[match(remove, tree$tip.label)])

write.tree(pruned, "pruned_tree.support.tree")


# Generate trees for separate assemblies ----------------------------------

# PACMX
pacmx <- read_tsv(here("data", "PACMX_metadata.txt"), col_names = TRUE) %>% 
  filter(dataset != "remove") %>% # 245 inds
  dplyr::select(Bioinformatics_ID)

all_samps <- as.data.frame(pruned$tip.label) %>% # 479 tips
  rename(Bioinformatics_ID = `pruned$tip.label`)

# Make a list of samples to EXCLUDE from the tree
remove <- anti_join(all_samps, pacmx, by = "Bioinformatics_ID")

pacmx_tree <- drop.tip(pruned, pruned$tip.label[match(remove$Bioinformatics_ID, pruned$tip.label)]) # 240 remain
write.tree(pacmx_tree, "PACMX.tree")

# ========= ATL_MXPL
atlmxpl <- read_tsv(here("data", "ATL_MXPL_metadata.txt"), col_names = TRUE) %>% 
  filter(dataset != "remove") %>% # 189 inds
  dplyr::select(Bioinformatics_ID)

all_samps <- as.data.frame(pruned$tip.label) %>% # 479 tips
  rename(Bioinformatics_ID = `pruned$tip.label`)

# Make a list of samples to EXCLUDE from the tree
remove <- anti_join(all_samps, atlmxpl, by = "Bioinformatics_ID") # 299; 9 in atlmxpl and not pruned tree

atlmxpl_tree <- drop.tip(pruned, pruned$tip.label[match(remove$Bioinformatics_ID, pruned$tip.label)]) # 180 remain
write.tree(atlmxpl_tree, "ATL_MXPL.tree")

# ========= forreri
forreri <- read_tsv(here("data", "forreri_metadata.txt"), col_names = TRUE) %>% 
  dplyr::select(Bioinformatics_ID)

remove <- anti_join(all_samps, forreri, by = "Bioinformatics_ID")

forreri_tree <- drop.tip(pruned, pruned$tip.label[match(remove$Bioinformatics_ID, pruned$tip.label)]) # 103 remain
write.tree(forreri_tree, "forreri.tree")
