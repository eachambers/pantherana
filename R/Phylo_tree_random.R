library(ape)
library(ggtree)
library(tidyverse)
setwd("~/Box Sync/Rana project/ddRADseq/ALL_RANA/Pooled_assembly/RAxML-ng")

tree <- read.tree("pruned_tree.support.tree")

bs <- tree$node.label
bs <- as.data.frame(bs)

bs$bs <- as.numeric(bs$bs)

min(bs$bs)
average(bs$bs)

bs %>% summarize(meanbs = mean(bs))

bs %>% filter(bs > 75) %>% count()

plot(tree)

ggplot(tree)

## Check node numbering
ggtree(tree) +
  geom_text2(aes(subset = !isTip, label = node)) +
  geom_tiplab(size = 2)

## Root the tree
tree_rooted <- root(tree, node = 925, edgelabel = TRUE)

ggtree(tree_rooted) +
  geom_text2(aes(subset = !isTip, label = node)) +
  geom_tiplab(size = 2)
