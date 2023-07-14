library(tidyverse)
library(here)

theme_set(theme_cowplot())

setwd("~/Box Sync/Rana project/ddRADseq/ALL_RANA/Separate_assemblies/ATL_MXPL/ADMIXTURE/")


# Look at CV error from ADMIXTURE ------------------------------------------

CV <- read_tsv(here("data", "ATL_MXPL_CV_error.txt"), col_names = TRUE)

CV %>% 
  ggplot(aes(x = K_value, y = CV_error)) +
  geom_point() +
  geom_line() +
  facet_grid(~dataset)


# Read in Q matrices for best Ks -------------------------------------------

# Get sample IDs in correct order
ids_relaxed <- read.table("ATL_MXPL_relaxed_0.25miss.fam") %>% 
  dplyr::select(V2) %>% 
  rename(Bioinformatics_ID = V2)
ids_stringent <- read.table("ATL_MXPL_stringent_0.25miss.fam") %>% 
  dplyr::select(V2) %>% 
  rename(Bioinformatics_ID = V2)

# Best K-values are: relaxed is 6-7, stringent is 7-8
# Get all four of these datasets
q_rel_k6 <- read_delim("ATL_MXPL_relaxed_0.25miss.6.Q", col_names = FALSE) %>% 
  mutate(K_value = 6)
q_rel_k6 <- bind_cols(q_rel_k6, ids_relaxed)

q_rel_k7 <- read_delim("ATL_MXPL_relaxed_0.25miss.7.Q", col_names = FALSE) %>% 
  mutate(K_value = 7)
q_rel_k7 <- bind_cols(q_rel_k7, ids_relaxed)
# Bind above results with samples
q_rel <- bind_rows(q_rel_k7, q_rel_k6)

q_str_k7 <- read_delim("ATL_MXPL_stringent_0.25miss.7.Q", col_names = FALSE) %>% 
  mutate(K_value = 7)
q_str_k7 <- bind_cols(q_str_k7, ids_stringent)

q_str_k8 <- read_delim("ATL_MXPL_stringent_0.25miss.8.Q", col_names = FALSE) %>% 
  mutate(K_value = 8)
q_str_k8 <- bind_cols(q_str_k8, ids_stringent)

# Bind above results with samples
q_str <- bind_rows(q_str_k8, q_str_k7)



# Build structure plot ----------------------------------------------------

results_rel_k7 <- q_rel %>%
  dplyr::filter(K_value == 7) %>% 
  # dplyr::select(!"X8") %>%
  tibble::rownames_to_column(var = "order") %>% 
  rowwise() %>% 
  mutate(max_K = which.max(c_across(X1:X7))) %>% 
  arrange(desc(max_K))

results_str_k8 <- q_str %>%
  dplyr::filter(K_value == 8) %>% 
  # dplyr::select(!"X8") %>%
  tibble::rownames_to_column(var = "order") %>% 
  rowwise() %>% 
  mutate(max_K = which.max(c_across(X1:X8))) %>% 
  arrange(desc(max_K))

# metadata <- read_tsv(here("data", "ATL_MXPL_metadata.txt"), col_names = TRUE)
final <- left_join(results_str_k8, metadata, by = "Bioinformatics_ID")
write_tsv(final, "stringent_k8_results.txt")





or <- results$order
results <- results %>% dplyr::arrange(factor(order, levels = or))
results$order <- factor(results$order, levels = results$order)

dat <-
  results %>% 
  tidyr::pivot_longer(names_to = "cluster", values_to = "proportion", -c(Bioinformatics_ID, K_value, order, max_K))

p_rel_k6 <-
  structure_plot_helper(dat, 
                        kcols = c("#73806d", "#ba94a6", "#6f82b7", "#e0895a", "#984625", "#9f9f9f"))
p_str_k8 <-
  structure_plot_helper(dat, 
                        kcols = c("#73806d", "#ba94a6", "#6f82b7", "#e0895a", "#984625", "#9f9f9f", "#eed2a7", "gray30"))


ids <- as.list(unique(dat$Bioinformatics_ID))

p_rel_k6 +
  scale_x_discrete(expand = c(0,0), labels = ids) +
  ggtitle("relaxed_k6 (n=199)") # 171 / 199

p_rel_k7 +
  scale_x_discrete(expand = c(0,0), labels = ids) +
  ggtitle("relaxed_k7 (n=199)") # 171 / 199

p_str_k8 +
  scale_x_discrete(expand = c(0,0), labels = ids) +
  ggtitle("stringent_k8 (n=171)") # 171 / 199
  


# Run PCA -----------------------------------------------------------------


