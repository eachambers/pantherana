library(here)
source(here("R", "Pop_gen.R"))

retrieve_hhsd_coding <- function(dataset_name, save_imap = FALSE) {
  if (dataset_name == "forreri") {
    forreri <- import_admix_data(path = here("data", "admixture"), prefix = "forreri_0.25miss_ldp", K_values = 5)
    dat <- forreri$dat$K5
    kcols <-
      dat %>% 
      dplyr::mutate(pop = case_when(max_K == "X1" ~ "forr",
                                    max_K == "X2" ~ "miad",
                                    max_K == "X3" ~ "hilli",
                                    max_K == "X4" ~ "arce",
                                    max_K == "X5" ~ "flor"),
                    kcols = case_when(max_K == "X1" ~ "gray74",
                                      max_K == "X2" ~ "#73806d",
                                      max_K == "X3" ~ "#428b9b",
                                      max_K == "X4" ~ "#d2a3a6",
                                      max_K == "X5" ~ "gray30"))
  }
  if (dataset_name == "mxpl") {
    atlmx <- import_admix_data(path = here("data", "admixture"), prefix = "ATL_MXPL_relaxed_0.25miss_ldp", K_values = 6)
    dat = atlmx$dat$K6
    kcols <-
      dat %>% 
      dplyr::filter(Bioinformatics_ID != "IRL57_LCA") %>% 
      dplyr::mutate(pop = case_when(max_K == "X1" ~ "macr2",
                                    max_K == "X2" ~ "macr",
                                    max_K == "X3" ~ "berl",
                                    max_K == "X4" ~ "spec",
                                    max_K == "X5" ~ "chic",
                                    max_K == "X6" ~ "neov"),
                    kcols = case_when(max_K == "X1" ~ "#e0b693",
                                      max_K == "X2" ~ "#e0895a",
                                      max_K == "X3" ~ "#984625",
                                      max_K == "X4" ~ "#80a4bc",
                                      max_K == "X5" ~ "#6f82b7",
                                      max_K == "X6" ~ "#e97490")) %>% 
      dplyr::filter(pop == "chic" | pop == "spec" | pop == "neov" | pop == "berl")
  }
    
  if (dataset_name == "foothills") {
    pacmx <- import_admix_data(path = here("data", "admixture"), prefix = "new_PACMX_relaxed_0.25miss_ldp", K_values = 6)
    dat = pacmx$dat$K6
    yava_samps <- c("T14298_Ryav_OUT", "T14438_OUT", "T3447_yava_OUT", "T3449_yava_OUT")
    aten_short <- c("T2020_Aten_PAC", "T2022_Aten_PAC", "T2024_PAC", "T2025_PAC")
    aten_long <- c("T282_Atoy_CMX", "T283_Atoy_CMX", "T290_Atoy_CMX")
    kcols1 <-
      dat %>% 
      dplyr::filter(!Bioinformatics_ID %in% c(yava_samps, aten_short, aten_long)) %>% 
      dplyr::mutate(pop = case_when(max_K == "X1" ~ "miad",
                                    max_K == "X2" ~ "omio",
                                    max_K == "X3" ~ "omig",
                                    max_K == "X4" ~ "macr",
                                    max_K == "X5" ~ "chic",
                                    max_K == "X6" ~ "magn"),
                    kcols = case_when(max_K == "X1" ~ "#73806d",
                                      max_K == "X2" ~ "#ba94a6",
                                      max_K == "X3" ~ "#a8a2ca",
                                      max_K == "X4" ~ "#e0895a",
                                      max_K == "X5" ~ "#6f82b7",
                                      max_K == "X6" ~ "#f7cd5e")) %>% 
      dplyr::filter(pop == "magn" | pop == "omig" | pop == "omio")
    kcols2 <-
      dat %>% 
      dplyr::filter(Bioinformatics_ID %in% c(yava_samps, aten_short, aten_long)) %>% 
      dplyr::mutate(pop = case_when(Bioinformatics_ID %in% yava_samps ~ "yava",
                                    Bioinformatics_ID %in% aten_short ~ "ates",
                                    Bioinformatics_ID %in% aten_long ~ "atel"),
                    kcols = case_when(Bioinformatics_ID %in% yava_samps ~ "#054051",
                                      Bioinformatics_ID %in% aten_short ~ "#c0a06f",
                                      Bioinformatics_ID %in% aten_long ~ "#428b9b"))
    kcols <- bind_rows(kcols1, kcols2)
  }
  if (save_imap) write_tsv(kcols %>% dplyr::select(Bioinformatics_ID, pop), paste0(here("data", "hhsd"), "/", dataset_name, "-Imap.txt"), col_names = FALSE)
  return(kcols)
}

