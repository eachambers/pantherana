# "#80a4bc" = spectabilis spectabilis
# "#6f82b7" = spectabilis chichicuahutla
# "#e0895a" = macroglossa
# "#e0b693" = macroglossa2
# "#a8a2ca" = omiltemana oaxaca
# "#ba94a6" = omiltemana guerrero
# "#984625" = berlandieri
# "#e97490" = berlandieri neovolcanica
# "#054051" = yavapaiensis / foothills
# "#f7cd5e" = magnaocularis
# "#428b9b" = atenquique long (sp 8)
# "#c0a06f" = atenquique short (sp 7)
# "#82ccc8" = sp nov CR/PM
# "#73806d" = forreri miadis / southern forreri
# "#eeb18f" = miadis
# "gray30" = forreri floresi
# "gray74" = forreri forreri
# "#647777" = forreri4
# "#3b808f" = forreri hillisi
# "#924c62" = lenca
# "#c8b656" = taylori


#' Retrieve colors for plotting Rana analyses
#'
#' @param dataset set colors specific to the analysis ("admixture", "tess" TODO)
#' @param dataset_name if dataset specified, specify dataset name ("ATL_MXPL", "PACMX", "CENTAM", "forreri", or "forreri_mvz")
#' @param save_imap whether to save Imap file for BPP
#' @param K_value if `dataset_name = "admixture"`, K value to retrieve colors for
#'
#' @return
#' @export
retrieve_kcols <- function(K_value, dataset, dataset_name = "ATL_MXPL", save_imap = FALSE){
  if (dataset == "admixture") {
    if (dataset_name == "ATL_MXPL") {
      if (K_value == 3) kcols = c(X1 = "#80a4bc", X3 = "#984625", X2 = "#e0895a")
      if (K_value == 4) kcols = c(X1 = "#e0895a", X3 = "#80a4bc", X2 = "#984625", X4 = "#e0b693")
      if (K_value == 5) kcols = c(X5 = "#80a4bc", X1 = "#6f82b7", X4 = "#e0b693", X3 = "#984625", X2 = "#e0895a")
      if (K_value == 6) kcols = c(X4 = "#80a4bc", X5 = "#6f82b7", X1 = "#e0b693", X3 = "#984625", 
                                  X2 = "#e0895a", X6 = "#e97490")
      if (K_value == 7) kcols = c(X1 = "#80a4bc", X7 = "#6f82b7", X4 = "#e0b693", X6 = "#e97490", 
                                  X2 = "#e0895a", X5 = "#984625", X3 = "#c8b656")
    }

    if (dataset_name == "CENTAM") {
      if (K_value == 3) kcols = c(X1 = "#73806d", X2 = "#924c62", X3 = "#e0895a")
      if (K_value == 4) kcols = c(X1 = "#e0895a", X2 = "#73806d", X3 = "#82ccc8", X4 = "#924c62")
      if (K_value == 5) kcols = c(X1 = "#c8b656", X2 = "#e0895a", X3 = "#e0b693", X4 = "#924c62", X5 = "#73806d")
      if (K_value == 6) kcols = c(X1 = "#e0b693", X2 = "#c8b656", X3 = "#924c62", X4 = "#73806d", X5 = "#82ccc8", X6 = "#e0895a")
      if (K_value == 7) kcols = c(X1 = "#c8b656", X2 = "#e0895a", X3 = "#924c62", X4 = "#73806d", X5 = "#82ccc8", X6 = "#ba94a6", X7 = "#e0b693")
    }

    if (dataset_name == "PACMX") {
      if (K_value == 3) kcols = c(X1 = "#73806d", X2 = "#ba94a6", X3 = "#e0895a")
      if (K_value == 4) kcols = c(X1 = "#ba94a6", X2 = "#73806d", X3 = "#6f82b7", X4 = "#e0895a")
      if (K_value == 5) kcols = c(X1 = "#73806d", X2 = "gray74", X3 = "#e0895a", X4 = "#ba94a6", X5 = "#6f82b7")
      if (K_value == 6) kcols = c(X1 = "#73806d", X2 = "#ba94a6", X3 = "#a8a2ca", X4 = "#e0895a", X5 = "#6f82b7", X6 = "#054051")
      if (K_value == 7) kcols = c(X1 = "#73806d", X2 = "#a8a2ca", X3 = "gray74", X4 = "#ba94a6", X5 = "#6f82b7", X6 = "#054051", X7 = "#e0895a")
      if (K_value == 8) kcols = c(X1 = "gray74", X2 = "#054051", X3 = "#e0895a", X4 = "gray30", X5 = "#73806d", X6 = "#428b9b", X7 = "#6f82b7", X8 = "#ba94a6")
      if (K_value == 9) kcols = c(X1 = "gray74", X2 = "#054051", X3 = "#ba94a6", X4 = "#6f82b7", X5 = "#e0895a", X6 = "#428b9b", X7 = "gray30", X8 = "#73806d", X9 = "#a8a2ca")
      if (K_value == 10) kcols = c(X1 = "#a8a2ca", X2 = "#ba94a6", X3 = "#73806d", X4 = "gray74", X5 = "#e0895a", X6 = "#e0b693", X7 = "#054051", X8 = "#428b9b", X9 = "gray30", X10 = "#6f82b7")
    }

    if (dataset_name == "forreri") {
      if (K_value == 2) kcols = c(X1 = "#73806d", X2 = "gray74")
      if (K_value == 3) kcols = c(X1 = "#73806d", X2 = "gray30", X3 = "gray74")
      if (K_value == 4) kcols = c(X1 = "gray30", X2 = "#73806d", X3 = "#eeb18f", X4 = "gray74")
      if (K_value == 5) kcols = c(X1 = "gray74", X2 = "#73806d", X3 = "#428b9b", X4 = "#d2a3a6", X5 = "gray30")
      if (K_value == 6) kcols = c(X1 = "gray30", X2 = "#eeb18f", X3 = "#d2a3a6", X4 = "gray74", X5 = "#428b9b", X6 = "#73806d")
    }
  }
  
  if (dataset == "spp_delim") {
    if (dataset_name == "forreri") {
      forreri <- import_admix_data(path = ".", file = "forreri/forreri_0.25miss_ldp", K_values = 5)
      dat = forreri[[1]]
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
      if (save_imap) write_tsv(kcols %>% dplyr::select(Bioinformatics_ID, pop), paste0("../../Landgen/spp_delim/", dataset_name, "-Imap.txt"), col_names = FALSE)
    }
    if (dataset_name == "spectabilis") {
      atlmx <- import_admix_data(path = ".", file = "ATL_MXPL/ATL_MXPL_relaxed_0.25miss_ldp", K_values = 5)
      dat = atlmx[[1]]
      kcols <-
        dat %>% 
        dplyr::mutate(pop = case_when(max_K == "X1" ~ "chic",
                                      max_K == "X2" ~ "macr",
                                      max_K == "X3" ~ "berl",
                                      max_K == "X4" ~ "macr2",
                                      max_K == "X5" ~ "spec"),
                      kcols = case_when(max_K == "X1" ~ "#6f82b7",
                                        max_K == "X2" ~ "#e0895a",
                                        max_K == "X3" ~ "#984625",
                                        max_K == "X4" ~ "#e0b693",
                                        max_K == "X5" ~ "#80a4bc")) %>% 
        dplyr::filter(pop == "chic" | pop == "spec")
      if (save_imap) write_tsv(kcols %>% dplyr::select(Bioinformatics_ID, pop), paste0("../../Landgen/spp_delim/", dataset_name, "-Imap.txt"), col_names = FALSE)
    }
    if (dataset_name == "berlandieri") {
      atlmx <- import_admix_data(path = ".", file = "ATL_MXPL/ATL_MXPL_relaxed_0.25miss_ldp", K_values = 7)
      dat = atlmx[[1]]
      kcols <-
        dat %>% 
        dplyr::filter(Bioinformatics_ID != "IRL57_LCA") %>% 
        dplyr::mutate(pop = case_when(max_K == "X1" ~ "spec",
                                      max_K == "X2" ~ "macr",
                                      max_K == "X3" ~ "tayl",
                                      max_K == "X4" ~ "macr2",
                                      max_K == "X5" ~ "berl",
                                      max_K == "X6" ~ "neov",
                                      max_K == "X7" ~ "chic"),
                      kcols = case_when(max_K == "X1" ~ "#80a4bc",
                                        max_K == "X2" ~ "#e0895a",
                                        max_K == "X3" ~ "#c8b656",
                                        max_K == "X4" ~ "#e0b693",
                                        max_K == "X5" ~ "#984625",
                                        max_K == "X6" ~ "#e97490",
                                        max_K == "X7" ~ "#6f82b7")) %>% 
        dplyr::filter(pop == "berl" | pop == "neov")
      if (save_imap) write_tsv(kcols %>% dplyr::select(Bioinformatics_ID, pop), paste0("../../Landgen/spp_delim/", dataset_name, "-Imap.txt"), col_names = FALSE)
    }
    if (dataset_name == "foothills") {
      pacmx <- import_admix_data(path = ".", file = "PACMX/new_PACMX_relaxed_0.25miss_ldp", K_values = 6)
      dat = pacmx[[1]]
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
      if (save_imap) write_tsv(kcols %>% dplyr::select(Bioinformatics_ID, pop), paste0("../../Landgen/spp_delim/", dataset_name, "-Imap.txt"), col_names = FALSE)
    }
  }
  return(kcols)
}
