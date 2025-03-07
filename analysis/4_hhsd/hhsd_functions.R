#' Retrieves how putative species are encoded from admixture results and generates
#' HHSD Imap file with population assignments
#'
#' @param dataset_name separate assembly name; options are "forreri", "mxpl", or "foothills"
#' @param save_imap whether to save an Imap file for running HHSD (defaults to FALSE)
#'
#' @returns
#' @export
retrieve_hhsd_coding <- function(dataset_name, save_imap = FALSE) {
  path_admix = here("data", "3_Analyses", "2_popgen")
  if (dataset_name == "forreri") {
    forreri <- import_admix_data(path = paste0(path_admix, "/forreri/"), 
                                 prefix = "forreri_0.25miss_ldp", K_values = 5)
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
    atlmx <- import_admix_data(path = paste0(path_admix, "/ATL_MXPL/"), 
                               prefix = "ATL_MXPL_relaxed_0.25miss_ldp", K_values = 6)
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
    pacmx <- import_admix_data(path = paste0(path_admix, "/PACMX/"), 
                               prefix = "new_PACMX_relaxed_0.25miss_ldp", K_values = 6)
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
  if (save_imap) write_tsv(kcols %>% dplyr::select(Bioinformatics_ID, pop), paste0(here("data", "3_Analyses", "3_hhsd"), "/", dataset_name, "/input_files/", dataset_name, "-Imap.txt"), col_names = FALSE)
  return(kcols)
}

#' Builds figure with HHSD results with gdi and HPDs, faceted on migprior params
#'
#' @param dat tidy HHSD "decision.csv" results, see `hhsd.R` for how to generate
#' @param migpriors number of migpriors provided; will be how plot is faceted. Options are "single" or "multiple"
#'
#' @returns
#' @export
gdi_plot <- function(dat, migpriors) {
  if (migpriors == "single") {
    p <- dat %>% 
      ggplot(aes(x = gdi, y = node)) + 
      geom_rect(aes(xmin = 0, xmax = 1, ymin = 0, ymax = 2.5),
                fill = "lightgrey", alpha = 0.5) +
      geom_vline(xintercept = c(0.2, 0.7), linetype = "dashed", color = "grey30") +
      geom_errorbar(aes(xmin = hpd_25, xmax = hpd_975), linewidth = 0.5,
                    position = position_dodge(0.05)) +
      geom_line() +
      geom_point(aes(fill = node), color = "black", pch = 21, size = 3) +
      scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), expand = c(0, 0)) +
      theme(axis.title.y = element_blank(),
            legend.position = "none",
            strip.background = element_blank()) +
      facet_grid(~algorithm)
  }
  if (migpriors == "multiple") {
    p <- dat %>% 
      ggplot(aes(x = gdi, y = node)) + 
      geom_rect(aes(xmin = 0, xmax = 1, ymin = 0, ymax = 2.5),
                fill = "lightgrey", alpha = 0.5) +
      geom_vline(xintercept = c(0.2, 0.7), linetype = "dashed", color = "grey30") +
      geom_errorbar(aes(xmin = hpd_25, xmax = hpd_975), linewidth = 0.5,
                    position = position_dodge(0.05)) +
      geom_line() +
      geom_point(aes(fill = node), color = "black", pch = 21, size = 3) +
      scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), expand = c(0, 0)) +
      theme(axis.title.y = element_blank(),
            legend.position = "none",
            strip.background = element_blank()) +
      facet_grid(~migprior)
    # Two-level facet grid on migprior and algorithm
    # facet_grid(rows = vars(migprior), cols = vars(algorithm), 
    #            labeller = plot_labeller, scales = "free_y", switch = "y")
  }
  return(p)
}
