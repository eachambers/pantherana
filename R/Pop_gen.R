
#' Plots CV error results for range of K-values
#'
#' @param dat path to tsv containing CV error results; can get this in bash by doing `grep -h CV prefix*.out`
#' @param faceted if FALSE (default), only a single dataset provided; if TRUE, will facet based on dataset name
#'
#' @return
#' @export
cv_error <- function(dat, faceted = FALSE){
  CV <- read_tsv(dat, col_names = TRUE)
  minK <- min(CV$K_value)
  maxK <- max(CV$K_value)
  
  p <-
    CV %>% 
    ggplot(aes(x = K_value, y = CV_error)) +
    geom_point() +
    geom_line() +
    scale_x_continuous(breaks = seq(minK, maxK, 1))
  
  if (faceted) p <- p + facet_grid(~dataset)
  
  return(p)
}

#' Imports Q matrices and IDs (from corresponding ADMIXTURE *.fam file)
#'
#' @param path path to Q matrices (defaults to current wd)
#' @param file file name (without suffix!)
#' @param K_values list of K-values to import
#'
#' @return
#' @export
import_admix_data <- function(path = ".", file, K_values){
  ids <- read_table(paste0(file, ".fam"), col_names = FALSE) %>% 
    dplyr::select(X2) %>%
    pull(X2)

  dat <- pmap(tibble(path = path, file = file, K_value = K_values), import_admix_helper, ids = ids)
  
  return(dat)
}

#' Helper function for importing ADMIXTURE results
#'
#' @param path path to Q matrices (defaults to current wd)
#' @param file file name (without suffix!)
#' @param K_value list of K-value(s) to import
#' @param ids vector of sample IDs
#'
#' @return
#' @export
import_admix_helper <- function(path, file, K_value, ids){
  dat <-
    read_delim(paste0(path, "/", file, ".", K_value, ".Q"), col_names = FALSE) %>%
    rowwise() %>% 
    mutate(max_K_val = max(c_across(everything()))) %>% 
    rowwise() %>% 
    mutate(max_K = colnames(.)[which.max(c_across())],
           K_val = K_value) %>% 
    cbind(., ids) %>% 
    rename(Bioinformatics_ID = ids) %>% 
    rownames_to_column(var = "order")
}

#' Main function to build structure-style bar plot
#'
#' @param dat output from `import_admix_data()` function (for one K-value)
#' @param K_value K-value you'd like plotted
#' @param order how to order samples in plot (defaults to maximum K-value "max_K"; other option is "manual")
#' @param man_order if order = "manual", manually provide vector of sample ordering (numbers, NOT sample IDs)
#' @param kcols colors to be used
#' @param xaxis_labels whether to label x axis with sample IDs (defaults to TRUE)
#' @param write_output whether to save output (and combine with metadata; defaults to TRUE)
#' @param metadata_path if write_output = TRUE, full path to metadata file (must have Bioinformatics_ID column)
#' @param output_name if write_output = TRUE, name to save output file a
#'
#' @return
#' @export
build_str_plot <- function(dat, K_value, order = "max_K", man_order = NULL, kcols, 
                           xaxis_labels = TRUE, write_output = TRUE, metadata_path = NULL, 
                           output_name = NULL, export_plot = NULL){
  if(order == "max_K") {
    dat <- dat %>% 
      dplyr::arrange(desc(max_K_val)) %>% 
      dplyr::arrange(max_K)
    or <- dat$order
    dat <- dat %>% 
      dplyr::arrange(factor(order, levels = or))
    dat$order <- factor(dat$order, levels = dat$order)
  }
  
  if(order == "manual") {
    dat <- dat %>% 
      dplyr::arrange(factor(order, levels = man_order))
    dat$order <- factor(dat$order, levels = dat$order)
  }
  
  dat <-
    dat %>% 
    tidyr::pivot_longer(names_to = "cluster", values_to = "proportion", 
                        -c(Bioinformatics_ID, K_val, order, max_K, max_K_val))

  p <- structure_plot_helper(dat, kcols)
  
  if (xaxis_labels) {
    ids <- as.list(unique(dat$Bioinformatics_ID))
    p <- p + 
      scale_x_discrete(expand = c(0,0), labels = ids) +
      theme(axis.text.x = element_text(angle = 90, size = 5))
  }
  
  if (write_output) {
    metadata <- read_tsv(metadata_path, col_names = TRUE)
    final <- left_join(dat, metadata, by = "Bioinformatics_ID") %>% 
      pivot_wider(names_from = cluster, values_from = proportion)
    write_tsv(final, paste0(output_name, "_K", K_value, ".txt"))
  }
  
  if (export_plot) {
    cowplot::save_plot(filename = paste0(output_name, "_K", K_value, ".pdf"), p, base_width = 14, base_height = 8)
  }
  return(p)
}

#' Helper function to actually build structure plot
#'
#' @param dat tidy df with order, proportion, and cluster cols
#' @param kcols colors to be used
#'
#' @return Structure-style plot
#' @export
structure_plot_helper <- function(dat, kcols){
  dat %>% 
    ggplot2::ggplot(aes(x = order, y = proportion, fill = cluster)) +
    ggplot2::geom_bar(stat = "identity") + 
    panel_border() + 
    scale_y_continuous(expand = c(0,0)) +
    # scale_x_discrete(expand = c(0,0)) +
    scale_fill_manual(values = kcols) +
    theme_cowplot() %+replace% theme(axis.line = element_line(colour = "black"),
                                     axis.text.x = element_blank(),
                                     axis.text.y = element_blank(),
                                     axis.ticks = element_blank(),
                                     axis.title.x = element_blank(),
                                     axis.title.y = element_blank(),
                                     legend.position = "none",
                                     panel.border = element_rect(fill = NA, colour = "black", linetype = "solid", linewidth = 1.5),
                                     strip.text.y = element_text(size = 30, face = "bold"),
                                     strip.background = element_rect(colour = "white", fill = "white"),
                                     panel.spacing = unit(-0.1, "lines"))
}


# -------------------------------------------------------------------------
# PIE CHARTS --------------------------------------------------------------

#' Make pie charts with data
#'
#' @param dat ADMIXTURE results for a single sample
#' @param kvalmaxcol column name corresponding to max K-value
#' @param kcols colors to assign, must be equal to K-value
#' @param alpha transparency of pies
#' @param legend.position whether to keep legend (default is "none")
#'
#' @return
#' @export
make_pies <- function(dat, kvalmaxcol, kcols, alpha = 0.5, legend.position = "none"){
  dat %>% 
    # TODO: need to automate the col names below somehow
    pivot_longer(names_to = "cluster", values_to = "proportion", X1:kvalmaxcol) %>% 
    ggplot() +
    geom_bar(aes(x = "", y = proportion, fill = cluster),
             stat = "identity", width = 1, alpha = alpha) +
    coord_polar("y") +
    theme_void() +
    theme(legend.position = legend.position) +
    scale_fill_manual(values = kcols)
}

#' Draws pie charts on to map
#'
#' @param plot map to plot pies on top of
#' @param lat latitude coord for pies
#' @param long longitude coord for pies
#' @param height height of pies
#' @param width width of pies
#'
#' @return
#' @export
draw_pies <- function(plot, lat, long, height = 1, width = 1) {
  draw_plot(plot = plot, 
            x = long, 
            y = lat,
            height = height,
            width = width, 
            hjust = 0.5, 
            vjust = 0.5)
}

#' Takes in coordinate data and makes new columns for repelled values
#'
#' @param dat dataframe containing "lat" and "long" cols for coordinates
#' @param padding_param repel_text::box.padding parameter setting
#' @param force_param repel_text::force parameter setting
#' @param ignore_repel if FALSE (default), no values to ignore with repelling; if TRUE, user must provide list of samples to ignore
#' @param ignore_repel_list if `ignore_repel = TRUE`, list of samples to not repel
#' @param ignore_repel_condition if `ignore_repel = TRUE`, condition to ignore repelling on (e.g., "long>-100")
#'
#' @return dataframe with appended columns "new_long" and "new_lat" that have repelled values
#' @export
repel_coords <- function(dat, padding_param = 0.5, force_param = 4, ignore_repel = FALSE, ignore_repel_list = NULL, ignore_repel_condition = NULL){
  dat_w_repel <- dat %>%
    dplyr::filter(!is.na(lat)) %>% 
    bind_cols(dat %>% 
                dplyr::filter(!is.na(lat)) %>% 
                dplyr::select(x = long, y = lat) %>% 
                mutate(label = 'XX') %>% 
                repel::repel_text(box.padding = padding_param, force = force_param) %>%
                dplyr::select(new_lat = y, new_long = x))
  
  if (ignore_repel) {
    if (!is.null(ignore_repel_condition)) {
      dat_w_repel <-
        dat_w_repel %>%
        mutate(to_repel = ignore_repel_condition) %>%
        mutate(new_lat = ifelse(to_repel, new_lat, lat),
               new_long = ifelse(to_repel, new_long, long))
    }
    
    if (!is.null(ignore_repel_list)) {
      dat_w_repel <-
        dat_w_repel %>%
        # TODO below is incorrect; fix later
        mutate(to_repel = ignore_repel_condition) %>%
        mutate(new_lat = ifelse(to_repel, new_lat, lat),
               new_long = ifelse(to_repel, new_long, long))
    }
  }
  return(dat_w_repel)
}

# -------------------------------------------------------------------------
# TESS3 FUNCTIONS ---------------------------------------------------------

# TODO potentially remove below?

tess_ggbarplot <- function(qmat, ggplot_fill = algatr_col_default("ggplot"), sort_by_Q = TRUE, legend = TRUE) {
  # Get K
  K <- ncol(qmat)
  dat <- as.data.frame(qmat) %>%
    tibble::rownames_to_column(var = "order")
  
  if (sort_by_Q) {
    gr <- apply(qmat, MARGIN = 1, which.max)
    gm <- max(gr)
    gr.o <- order(sapply(1:gm, FUN = function(g) mean(qmat[, g])))
    gr <- sapply(gr, FUN = function(i) gr.o[i])
    or <- order(gr)

    dat <- dat %>%
      dplyr::arrange(factor(order, levels = or))
    dat$order <- factor(dat$order, levels = dat$order)
  }
  
  # Make into tidy df
  gg_df <-
    dat %>%
    tidyr::pivot_longer(names_to = "K_value", values_to = "Q_value",
                        -c(order))
  
  # Build plot using helper function
  plt <- ggbarplot_helper(gg_df) + ggplot_fill
  
  # Remove legend
  if (!legend) plt <- plt + ggplot2::theme(legend.position = "none")
  
  return(plt)
}

#' Helper function for TESS barplots using ggplot
#'
#' @param dat Q matrix
#'
#' @return barplot with Q-values and individuals, colorized by K-value
#' @export
ggbarplot_helper <- function(dat) {
  dat %>%
    ggplot2::ggplot(aes(x = order, y = Q_value, fill = K_value)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_y_continuous(expand = c(0,0)) +
    ggplot2::scale_x_discrete(expand = c(0,0)) +
    ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"),
                   axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank(),
                   panel.border = ggplot2::element_rect(fill = NA, colour = "black", linetype = "solid", linewidth = 1.5),
                   strip.text.y = ggplot2::element_text(size = 30, face = "bold"),
                   strip.background = ggplot2::element_rect(colour = "white", fill = "white"),
                   panel.spacing = unit(-0.1, "lines"))
}

#' Save Q-matrices resulting from running TESS
#'
#' @param tess3_obj tess3 object returned from running TESS
#' @param dos dosage matrix used as input into TESS
#' @param save_path path to save file to; saves in pwd if not specified
#' @param dataset_name name of dataset, for file naming
#'
#' @return saves tsv file with results
#' @export
save_tess_results <- function(tess3_obj, dos, save_path = ".", dataset_name) {
  # Get best K
  bestK <- tess3_obj[["K"]]
  
  # Get Q values for best K
  qmat <- tess3r::qmatrix(tess3_obj, K = bestK)
  
  # Retrieve sample IDs from dosage matrix
  dos_names <- as.data.frame(rownames(dos)) %>% 
    rename(Bioinformatics_ID = `rownames(dos)`)
  
  dat <- cbind(dos_names, qmat) %>% 
    dplyr::mutate(K_value = bestK)
  
  write_tsv(dat, paste0(save_path, "/", dataset_name, "_TESS3_results_K", bestK, ".txt"))
  
  return(dat)
}
