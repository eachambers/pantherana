#' Imports Q matrices, output files, and IDs (from corresponding ADMIXTURE *.fam file)
#'
#' @param path path to Q matrices, .fam file, and .out files (defaults to current wd)
#' @param prefix file name (without suffix!)
#' @param K_values list of K-values to import Q matrices for; CV error will be imported for all K values that are in dir
#'
#' @return
#' @export
import_admix_data <- function(path = ".", prefix, K_values){
  ids <- read_table(paste0(path, "/", prefix, ".fam"), col_names = FALSE) %>% 
    dplyr::select(X2) %>%
    pull(X2)
  
  files <- intersect(list.files(path = path, pattern = ".out", full.names = TRUE),
                     list.files(path = path, pattern = prefix, full.names = TRUE))
  shortfiles <- intersect(list.files(path = path, pattern = ".out", full.names = FALSE),
                          list.files(path = path, pattern = prefix, full.names = FALSE))

  cv_scores <-
    1:length(files) %>% 
    lapply(function(x) {
      cv <- readLines(con = files[[x]])
      cv <- cv[grepl("^CV error", cv)]
      cv <- as.data.frame(cv) %>% dplyr::mutate(filename = shortfiles[[x]])
      return(cv)
    }) %>% 
    dplyr::bind_rows() %>% 
    tidyr::separate(cv, sep = ": ", into = c("temp", "cv")) %>% 
    tidyr::separate(filename, sep = "\\.", into = c("temp2", "temp3", "K", "rep", "temp4")) %>% 
    dplyr::mutate(kval = readr::parse_number(temp),
                  rep = readr::parse_number(rep)) %>%
    dplyr::select(-c(temp, temp2, temp3, temp4, kval)) %>% 
    dplyr::mutate(cv = as.numeric(cv),
                  rep = as.numeric(rep),
                  K = as.numeric(K))
  
  # Extract rep with minimum CV error for each K
  rep_min_cv <-
    cv_scores %>% 
    dplyr::group_by(K) %>% 
    dplyr::mutate(rep_min_cv = rep[which.min(cv)]) %>% 
    dplyr::filter(K %in% K_values)
  
  # Import Q matrices for specified K values and rep with min CV error
  dat <- pmap(tibble(path = path, prefix = prefix, rep_min_cv %>% dplyr::select(K, rep_min_cv) %>% rename(K_value = K, rep = rep_min_cv) %>% distinct() %>% arrange(K_value)), 
              import_admix_helper, 
              ids = ids)
  names(dat) <- paste0("K", K_values)
  
  return(list(cv_scores = cv_scores, dat = dat))
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
import_admix_helper <- function(path, prefix, K_value, rep, ids){
  dat <-
    read_delim(paste0(path, "/", prefix, ".", K_value, ".rep", rep, ".Q"), col_names = FALSE) %>%
    rowwise() %>% 
    mutate(max_K_val = max(c_across(everything()))) %>% 
    rowwise() %>% 
    mutate(max_K = colnames(.)[which.max(c_across())],
           K_val = K_value) %>% 
    cbind(., ids) %>% 
    rename(Bioinformatics_ID = ids) %>% 
    rownames_to_column(var = "order")
}

#' Plots CV error results for range of K-values
#'
#' @param cv_scores cv_scores object from `import_admix_data()` function output
#' @param faceted if FALSE (default), only a single dataset provided; if TRUE, will facet based on dataset name
#' @param bestk best K value (will add vertical red line)
#' @param hilite range of K values to highlight
#'
#' @return
#' @export
plot_cv_error <- function(cv_scores, bestk, hilite, faceted = FALSE){
  cv_scores %>% 
    arrange(K, rep) %>% 
    group_by(K) %>% 
    dplyr::mutate(mean = mean(cv)) %>% 
    dplyr::mutate(sd = sd(cv)) %>% 
    dplyr::ungroup() %>% 
    ggplot2::ggplot(aes(x = K, y = mean)) +
    annotate("rect", xmin = min(hilite), xmax = max(hilite), ymin = 0, ymax = max(cv_scores$cv),
             alpha = 0.5, fill = "#88a1c6") +
    ggplot2::geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), width = 0.2, color = "darkgrey") +
    ggplot2::geom_point() +
    ggplot2::geom_line(aes(group = 1)) +
    scale_x_continuous(breaks = unique(cv_scores$K)) +
    ylab("Cross-validation error") +
    xlab("K value") +
    theme(panel.grid.major.x = element_line(color = "gray95", size = 0.5)) +
    geom_vline(color = "red", xintercept = bestk)
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
#' @param output_path if write_output = TRUE, name to save output file a
#'
#' @return
#' @export
build_str_plot <- function(dat, K_value, order = "max_K", man_order = NULL, kcols, 
                           xaxis_labels = TRUE, write_output = TRUE, metadata_path = NULL, 
                           output_path = NULL, export_plot = NULL){
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
    # dat <- dat %>% 
    #   dplyr::arrange(factor(order, levels = man_order))
    dat$order <- factor(dat$order, levels = man_order)
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
    cowplot::save_plot(filename = paste0(output_path, "/", dataset_name, "_K", K_value, ".pdf"), p, base_width = 14, base_height = 8)
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
