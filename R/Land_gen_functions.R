#' unfold converts the lower diagonal elements of a matrix into a vector
#'
#' @param X is a distance matrix
#' @param scale if TRUE then matrices will be standardized (defaults to TRUE)
#' @family MMRR functions
unfold <- function(X, scale = TRUE) {
  x <- vector()
  for (i in 2:nrow(X)) x <- c(x, X[i, 1:i - 1])
  if (scale == TRUE) x <- scale(x, center = TRUE, scale = TRUE)
  return(x)
}

# convert from matrix/data.frame/sf to formatted df
coords_to_df <- function(coords) {
  if (inherits(coords, "sf")) coords <- sf::st_coordinates(coords)
  if (is.matrix(coords)) coords <- data.frame(coords)
  colnames(coords) <- c("x", "y")
  return(coords)
}

#' Make map from model
#'
#' @param gdm_model GDM model
#' @param envlayers SpatRaster or Raster* object (LAYER NAMES MUST CORRESPOND WITH GDM MODEL)
#' @param coords data frame with x and y coordinates
#' @param scl constant for rescaling variable vectors for plotting (defaults to 1)
#' @param display_axes display PC axes text, labels, and ticks (defaults to FALSE)
#' @inheritParams gdm_do_everything
#'
#' @return GDM RGB map
#'
#' @family GDM functions
#'
#' @export
build_gdm_map <- function(gdm_model, envlayers, coords, plot_vars = TRUE, scl = 1, display_axes = FALSE, quiet = FALSE) {
  # convert envlayers to SpatRaster
  if (!inherits(envlayers, "SpatRaster")) envlayers <- terra::rast(envlayers)
  
  # convert coords to df
  coords <- coords_to_df(coords)
  
  # CHECK that all of the model variables are included in the stack of environmental layers
  # Create list of environmental predictors (everything but Geographic)
  check_geo <- gdm_model$predictors == "Geographic"
  if (any(check_geo)) {
    model_vars <- gdm_model$predictors[-which(check_geo)]
  } else {
    model_vars <- gdm_model$predictors
  }
  
  # Check that model variables are included in names of envlayers
  var_check <- model_vars %in% names(envlayers)
  
  # Print error with missing layers
  if (!all(var_check)) {
    stop(paste("missing model variable(s) from raster stack:", model_vars[!var_check]))
  }
  
  # Subset envlayers to only include variables in final model
  envlayers_sub <- terra::subset(envlayers, model_vars)
  
  # CREATE MAP ----------------------------------------------------------------------------------------------------
  
  # Transform GIS layers
  # Convert envlayers to raster
  envlayers_sub <- raster::stack(envlayers_sub)
  rastTrans <- gdm::gdm.transform(gdm_model, envlayers_sub)
  rastTrans <- terra::rast(rastTrans)
  
  # Remove NA values
  rastDat <- na.omit(terra::values(rastTrans))
  
  # Run PCA
  pcaSamp <- stats::prcomp(rastDat)
  
  # Count number of layers
  n_layers <- terra::nlyr(rastTrans)
  # Max number of layers to plot is 3, so adjust n_layers accordingly
  if (n_layers > 3) {
    n_layers <- 3
  }
  
  # Check if there are only coordinate layers
  # If there are only coordinate layers (i.e., no env layers) than you cannot create the map
  if (all(names(rastTrans) %in% c("xCoord", "yCoord"))) stop("All model splines for environmental variables are zero")
  
  # Make PCA raster
  pcaRast <- terra::predict(rastTrans, pcaSamp, index = 1:n_layers)
  
  # Scale rasters to get colors (each layer will correspond with R, G, or B in the final plot)
  pcaRastRGB <- stack_to_rgb(pcaRast)
  
  # If there are fewer than 3 n_layers (e.g., <3 variables), the RGB plot won't work (because there isn't an R, G, and B)
  # To get around this, create a blank raster (i.e., a white raster), and add it to the stack
  if (n_layers < 3) {
    warning("Fewer than three non-zero coefficients provided, adding white substitute layers to RGB plot")
    # Create white raster by multiplying a layer of pcaRast by 0 and adding 255
    white_raster <- pcaRastRGB[[1]] * 0 + 255
  }
  
  # If n_layers = 2, you end up making a bivariate map
  if (n_layers == 2) {
    pcaRastRGB <- c(pcaRastRGB, white_raster)
  }
  
  # If n_layers = 1, you end up making a univariate map
  if (n_layers == 1) {
    pcaRastRGB <- c(pcaRastRGB, white_raster, white_raster)
  }
  
  # Plot raster if quiet = FALSE
  if (!quiet) terra::plotRGB(pcaRastRGB, r = 1, g = 2, b = 3)
  if (!is.null(coords) & !quiet) points(coords, cex = 1.5)
  
  # Plot variable vectors
  if (!quiet & plot_vars & (n_layers == 3)) {
    gdm_plot_vars(pcaSamp, pcaRast, pcaRastRGB, coords, x = "PC1", y = "PC2", scl = scl, display_axes = display_axes)
  }
  
  if (!quiet & plot_vars & (n_layers != 3)) {
    warning("variable vector plot is not available for model with fewer than 3 final variables, skipping...")
  }
  
  s <- list(rastTrans, pcaRastRGB, pcaSamp, pcaRast)
  names(s) <- c("rastTrans", "pcaRastRGB", "pcaSamp", "pcaRast")
  return(s)
}