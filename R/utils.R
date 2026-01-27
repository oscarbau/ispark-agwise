#' Check Input Data Schema
#'
#' @description
#' Validates that the input dataset is not \code{NULL} and contains required columns.
#'
#' @param input \code{data.frame}, \code{tibble}, or \code{NULL}. The input dataset to validate.
#'
#' @param vars \code{character}. Vector of required column names.
#'
#' @return Invisibly returns \code{NULL}. Raises an error if validation fails.
#'
#' @details
#' This is an internal utility function used to validate input datasets before 
#' processing. It checks:
#' \enumerate{
#'   \item Input is not \code{NULL}
#'   \item All required columns are present (case-sensitive)
#' }
#'
#' @section Error messages:
#' \itemize{
#'   \item "Please, provide an input in CSV." if \code{input} is \code{NULL}
#'   \item "Input format not supported. Please check..." if required columns are missing
#' }
#'
#' @keywords internal
#'
.check_schema <- function(input = NULL, vars = NULL){

  # Check inputs
  if(is.null(input)){
    stop("Please, provide an input in CSV.\n")
  }
  if(!all(vars %in% colnames(input))){
    stop("Input format not supported. Please check https://github.com/egbendito/ispark/blob/main/README.md for more information.\n")
  }

}


#' Set Data Path for Processing Levels
#'
#' @description
#' Creates and returns the standardized directory path for intermediate outputs 
#' organized by processing level (raw, input, intermediate, output).
#'
#' @param level \code{character}. Required. The processing level. 
#'   Must be one of: \code{"raw"}, \code{"inp"}, \code{"int"}, or \code{"out"}.
#'
#' @param data.path \code{character} or \code{NULL}. Optional. Base data directory.
#'   If \code{NULL}, uses \code{tempdir()} as the base.
#'
#' @return \code{character}. The full path to the data directory for the specified 
#'   level (creating it if necessary).
#'
#' @details
#' This function creates a standardized directory structure:
#' \itemize{
#'   \item \code{"raw"}: For downloaded raw data (e.g., DEM, soil, climate rasters)
#'   \item \code{"inp"}: For processed input data (e.g., cleaned field data)
#'   \item \code{"int"}: For intermediate outputs (e.g., BLUPs, covariates, ML data)
#'   \item \code{"out"}: For final outputs and results
#' }
#'
#' @section Side effects:
#' Creates the directory recursively if it does not exist.
#'
#' @keywords internal
#'
.setpath <- function(level = NULL, data.path = NULL){

  # Set path
  if(is.null(level)){
    stop("\nIndicate the level of the output. Can be one of: 'raw', 'inp', 'int' or 'out'.\n")
  }
  if(is.null(data.path)){
    p <- file.path(tempdir(), level)
  } else {
    p <- file.path(data.path, level)
  }

  dir.create(p, showWarnings = FALSE, recursive = TRUE)

  return(p)

}


#' Build Spatial Field Clusters
#'
#' @description
#' Groups observations into spatial clusters based on hierarchical clustering 
#' of geographic distances. Used to identify field boundaries from point data.
#'
#' @param xy_sf \code{sf} object. An \code{sf} spatial points object with coordinates.
#'
#' @param hs \code{numeric}. Vector of hierarchical clustering distance threshold(s) 
#'   (in the same units as coordinates, typically meters). The function attempts 
#'   each threshold until it produces fewer clusters than observations.
#'
#' @return \code{character}. Vector of field cluster IDs in the format \code{"field<N>"}.
#'
#' @details
#' This function uses hierarchical clustering (complete linkage) on pairwise distances 
#' to group nearby observations. It is useful for identifying distinct field locations 
#' from plot-level data in agronomic studies.
#'
#' The distance threshold \code{hs} should be based on expected field-to-field distances. 
#' For example, with a 100 m threshold, observations within 100 m of each other will 
#' tend to be clustered together.
#'
#' @section Dependencies:
#' Requires package: \pkg{sf}
#'
#' @keywords internal
#'
#' @seealso \code{\link[sf]{st_distance}} for distance calculation
#'
#' @examples
#' \dontrun{
#' library(sf)
#' coords <- data.frame(
#'   x = c(36.80, 36.81, 37.00),
#'   y = c(-1.30, -1.31, -1.50)
#' )
#' xy_sf <- sf::st_as_sf(coords, coords = c("x", "y"), crs = 4326)
#' 
#' clusters <- .mk_field_clusters(xy_sf, hs = c(1000, 5000, 10000))
#' # Returns field IDs: "field1", "field1", "field2" (example)
#' }
#'
# Helper to build robust field clusters at multiple radii
.mk_field_clusters <- function(xy_sf, hs = distance) {
  md <- units::drop_units(sf::st_distance(xy_sf))
  hc <- hclust(as.dist(md), method = "complete")
  cl <- cutree(hc, h = hs[1])
  for (h in hs) {
    cl <- cutree(hc, h = h)
    if (length(unique(cl)) < nrow(xy_sf)) return(paste0("field", cl))
  }
  paste0("field", cl)
}


#' Calculate Distribution Skewness
#'
#' @description
#' Computes the coefficient of skewness for a numeric vector, 
#' measuring asymmetry of the distribution.
#'
#' @param x \code{numeric}. A numeric vector.
#'
#' @return \code{numeric}. Skewness coefficient. Values > 0 indicate right skew,
#'   < 0 indicate left skew, and ≈ 0 indicate symmetry.
#'
#' @details
#' Calculates the mean of standardized deviations cubed. Missing values 
#' (non-finite numbers) are removed before calculation.
#'
#' @keywords internal
#'
.skewness <- function(x) {
  x <- x[is.finite(x)]
  m <- mean(x)
  s <- sd(x)
  mean(((x - m) / s)^3)
}


#' Calculate Outlier Proportion
#'
#' @description
#' Estimates the proportion of outliers in a numeric vector using 
#' the interquartile range (IQR) method.
#'
#' @param x \code{numeric}. A numeric vector.
#'
#' @return \code{numeric}. Proportion of values that are outliers 
#'   (range 0-1). Uses threshold: Q1 - 1.5×IQR to Q3 + 1.5×IQR.
#'
#' @details
#' Applies the standard boxplot outlier rule. Missing values are removed 
#' before calculation.
#'
#' @keywords internal
#'
.outlier <- function(x) {
  x <- x[is.finite(x)]
  q <- quantile(x, c(0.25, 0.75))
  iqr <- q[2] - q[1]
  mean(x < (q[1] - 1.5 * iqr) | x > (q[2] + 1.5 * iqr))
}


#' Calculate Log-Scale Span
#'
#' @description
#' Computes the span of positive numeric values on a log10 scale.
#'
#' @param x \code{numeric}. A numeric vector (only positive values used).
#'
#' @return \code{numeric}. Difference between log10(max) and log10(min). 
#'   Roughly indicates the number of orders of magnitude.
#'
#' @details
#' Useful for assessing whether data spans multiple orders of magnitude 
#' and might benefit from log transformation.
#'
#' @keywords internal
#'
.sspan <- function(x) {
  x <- x[x > 0 & is.finite(x)]
  log10(max(x)) - log10(min(x))
}


#' Automatic Variable Transformation and Scaling
#'
#' @description
#' Applies data-driven transformation and scaling based on skewness, 
#' outlier prevalence, and span. Useful for preparing covariates for 
#' regression or machine learning.
#'
#' @param x \code{numeric}. A numeric vector to transform.
#'
#' @return \code{list} with two elements:
#'   \describe{
#'     \item{\code{transformed}}{numeric; transformed and scaled values}
#'     \item{\code{params}}{list; transformation and scaling parameters for reversal}
#'   }
#'   The \code{params} list contains:
#'   \describe{
#'     \item{\code{transform}}{character; transformation applied: "log1p", "signed_sqrt", or "none"}
#'     \item{\code{scaling}}{character; scaling method: "robust" (median/IQR) or "zscore"}
#'     \item{\code{center}}{numeric; center value used in scaling}
#'     \item{\code{scale}}{numeric; scale value used in scaling}
#'   }
#'
#' @details
#' \strong{Transformation decision tree:}
#' \enumerate{
#'   \item If all positive AND span > 2 AND skewness > 1 → log1p transform
#'   \item Else if |skewness| > 1 → signed square root transform
#'   \item Else → no transformation
#' }
#'
#' \strong{Scaling decision:}
#' \itemize{
#'   \item If outlier proportion > 2\% → robust scaling (median/IQR)
#'   \item Else → z-score scaling (mean/SD)
#' }
#'
#' Use \code{\link{.autobacktransform}} to reverse the transformation.
#'
#' @keywords internal
#'
#' @seealso \code{\link{.autobacktransform}} for reversing transformations
#'
#' @examples
#' \dontrun{
#' # Right-skewed variable
#' x <- rlnorm(100, meanlog = 0, sdlog = 1)
#' result <- .autotransform(x)
#' # result$transformed will be log1p transformed and scaled
#' # Use result$params to transform new data identically
#' }
#'
.autotransform <- function(x) {
  sk <- .skewness(x)
  of <- .outlier(x)
  span <- if (all(x > 0, na.rm = TRUE)) .sspan(x) else NA

  # Step 1: transform
  if (!is.na(span) && span > 2 && sk > 1) {
    xt <- log1p(x)
    transform <- "log1p"
  } else if (abs(sk) > 1) {
    xt <- sign(x) * sqrt(abs(x))
    transform <- "signed_sqrt"
  } else {
    xt <- x
    transform <- "none"
  }

  # Step 2: scale
  if (of > 0.02) {
    xs <- (xt - median(xt, na.rm = TRUE)) / IQR(xt, na.rm = TRUE)
    scaling <- "robust"
    center <- median(xt, na.rm = TRUE)
    scale <- IQR(xt, na.rm = TRUE)
  } else {
    xs <- (xt - mean(xt, na.rm = TRUE)) / sd(xt, na.rm = TRUE)
    scaling <- "zscore"
    center <- mean(xt, na.rm = TRUE)
    scale <- sd(xt, na.rm = TRUE)
  }

  vals <- list(
    transformed = xs,
    params = list(
      transform = transform,
      scaling = scaling,
      center = center,
      scale = scale
    )
  )
  return(vals)
}


#' Reverse Automatic Transformation and Scaling
#'
#' @description
#' Reverses a transformation and scaling applied by \code{\link{.autotransform}}.
#'
#' @param x \code{numeric}. A transformed and scaled vector.
#'
#' @param params \code{list}. A parameters list from \code{\link{.autotransform}} 
#'   output containing \code{transform}, \code{scaling}, \code{center}, and \code{scale}.
#'
#' @return \code{numeric}. The original-scale values (approximately, 
#'   given rounding and numeric precision).
#'
#' @details
#' \strong{Process:}
#' \enumerate{
#'   \item Reverse scaling: \code{x * scale + center}
#'   \item Reverse transformation: Invert "log1p", "signed_sqrt", or apply none
#' }
#'
#' This is the inverse of \code{\link{.autotransform}}.
#'
#' @keywords internal
#'
#' @seealso \code{\link{.autotransform}} for forward transformation
#'
#' @examples
#' \dontrun{
#' x_orig <- rlnorm(100, meanlog = 0, sdlog = 1)
#' result <- .autotransform(x_orig)
#' x_recovered <- .autobacktransform(result$transformed, result$params)
#' # x_recovered ≈ x_orig
#' }
#'
.autobacktransform <- function(x, params) {

  ## ---- undo scaling ----
  xt <- x * params$scale + params$center

  ## ---- undo transform ----
  if (params$transform == "log1p") {
    x <- expm1(xt)
  } else if (params$transform == "signed_sqrt") {
    x <- sign(xt) * (xt^2)
  } else {
    x <- xt
  }

  return(x)

}

# .autobacktransform(x = l[[3]][["prec"]], params = meta[meta$variable == "prec", ])

