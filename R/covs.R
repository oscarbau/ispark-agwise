#' Extract and Process Covariate Data
#'
#' @description
#' Downloads and processes covariate data including iSDA soil properties, digital 
#' elevation model (DEM) derivatives, and WorldClim climate data. Combines these sources
#' into a single covariate dataset with metadata.
#'
#' @param data.path \code{character} or \code{NULL}. Optional. Target directory where data 
#'   will be written to. If left empty, will write to \code{tempdir()}.
#'
#' @param data.input \code{data.frame} or \code{tibble}. Required. Must contain at least 
#'   the following columns:
#'   \describe{
#'     \item{\code{id}}{unique row identifier}
#'     \item{\code{longitude}}{numeric; geographic longitude (WGS84 or CRS matching the downloaded \pkg{geodata} GADM layer)}
#'     \item{\code{latitude}}{numeric; geographic latitude}
#'   }
#'
#' @param out.filename \code{character} or \code{NULL}. Optional. Output filename 
#'   (without extension). Defaults to \code{"covs"}. The function will create two files:
#'   \code{<out.filename>.csv} for covariate data and 
#'   \code{<out.filename>-metadata.csv} for metadata.
#'
#' @param iso \code{character} or \code{NULL}. Optional. ISO3 country codes. 
#'   If provided, limits data extraction to specified countries. If \code{NULL}, 
#'   will be inferred from coordinates.
#'
#' @param years \code{numeric} or \code{NULL}. Optional. Year(s) for climate data. 
#'   Currently not used in the implementation.
#'
#' @param months \code{numeric} or \code{NULL}. Optional. Month(s) for WorldClim data.
#'   Can be used to subset climate variables by month.
#'
#' @return Invisibly returns \code{NULL}. Creates two CSV files in the 
#'   \code{<data.path>/int/covs/} directory:
#'   \describe{
#'     \item{\code{<out.filename>.csv}}{Covariate data with \code{id}, \code{longitude}, 
#'       \code{latitude}, and extracted covariate variables}
#'     \item{\code{<out.filename>-metadata.csv}}{Metadata about covariate sources and 
#'       transformation parameters}
#'   }
#'
#' @details
#' This function integrates three main covariate sources:
#' \enumerate{
#'   \item \strong{iSDA Soil Data}: 17+ soil properties at two depths (20cm, 50cm)
#'   \item \strong{DEM Derivatives}: Elevation, slope, TRI, and TPI from SRTM30
#'   \item \strong{WorldClim Climate}: Bioclimatic variables and monthly data
#' }
#' 
#' Variables are automatically transformed (log, sqrt) and scaled based on skewness 
#' and outlier prevalence using \code{.autotransform()}.
#'
#' @section Side effects:
#' - Creates directory structure in \code{<data.path>/int/covs/}
#' - Downloads covariate rasters if not cached
#' - Writes two CSV files to disk
#'
#' @section Dependencies:
#' Requires packages: \pkg{dplyr}, \pkg{purrr}, and internal utility functions 
#' (\code{.check_schema()}, \code{.setpath()}, \code{.gadms()}, \code{.isda()}, 
#' \code{.dem()}, \code{.worldclim()}).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' df <- data.frame(
#'   id = 1:3,
#'   longitude = c(36.8, 36.81, 36.79),
#'   latitude = c(-1.3, -1.31, -1.29)
#' )
#'
#' covs(
#'   data.path = "./my/data",
#'   data.input = df,
#'   iso = "KEN",
#'   out.filename = "my_covs"
#' )
#' }
#'
#' @author
#' Eduardo Garcia Bendito
#'
#' @seealso \code{\link{blups}} for field-level random effects estimation
#'
covs <- function(data.path = NULL, data.input = NULL, out.filename = NULL,
                 iso = NULL, years = NULL, months = NULL){

  require(dplyr)
  require(purrr)

  # Check inputs
  .check_schema(data.input, vars = c("id", "longitude", "latitude"))

  # Set path
  covs.path <- .setpath(level = "int", data.path = data.path)

  # Get GADM
  iso <- .gadms(iso, data.path)[, 1, drop = TRUE]

  locs <- unique(data.input[,c("id", "longitude", "latitude")])

  # Data sources
  # iSDA Soil
  i <- .isda(data.input = data.input, data.path = data.path)

  # DEM
  d <- .dem(data.input = data.input, data.path = data.path, iso = iso)

  # WorldClim
  w <- .worldclim(data.input = data.input, data.path = data.path, iso = iso, months = months)

  l <- list(i, d, w)

  data_list <- compact(map(l, "data"))
  meta_list <- compact(map(l, "meta"))

  covs <- reduce(
    data_list,
    ~ left_join(.x, select(.y, -longitude, -latitude), by = "id"),
    .init = locs
  )

  mcovs <- bind_rows(meta_list)

  dir.create(file.path(covs.path, "covs"), showWarnings = FALSE, recursive = TRUE)
  if(is.null(out.filename)){
    out.filename <- "covs"
  }

  out.filename <- gsub(".csv", "", out.filename)
  write.csv(covs, file.path(covs.path, "covs", paste0(out.filename, ".csv")), row.names = FALSE)
  write.csv(mcovs, file.path(covs.path, "covs", paste0(out.filename, "-metadata.csv")), row.names = FALSE)

}
