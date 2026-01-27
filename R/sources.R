#' Extract Digital Elevation Model (DEM) Data
#'
#' @description
#' Downloads SRTM30 elevation data and derives terrain variables (slope, TRI, TPI) 
#' for given coordinates. Applies automatic transformation and scaling.
#'
#' @param data.input \code{data.frame} or \code{tibble}. Required. Must contain 
#'   at least \code{longitude} and \code{latitude} columns.
#'
#' @param data.path \code{character} or \code{NULL}. Optional. Root data directory. 
#'   If \code{NULL}, uses \code{tempdir()}.
#'
#' @param iso \code{character} or \code{NULL}. Optional. ISO3 country code(s). 
#'   If \code{NULL}, will be inferred from coordinates using GADM layer.
#'
#' @return \code{list} with two elements:
#'   \describe{
#'     \item{\code{data}}{data.frame with location identifiers and extracted DEM variables: 
#'       \code{elevation}, \code{slope}, \code{tri} (Terrain Roughness Index), 
#'       \code{tpi} (Topographic Position Index)}
#'     \item{\code{meta}}{data.frame with transformation and scaling parameters 
#'       for each variable (transform type, scaling method, center, scale)}
#'   }
#'
#' @details
#' Extracts four DEM-derived variables:
#' \enumerate{
#'   \item \code{elevation}: Raw elevation in meters
#'   \item \code{slope}: Terrain slope in degrees (calculated from 8 neighbors)
#'   \item \code{tri}: Terrain Roughness Index
#'   \item \code{tpi}: Topographic Position Index (position relative to neighborhood)
#' }
#'
#' All values are automatically transformed and scaled using \code{\link{.autotransform}}.
#'
#' @section Dependencies:
#' Requires packages: \pkg{geodata}, \pkg{terra}
#'
#' @keywords internal
#'
#' @seealso \code{\link{covs}}, \code{\link{.gadms}}, \code{\link{.autotransform}}
#'
.dem <- function(data.input = NULL, data.path = NULL, iso = NULL){

  require(geodata)

  # Check inputs
  .check_schema(data.input, vars = c("longitude", "latitude"))

  df <- unique(data.input[,c("longitude", "latitude")])

  # Check ISO3
  if(is.null(iso)){
    files <- list.files(file.path(data.path, "raw", "gadm"))
    iso <- sub(".*?([A-Z]+).*", "\\1", files)
    if(!length(files)){
      i <- .gadms(iso, data.path)
      iso <- terra::extract(i, df)[,2]
    }
  }

  # Set path
  dem.path <- .setpath(level = "raw", data.path = data.path)

  out <- df
  out$location <- paste(out[["longitude"]], out[["latitude"]], sep = "_")

  mout <- data.frame()

  dem <- geodata::elevation_30s(country = iso, path = dem.path)
  slp <- terra::terrain(dem, v = "slope", unit = "degrees", neighbors = 8)
  tri <- terra::terrain(dem, v = "TRI", unit = "degrees", neighbors = 8)
  tpi <- terra::terrain(dem, v = "TPI", unit = "degrees", neighbors = 8)
  d <- c(dem, slp, tri, tpi)
  names(d) <- c("elevation", "slope", "tri", "tpi")

  for (lyr in 1:terra::nlyr(d)) {
    vals <- NULL

    val <- terra::extract(d[[lyr]], df)[,2]
    vals <- c(vals, val)
    # Transform (if needed) and scale values
    vals <- .autotransform(vals)
    vals$params$variable <- names(d[[lyr]])

    out <- cbind(out, vals$transformed)
    colnames(out)[ncol(out)] <- names(d[[lyr]])

    mout <- rbind(mout, vals$params)

  }

  r <- list(data = out,
            meta = mout)

  return(r)

}


#' Download GADM Administrative Boundaries
#'
#' @description
#' Downloads Global Administrative Areas (GADM) boundaries for specified country 
#' or all countries. Used to define geographic extent and coordinate reference 
#' system for data extraction.
#'
#' @param iso \code{character} or \code{NULL}. Optional. ISO3 country code(s). 
#'   If \code{NULL}, downloads all countries (slow).
#'
#' @param data.path \code{character} or \code{NULL}. Optional. Directory to cache 
#'   GADM downloads. If \code{NULL}, uses \code{tempdir()}.
#'
#' @param level \code{numeric}. Optional. Administrative boundary level 
#'   (0 = country, 1 = state/province, 2 = county). Defaults to 0.
#'
#' @return \code{SpatVect} or \code{sf} object with GADM boundaries (typically 
#'   returned as \pkg{terra} \code{SpatVect}).
#'
#' @details
#' Downloads from the latest GADM dataset (resolution = 2) via \pkg{geodata}. 
#' Useful for:
#' \itemize{
#'   \item Setting geographic extent for data extraction
#'   \item Inferring country from coordinates
#'   \item Projecting data to consistent CRS
#' }
#'
#' @section Dependencies:
#' Requires package: \pkg{geodata}
#'
#' @keywords internal
#'
#' @seealso \code{\link{.dem}}, \code{\link{.isda}}, \code{\link{.worldclim}}
#'
#' @examples
#' \dontrun{
#' # Download Kenya boundaries
#' gadm_ken <- .gadms(iso = "KEN", data.path = "./data")
#'
#' # Download all countries (slow!)
#' gadm_all <- .gadms(data.path = "./data")
#' }
#'
.gadms <- function(iso = NULL, data.path = NULL, level = 0){

  require(geodata)

  gadm.path <- .setpath(data.path = data.path, level = "raw")

  # Get GADM
  if(is.null(iso)){
    i <- geodata::gadm(country = geodata::country_codes()[[2]], level = level, version = "latest", resolution = 2, path = gadm.path)
  } else {
    i <- geodata::gadm(country = iso, level = level, version = "latest", resolution = 2, path = gadm.path)
  }

  return(i)

}


#' Extract iSDA Soil Properties
#'
#' @description
#' Downloads and extracts 17 soil properties from the iSDA (Innovative Solutions 
#' for Decision Agriculture) database for Africa at two soil depths. Applies 
#' automatic transformation and scaling.
#'
#' @param data.input \code{data.frame} or \code{tibble}. Required. Must contain 
#'   at least \code{longitude} and \code{latitude} columns.
#'
#' @param data.path \code{character} or \code{NULL}. Optional. Root data directory. 
#'   If \code{NULL}, uses \code{tempdir()}.
#'
#' @return \code{list} with two elements:
#'   \describe{
#'     \item{\code{data}}{data.frame with location identifiers and extracted soil 
#'       properties (e.g., \code{al_20}, \code{al_50}, \code{clay_20}, \code{clay_50}, ...)}
#'     \item{\code{meta}}{data.frame with transformation and scaling parameters 
#'       for each variable}
#'   }
#'
#' @details
#' Extracts 17 soil variables at 2 soil depths (20 cm, 50 cm), yielding 34 features:
#' \itemize{
#'   \item Aluminum (al), Clay, Total Carbon (c.tot), Calcium (ca)
#'   \item Bulk Density (db.od), ECEC (ecec.f), Iron (fe), Potassium (k), Magnesium (mg)
#'   \item Total Nitrogen (n.tot), Organic Carbon (oc), Phosphorus (p), pH (ph.h2o)
#'   \item Sand, Silt, Sulfur (s), Water Holding Capacity (wpg2), Zinc (zn)
#' }
#'
#' Requires internet access to download iSDA rasters on first use. Data are cached.
#'
#' @section Dependencies:
#' Requires packages: \pkg{geodata}, \pkg{terra}
#'
#' @keywords internal
#'
#' @seealso \code{\link{covs}}, \code{\link{.dem}}, \code{\link{.worldclim}}, \code{\link{.autotransform}}
#'
.isda <- function(data.input = NULL, data.path = NULL){

  require(geodata)

  vars = c("al", "clay", "c.tot", "ca", "db.od", "ecec.f", "fe", "k", "mg",
           "n.tot", "oc", "p", "ph.h2o", "sand", "silt", "s", "wpg2", "zn")
  depths = c(20, 50)

  # Check inputs
  .check_schema(input = data.input, vars = c("longitude", "latitude"))

  df <- unique(data.input[,c("longitude", "latitude")])

  # Set path
  isda.path <- .setpath(level = "raw", data.path = data.path)

  out <- df
  out$location <- paste(out[["longitude"]], out[["latitude"]], sep = "_")

  mout <- data.frame()

  for (v in vars) {
    for (d in depths) {
      data <- geodata::soil_af_isda(var = v, depth = d, path = isda.path)

      vals <- NULL

      val <- terra::extract(data, df)[,2]
      vals <- c(vals, val)
      # Transform (if needed) and scale values
      vals <- .autotransform(vals)
      vals$params$variable <- paste(v, d, sep = "_")

      out <- cbind(out, vals$transformed)
      colnames(out)[ncol(out)] <- paste(v, d, sep = "_")

      mout <- rbind(mout, vals$params)

    }
  }

  r <- list(data = out,
            meta = mout)

  return(r)

}


#' Extract WorldClim Climate Data
#'
#' @description
#' Downloads WorldClim bioclimatic and monthly variables (temperature, average 
#' temperature, precipitation) for given coordinates. Applies automatic 
#' transformation and scaling.
#'
#' @param data.input \code{data.frame} or \code{tibble}. Required. Must contain 
#'   at least \code{longitude} and \code{latitude} columns.
#'
#' @param iso \code{character} or \code{NULL}. Optional. ISO3 country code(s). 
#'   If \code{NULL}, will be inferred from coordinates using GADM layer.
#'
#' @param months \code{numeric} or \code{NULL}. Optional. Months to include (1-12). 
#'   If \code{NULL}, uses all 12 months.
#'
#' @param data.path \code{character} or \code{NULL}. Optional. Root data directory. 
#'   If \code{NULL}, uses \code{tempdir()}.
#'
#' @return \code{list} with two elements:
#'   \describe{
#'     \item{\code{data}}{data.frame with location identifiers and extracted climate 
#'       variables (\code{tmax}, \code{tavg}, \code{prec})}
#'     \item{\code{meta}}{data.frame with transformation and scaling parameters}
#'   }
#'
#' @details
#' Extracts three climate variables:
#' \enumerate{
#'   \item \code{tmax}: Maximum temperature (mean across selected months, in °C)
#'   \item \code{tavg}: Average temperature (mean across selected months, in °C)
#'   \item \code{prec}: Precipitation (sum across selected months, in mm)
#' }
#'
#' WorldClim data are at 30-second resolution. All values are automatically 
#' transformed and scaled using \code{\link{.autotransform}}.
#'
#' @section Dependencies:
#' Requires packages: \pkg{geodata}, \pkg{terra}
#'
#' @keywords internal
#'
#' @seealso \code{\link{covs}}, \code{\link{.dem}}, \code{\link{.isda}}, \code{\link{.autotransform}}
#'
.worldclim <- function(data.input = NULL, iso = NULL, months = NULL, data.path = NULL){

  require(geodata)

  vars = c("tmax", "tavg", "prec")

  # Check inputs
  .check_schema(data.input, vars = c("longitude", "latitude"))

  df <- unique(data.input[,c("longitude", "latitude")])

  # Check ISO3
  if(is.null(iso)){
    files <- list.files(file.path(data.path, "raw", "gadm"))
    iso <- sub(".*?([A-Z]+).*", "\\1", files)
    if(!length(files)){
      i <- .gadms(iso, data.path)
      iso <- terra::extract(i, df)[,2]
    }
  }

  # Set path
  worldclim.path <- .setpath(level = "raw", data.path = data.path)

  out <- df
  out$location <- paste(out[["longitude"]], out[["latitude"]], sep = "_")

  mout <- data.frame()

  months <- if (is.null(months)) 1:12 else months

  for (v in vars) {
    data <- geodata::worldclim_country(var = v, country = iso, res = 0.5, path = worldclim.path)[[months]]

    if(v != "prec"){
      data <- terra::mean(data, na.rm = T)
    } else {
      data <- sum(data, na.rm = T)
    }

    vals <- NULL

    val <- terra::extract(data, df)[,2]
    vals <- c(vals, val)
    # Transform (if needed) and scale values
    vals <- .autotransform(vals)
    vals$params$variable <- v

    out <- cbind(out, vals$transformed)
    colnames(out)[ncol(out)] <- v

    mout <- rbind(mout, vals$params)

  }

  r <- list(data = out,
            meta = mout)

  return(r)

}
