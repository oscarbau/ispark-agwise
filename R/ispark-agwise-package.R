#' isparkagwise: Agwise Modelling Framework for iSPARK Project
#'
#' @description
#' The `ispark-agwise` package provides a reproducible end-to-end pipeline to generate
#' fertilizer recommendations for African countries. The pipeline can be implemented
#' for any country in Africa as it relies on Africa-only data sources (iSDA Soil).
#'
#' @details
#' ## Main Features
#'
#' - **Field-level analysis**: Estimate best linear unbiased predictions (BLUPs) of
#'   field-level response to fertilizer inputs
#' - **Environmental covariates**: Integrate soil properties (iSDA), topography (DEM),
#'   and climate (WorldClim) data
#' - **Machine learning**: Train predictive models for site-specific fertilizer recommendations
#' - **Reproducible pipeline**: Fully documented workflow from raw data to recommendations
#'
#' ## Workflow
#'
#' The typical data pipeline consists of:
#'
#' 1. **Field-level modeling** (`blups()`): Analyze yield response to N/P/K fertilizers
#' 2. **Covariate extraction** (`covs()`): Download and process environmental data
#' 3. **Data preparation** (`data.ml()`): Merge field and environmental data
#' 4. **Machine learning** (`ml.execute()`): Train site-specific recommendation models
#'
#' ## Data Sources
#'
#' - **Soil**: iSDA (Innovative Solutions for Decision Agriculture) - 17 soil properties
#'   at two depths (20cm, 50cm)
#' - **Topography**: SRTM30 digital elevation model - elevation, slope, TRI, TPI
#' - **Climate**: WorldClim - temperature (max, mean) and precipitation
#' - **Administrative**: GADM - country and regional boundaries
#'
#' ## Example
#'
#' ```r
#' library(ispark.agwise)
#'
#' # Prepare field-level data with yield and fertilizer rates
#' field_data <- read.csv("field_observations.csv")
#'
#' # Step 1: Estimate field BLUPs
#' blups_output <- blups(
#'   data.path = "./agwise_data",
#'   data.input = field_data,
#'   iso = "KEN"
#' )
#'
#' # Step 2: Extract environmental covariates
#' covariates <- covs(
#'   data.path = "./agwise_data",
#'   data.input = field_data,
#'   iso = "KEN"
#' )
#'
#' # Step 3: Prepare ML dataset
#' data.ml(data.path = "./agwise_data")
#'
#' # Step 4: Train models
#' ml.execute(
#'   data.input = read.csv("./agwise_data/int/ml/data/ml-data.csv")
#' )
#' ```
#'
#' @section Data Organization:
#' The package uses a standardized directory structure:
#' \itemize{
#'   \item `raw/`: Downloaded rasters and spatial data
#'   \item `inp/`: Cleaned input datasets
#'   \item `int/`: Intermediate outputs (BLUPs, covariates, ML data)
#'   \item `out/`: Final recommendations and model objects
#' }
#'
#' @section Dependencies:
#' Key packages used throughout the pipeline:
#' \itemize{
#'   \item \pkg{dplyr}: Data manipulation
#'   \item \pkg{lme4}: Linear mixed-effects models
#'   \item \pkg{sf}: Spatial data handling
#'   \item \pkg{terra}: Raster data processing
#'   \item \pkg{geodata}: Download elevation, soil, and climate data
#'   \item \pkg{h2o}: Machine learning models
#' }
#'
#' @docType package
#' @name ispark-agwise
#' @keywords internal
#'
#' @importFrom dplyr bind_rows compact left_join map reduce select
#' @importFrom purrr map_df
#' @importFrom lme4 lmer
#' @importFrom sf st_as_sf st_distance
#' @importFrom terra extract terrain nlyr c
#' @importFrom geodata gadm elevation_30s soil_af_isda worldclim_country country_codes
#' @importFrom h2o h2o.init h2o.as.h2o h2o.splitFrame h2o.grid h2o.gbm h2o.saveModel h2o.shutdown
#'
NULL
