#' Iterative BLUPs of Field Clusters with Fertilizer Fixed Effects
#'
#' @description
#' Fits a linear mixed-effects model of yield on fertilizer rates with a
#' random intercept for spatially-clustered fields, then iteratively removes
#' outliers (by standardized residual) to increase model fit. Returns the input
#' data augmented with predicted values, residuals, and per-field BLUPs.
#'
#' blups:
#' \itemize{
#'   \item cleans and coerces coordinates and yield units
#'   \item replaces missing fertilizer rates with zeros
#'   \item constructs treatment labels from fertilizer quantile bins
#'   \item creates a "NOT" factor summarizing N/P/K presence (Zero, N, P, K, PK, NP, NK, NPK)
#'   \item builds spatial points, forms field clusters via \code{mk_field_clusters()}
#'   \item fits \code{lme4::lmer(yield ~ N + P + K + (1|fieldNr))}
#'   \item computes field BLUPs and predicted yields
#'   \item iteratively filters observations with |z-residual| > threshold until
#'         target R^2 or max iterations is reached
#'   \item reports RMSE and final R^2
#' }
#'
#' @param data.path \code{character} or \code{NULL}. Optional. Target directory where data will be written to.
#'   If left empty, will write to \code{tmp()}.
#'
#' @param data.input \code{data.frame} or \code{tibble}. Required. Must contain at least
#'   the following columns:
#'   \describe{
#'     \item{\code{longitude}}{numeric; geographic longitude (WGS84 or CRS matching the downloaded \pkg{geodata} GADM layer)}
#'     \item{\code{latitude}}{numeric; geographic latitude}
#'     \item{\code{N_fertilizer}}{numeric; nitrogen rate (kg/ha)}
#'     \item{\code{P_fertilizer}}{numeric; phosphorus rate (kg/ha)}
#'     \item{\code{K_fertilizer}}{numeric; potassium rate (kg/ha)}
#'     \item{\code{yield}}{numeric; yield (kg/ha)}
#'   }
#'   The function removes rows with missing coordinates and coerces \code{longitude}/\code{latitude} to numeric.
#'
#' @param distance \code{numeric} or \code{NULL}. Optional. Intended clustering distance
#'   parameter for \code{mk_field_clusters()}. Not used in the provided implementation.
#'
#' @return A \code{data.frame}/\code{tibble} with the original columns plus:
#'   \describe{
#'     \item{\code{id}}{row identifier}
#'     \item{\code{treatment}}{factor-like character label formed from quantile bins of N/P/K; "N0P0K0" for zeros}
#'     \item{\code{NOT}}{character summary of nutrient presence: Zero, N, P, K, PK, NP, NK, NPK}
#'     \item{\code{fieldNr}}{factor of spatial field cluster IDs}
#'     \item{\code{field_blup}}{numeric random intercept (BLUP) for \code{fieldNr}}
#'     \item{\code{predicted_yield_blue}}{numeric; fixed-effects prediction}
#'     \item{\code{predicted_yield_blup}}{numeric; fixed + random (BLUP) prediction}
#'     \item{\code{residual}}{numeric; \code{yield - predicted_yield_blup}}
#'     \item{\code{resid_z}}{numeric; standardized residual (only present during/after filtering iterations)}
#'   }
#'
#' @details
#' \strong{Model and iteration:}
#' The function fits:
#' \preformatted{
#'   yield ~ N_fertilizer + P_fertilizer + K_fertilizer + (1 | fieldNr)
#' }
#' then extracts fixed effects and per-cluster BLUPs. It iteratively removes points
#' with standardized residual magnitude greater than \code{threshold_z = 1} until
#' \code{R^2 >= r2_target = 0.99} or \code{max_iter = 100}.
#'
#' \strong{Spatial data:}
#' It uses GADM administrative layer from \pkg{geodata} to set CRS and to define the spatial object for clustering.
#'
#' \strong{Treatment construction:}
#' N, P, K are binned by quantiles at \code{probs = c(0, 0.25, 0.5, 0.75, 1)} with \code{include.lowest = TRUE}
#' and labels \code{1:4}. The "all-zero" case is labeled \code{"N0P0K0"}.
#'
#' @section Side effects:
#' - Prints iteration logs with R² and cluster counts to the console.
#'
#' @section Dependencies:
#' Requires packages: \pkg{dplyr}, \pkg{lme4}, \pkg{Metrics}, \pkg{sf}, \pkg{terra}.
#' Also requires a user-supplied function \code{mk_field_clusters()} to produce cluster IDs from an \code{sf} point object.
#'
#' @examples
#' \dontrun{
#' # Minimal example (assuming mk_field_clusters and spatial files exist)
#' library(dplyr)
#' df <- tibble::tibble(
#'   longitude    = c(36.8, 36.81, 36.79, 36.8),
#'   latitude     = c(-1.3, -1.31, -1.29, -1.305),
#'   N_fertilizer = c(0, 30, 60, 90),
#'   P_fertilizer = c(0, 10, 0, 20),
#'   K_fertilizer = c(0, 0, 15, 15),
#'   yield        = c(1.8, 2.4, 2.1, 2.7)
#' )
#'
#' out <- blups(data.path = "./my/data/folder", data.input = df, distance = 5000, iso = "KEN")
#' dplyr::glimpse(out)
#'
#' # Inspect model fit columns
#' summary(out$predicted_yield_blup)
#' mean((out$yield - out$predicted_yield_blup)^2)^0.5  # RMSE (equivalent to Metrics::rmse)
#' }
#'
#' @author
#' Eduardo Garcia Bendito
#'
#' @seealso \code{\link[lme4]{lmer}}, \code{\link[sf]{st_as_sf}}, \code{\link[Metrics]{rmse}}
#'
#' @export

blups <- function(data.path = NULL, data.input = NULL, distance = NULL, iso = NULL){
  require(dplyr)
  require(tibble)
  require(lme4)
  require(emmeans)
  require(Metrics)
  require(sf)

  # Check inputs
  .check_schema(data.input, vars = c("longitude", "latitude", "N_fertilizer", "P_fertilizer", "K_fertilizer", "yield"))

  # Set path
  blups.path <- .setpath(level = "int", data.path = data.path)

  # Get GADM
  i <- .gadms(iso, data.path)

  data <- data.input

  # Clean inputs
  data$longitude <- as.numeric(data$longitude)
  data$latitude <- as.numeric(data$latitude)
  data <- data[!is.na(data$longitude) & !is.na(data$latitude),]
  data$id <- 1:nrow(data)
  # Correcting some wrong (?) observations
  data$yield <- ifelse(data$yield < 10 & data$yield > 0, data$yield * 1000, data$yield)
  data$yield <- data$yield/1000

  # Making NA into 0's
  data$N_fertilizer <- ifelse(is.na(data$N_fertilizer), 0, data$N_fertilizer)
  data$P_fertilizer <- ifelse(is.na(data$P_fertilizer), 0, data$P_fertilizer)
  data$K_fertilizer <- ifelse(is.na(data$K_fertilizer), 0, data$K_fertilizer)

  # Formatting "treatment" variable
  probs <- c(0, 0.25, 0.5, 0.75, 1)
  data$treatment <- paste0("N", cut(data$N_fertilizer, right = TRUE, include.lowest	= TRUE,
                                    breaks = unique(as.vector(quantile(data$N_fertilizer, probs))),
                                    labels = 1:(length(unique(as.vector(quantile(data$N_fertilizer, probs))))-1)),
                           "P", cut(data$P_fertilizer, right = TRUE, include.lowest	= TRUE,
                                    breaks = unique(as.vector(quantile(data$P_fertilizer, probs))),
                                    labels = 1:(length(unique(as.vector(quantile(data$P_fertilizer, probs))))-1)),
                           "K", cut(data$K_fertilizer, right = TRUE, include.lowest	= TRUE,
                                    breaks = unique(as.vector(quantile(data$K_fertilizer, probs))),
                                    labels = 1:(length(unique(as.vector(quantile(data$K_fertilizer, probs))))-1)))
  data$treatment <- ifelse(data$N_fertilizer == 0 & data$P_fertilizer == 0 & data$K_fertilizer == 0, "N0P0K0", data$treatment)
  data$NOT <- NA
  data$NOT <- ifelse(data$N_fertilizer == 0 & data$P_fertilizer == 0 & data$K_fertilizer == 0, "Zero", data$NOT)
  data$NOT <- ifelse(data$N_fertilizer != 0 & data$P_fertilizer == 0 & data$K_fertilizer == 0, "N", data$NOT)
  data$NOT <- ifelse(data$N_fertilizer == 0 & data$P_fertilizer != 0 & data$K_fertilizer == 0, "P", data$NOT)
  data$NOT <- ifelse(data$N_fertilizer == 0 & data$P_fertilizer == 0 & data$K_fertilizer != 0, "K", data$NOT)
  data$NOT <- ifelse(data$N_fertilizer == 0 & data$P_fertilizer != 0 & data$K_fertilizer != 0, "PK", data$NOT)
  data$NOT <- ifelse(data$N_fertilizer != 0 & data$P_fertilizer != 0 & data$K_fertilizer == 0, "NP", data$NOT)
  data$NOT <- ifelse(data$N_fertilizer != 0 & data$P_fertilizer == 0 & data$K_fertilizer != 0, "NK", data$NOT)
  data$NOT <- ifelse(data$N_fertilizer != 0 & data$P_fertilizer != 0 & data$K_fertilizer != 0, "NPK", data$NOT)

  # Create the necessary farm_id aggregating by proximity
  # Build xy sf and clusters
  xy <- data[, c("longitude", "latitude", "id")] |>
    sf::st_as_sf(coords = c("longitude", "latitude"), crs = sf::st_crs(i))
  fieldNr <- .mk_field_clusters(xy_sf = xy, hs = distance)

  # Join cluster id
  id_map <- tibble(id = xy$id, fieldNr = fieldNr)
  data <- left_join(data, id_map, by = "id")

  # Factorize the field cluster
  data <- data %>%
    mutate(fieldNr = factor(fieldNr))

  ########################################################################################################

  cat("n obs:", nrow(data),
      "| n fieldNr:", nlevels(data$fieldNr), "\n")

  # Iterative cleaning and BLUPs with random intercept for fieldNr only
  r2_target   <- 0.99
  r2_current  <- 0
  iteration   <- 1
  max_iter    <- 100
  threshold_z <- 1

  # Loop over the model fitting to maximize fittness
  repeat {
    model <- lmer(
      yield ~ N_fertilizer + P_fertilizer + K_fertilizer + (1 | fieldNr),
      data = data
    )
    fixed_eff <- fixef(model)

    field_blups <- ranef(model)$fieldNr %>%
      as.data.frame() %>%
      rownames_to_column("fieldNr") %>%
      rename(field_blup = `(Intercept)`) %>%
      mutate(fieldNr = factor(fieldNr, levels = levels(data$fieldNr)))

    data <- data %>%
      select(-any_of(c("field_blup"))) %>%
      left_join(field_blups, by = "fieldNr") %>%
      mutate(
        predicted_yield_blue = fixed_eff[1] +
          fixed_eff["N_fertilizer"] * N_fertilizer +
          fixed_eff["P_fertilizer"] * P_fertilizer +
          fixed_eff["K_fertilizer"] * K_fertilizer,
        predicted_yield_blup = predicted_yield_blue + field_blup,
        residual = yield - predicted_yield_blup
      )

    r2_current <- cor(data$yield, data$predicted_yield_blup, use = "complete.obs")^2
    cat(sprintf("Iteration %d: R² = %.3f | n = %d | clusters = %d\n",
                iteration, r2_current, nrow(data), nlevels(data$fieldNr)))

    if (r2_current >= r2_target | iteration >= max_iter) break

    data <- data %>%
      mutate(resid_z = as.numeric(scale(residual))) %>%
      filter(abs(resid_z) <= threshold_z)

    iteration <- iteration + 1
  }

  rmse_blup <- rmse(data$yield, data$predicted_yield_blup)
  cat(sprintf("Final RMSE: %.2f | Final R²: %.3f\n", rmse_blup, r2_current), "\n")

  # Back to kg/ha
  data <- data %>%
    mutate(yield = yield * 1000,
           predicted_yield_blue = predicted_yield_blue * 1000,
           predicted_yield_blup = predicted_yield_blup * 1000,
           field_blup = field_blup * 1000)
  data$location <- paste(data[["longitude"]], data[["latitude"]], sep = "_")

  dir.create(file.path(blups.path, "not"), showWarnings = FALSE, recursive = TRUE)
  write.csv(data, file.path(blups.path, "not", paste0("blups.csv")), row.names = FALSE)

  return(data)

}

# x <- blups(data.path = "../../data",
#            data.input = read.csv("../../data/inp/not/not.csv"),
#            distance = 1000, iso = "KE")
# write.csv(x, "~/ispark-agwise/ispark-data/intermediate/not/blups.csv", row.names = FALSE)
