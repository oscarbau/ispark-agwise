#' Prepare Machine Learning Dataset
#'
#' @description
#' Merges BLUPs (Best Linear Unbiased Predictions) and covariate data into a 
#' single dataset suitable for machine learning modeling.
#'
#' @param data.path \code{character} or \code{NULL}. Optional. Root data directory 
#'   where intermediate outputs are stored. If \code{NULL}, reads from \code{tempdir()}.
#'   Expected to contain \code{int/not/blups.csv} and \code{int/covs/covs.csv}.
#'
#' @return Invisibly returns \code{NULL}. Creates a merged dataset written to 
#'   \code{<data.path>/int/ml/data/ml-data.csv} containing the original columns 
#'   from both BLUPs and covariates, joined by \code{id}.
#'
#' @details
#' This function performs a left join of BLUPs and covariate data, keeping all 
#' rows from BLUPs and matching covariates. Geographic coordinates (\code{longitude}, 
#' \code{latitude}) are kept from the BLUPs dataset and duplicated columns from 
#' covariates are dropped.
#'
#' \strong{Required files:}
#' \enumerate{
#'   \item \code{int/not/blups.csv} from \code{\link{blups}}
#'   \item \code{int/covs/covs.csv} from \code{\link{covs}}
#' }
#'
#' @section Side effects:
#' - Creates directory \code{<data.path>/int/ml/data/} if it does not exist
#' - Writes \code{ml-data.csv} to disk
#'
#' @section Dependencies:
#' Requires package: \pkg{dplyr}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data.ml(data.path = "./my/data")
#' # Reads: ./my/data/int/not/blups.csv and ./my/data/int/covs/covs.csv
#' # Writes: ./my/data/int/ml/data/ml-data.csv
#' }
#'
#' @author
#' Eduardo Garcia Bendito
#'
#' @seealso \code{\link{blups}}, \code{\link{covs}}, \code{\link{ml.execute}}
#'
data.ml <- function(data.path = NULL){

  require(dplyr)

  # Check required data exists
  need <- all(file.exists(file.path(data.path, "int", "not", paste0("blups.csv"))),
              file.exists(file.path(data.path, "int", "covs", paste0("covs.csv"))))
  if(!need){
    stop("The necessary files (BLUPS and Covariates) do not exists.\n")
  }

  # Set path
  ml.path <- .setpath(level = "int", data.path = data.path)

  blups <- read.csv(file.path(data.path, "int", "not", paste0("blups.csv")))
  covs <- read.csv(file.path(data.path, "int", "covs", paste0("covs.csv")))

  d <- blups %>%
    left_join(
      covs %>% select(-longitude, -latitude),
      by = "id"
    )

  dir.create(file.path(ml.path, "ml", "data"), showWarnings = FALSE, recursive = TRUE)
  write.csv(d, file.path(ml.path, "ml", "data", "ml-data.csv"), row.names = FALSE)

}

#' Train Machine Learning Models with Hyperparameter Tuning
#'
#' @description
#' Trains a gradient boosting machine (GBM) model using H2O with automated 
#' hyperparameter tuning via grid search. Splits data into training and test sets,
#' performs hyperparameter optimization, and saves the best model.
#'
#' @param data.input \code{data.frame} or \code{tibble}. Required. The ML dataset 
#'   (typically from \code{\link{data.ml}}). Must contain the response variable 
#'   and all predictor variables.
#'
#' @param response \code{character} or \code{NULL}. Optional. Name of the response 
#'   variable. Defaults to \code{"predicted_yield_blup"} if not specified.
#'
#' @param predictors \code{character} or \code{NULL}. Optional. Vector of predictor 
#'   variable names. If \code{NULL}, automatically selects all columns except 
#'   a predefined set of non-predictor columns (\code{fieldNr}, \code{id}, 
#'   \code{treatment}, \code{NOT}, \code{predicted_yield_blue}, \code{predicted_yield_blup}, 
#'   \code{field_blup}, \code{yield}, \code{location}, \code{latitude}, \code{longitude}, 
#'   \code{source}, \code{VARIETY}, \code{n}, \code{top.90}, \code{median}, \code{se}, 
#'   \code{bottom.10}, \code{geometry}).
#'
#' @details
#' \strong{Model configuration:}
#' - Algorithm: Gradient Boosting Machine (GBM)
#' - Train/test split: 70/30
#' - Hyperparameter grid:
#'   \itemize{
#'     \item \code{ntrees}: 500, 600, 700, 800, 900, 1000
#'     \item \code{max_depth}: 4, 6, 8
#'   }
#' - Cross-validation: 5-fold
#' - Random seed: 444 (for reproducibility)
#'
#' \strong{Output files:}
#' Grid search results are saved to \code{<data.path>/int/ml/models/gridGB.RDS}
#' The best model is saved to \code{<data.path>/int/ml/models/} as an H2O model object.
#'
#' @section Side effects:
#' - Initializes H2O cluster and shuts it down at the end
#' - Creates directory \code{<data.path>/int/ml/models/} if it does not exist
#' - Saves RDS and H2O model files to disk
#' - Prints grid search and model training progress to console
#'
#' @section Dependencies:
#' Requires package: \pkg{h2o}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Load prepared ML data
#' ml_data <- read.csv("./my/data/int/ml/data/ml-data.csv")
#'
#' # Train model with default response and auto-selected predictors
#' ml.execute(data.input = ml_data)
#'
#' # Train with custom response and predictors
#' ml.execute(
#'   data.input = ml_data,
#'   response = "yield",
#'   predictors = c("elevation", "slope", "temperature", "rainfall")
#' )
#' }
#'
#' @author
#' Eduardo Garcia Bendito
#'
#' @seealso \code{\link{data.ml}} for preparing input data
#'
ml.execute <- function(data.input = NULL, response = NULL, predictors = NULL){

  require(h2o)

  if(is.null(data.input)){
    stop("`data` can't be NULL. Need to provide data inputs.\n")
  }

  # Set path
  ml.path <- .setpath(level = "int", data.path = data.path)
  dir.create(file.path(ml.path, "ml", "models"), showWarnings = FALSE, recursive = TRUE)

  h2o::h2o.init()
  d.h2o <- h2o::as.h2o(data.input)

  if(is.null(response)){
    response <- "predicted_yield_blup"
  }

  if(is.null(predictors)){
    predictors <- data.input |> names()
    predictors <- predictors[!predictors %in% c("fieldNr", "id", "treatment", "NOT", "predicted_yield_blue", "predicted_yield_blup", "field_blup", "yield", "location", "latitude", "longitude", "source", "VARIETY", "n", "top.90", "median", "se", "bottom.10", "geometry")]
  }

  #create a random training-test split of our data ## should be possible to do it by missing one
  d.h2o.split <- h2o::h2o.splitFrame(data = d.h2o, ratios = 0.7, seed = 444)
  train <- d.h2o.split[[1]]
  test <- d.h2o.split[[2]]

  #  grid search to tune hyper parameters
  hyperparams.GB <- list(
    ntrees = seq(500, 1000, 100),
    max_depth = seq(4, 8, 2)
  )

  grid.GB <- h2o.grid(
    algorithm = "gbm",
    x = predictors,
    y = response,
    grid_id = "hyperparams_gbm",
    hyper_params = hyperparams.GB,
    training_frame = train,
    validation_frame = test,
    seed = 444
  )

  # Get the best hyper parameters
  # Write data to the file
  saveRDS(grid.GB, file.path(ml.path, "ml", "models", "gridGB.RDS"))
  hyperparams.GB <- h2o.getModel(grid.GB@model_ids[[1]])

  ## fit the model with the tuned hyper parameters:
  ML.GB <- h2o::h2o.gbm(x = predictors,
                        y = response,
                        ntrees = hyperparams.GB@parameters$ntrees,
                        max_depth = hyperparams.GB@parameters$max_depth,
                        training_frame = train,
                        validation_frame = test,
                        keep_cross_validation_predictions = TRUE,
                        nfolds = 5,
                        seed = 444)
  h2o::h2o.saveModel(object = ML.GB, path = file.path(ml.path, "ml", "models"), force = TRUE)

  h2o::h2o.shutdown(prompt = FALSE)

}
