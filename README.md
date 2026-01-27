# iSPARK - Agwise Package Documentation Summary

This R package `ispark-agwise` has been fully documented using roxygen2 style comments. Documentation includes:

## Main User-Facing Functions

### Agronomic Analysis
- **`blups()`** - Fits iterative BLUPs of field clusters with fertilizer fixed effects
  - Creates spatial field clusters from point data
  - Fits linear mixed-effects model with random field intercepts
  - Performs iterative outlier removal for improved model fit
  - Returns augmented data with predictions and residuals

### Covariate Data Processing
- **`covs()`** - Extracts and processes covariate data
  - Integrates iSDA soil properties, DEM derivatives, and WorldClim climate data
  - Applies automatic transformation and scaling to all covariates
  - Outputs data and metadata files organized by processing level

### Machine Learning Pipeline
- **`data.ml()`** - Prepares machine learning dataset
  - Merges BLUPs and covariate data
  - Creates unified dataset for model training
  
- **`ml.execute()`** - Trains machine learning models with hyperparameter tuning
  - Implements gradient boosting (GBM) via H2O
  - Performs automated grid search for hyperparameter optimization
  - Supports custom response and predictor specification

## Internal Utility Functions

### Data Validation
- **`.check_schema()`** - Validates input data format and required columns
- **`.check_gee_projectid()`** - Validates Google Earth Engine project ID
- **`.check_gee_inputs()`** - Validates GEE data request format

### Path Management
- **`.setpath()`** - Creates standardized directory structure for processing levels
  - Supports: raw, inp (input), int (intermediate), out (output)

### Spatial Operations
- **`.mk_field_clusters()`** - Builds spatial field clusters via hierarchical clustering
  - Groups observations based on geographic distance
  - Useful for identifying field boundaries from plot data

### Data Transformation
- **`.skewness()`** - Calculates distribution skewness coefficient
- **`.outlier()`** - Estimates proportion of outliers (IQR method)
- **`.sspan()`** - Calculates log-scale span for identifying order-of-magnitude ranges
- **`.autotransform()`** - Applies automatic transformation and scaling
  - Intelligently chooses: log1p, signed sqrt, or no transformation
  - Applies robust or z-score scaling based on outlier prevalence
  - Returns transformation parameters for reproducible reversal
- **`.autobacktransform()`** - Reverses automatic transformations
  - Uses stored parameters to convert back to original scale

### External Data Source Functions
- **`.gadms()`** - Downloads Global Administrative Areas (GADM) boundaries
  - Defines geographic extent and coordinate reference systems
  - Can download specific countries or all countries
  
- **`.dem()`** - Extracts Digital Elevation Model (DEM) data
  - Elevation, slope, TRI (Terrain Roughness Index), TPI (Topographic Position Index)
  - Downloads from SRTM30 via geodata package
  
- **`.isda()`** - Extracts iSDA soil properties
  - 17 soil variables at 2 depths (20cm, 50cm) = 34 features
  - Includes: Al, clay, carbon, Ca, Mg, Fe, K, P, N, pH, sand, silt, etc.
  
- **`.worldclim()`** - Extracts WorldClim climate data
  - Temperature (max, average) and precipitation
  - Aggregates monthly data into annual or custom month subsets

### Configuration
- **`.max_workers()`** - Validates and limits worker thread count (max 3)

## Documentation Standard

All functions include roxygen2 documentation with:
- **`@description`** - Clear description of function purpose and workflow
- **`@param`** - Documented parameters with types and descriptions
- **`@return`** - Specification of return values and side effects
- **`@details`** - Implementation details and decision logic
- **`@section`** - Special sections for dependencies and side effects
- **`@keywords internal`** - Marking of internal utility functions
- **`@export`** - Specification of exported vs. internal functions
- **`@seealso`** - Cross-references to related functions
- **`@examples`** - Usage examples (wrapped in `\dontrun{}` where appropriate)
- **`@author`** - Author attribution

## Next Steps

To generate the `.Rd` (documentation) files and update NAMESPACE:

```R
# In R, within the package directory:
devtools::document()
```

Or use RStudio's keyboard shortcut:
- **Ctrl+Shift+D** (Windows/Linux)
- **Cmd+Shift+D** (macOS)

This will create `.Rd` files in the `man/` directory and update the `NAMESPACE` file to reflect exported functions.

## Data Pipeline

The typical workflow is:

1. **`blups()`** - Analyze field-level yield response to fertilizer inputs
2. **`covs()`** - Extract environmental covariates at field locations
3. **`data.ml()`** - Merge results into unified ML dataset
4. **`ml.execute()`** - Train predictive models

Output structure:
```
data/
├── raw/        # Downloaded rasters and administrative boundaries
├── inp/        # Cleaned input data
├── int/        # Intermediate outputs (BLUPs, covariates, ML data)
└── out/        # Final models and results
```
