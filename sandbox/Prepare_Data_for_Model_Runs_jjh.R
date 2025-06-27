##Prepare observation data (DSCS presence/absence) for modeling

## Load libraries
library(terra)
library(fuzzySim)

## Define data directories
root_dir <- file.path('D:', 'dscs_nerto2025', 'data')
#root_dir <- file.path('C:', '_BioGeoProjects', 'NortheastCanyons_DSCS_Modeling')

## Read merged observation data into data frame
df <- read.csv(file = file.path(root_dir, 
                                #'Analysis', 
                                #'Merged_Samples', 
                                #'prediction_predictor_rasters',
                                #'data',
                                'dat_mod_final.csv'), 
               header = TRUE)

## Define set of potential environmental covariates
predictors <- list.files(path = file.path(root_dir, 
                                          #'Analysis', 
                                          #'Predictor_Development', 
                                          #'Predictor_Geotiffs',
                                          'predictor_rasters',
                                          'adjsd_bpi'), 
                         pattern = '.tif')

## Extract covariate values at each sample location and add to data frame
for (i in 1:length(predictors)) {
  r <- rast(x = file.path(root_dir, 
                          #'Analysis', 
                          #'Predictor_Development', 
                          #'Predictor_Geotiffs', 
                          'predictor_rasters',
                          'adjsd_bpi',
                          predictors[i]))
  
  df[, gsub(pattern = '.tif', replacement = '', 
            x = predictors[i])] <- extract(x = r, y = df[, c('x', 'y')], 
                                           method = 'simple', ID = FALSE)
}

## Select among correlated variables based on VIF
corr_analysis <- corSelect(data = df, 
                           var.cols = which(names(df) %in% gsub(pattern = '.tif', replacement = '', x = predictors)), 
                           cor.thresh = 0.7, method = 'spearman', select = 'VIF')

predictors <- paste(corr_analysis$selected.vars, '.tif', sep = '')

write.csv(df, "ne_canyons_dscs_adjsd_bpi_20250627.csv")

