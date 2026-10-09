##############################################################################################################
# Run machine-learning training (NNET or XGBoost) based on psbo data (without heterotrophs)
##############################################################################################################

# Load required libraries
library(adiv); library(boot); library(DescTools); library(dendextend)
library(dplyr); library(doFuture); library(doParallel); library(doRNG)
library(foreach); library(forecast); library(future); library(ggdendro)
library(ggplot2); library(ggpubr); library(ggtree); library(ggtreeDendro)
library(ggtreeExtra); library(hydroTSM); library(iNEXT); library(lubridate)
library(maps); library(Metrics); library(microseq); library(ModelMetrics)
library(neuralnet); library(NeuralNetTools); library(nnet); library(parallelDist)
library(patchwork); library(phyloseq); library(phyloseqCompanion); library(phytoclass)
library(progressr); library(readxl); library(reshape2); library(scales)
library(stringr); library(tidyr); library(tidyverse); library(vegan)
library(xgboost); library(tidymodels); library(bonsai); library(purrr); library(caret)

# Transformation function
replace_zero_halfmin <- function(x){
  min_nonzero <- min(x[x > 0], na.rm = TRUE)
  x[x == 0] <- min_nonzero / 2
  return(x)
}

# Data Import and formatting
###########################################################################

data.pft.pig <- read.csv("~/Desktop/PETRIMED/DATA/PIGMENTS/sola_pigments_pst_ptot.csv")
pal.pft <- readRDS("~/Desktop/PETRIMED/DATA/PIGMENTS/pst_palette.RData")
pal.pft <- pal.pft[names(pal.pft) %in% colnames(data.pft.pig)]

info      <- colnames(data.pft.pig)[1:7]
pigments  <- colnames(data.pft.pig)[30:ncol(data.pft.pig)]
pfts      <- names(pal.pft)
proks     <- pfts[19:length(pfts)]

# Data preparation: transformations applied to full dataset before splitting
###########################################################################

set.seed(123456)
ds  <- data.pft.pig[data.pft.pig$SF == "PICO" & data.pft.pig$ST == "SOLA", ];dim(ds)
dss <- ds[sample(nrow(ds)),]

predicted_var <- pfts
predictors    <- pigments

# Step 1: zero -> half-min on pigments and PFT abundances
dss[predictors]    <- lapply(dss[predictors],    replace_zero_halfmin)
dss[predicted_var] <- lapply(dss[predicted_var], replace_zero_halfmin)

# Step 2: relative abundances on pigments (row-wise, after zero replacement)
# This is very important to export the models in more oligotrophic or eutrophic ecosystems (BBMO)
dss[predictors] <- decostand(dss[predictors], MARGIN = 1, method = "total")

# Step 3: log10 transform on PFT abundances only
# NO log10 TRANSFORMATION
#dss[predicted_var] <- log10(dss[predicted_var])

# XGBoost training loop
###########################################################################

# Loop over all PFTs
list_top_models       <- list()
list_top_models_info  <- list()
list_top_models_olden <- list()
list_ds               <- list()

for (i in predicted_var) {
  message("Training models for ", i)
  if (all(dss[[i]] == 0 | dss[[i]] == Inf , na.rm = TRUE)) { list_top_models[[i]] <- NULL; list_top_models_info[[i]] <- NA; list_top_models_olden[[i]] <- NA; list_ds[[i]] <- list(train_df = NA, valid_df = NA); next }
  
  ## 1/ Data & model prep + parametrization
  # ---- Train/Validation split
  train_idx <- createDataPartition(dss[[i]], p = 0.8, list = FALSE)
  train_df <- dss[train_idx, c(pigments, i)]
  valid_df <- dss[-train_idx, c(pigments, i)]
  
  ## 2/ Recipe and model specification
  rec <- recipe(as.formula(paste(i, "~ .")), data = train_df) %>% step_zv(all_predictors())
  n_pred <- length(predictors)
  
  xgb_spec <- boost_tree(
    trees          = tune(),
    tree_depth     = tune(),
    learn_rate     = tune(),
    loss_reduction = tune(),
    min_n          = tune(),
    sample_size    = tune(),
    mtry           = tune()
  ) |>
    set_engine("xgboost", objective = "reg:squarederror", verbose = 0) |>
    set_mode("regression")
  
  wf    <- workflow() |> add_recipe(rec) |> add_model(xgb_spec)
  folds <- vfold_cv(train_df, v = 5)
  
  ## 3/ Reduced hyperparameter grid
  n_pred <- ncol(train_df) - 1  # exclude response
  grid_tm <- expand_grid(
    trees = c(10,20),
    tree_depth = c(1,2,3),
    learn_rate = c(0.5,1,1.5,2),
    loss_reduction = c(0.1,0.5,1,2),
    sample_size = c(0.6,0.8),
    min_n = c(1,3,5),
    mtry = c(floor(n_pred / 2), n_pred)
  )
  
  ## 4/ Hyperparameter tuning via cross-validation
  tuned <- tune_grid(
    wf,
    resamples = folds,
    grid      = grid_tm,
    metrics   = metric_set(rmse, rsq),
    control   = control_grid(save_pred = TRUE)
  )
  
  # Collect CV metrics (on log10 scale, as trained)
  results <- collect_metrics(tuned) |>
    select(-.estimator) |>
    pivot_wider(
      names_from  = .metric,
      values_from = c(mean, std_err),
      values_fn   = mean   # collapse duplicates if any
    ) |>
    distinct(across(all_of(names(grid_tm))), .keep_all = TRUE) |>  # safety dedup
    mutate(model_id = paste0("model_", row_number()))
  
  if (nrow(results) == 0) {
    warning("No valid models for ", i)
    next
  }
  
  ## 5/ Refit all candidate models, evaluate on train + validation sets
  models  <- list()
  info_tr <- data.frame()
  info_cv <- data.frame()
  olden   <- list()
  
  for (k in seq_len(nrow(results))) {
    param_names <- names(grid_tm)
    params      <- results[k, param_names]
    final_wf    <- finalize_workflow(wf, params)
    fit         <- fit(final_wf, data = train_df)
    current_id  <- results$model_id[k]
    
    # Predictions on log10 scale
    pred_tr <- predict(fit, train_df)$.pred
    pred_cv <- predict(fit, valid_df)$.pred
    obs_tr  <- train_df[[i]]
    obs_cv  <- valid_df[[i]]
    
    # Back-transform to original scale for interpretable metrics
    #pred_tr <- 10^pred_tr_log 
    #pred_cv <- 10^pred_cv_log 
    #obs_tr  <- 10^obs_tr_log
    #obs_cv  <- 10^obs_cv_log  
    
    info_tr <- rbind(info_tr, data.frame(
      model_id   = current_id,
      rmse_train = Metrics::rmse(obs_tr, pred_tr),
      cor_train  = cor(obs_tr, pred_tr, use = "complete.obs")
    ))
    
    info_cv <- rbind(info_cv, data.frame(
      model_id   = current_id,
      rmse_valid = Metrics::rmse(obs_cv, pred_cv),
      cor_valid  = cor(obs_cv, pred_cv, use = "complete.obs")
    ))
    
    # Store raw XGBoost model
    xgb_fit <- extract_fit_parsnip(fit)$fit
    models[[current_id]] <- xgb_fit
    
    # Feature importance
    imp <- xgb.importance(model = xgb_fit)
    if (nrow(imp) > 0) {
      imp$Model <- current_id
      olden[[current_id]] <- imp
    }
  }
  
  # Merge CV tuning results with train/validation evaluation metrics
  param.cv <- results |>
    left_join(info_tr, by = "model_id") |>
    left_join(info_cv, by = "model_id")
  
  list_top_models[[i]]       <- models
  list_top_models_info[[i]]  <- param.cv
  list_top_models_olden[[i]] <- bind_rows(olden)
  list_ds[[i]]$train_df      <- train_df
  list_ds[[i]]$valid_df      <- valid_df
}

# Export outputs
saveRDS(list_top_models,"Desktop/PETRIMED/ANALYSES/ML/FIN/metab_pico_models.RData")
saveRDS(list_top_models_info,"Desktop/PETRIMED/ANALYSES/ML/FIN/metab_pico_models_info.RData")
saveRDS(list_top_models_olden,"Desktop/PETRIMED/ANALYSES/ML/FIN/metab_pico_models_importance.RData")
saveRDS(list_ds,"Desktop/PETRIMED/ANALYSES/ML/FIN/metab_pico_ds.RData")


