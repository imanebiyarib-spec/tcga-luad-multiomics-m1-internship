#' =============================================================================
#' Project:  Genomics Project - Phase 1 (LUAD)
#' Task:     Supervised Driver-Gene Classification & Model Diagnostics
#' Author:   Imane BIYAR
#' Date:     July 2026
#'
#' Objective:
#' Train three independent classifiers (Random Forest, Elastic Net GLM, 
#' and Native XGBoost) on the Tier 2 candidate pool using a strict 
#' LUAD ground-truth compendium. Generate feature importance, global/local 
#' SHAP analysis, and export the top 'novel' candidate driver panels.
#'
#' Inputs:
#'   - Results/Integration/master_table_pathways.rds
#'   - Data/Labels/Census_allFri Jul 31 23_43_35_2026.csv (COSMIC)
#'   - Data/Labels/oncokb_biomarker.tsv (OncoKB)
#'
#' Outputs:
#'   - Results/ML/*.rds (Saved ML models to avoid retraining)
#'   - Results/ML/*.csv (Metrics, coefficients, SHAP importance)
#'   - Results/ML/*.png (Diagnostic plots, PCA biplots, SHAP summaries)
#'   - Results/Integration/*.csv (Ranked candidate panels)
#' =============================================================================

library(readxl)  
library(caret)
library(corrplot)
library(ggrepel)
library(ranger)
library(randomForest)
library(glmnet)
library(xgboost)
library(PRROC)
library(pROC)
library(dplyr)
library(tidyr)
library(ggplot2)
library(car)               
library(ResourceSelection) 
library(SHAPforxgboost)    
library(DiagrammeR)

dir.create("Results/ML", showWarnings = FALSE, recursive = TRUE)
dir.create("Results/Integration", showWarnings = FALSE, recursive = TRUE)

set.seed(1) 

# =============================================================================
# 7.1 - Load Master Table & Ground-Truth Labels 
# =============================================================================
master_table <- readRDS("Results/Integration/master_table_pathways.rds")

# 1. COSMIC Cancer Gene Census: Pan-lung somatic drivers
cgc_all <- read.csv("Data/Labels/Census_allFri Jul 31 23_43_35_2026.csv", stringsAsFactors = FALSE)
positive_cgc <- cgc_all %>%
  filter(grepl("lung adenocarcinoma|NSCLC|non-small cell lung", Tumour.Types.Somatic., ignore.case = TRUE)) %>%
  pull(Gene.Symbol) %>%
  unique()

# 2. OncoKB: FDA-recognized Precision Oncology Knowledge Base
oncokb_path <- list.files("Data/Labels", pattern = "oncokb.*\\.tsv|.*biomarker.*\\.tsv", full.names = TRUE)[1]

if (!is.na(oncokb_path) && file.exists(oncokb_path)) {
  oncokb_all <- read.delim(oncokb_path, sep = "\t", check.names = TRUE)
  positive_oncokb <- oncokb_all %>%
    filter(grepl("Lung|LUAD|NSCLC|Non-Small Cell Lung Cancer|All Solid Tumors", Cancer.Types, ignore.case = TRUE)) %>%
    pull(Gene) %>%
    unique()
} else {
  positive_oncokb <- character(0)
}
saveRDS(positive_oncokb, file = "Results/ML/positive_oncokb.rds")

# COMBINE STRICT LUAD COMPENDIA INTO THE GOLD STANDARD
positives <- Reduce(union, list(positive_cgc, positive_oncokb))
saveRDS(positives, file = "Results/ML/positive.rds")

# =============================================================================
# 7.2 - Define Explicit Named Feature Vectors & Build Labeled Matrix
# =============================================================================
feature_cols <- c("log2FC", "padj", "mut_freq_norm", "mut_q",
                  "Pct_Amplified", "Pct_HomDel",
                  "in_CNA_gistic", "in_CNA_cnv_threshold",
                  "in_GO_pathway", "in_KEGG_pathway")

master_table <- master_table %>%
  mutate(
    in_CNA_gistic        = as.numeric(in_CNA_gistic),
    in_CNA_cnv_threshold = as.numeric(in_CNA_cnv_threshold)
  )

stopifnot(all(feature_cols %in% colnames(master_table)))

candidate_pool <- master_table$gene[master_table$in_tier2]
negative_pool  <- setdiff(master_table$gene, union(candidate_pool, positives))

set.seed(1)
neg_sample <- sample(negative_pool, size = length(positives) * 5)

feature_genes <- unique(c(intersect(candidate_pool, master_table$gene),
                          intersect(positives, master_table$gene),
                          neg_sample))

feat <- master_table[match(feature_genes, master_table$gene), feature_cols]
rownames(feat) <- feature_genes
feat$label <- factor(ifelse(feature_genes %in% positives, "driver", "non_driver"),
                     levels = c("non_driver", "driver"))

# =============================================================================
# 7.3 - Stratified Train/Test Partition & PCA Biplot Data-Shift Control
# =============================================================================
set.seed(1)
train_idx <- caret::createDataPartition(feat$label, p = 0.80, list = FALSE)
train_set <- feat[train_idx, ]
test_set  <- feat[-train_idx, ]

pca_input <- bind_rows(
  train_set[, feature_cols] %>% mutate(Split = "Train (80%)"),
  test_set[, feature_cols]  %>% mutate(Split = "Test (20%)")
)

pca_res <- prcomp(pca_input[, feature_cols], center = TRUE, scale. = TRUE)
var_exp <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

pca_df <- data.frame(
  PC1   = pca_res$x[, 1],
  PC2   = pca_res$x[, 2],
  Split = factor(pca_input$Split, levels = c("Train (80%)", "Test (20%)")),
  Label = c(as.character(train_set$label), as.character(test_set$label))
)

loadings_df <- as.data.frame(pca_res$rotation[, 1:2])
colnames(loadings_df) <- c("PC1", "PC2")
loadings_df$Feature <- rownames(loadings_df)

arrow_mult <- min(max(abs(pca_df$PC1)) / max(abs(loadings_df$PC1)), 
                  max(abs(pca_df$PC2)) / max(abs(loadings_df$PC2))) * 0.85

loadings_df$PC1_scaled <- loadings_df$PC1 * arrow_mult
loadings_df$PC2_scaled <- loadings_df$PC2 * arrow_mult

p_pca_split <- ggplot() +
  geom_point(data = pca_df, aes(x = PC1, y = PC2, color = Split, shape = Label), 
             alpha = 0.65, size = 2.5) +
  geom_segment(data = loadings_df, 
               aes(x = 0, y = 0, xend = PC1_scaled, yend = PC2_scaled),
               arrow = arrow(length = unit(0.2, "cm"), type = "closed"), 
               color = "#37474F", linewidth = 0.8, alpha = 0.9) +
  geom_text_repel(data = loadings_df, 
                  aes(x = PC1_scaled, y = PC2_scaled, label = Feature),
                  color = "#212121", fontface = "bold", size = 3.8,
                  bg.color = "white", bg.r = 0.15, box.padding = 0.4,
                  point.padding = 0.2, max.overlaps = 20) +
  scale_color_manual(values = c("Train (80%)" = "#1976D2", "Test (20%)" = "#E65100"),
                     name = "Dataset Split") +
  scale_shape_manual(values = c("driver" = 17, "non_driver" = 16),
                     name = "Biological Class") +
  theme_bw(base_size = 12) +
  labs(title = "Multi-Omics PCA Biplot: Representativeness Control & Feature Loadings",
       subtitle = "Verifying absence of dataset shift and visualizing variable influence across PC1/PC2",
       x = paste0("PC1 (", var_exp[1], "% variance)"),
       y = paste0("PC2 (", var_exp[2], "% variance)")) +
  theme(plot.title = element_text(face = "bold", size = 13),
        legend.position = "right",
        panel.grid.minor = element_blank())

ggsave("Results/ML/pca_diagnostic_train_test_split.png", plot = p_pca_split, width = 11, height = 7.5, dpi = 300)

# =============================================================================
# 7.4 - STATISTICAL DIAGNOSTICS: Collinearity, VIF & Companion GLM
# =============================================================================
cor_matrix <- cor(train_set[, feature_cols], use = "pairwise.complete.obs")
write.csv(round(cor_matrix, 3), "Results/ML/feature_correlation_matrix.csv")

png("Results/ML/feature_correlation_circleplot.png", width = 1800, height = 1800, res = 200)
corrplot(cor_matrix, 
         method      = "color",           
         type        = "upper",           
         order       = "hclust",          
         addCoef.col = "black",           
         number.cex  = 0.8,               
         tl.col      = "black",          
         tl.srt      = 45,                 
         col         = colorRampPalette(c("#1976D2", "white", "#D32F2F"))(200),
         title       = "Multi-Omics Feature Collinearity Diagnostic",
         mar         = c(0, 0, 1, 0))
dev.off()

preproc <- preProcess(train_set[, feature_cols], method = c("center", "scale"))
train_scaled <- predict(preproc, train_set[, feature_cols])
train_scaled$label <- train_set$label
test_scaled  <- predict(preproc, test_set[, feature_cols])
test_scaled$label <- test_set$label

diagnostic_glm <- glm(label ~ ., data = train_scaled, family = binomial())

vif_values <- car::vif(diagnostic_glm)
write.csv(data.frame(Feature = names(vif_values), VIF = vif_values),
          "Results/ML/vif_diagnostic.csv", row.names = FALSE)

null_glm <- glm(label ~ 1, data = train_scaled, family = binomial())
mcfadden_r2 <- 1 - (as.numeric(logLik(diagnostic_glm)) / as.numeric(logLik(null_glm)))

png("Results/ML/glm_diagnostic_residual_plots.png", width = 2000, height = 1200, res = 200)
par(mfrow = c(2, 2)); plot(diagnostic_glm, which = 1:4); par(mfrow = c(1, 1))
dev.off()

hl_test <- ResourceSelection::hoslem.test(
  x = as.numeric(train_scaled$label == "driver"),
  y = fitted(diagnostic_glm), g = 10
)

glm_main <- glm(label ~ ., data = train_scaled, family = binomial())
glm_int1 <- glm(label ~ . + mut_freq_norm:in_CNA_gistic, data = train_scaled, family = binomial())
glm_int2 <- glm(label ~ . + log2FC:in_CNA_gistic, data = train_scaled, family = binomial())

# =============================================================================
# 7.5 - Cross-Validation Architecture & Custom Summary Function (with MCC)
# =============================================================================
compute_mcc <- function(pred, obs, positive = "driver") {
  cm <- table(Predicted = pred, Actual = obs)
  if (!all(c("driver", "non_driver") %in% rownames(cm))) return(NA)
  
  TP <- cm[positive, positive]
  TN <- sum(diag(cm)) - TP
  FP <- sum(cm[positive, ]) - TP
  FN <- sum(cm[, positive]) - TP
  
  num <- (TP * TN) - (FP * FN)
  den <- sqrt(as.numeric(TP + FP) * as.numeric(TP + FN) * 
                as.numeric(TN + FP) * as.numeric(TN + FN))
  
  if (den == 0) return(NA)
  return(num / den)
}

full_summary <- function(data, lev = NULL, model = NULL) {
  cm <- confusionMatrix(data$pred, data$obs, positive = "driver")
  roc_obj <- pROC::roc(response = data$obs, 
                       predictor = as.numeric(data$driver), 
                       levels = rev(lev), 
                       quiet = TRUE)
  
  mcc_val <- compute_mcc(data$pred, data$obs, positive = "driver")
  
  c(ROC              = as.numeric(pROC::auc(roc_obj)),
    Accuracy         = as.numeric(cm$overall["Accuracy"]),
    BalancedAccuracy = as.numeric(cm$byClass["Balanced Accuracy"]),
    Sensitivity      = as.numeric(cm$byClass["Sensitivity"]),
    Specificity      = as.numeric(cm$byClass["Specificity"]),
    F1               = as.numeric(cm$byClass["F1"]),
    MCC              = as.numeric(mcc_val))
}

ctrl <- trainControl(method = "repeatedcv", 
                     number = 5, 
                     repeats = 5,
                     classProbs = TRUE, 
                     summaryFunction = full_summary,
                     savePredictions = "final")

# =============================================================================
# 7.6 - Train Supervised Models: Random Forest, Elastic Net & XGBoost
# =============================================================================
imbalance_ratio <- sum(train_set$label == "non_driver") / sum(train_set$label == "driver")

if(!file.exists("Results/ML/rf_fit.rds")) {
  set.seed(1)
  rf_fit <- train(label ~ ., data = train_set, method = "rf",
                  trControl = ctrl, metric = "F1",
                  tuneGrid = expand.grid(mtry = seq(1, 11, 1)),
                  ntree = 500, 
                  classwt = c("non_driver" = 1, "driver" = imbalance_ratio))
  saveRDS(rf_fit, "Results/ML/rf_fit.rds")
}
rf_fit <- readRDS("Results/ML/rf_fit.rds")

png("Results/ML/rf_variable_importance.png", width = 1600, height = 1200, res = 200)
varImpPlot(rf_fit$finalModel, main = "Random Forest - Variable Importance")
dev.off()

train_weights <- ifelse(train_scaled$label == "driver", imbalance_ratio, 1)

if(!file.exists("Results/ML/glmnet_fit.rds")) {
  set.seed(1)
  glmnet_fit <- train(label ~ ., data = train_scaled, method = "glmnet",
                      trControl = ctrl, metric = "F1",
                      weights = train_weights, 
                      tuneGrid = expand.grid(alpha = seq(0, 1, 0.1),
                                             lambda = 10^seq(-3, 0, length = 20)))
  saveRDS(glmnet_fit, "Results/ML/glmnet_fit.rds")
}
glmnet_fit <- readRDS("Results/ML/glmnet_fit.rds")

glmnet_coefs <- as.matrix(coef(glmnet_fit$finalModel, glmnet_fit$bestTune$lambda))
coef_df <- data.frame(Variable = rownames(glmnet_coefs), Coefficient = as.numeric(glmnet_coefs)) %>%
  arrange(desc(abs(Coefficient)))
write.csv(coef_df, "Results/ML/glmnet_coefficients.csv", row.names = FALSE)

dtrain <- xgb.DMatrix(data = as.matrix(train_set[, feature_cols]), 
                      label = ifelse(train_set$label == "driver", 1, 0))
dtest  <- xgb.DMatrix(data = as.matrix(test_set[, feature_cols]), 
                      label = ifelse(test_set$label == "driver", 1, 0))

# Multi-Parameter Grid Search: Strategic Optimization (144 Combinations)
if(!file.exists("Results/ML/xgb_fit_native.rds")) {
  tune_grid <- expand.grid(
    max_depth        = c(2, 4, 6, 8),            
    eta              = c(0.01, 0.05, 0.10, 0.20),
    subsample        = c(0.6, 0.8, 1.0),         
    colsample_bytree = c(0.6, 0.8, 1.0)          
  )                                              
  
  grid_summary <- list()
  tuning_log   <- list()
  
  set.seed(1)
  for(i in 1:nrow(tune_grid)) {
    params_cv <- list(
      booster          = "gbtree",
      objective        = "binary:logistic",
      eval_metric      = "auc",
      eta              = tune_grid$eta[i],
      max_depth        = tune_grid$max_depth[i],
      subsample        = tune_grid$subsample[i],
      colsample_bytree = tune_grid$colsample_bytree[i],
      scale_pos_weight = imbalance_ratio           
    )
    
    cv_fit <- xgb.cv(
      params                = params_cv,
      data                  = dtrain,
      nrounds               = 500,
      nfold                 = 5,
      showsd                = TRUE,
      stratified            = TRUE,
      verbose               = 0,
      early_stopping_rounds = 20
    )
    
    best_iter <- cv_fit$best_iteration
    if(is.null(best_iter) || length(best_iter) == 0) {
      best_iter <- which.max(cv_fit$evaluation_log$test_auc_mean)
    }
    best_auc <- max(cv_fit$evaluation_log$test_auc_mean, na.rm = TRUE)
    
    grid_summary[[i]] <- data.frame(
      grid_id          = i,
      max_depth        = tune_grid$max_depth[i],
      eta              = tune_grid$eta[i],
      subsample        = tune_grid$subsample[i],
      colsample_bytree = tune_grid$colsample_bytree[i],
      best_iter        = best_iter,
      val_auc          = best_auc
    )
    
    tuning_log[[i]] <- cv_fit$evaluation_log %>%
      dplyr::select(iter, val_auc = test_auc_mean) %>%
      mutate(grid_id = i)
  }
  
  df_grid_results <- bind_rows(grid_summary) %>% 
    arrange(desc(val_auc)) %>%
    mutate(rank   = row_number(),
           Config = paste0("Rank #", rank, ": Depth=", max_depth, 
                           " | Eta=", eta, " | Sub=", subsample, 
                           " | Col=", colsample_bytree))
  
  best_config <- df_grid_results[1, ]
  
  xgb_params <- list(
    booster          = "gbtree",
    objective        = "binary:logistic",
    eval_metric      = "auc",
    eta              = best_config$eta,
    max_depth        = best_config$max_depth,
    subsample        = best_config$subsample,
    colsample_bytree = best_config$colsample_bytree,
    scale_pos_weight = imbalance_ratio           
  )
  
  set.seed(1)
  xgb_fit_native <- xgb.train(params = xgb_params, data = dtrain, nrounds = best_config$best_iter, verbose = 0)
  
  saveRDS(df_grid_results, "Results/ML/xgb_df_grid_results.rds")
  saveRDS(best_config, "Results/ML/xgb_best_config.rds")
  saveRDS(tuning_log, "Results/ML/xgb_tuning_log.rds")
  saveRDS(xgb_fit_native, "Results/ML/xgb_fit_native.rds")
}

df_grid_results <- readRDS("Results/ML/xgb_df_grid_results.rds")
best_config     <- readRDS("Results/ML/xgb_best_config.rds")
tuning_log      <- readRDS("Results/ML/xgb_tuning_log.rds")
xgb_fit_native  <- readRDS("Results/ML/xgb_fit_native.rds")

top_5_ids <- df_grid_results$grid_id[1:5]
plot_data <- bind_rows(tuning_log) %>% 
  filter(grid_id %in% top_5_ids) %>%
  left_join(df_grid_results %>% dplyr::select(grid_id, Config), by = "grid_id")

p_tuning <- ggplot(plot_data, aes(x = iter, y = val_auc, color = Config)) +
  geom_line(linewidth = 1) +
  labs(
    title    = "XGBoost Grid Search: Validation AUC Trajectories for Top 5 Configurations",
    subtitle = paste0("Global Supremum: AUC = ", round(best_config$val_auc, 4), 
                      " at nrounds = ", best_config$best_iter, 
                      " (Depth=", best_config$max_depth, ", Eta=", best_config$eta, ")"),
    x        = "Number of Boosting Iterations (nrounds)",
    y        = "Cross-Validated AUC (5-Fold)",
    color    = "Hyperparameter Configuration"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", legend.direction = "vertical",
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

ggsave("Results/ML/XGBoost_tuning_curve.png", plot = p_tuning, width = 11, height = 7, dpi = 300)

# Native Computation of the 25 CV Resamples 
xgb_params <- list(
  booster          = "gbtree",
  objective        = "binary:logistic",
  eval_metric      = "auc",
  eta              = best_config$eta,
  max_depth        = best_config$max_depth,
  subsample        = best_config$subsample,
  colsample_bytree = best_config$colsample_bytree,
  scale_pos_weight = imbalance_ratio           
)

if(!file.exists("Results/ML/xgb_cv_resamples.rds")) {
  set.seed(1)
  cv_folds <- caret::createMultiFolds(train_set$label, k = 5, times = 5)
  
  xgb_metrics_list <- lapply(names(cv_folds), function(fold_name) {
    train_idx_fold <- cv_folds[[fold_name]]
    
    dtrain_fold <- xgb.DMatrix(as.matrix(train_set[train_idx_fold, feature_cols]), 
                               label = ifelse(train_set$label[train_idx_fold] == "driver", 1, 0))
    dval_fold   <- xgb.DMatrix(as.matrix(train_set[-train_idx_fold, feature_cols]), 
                               label = ifelse(train_set$label[-train_idx_fold] == "driver", 1, 0))
    
    bst_fold <- xgb.train(params = xgb_params, data = dtrain_fold, nrounds = best_config$best_iter, verbose = 0)
    
    prob_val <- predict(bst_fold, dval_fold)
    pred_val <- factor(ifelse(prob_val > 0.5, "driver", "non_driver"), levels = c("non_driver", "driver"))
    
    eval_data <- data.frame(
      pred       = pred_val,
      obs        = train_set$label[-train_idx_fold],
      driver     = prob_val,
      non_driver = 1 - prob_val
    )
    
    metrics <- full_summary(eval_data, lev = levels(train_set$label))
    df_metrics <- as.data.frame(as.list(metrics))
    df_metrics$Resample <- fold_name
    return(df_metrics)
  })
  saveRDS(xgb_metrics_list, "Results/ML/xgb_cv_resamples.rds")
}
xgb_metrics_list <- readRDS("Results/ML/xgb_cv_resamples.rds")
xgb_fit <- list(resample = bind_rows(xgb_metrics_list))

# XGBoost Diagnostics: Feature Importance
xgb_importance <- xgb.importance(feature_names = feature_cols, model = xgb_fit_native)

png("Results/ML/xgboost_feature_importance.png", width = 1600, height = 1200, res = 200)
xgb.plot.importance(xgb_importance, 
                    top_n = 10, 
                    measure = "Gain", 
                    main = "XGBoost - Top Features by Gain (Predictive Contribution)")
dev.off()

# =============================================================================
# 7.7 - Cross-Validation Synthesis & Visual Comparison
# =============================================================================
summarise_resamples <- function(resample_df, model_name) {
  resample_df %>%
    dplyr::select(-Resample) %>%
    summarise(across(everything(), list(mean = mean, sd = sd, var = var), na.rm = TRUE)) %>%
    pivot_longer(everything(), names_to = c("Metric", "Stat"), names_sep = "_(?=[^_]+$)") %>%
    pivot_wider(names_from = Stat, values_from = value) %>%
    mutate(Model = model_name, .before = 1)
}

rf_summary     <- summarise_resamples(rf_fit$resample,     "Random Forest")
glmnet_summary <- summarise_resamples(glmnet_fit$resample, "Elastic Net")
xgb_summary    <- summarise_resamples(xgb_fit$resample,    "XGBoost")

combined_cv_summary <- bind_rows(rf_summary, glmnet_summary, xgb_summary)
write.csv(combined_cv_summary, "Results/ML/cv_summary_all_models.csv", row.names = FALSE)

per_fold_long <- bind_rows(
  rf_fit$resample     %>% mutate(Model = "Random Forest"),
  glmnet_fit$resample %>% mutate(Model = "Elastic Net"),
  xgb_fit$resample    %>% mutate(Model = "XGBoost")
) %>% pivot_longer(cols = c(ROC, Accuracy, BalancedAccuracy, Sensitivity, Specificity, F1, MCC),
                   names_to = "Metric", values_to = "Value")

p_cv_box <- ggplot(per_fold_long, aes(x = Metric, y = Value, fill = Model)) +
  geom_boxplot(position = position_dodge(0.8), outlier.size = 0.8, alpha = 0.85) +
  scale_fill_manual(values = c("Random Forest" = "#1976D2", "Elastic Net" = "#D32F2F", "XGBoost" = "#388E3C")) +
  theme_bw(base_size = 12) +
  labs(title = "Cross-Validated Multi-Omics Classification Performance",
       subtitle = "Comparison across 25 repeated cross-validation resamples (5-fold x 5-repeat)",
       y = "Metric Score", x = NULL) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, face = "bold"),
        legend.position = "top")
ggsave("Results/ML/cv_metrics_boxplot_all_models.png", plot = p_cv_box, width = 11, height = 6, dpi = 300)

# =============================================================================
# 7.8 - Held-Out Test Set Evaluation (Confirmatory Check)
# =============================================================================
eval_test <- function(fit, test_data, model_name, prob_data = NULL) {
  if (inherits(fit, "xgb.Booster")) {
    prob <- predict(fit, test_data)
    pred <- factor(ifelse(prob > 0.5, "driver", "non_driver"), levels = c("non_driver", "driver"))
  } else {
    pred <- predict(fit, test_data)
    prob <- if(is.null(prob_data)) predict(fit, test_data, type = "prob")[, "driver"] else prob_data
  }
  
  cm   <- confusionMatrix(pred, test_set$label, positive = "driver")
  roc_val <- as.numeric(pROC::roc(test_set$label, prob, levels = rev(levels(test_set$label)), quiet = TRUE)$auc)
  mcc_val <- compute_mcc(pred, test_set$label, positive = "driver")
  
  data.frame(
    Model            = model_name,
    ROC              = roc_val,
    Accuracy         = cm$overall["Accuracy"],
    BalancedAccuracy = cm$byClass["Balanced Accuracy"],
    Sensitivity      = cm$byClass["Sensitivity"],
    Specificity      = cm$byClass["Specificity"],
    F1               = cm$byClass["F1"],
    MCC              = mcc_val,
    row.names        = NULL
  )
}

test_summary_all <- bind_rows(
  eval_test(rf_fit,         test_set,    "Random Forest"),
  eval_test(glmnet_fit,     test_scaled, "Elastic Net"),
  eval_test(xgb_fit_native, dtest,       "XGBoost")
)
write.csv(test_summary_all, "Results/ML/held_out_test_summary_all_models.csv", row.names = FALSE)

png("Results/ML/roc_curves_all_models.png", width = 1400, height = 1400, res = 200)
roc_rf  <- pROC::roc(test_set$label, predict(rf_fit, test_set, type = "prob")[, "driver"], quiet = TRUE)
roc_en  <- pROC::roc(test_set$label, predict(glmnet_fit, test_scaled, type = "prob")[, "driver"], quiet = TRUE)
roc_xgb <- pROC::roc(test_set$label, predict(xgb_fit_native, dtest), quiet = TRUE)

plot(roc_rf, col = "#1976D2", lwd = 2.5, main = "ROC Curve Comparison - Held-Out Test Set")
plot(roc_en, col = "#D32F2F", lwd = 2.5, add = TRUE)
plot(roc_xgb, col = "#388E3C", lwd = 2.5, add = TRUE)
legend("bottomright", legend = c(paste0("Random Forest (AUC = ", round(roc_rf$auc, 3), ")"),
                                 paste0("XGBoost (AUC = ", round(roc_xgb$auc, 3), ")"),
                                 paste0("Elastic Net (AUC = ", round(roc_en$auc, 3), ")")),
       col = c("#1976D2", "#388E3C", "#D32F2F"), lwd = 2.5, bty = "n")
dev.off()

# =============================================================================
# 7.9 - Final Driver Ranking: Rescuing Candidates from Tier 2 Pool - RF 
# =============================================================================
tier2_genes <- master_table$gene[master_table$in_tier2]
tier2_feat  <- master_table[match(tier2_genes, master_table$gene), c("gene", feature_cols)]

tier2_rf_prob <- predict(rf_fit, tier2_feat[, feature_cols], type = "prob")[, "driver"]
names(tier2_rf_prob) <- tier2_genes

ranked_panel_df <- data.frame(
  gene                 = names(tier2_rf_prob),
  rf_probability       = as.numeric(tier2_rf_prob),
  already_known_driver = names(tier2_rf_prob) %in% positives
) %>% arrange(desc(rf_probability))

write.csv(ranked_panel_df, "Results/Integration/tier2_ranked_by_RF.csv", row.names = FALSE)
write.csv(head(ranked_panel_df, 10),  "Results/Integration/top_10_drivers.csv",  row.names = FALSE)
write.csv(head(ranked_panel_df, 50),  "Results/Integration/top_50_drivers.csv",  row.names = FALSE)
write.csv(head(ranked_panel_df, 100), "Results/Integration/top_100_drivers.csv", row.names = FALSE)

top_10_novel <- ranked_panel_df %>% filter(!already_known_driver) %>% head(10)
write.csv(top_10_novel, "Results/Integration/top_10_novel_for_literature_validation.csv", row.names = FALSE)
saveRDS(head(ranked_panel_df, 50)$gene, "Results/Integration/top_50_drivers.rds")

# =============================================================================
# 7.10 - Final Driver Ranking: Rescuing Candidates from Tier 2 Pool - XGBoost
# =============================================================================
dmatrix_tier2 <- xgb.DMatrix(data = as.matrix(tier2_feat[, feature_cols]))

tier2_xgb_prob <- predict(xgb_fit_native, dmatrix_tier2)
names(tier2_xgb_prob) <- tier2_genes

ranked_panel_xgb_df <- data.frame(
  gene                 = names(tier2_xgb_prob),
  xgb_probability      = as.numeric(tier2_xgb_prob),
  already_known_driver = names(tier2_xgb_prob) %in% positives
) %>% arrange(desc(xgb_probability))

write.csv(ranked_panel_xgb_df, "Results/Integration/tier2_ranked_by_XGBoost.csv", row.names = FALSE)
write.csv(head(ranked_panel_xgb_df, 10),  "Results/Integration/top_10_drivers_XGBoost.csv",  row.names = FALSE)
write.csv(head(ranked_panel_xgb_df, 50),  "Results/Integration/top_50_drivers_XGBoost.csv",  row.names = FALSE)
write.csv(head(ranked_panel_xgb_df, 100), "Results/Integration/top_100_drivers_XGBoost.csv", row.names = FALSE)

top_10_novel_xgb <- ranked_panel_xgb_df %>% filter(!already_known_driver) %>% head(10)
write.csv(top_10_novel_xgb, "Results/Integration/top_10_novel_for_literature_validation_XGBoost.csv", row.names = FALSE)
saveRDS(head(ranked_panel_xgb_df, 50)$gene, "Results/Integration/top_50_drivers_XGBoost.rds")

# =============================================================================
# 7.11 - Explainable AI (XAI): SHAP Analysis for Native XGBoost (Global)
# =============================================================================
shap_data <- as.matrix(train_set[, feature_cols])

shap_values_xgb <- shap.values(xgb_model = xgb_fit_native, X_train = shap_data)

shap_importance_xgb_df <- data.frame(
  Feature = names(shap_values_xgb$mean_shap_score), 
  Mean_Absolute_SHAP = as.numeric(shap_values_xgb$mean_shap_score)
)
write.csv(shap_importance_xgb_df, "Results/ML/xgboost_shap_importance_ranking.csv", row.names = FALSE)

shap_long <- shap.prep(xgb_model = xgb_fit_native, X_train = shap_data)

p_shap_summary <- shap.plot.summary(shap_long) +
  scale_color_gradient(low = "#1976D2", high = "#D32F2F", name = "Feature value") +
  theme_bw(base_size = 12) +
  labs(title = "SHAP Summary Plot: XGBoost Driver Classification",
       subtitle = "Visualizing omics feature impact and directionality on predictions",
       x = "Global Importance (Mean Absolute Impact)") +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "right")

ggsave("Results/ML/xgboost_shap_summary.png", plot = p_shap_summary, width = 10, height = 7, dpi = 300)

# =============================================================================
# 7.12 - SHAP Local Explanations: Top 10 Novel & Top 10 Overall XGBoost Panels
# =============================================================================
shap_data_tier2 <- as.matrix(tier2_feat[, feature_cols])
shap_values_tier2 <- shap.values(xgb_model = xgb_fit_native, X_train = shap_data_tier2)

shap_contributions_tier2 <- as.data.frame(shap_values_tier2$shap_score)
shap_contributions_tier2$gene <- tier2_genes

plot_local_shap_xgb <- function(gene_name, subfolder_prefix) {
  gene_explanation <- shap_contributions_tier2 %>%
    filter(gene == gene_name) %>%
    dplyr::select(-gene) %>%
    pivot_longer(cols = everything(), names_to = "Feature", values_to = "SHAP_Value") %>%
    arrange(SHAP_Value) %>%
    mutate(Direction = ifelse(SHAP_Value > 0, "Pro-Driver (+)", "Anti-Driver (-)"))
  
  p <- ggplot(gene_explanation, aes(x = reorder(Feature, SHAP_Value), y = SHAP_Value, fill = Direction)) +
    geom_col(color = "black", alpha = 0.85) +
    coord_flip() +
    scale_fill_manual(values = c("Pro-Driver (+)" = "#D32F2F", "Anti-Driver (-)" = "#1976D2")) +
    theme_bw(base_size = 12) +
    labs(title = paste("XGBoost Local SHAP Profile:", gene_name),
         subtitle = "Feature contributions driving the model's prediction",
         x = "Omics Feature", y = "SHAP Value (Contribution)") +
    theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
  
  dir.create(paste0("Results/ML/", subfolder_prefix), showWarnings = FALSE, recursive = TRUE)
  ggsave(paste0("Results/ML/", subfolder_prefix, "/", gene_name, "_SHAP_profile.png"), plot = p, width = 8, height = 5, dpi = 300)
}

top_10_xgb_novel_exps <- shap_contributions_tier2 %>%
  filter(gene %in% top_10_novel_xgb$gene)
write.csv(top_10_xgb_novel_exps, "Results/ML/xgboost_shap_top10_novel_explanations.csv", row.names = FALSE)

for (g in top_10_novel_xgb$gene) {
  plot_local_shap_xgb(gene_name = g, subfolder_prefix = "Top10_Novel_Shap_Plots")
}

top_10_overall_genes <- head(ranked_panel_xgb_df, 10)$gene
top_10_xgb_overall_exps <- shap_contributions_tier2 %>%
  filter(gene %in% top_10_overall_genes)
write.csv(top_10_xgb_overall_exps, "Results/ML/xgboost_shap_top10_overall_explanations.csv", row.names = FALSE)

for (g in top_10_overall_genes) {
  plot_local_shap_xgb(gene_name = g, subfolder_prefix = "Top10_Overall_Shap_Plots")
}