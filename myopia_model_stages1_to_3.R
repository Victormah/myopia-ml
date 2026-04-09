# ==============================================================================
# Machine Learning Analysis Framework: Myopia Chronobiological Signature
# ==============================================================================
# Stage 1: Primary Discovery Model (Retina Onset Dataset)
# Stage 2: Cross-Tissue Validation (Choroid Onset Dataset)
# Stage 3: Cross-Stage Validation (Independent Progression Dataset)
# Stage 4: External Validation (another file)
# ==============================================================================

# === 0. Load Libraries ===
suppressPackageStartupMessages({
  library(edgeR)
  library(tidyverse)
  library(caret)
  library(randomForest)
  library(pROC)
  library(iml)         # For SHAP interpretation
  library(FactoMineR)  # Optional for PCA
  library(glmnet)      # For LASSO and Elastic Net
  library(Boruta)      # For Boruta feature selection
  library(e1071)       # For SVM and Naive Bayes
})

# ==============================================================================
# === 1. Helper Functions ===
# ==============================================================================

# Feature Selection Function (Boruta / LASSO)
repeat_feature_selection <- function(method = c("LASSO", "Boruta"), 
                                     train_data, pheno_data, logCPM_data = NULL,
                                     train_idx, ntime = 10, verbose = TRUE) {
  method <- match.arg(method)
  all_selected <- list()
  
  for (i in 1:ntime) {
    set.seed(i)
    # Add circular time variables
    train_data_mod <- train_data %>%
      dplyr::select(-time_class) %>%
      dplyr::mutate(
        sin_time = pheno_data$sin_time[train_idx],
        cos_time = pheno_data$cos_time[train_idx]
      )
    
    if (method == "LASSO") {
      X <- as.matrix(train_data_mod)
      y <- as.numeric(factor(pheno_data$time_class[train_idx])) 
      lasso_model <- cv.glmnet(X, y, alpha = 1)
      coef_mat <- as.matrix(coef(lasso_model, s = "lambda.min"))
      selected_genes <- rownames(coef_mat)[coef_mat[,1] != 0]
      selected_genes <- selected_genes[selected_genes != "(Intercept)"]
      
    } else if (method == "Boruta") {
      train_boruta <- train_data_mod %>%
        dplyr::mutate(time_class = as.factor(pheno_data$time_class[train_idx]))
      boruta_model <- Boruta(time_class ~ ., data = train_boruta, maxRuns = 100, doTrace = 0)
      selected_genes <- getSelectedAttributes(boruta_model, withTentative = FALSE)
    }
    
    if (verbose) cat(sprintf("[%s] Iteration %d: %d genes selected\n", method, i, length(selected_genes)))
    all_selected[[i]] <- selected_genes
  }
  
  # Flatten and tabulate frequencies
  selected_all <- unlist(all_selected)
  freq_table <- sort(table(selected_all), decreasing = TRUE)
  freq_df <- as.data.frame(freq_table)
  colnames(freq_df) <- c("Gene", "Frequency")
  
  return(freq_df)
}

# Classification Function (Nested CV Evaluation)
repeat_classification <- function(data, model_type = "RF", nround = 50, label_col = "time_class", train_ratio = 0.8) {
  Accuracy <- numeric(nround); Precision <- numeric(nround)
  Recall <- numeric(nround); F1 <- numeric(nround)
  AUC <- numeric(nround); MCC <- numeric(nround)
  
  data[[label_col]] <- as.factor(data[[label_col]])
  pos <- levels(data[[label_col]])[2] # Target: ZT_812
  neg <- levels(data[[label_col]])[1] # Target: ZT_other
  
  # Inner Loop: 5-fold CV for Hyperparameter Tuning
  inner_ctrl <- caret::trainControl(method = "cv", number = 5, classProbs = TRUE, summaryFunction = caret::twoClassSummary)
  
  for (i in 1:nround) {
    set.seed(i)
    cat(sprintf("Training %s - Round %d/%d\n", model_type, i, nround))
    
    # Outer Loop: Stratified 80/20 Resampling
    train_idx <- caret::createDataPartition(data[[label_col]], p = train_ratio, list = FALSE)
    train_data <- data[train_idx, ]
    test_data <- data[-train_idx, ]
    
    train_data[[label_col]] <- as.factor(train_data[[label_col]])
    test_data[[label_col]] <- factor(test_data[[label_col]], levels = levels(train_data[[label_col]]))
    
    # Train Models & Tune Hyperparameters via caret
    if (model_type == "RF") {
      model <- caret::train(as.formula(paste(label_col, "~ .")), data = train_data, method = "rf", metric = "ROC", trControl = inner_ctrl, tuneLength = 3, ntree = 500)
    } else if (model_type == "SVM") {
      model <- caret::train(as.formula(paste(label_col, "~ .")), data = train_data, method = "svmRadial", metric = "ROC", trControl = inner_ctrl, preProcess = c("center", "scale"), tuneLength = 5)
    } else if (model_type == "NB") {
      model <- caret::train(as.formula(paste(label_col, "~ .")), data = train_data, method = "nb", metric = "ROC", trControl = inner_ctrl, tuneLength = 3)
    } else if (model_type == "LR") {
      model <- caret::train(as.formula(paste(label_col, "~ .")), data = train_data, method = "glmnet", metric = "ROC", trControl = inner_ctrl, tuneLength = 5)
    } else {
      stop("Unsupported model type.")
    }
    
    # Predict on UNSEEN test_data
    pred <- predict(model, test_data)
    prob <- predict(model, test_data, type = "prob")[, pos]
    
    # Evaluate Metrics
    cm <- table(Predicted = pred, Actual = test_data[[label_col]])
    TP <- ifelse(pos %in% rownames(cm) & pos %in% colnames(cm), cm[pos, pos], 0)
    TN <- ifelse(neg %in% rownames(cm) & neg %in% colnames(cm), cm[neg, neg], 0)
    FP <- ifelse(pos %in% rownames(cm) & neg %in% colnames(cm), cm[pos, neg], 0)
    FN <- ifelse(neg %in% rownames(cm) & pos %in% colnames(cm), cm[neg, pos], 0)
    
    Accuracy[i] <- (TP + TN) / sum(cm)
    Precision[i] <- ifelse((TP + FP) == 0, 0, TP / (TP + FP))
    Recall[i] <- ifelse((TP + FN) == 0, 0, TP / (TP + FN))
    F1[i] <- ifelse((Precision[i] + Recall[i]) == 0, 0, 2 * (Precision[i] * Recall[i]) / (Precision[i] + Recall[i]))
    
    denom <- sqrt(as.numeric(TP + FP) * as.numeric(TP + FN) * as.numeric(TN + FP) * as.numeric(TN + FN))
    MCC[i] <- ifelse(denom == 0, 0, (TP * TN - FP * FN) / denom)
    AUC[i] <- tryCatch({ pROC::roc(test_data[[label_col]], prob, quiet = TRUE)$auc }, error = function(e) NA)
  }
  
  # Summarize Results
  df_res <- data.frame(Accuracy, Precision, Recall, F1, AUC, MCC)
  summary_stats <- df_res %>%
    summarise(
      Accuracy_Mean = mean(Accuracy, na.rm = TRUE), Accuracy_SD = sd(Accuracy, na.rm = TRUE),
      Precision_Mean = mean(Precision, na.rm = TRUE), Precision_SD = sd(Precision, na.rm = TRUE),
      Recall_Mean = mean(Recall, na.rm = TRUE), Recall_SD = sd(Recall, na.rm = TRUE),
      F1_Mean = mean(F1, na.rm = TRUE), F1_SD = sd(F1, na.rm = TRUE),
      AUC_Mean = mean(AUC, na.rm = TRUE), AUC_SD = sd(AUC, na.rm = TRUE),
      MCC_Mean = mean(MCC, na.rm = TRUE), MCC_SD = sd(MCC, na.rm = TRUE)
    ) %>%
    pivot_longer(cols = everything(), names_to = c("Metric", ".value"), names_sep = "_") %>%
    mutate(CI95 = 2 * sqrt((Mean * (1 - Mean)) / nround))
  
  return(list(Raw = df_res, Summary = summary_stats))
}

# Results Formatting Functions
combine_results <- function(...) {
  result_list <- list(...)
  method_names <- names(result_list)
  df_list <- lapply(seq_along(result_list), function(i) {
    df <- result_list[[i]]$Summary
    df$Model <- method_names[i]
    return(df)
  })
  df_combined <- do.call(rbind, df_list)
  return(df_combined[, c("Model", "Metric", "Mean", "SD", "CI95")])
}

combine_means_wide <- function(...) {
  result_list <- list(...)
  method_names <- names(result_list)
  df_list <- lapply(seq_along(result_list), function(i) {
    df <- result_list[[i]]$Summary
    df$Model <- method_names[i]
    return(df[, c("Metric", "Mean", "Model")])
  })
  df_combined <- do.call(rbind, df_list)
  df_wide <- df_combined %>%
    pivot_wider(names_from = Model, values_from = Mean) %>%
    arrange(factor(Metric, levels = c("Accuracy", "Precision", "Recall", "F1", "AUC", "MCC")))
  return(df_wide)
}


# ==============================================================================
# === 2. Load & Preprocess Primary Data (Onset) ===
# ==============================================================================
pheno <- read.csv("GSE227724_pheno_onset.csv")
counts_mat <- read.csv("GSE227724_counts_onset.csv", row.names = 1)

# Clean Column Names and Align
pheno <- pheno %>% mutate(across(everything(), ~str_trim(.)))
colnames(counts_mat) <- str_trim(colnames(counts_mat))
pheno$Sample_ID <- str_trim(pheno$Sample_ID)

common_ids <- intersect(pheno$Sample_ID, colnames(counts_mat))
pheno <- pheno %>% filter(Sample_ID %in% common_ids)
counts_mat <- counts_mat[, pheno$Sample_ID]

# Time Encoding (Circadian periodicity)
pheno <- pheno %>%
  mutate(
    Time = str_trim(time_point),
    Time_num = as.numeric(gsub("Z", "", Time)),
    Time_num = ifelse(is.na(Time_num), 0, Time_num),
    Time_rad = 2 * pi * Time_num / 24,
    sin_time = sin(Time_rad),
    cos_time = cos(Time_rad)
  )

# Create DGE, Normalize, and Log-CPM Transformation
dge <- DGEList(counts = counts_mat)
keep <- filterByExpr(dge, group = pheno$eye_treatment)
dge <- dge[keep, , keep.lib.sizes = TRUE]
dge <- calcNormFactors(dge)

logCPM <- cpm(dge, log = TRUE, prior.count = 1)
logCPM <- t(logCPM)
logCPM_df <- as.data.frame(logCPM)

logCPM_df$eye_treatment <- factor(pheno$eye_treatment)
logCPM_df$sin_time <- pheno$sin_time
logCPM_df$cos_time <- pheno$cos_time

# Separate Tissues based on preliminary PCA
pheno_retina <- pheno[pheno$Sample == "retina", ]
pheno_choroid <- pheno[pheno$Sample == "choroid", ]

logCPM_retina <- logCPM_df[pheno_retina$Sample_ID, ]
logCPM_choroid <- logCPM_df[pheno_choroid$Sample_ID, ]

# Assign Time Classes (Target Formulation)
pheno_retina$time_class <- ifelse(pheno_retina$time_point %in% c("Z08", "Z12"), "ZT_812", "ZT_other")
pheno_choroid$time_class <- ifelse(pheno_choroid$time_point %in% c("Z08", "Z12"), "ZT_812", "ZT_other")

# Finalize Retina Data 
logCPM_retina$eye_treatment <- NULL
logCPM_retina <- data.frame(logCPM_retina, Sex = pheno_retina$Sex, time_class = pheno_retina$time_class)

# Finalize Choroid Data
logCPM_choroid$eye_treatment <- NULL
logCPM_choroid <- data.frame(logCPM_choroid, Sex = pheno_choroid$Sex, time_class = pheno_choroid$time_class)


# ==============================================================================
# === 3. STAGE 1: Primary Discovery Model (Retina) ===
# ==============================================================================
set.seed(42)
# Strictly limited to training partition to avoid data leakage
train_idx <- createDataPartition(logCPM_retina$time_class, p = 0.8, list = FALSE)
train_data <- logCPM_retina[train_idx, ]
test_data  <- logCPM_retina[-train_idx, ]

# Feature Selection
boruta_summary <- repeat_feature_selection(method = "Boruta", train_data = train_data, pheno_data = pheno_retina, train_idx = train_idx, ntime = 50)
lasso_summary <- repeat_feature_selection(method = "LASSO", train_data = train_data, pheno_data = pheno_retina, train_idx = train_idx, ntime = 50)

# Select Consensus Genes (Found 50/50 times)
boruta_selected_genes <- boruta_summary$Gene[boruta_summary$Frequency >= 50]
lasso_selected_genes <- lasso_summary$Gene[lasso_summary$Frequency >= 50]

# Prepare Final Retina Datasets
train_data_boruta_selected <- data.frame(train_data[, boruta_selected_genes, drop=FALSE], Sex = train_data$Sex, time_class = train_data$time_class)
test_data_boruta_selected <- data.frame(test_data[, boruta_selected_genes, drop=FALSE], Sex = test_data$Sex, time_class = test_data$time_class)
full_data_boruta_selected <- rbind(train_data_boruta_selected, test_data_boruta_selected)

train_data_lasso_selected <- data.frame(train_data[, lasso_selected_genes, drop=FALSE], Sex = train_data$Sex, time_class = train_data$time_class)
test_data_lasso_selected <- data.frame(test_data[, lasso_selected_genes, drop=FALSE], Sex = test_data$Sex, time_class = test_data$time_class)
full_data_lasso_selected <- rbind(train_data_lasso_selected, test_data_lasso_selected)

# Run Retina Classification (Nested CV)
cat("\n--- Running Stage 1: Retina Classification (Boruta Signatures) ---\n")
result_boruta_rf <- repeat_classification(data = full_data_boruta_selected, model_type = "RF", nround = 50)
result_boruta_svm <- repeat_classification(data = full_data_boruta_selected, model_type = "SVM", nround = 50)
result_boruta_nb <- repeat_classification(data = full_data_boruta_selected, model_type = "NB", nround = 50)
result_boruta_lr <- repeat_classification(data = full_data_boruta_selected, model_type = "LR", nround = 50)

cat("\n--- Running Stage 1: Retina Classification (LASSO Signatures) ---\n")
result_lasso_rf <- repeat_classification(data = full_data_lasso_selected, model_type = "RF", nround = 50)
result_lasso_svm <- repeat_classification(data = full_data_lasso_selected, model_type = "SVM", nround = 50)
result_lasso_nb <- repeat_classification(data = full_data_lasso_selected, model_type = "NB", nround = 50)
result_lasso_lr <- repeat_classification(data = full_data_lasso_selected, model_type = "LR", nround = 50)


# ==============================================================================
# === 4. STAGE 2: Cross-Tissue Validation (Choroid) ===
# ==============================================================================
cat("\n--- Running Stage 2: Cross-Tissue Validation (Choroid) ---\n")
# Filter choroid data to use ONLY the Retina-derived signatures
full_data_choroid_boruta <- data.frame(
  logCPM_choroid[, intersect(boruta_selected_genes, colnames(logCPM_choroid)), drop=FALSE],
  Sex = logCPM_choroid$Sex,
  time_class = logCPM_choroid$time_class
)

# Run Choroid Classification using Retina signatures
result_choroid_rf <- repeat_classification(data = full_data_choroid_boruta, model_type = "RF", nround = 50)
result_choroid_svm <- repeat_classification(data = full_data_choroid_boruta, model_type = "SVM", nround = 50)
result_choroid_nb <- repeat_classification(data = full_data_choroid_boruta, model_type = "NB", nround = 50)
result_choroid_lr <- repeat_classification(data = full_data_choroid_boruta, model_type = "LR", nround = 50)


# ==============================================================================
# === 5. STAGE 3: Cross-Stage Validation (Progression Dataset) ===
# ==============================================================================
cat("\n--- Running Stage 3: Cross-Stage Validation (Progression Dataset) ---\n")
# NOTE FOR GITHUB: Replace 'GSE_Progression_pheno.csv' and 'GSE_Progression_counts.csv' with actual filenames.
pheno_prog <- read.csv("GSE261232_counts_progress.csv")
counts_prog <- read.csv("GSE261232_pheno_progress.csv", row.names = 1)

# Align Data
common_ids_prog <- intersect(str_trim(pheno_prog$Sample_ID), str_trim(colnames(counts_prog)))
pheno_prog <- pheno_prog %>% filter(str_trim(Sample_ID) %in% common_ids_prog)
counts_prog <- counts_prog[, pheno_prog$Sample_ID]

# Normalize Progression Data via edgeR
dge_prog <- DGEList(counts = counts_prog)
keep_prog <- filterByExpr(dge_prog) 
dge_prog <- dge_prog[keep_prog, , keep.lib.sizes = TRUE]
dge_prog <- calcNormFactors(dge_prog)

logCPM_prog <- t(cpm(dge_prog, log = TRUE, prior.count = 1))
logCPM_prog_df <- as.data.frame(logCPM_prog)

# Assign Time Classes based on established criteria
pheno_prog$time_class <- ifelse(pheno_prog$time_point %in% c("Z08", "Z12"), "ZT_812", "ZT_other")

# Filter progression data to use ONLY the Retina-derived signatures
full_data_progression_boruta <- data.frame(
  logCPM_prog_df[, intersect(boruta_selected_genes, colnames(logCPM_prog_df)), drop=FALSE],
  Sex = pheno_prog$Sex,
  time_class = pheno_prog$time_class
)

# Run Progression Classification using Retina signatures
result_progression_rf <- repeat_classification(data = full_data_progression_boruta, model_type = "RF", nround = 50)
result_progression_svm <- repeat_classification(data = full_data_progression_boruta, model_type = "SVM", nround = 50)
result_progression_nb <- repeat_classification(data = full_data_progression_boruta, model_type = "NB", nround = 50)
result_progression_lr <- repeat_classification(data = full_data_progression_boruta, model_type = "LR", nround = 50)


# ==============================================================================
# === 6. Combine & Export Results ===
# ==============================================================================

# Compile Stage 1 Results
all_results_retina <- combine_results(
  lasso_rf = result_lasso_rf, lasso_svm = result_lasso_svm,
  lasso_nb = result_lasso_nb, lasso_lr = result_lasso_lr,
  boruta_rf = result_boruta_rf, boruta_svm = result_boruta_svm,
  boruta_nb = result_boruta_nb, boruta_lr = result_boruta_lr
)
mean_wide_retina <- combine_means_wide(
  lasso_rf = result_lasso_rf, lasso_svm = result_lasso_svm,
  lasso_nb = result_lasso_nb, lasso_lr = result_lasso_lr,
  boruta_rf = result_boruta_rf, boruta_svm = result_boruta_svm,
  boruta_nb = result_boruta_nb, boruta_lr = result_boruta_lr
)

# Compile Stage 2 Results (Choroid)
all_results_choroid <- combine_results(
  choroid_rf = result_choroid_rf, choroid_svm = result_choroid_svm,
  choroid_nb = result_choroid_nb, choroid_lr = result_choroid_lr
)

# Compile Stage 3 Results (Progression)
all_results_progression <- combine_results(
  progression_rf = result_progression_rf, progression_svm = result_progression_svm,
  progression_nb = result_progression_nb, progression_lr = result_progression_lr
)

# Export All Tables
write.csv(lasso_summary, "lasso_summary_retina_onset_timeclass.csv", row.names = FALSE)
write.csv(boruta_summary, "boruta_summary_retina_onset_timeclass.csv", row.names = FALSE)
write.csv(all_results_retina, "all_results_retina_onset_timeclass.csv", row.names = FALSE)
write.csv(mean_wide_retina, "mean_retina_onset_timeclass.csv", row.names = FALSE)
write.csv(all_results_choroid, "all_results_choroid_validation.csv", row.names = FALSE)
write.csv(all_results_progression, "all_results_progression_validation.csv", row.names = FALSE)

# Save Workspace
save.image("myopia_model_stages1_to_3.rdata")
cat("\nPipeline Execution Complete. Environment Saved.\n")