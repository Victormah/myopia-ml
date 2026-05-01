
# ------------------------------------------------------------------------------
# 0. Load libraries
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(edgeR)
  library(tidyverse)
  library(caret)
  library(randomForest)
  library(pROC)
  library(glmnet)
  library(Boruta)
  library(e1071)
  library(ggplot2)
  library(iml)       # optional SHAP/model interpretation
})

set.seed(123)

# ------------------------------------------------------------------------------
# 1. User parameters
# ------------------------------------------------------------------------------
params <- list(
  onset_pheno_file       = "GSE227724_pheno_onset.csv",
  onset_counts_file      = "GSE227724_counts_onset.csv",
  progression_pheno_file = "GSE261232_pheno_progress.csv",
  progression_counts_file= "GSE261232_counts_progress.csv",
  sample_id_col          = "Sample_ID",
  tissue_col             = "Sample",
  time_col               = "time_point",
  treatment_col          = "eye_treatment",
  sex_col                = "Sex",
  zt_window              = c(8, 12),
  n_feature_runs         = 50,
  boruta_max_runs        = 100,
  n_outer_rounds         = 50,
  train_ratio            = 0.80,
  n_permutations         = 100,   # increase to 500-1000 for final submission
  n_random_gene_sets     = 100,   # increase for final submission
  control_nround         = 20     # lower than final model to reduce runtime
)

# Candidate reference-gene symbols for optional negative-control analysis.
# These are common housekeeping/reference genes used broadly in expression studies;
# they are NOT assumed to be stable in chick retina/choroid. The code below converts
# them to Gallus gallus Ensembl IDs when the expression matrix uses Ensembl IDs.
housekeeping_6 <- c(
  "ENSGALG00000016462",
  "ENSGALG00000004505",
  "ENSGALG00000008103",
  "ENSGALG00000016259",
  "ENSGALG00000003099",
  "ENSGALG00000008151"
)

housekeeping_6 <- intersect(housekeeping_6, colnames(retina_onset$data))
length(housekeeping_6)
housekeeping_6


# ============================================================================== 
# 2. Helper functions
# ============================================================================== 

parse_zt <- function(x) {
  # Converts labels such as Z08, Z8, ZT08, ZT8, 8 into numeric Zeitgeber time.
  z <- gsub("ZT", "", x, ignore.case = TRUE)
  z <- gsub("Z", "", z, ignore.case = TRUE)
  z <- gsub("[^0-9.]", "", z)
  as.numeric(z)
}

add_circadian_metadata <- function(pheno, time_col = "time_point") {
  pheno %>%
    mutate(
      Time_raw = .data[[time_col]],
      Time_num = parse_zt(.data[[time_col]]),
      Time_rad = 2 * pi * Time_num / 24,
      sin_time = sin(Time_rad),
      cos_time = cos(Time_rad)
    )
}

add_zt_window_label <- function(pheno, zt_window = c(8, 12), label_col = "time_class") {
  pheno[[label_col]] <- ifelse(pheno$Time_num %in% zt_window, "ZT8_12", "ZT_other")
  pheno[[label_col]] <- factor(pheno[[label_col]], levels = c("ZT8_12", "ZT_other"))
  pheno
}

mcc_binary <- function(pred, truth, positive = "ZT8_12") {
  pred <- factor(pred, levels = levels(truth))
  truth <- factor(truth, levels = levels(truth))
  negative <- setdiff(levels(truth), positive)[1]
  cm <- table(Predicted = pred, Actual = truth)
  TP <- ifelse(positive %in% rownames(cm) && positive %in% colnames(cm), cm[positive, positive], 0)
  TN <- ifelse(negative %in% rownames(cm) && negative %in% colnames(cm), cm[negative, negative], 0)
  FP <- ifelse(positive %in% rownames(cm) && negative %in% colnames(cm), cm[positive, negative], 0)
  FN <- ifelse(negative %in% rownames(cm) && positive %in% colnames(cm), cm[negative, positive], 0)
  denom <- sqrt(as.numeric(TP + FP) * as.numeric(TP + FN) * as.numeric(TN + FP) * as.numeric(TN + FN))
  ifelse(denom == 0, NA_real_, (TP * TN - FP * FN) / denom)
}

metric_summary <- function(df_res) {
  df_res %>%
    pivot_longer(cols = everything(), names_to = "Metric", values_to = "Value") %>%
    group_by(Metric) %>%
    summarise(
      Mean = mean(Value, na.rm = TRUE),
      SD = sd(Value, na.rm = TRUE),
      CI95_low = quantile(Value, 0.025, na.rm = TRUE),
      CI95_high = quantile(Value, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
}

# Select gene-expression predictors only. This is the main leakage-prevention step.
get_expression_features <- function(data, label_col = "time_class") {
  forbidden <- c(
    label_col, "time_class", "Time_raw", "Time_num", "Time_rad", "sin_time", "cos_time",
    "time_point", "ZT", "day", "day_num", "Sample_ID", "Sample", "eye_treatment",
    "Sex", "sex", "tissue", "group"
  )
  candidate <- setdiff(colnames(data), forbidden)
  candidate[vapply(data[candidate], is.numeric, logical(1))]
}

# ------------------------------------------------------------------------------
# Optional housekeeping/reference-gene control: symbol-to-Ensembl conversion
# ------------------------------------------------------------------------------
clean_ensembl_id <- function(x) sub("\\..*$", "", x)

map_gallus_symbols_to_ensembl <- function(symbols,
                                          cache_file = "housekeeping_symbol_to_ensembl_gallus.csv",
                                          ensembl_version = NULL) {
  # Uses Ensembl BioMart to convert Gallus gallus external gene symbols to
  # Ensembl gene IDs. A cache file is written so the final analysis is reproducible.
  if (file.exists(cache_file)) {
    return(read.csv(cache_file, stringsAsFactors = FALSE))
  }

  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    warning("Package 'biomaRt' is not installed. Install it or provide a mapping CSV.")
    return(data.frame())
  }

  mart <- tryCatch({
    if (is.null(ensembl_version)) {
      biomaRt::useEnsembl(biomart = "genes", dataset = "ggallus_gene_ensembl")
    } else {
      biomaRt::useEnsembl(biomart = "genes", dataset = "ggallus_gene_ensembl",
                          version = ensembl_version)
    }
  }, error = function(e) {
    warning("Could not connect to Ensembl BioMart: ", conditionMessage(e))
    NULL
  })

  if (is.null(mart)) return(data.frame())

  mapping <- tryCatch({
    biomaRt::getBM(
      attributes = c("external_gene_name", "ensembl_gene_id", "gene_biotype", "description"),
      filters = "external_gene_name",
      values = unique(symbols),
      mart = mart
    )
  }, error = function(e) {
    warning("BioMart query failed: ", conditionMessage(e))
    data.frame()
  })

  if (nrow(mapping) > 0) {
    mapping <- mapping %>%
      mutate(requested_symbol = external_gene_name) %>%
      arrange(requested_symbol, desc(gene_biotype == "protein_coding")) %>%
      distinct(requested_symbol, ensembl_gene_id, .keep_all = TRUE)
    write.csv(mapping, cache_file, row.names = FALSE)
  }
  mapping
}

prepare_housekeeping_control_features <- function(data,
                                                  symbols = housekeeping_symbols,
                                                  mapping_file = "housekeeping_symbol_to_ensembl_gallus.csv",
                                                  ensembl_version = NULL) {
  # Returns features matching the expression matrix column identifiers.
  # If columns are gene symbols, it uses symbols directly. If columns are Ensembl
  # IDs, it converts symbols to Gallus gallus Ensembl IDs and matches those IDs.
  features <- get_expression_features(data)

  # Case 1: expression matrix already uses gene symbols.
  symbol_hits <- intersect(symbols, features)
  if (length(symbol_hits) >= 2) {
    write.csv(data.frame(feature_used = symbol_hits, input_type = "gene_symbol"),
              "housekeeping_control_features_used.csv", row.names = FALSE)
    return(symbol_hits)
  }

  # Case 2: expression matrix uses Ensembl IDs, with or without version suffix.
  mapping <- map_gallus_symbols_to_ensembl(symbols, cache_file = mapping_file,
                                           ensembl_version = ensembl_version)
  if (nrow(mapping) == 0 || !("ensembl_gene_id" %in% colnames(mapping))) {
    warning("No housekeeping/reference-gene mapping available; skipping housekeeping control.")
    return(character(0))
  }

  feature_tbl <- data.frame(
    original_feature = features,
    clean_feature = clean_ensembl_id(features),
    stringsAsFactors = FALSE
  )
  hits <- feature_tbl$original_feature[feature_tbl$clean_feature %in% mapping$ensembl_gene_id]

  used <- mapping %>%
    filter(ensembl_gene_id %in% clean_ensembl_id(hits)) %>%
    select(requested_symbol, external_gene_name, ensembl_gene_id, gene_biotype, description)
  write.csv(used, "housekeeping_symbol_to_ensembl_used.csv", row.names = FALSE)
  write.csv(data.frame(feature_used = hits, input_type = "ensembl_gene_id"),
            "housekeeping_control_features_used.csv", row.names = FALSE)

  if (length(hits) < 2) {
    warning("Fewer than two housekeeping/reference genes matched the expression matrix after Ensembl conversion.")
  }
  hits
}

preprocess_rnaseq <- function(counts_file, pheno_file, sample_id_col = "Sample_ID",
                              time_col = "time_point", group_col = NULL) {
  pheno <- read.csv(pheno_file, stringsAsFactors = FALSE, check.names = FALSE) %>%
    mutate(across(where(is.character), stringr::str_trim))
  counts <- read.csv(counts_file, row.names = 1, check.names = FALSE)
  colnames(counts) <- stringr::str_trim(colnames(counts))
  pheno[[sample_id_col]] <- stringr::str_trim(pheno[[sample_id_col]])

  common_ids <- intersect(pheno[[sample_id_col]], colnames(counts))
  if (length(common_ids) == 0) stop("No overlapping sample IDs between metadata and count matrix.")
  pheno <- pheno %>% filter(.data[[sample_id_col]] %in% common_ids)
  counts <- counts[, pheno[[sample_id_col]], drop = FALSE]

  pheno <- add_circadian_metadata(pheno, time_col = time_col)
  pheno <- add_zt_window_label(pheno, zt_window = params$zt_window)

  dge <- DGEList(counts = counts)
  if (!is.null(group_col) && group_col %in% colnames(pheno)) {
    keep <- filterByExpr(dge, group = pheno[[group_col]])
  } else {
    keep <- filterByExpr(dge)
  }
  dge <- dge[keep, , keep.lib.sizes = TRUE]
  dge <- calcNormFactors(dge)
  logCPM <- t(cpm(dge, log = TRUE, prior.count = 1))
  logCPM <- as.data.frame(logCPM, check.names = FALSE)
  list(pheno = pheno, logCPM = logCPM, dge = dge)
}

make_tissue_dataset <- function(logCPM, pheno, tissue = "retina",
                                sample_id_col = "Sample_ID", tissue_col = "Sample") {
  p <- pheno %>% filter(.data[[tissue_col]] == tissue)
  x <- logCPM[p[[sample_id_col]], , drop = FALSE]
  stopifnot(nrow(x) == nrow(p))
  # Metadata are retained for plotting/sensitivity analyses, but predictors are later restricted to genes only.
  out <- cbind(
    x,
    Time_num = p$Time_num,
    Time_rad = p$Time_rad,
    sin_time = p$sin_time,
    cos_time = p$cos_time,
    time_class = p$time_class
  )
  if (params$sex_col %in% colnames(p)) out$Sex <- p[[params$sex_col]]
  if (params$treatment_col %in% colnames(p)) out$eye_treatment <- p[[params$treatment_col]]
  rownames(out) <- p[[sample_id_col]]
  list(data = out, pheno = p)
}

# ------------------------------------------------------------------------------
# Harmonic gene-level model for periodicity
# ------------------------------------------------------------------------------
fit_harmonic_models <- function(expr_df, pheno, covariates = c("eye_treatment", "Sex"),
                                label = "retina_onset") {
  # expr_df: samples x genes, matching row order of pheno.
  stopifnot(nrow(expr_df) == nrow(pheno))
  covariates <- covariates[covariates %in% colnames(pheno)]
  meta <- pheno %>% select(any_of(c("sin_time", "cos_time", covariates)))

  results <- lapply(colnames(expr_df), function(g) {
    df <- data.frame(y = expr_df[[g]], meta, check.names = FALSE)

    full_terms <- c("sin_time", "cos_time", covariates)
    null_terms <- covariates
    full_formula <- as.formula(paste("y ~", paste(full_terms, collapse = " + ")))
    null_formula <- if (length(null_terms) > 0) {
      as.formula(paste("y ~", paste(null_terms, collapse = " + ")))
    } else {
      y ~ 1
    }

    fit_full <- tryCatch(lm(full_formula, data = df), error = function(e) NULL)
    fit_null <- tryCatch(lm(null_formula, data = df), error = function(e) NULL)
    if (is.null(fit_full) || is.null(fit_null)) return(NULL)

    an <- tryCatch(anova(fit_null, fit_full), error = function(e) NULL)
    p_circ <- if (!is.null(an) && nrow(an) >= 2) an$`Pr(>F)`[2] else NA_real_
    co <- coef(fit_full)
    b_sin <- ifelse("sin_time" %in% names(co), co[["sin_time"]], NA_real_)
    b_cos <- ifelse("cos_time" %in% names(co), co[["cos_time"]], NA_real_)
    amplitude <- sqrt(b_sin^2 + b_cos^2)
    phase_rad <- atan2(b_sin, b_cos)
    peak_time <- (phase_rad %% (2 * pi)) * 24 / (2 * pi)

    data.frame(
      Gene = g,
      p_circadian = p_circ,
      beta_sin = b_sin,
      beta_cos = b_cos,
      amplitude = amplitude,
      peak_time_ZT = peak_time,
      stringsAsFactors = FALSE
    )
  })

  res <- bind_rows(results) %>%
    mutate(FDR = p.adjust(p_circadian, method = "BH")) %>%
    arrange(FDR, desc(amplitude))
  write.csv(res, paste0("harmonic_periodicity_", label, ".csv"), row.names = FALSE)
  res
}

plot_harmonic_summary <- function(harmonic_results, output_file = "harmonic_periodicity_volcano.pdf") {
  p <- ggplot(harmonic_results, aes(x = peak_time_ZT, y = -log10(FDR), size = amplitude)) +
    geom_point(alpha = 0.6) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    scale_x_continuous(breaks = seq(0, 24, 4), limits = c(0, 24)) +
    labs(x = "Estimated peak time (ZT)", y = "-log10(FDR)", size = "Amplitude",
         title = "Harmonic model of circadian-associated gene expression") +
    theme_bw()
  ggsave(output_file, p, width = 8, height = 5)
  p
}

# ------------------------------------------------------------------------------
# Leakage-free feature selection
# ------------------------------------------------------------------------------
repeat_feature_selection <- function(method = c("LASSO", "Boruta"), data, label_col = "time_class",
                                     n_runs = 50, boruta_max_runs = 100, verbose = TRUE) {
  method <- match.arg(method)
  features <- get_expression_features(data, label_col = label_col)
  x <- data[, features, drop = FALSE]
  y <- factor(data[[label_col]], levels = c("ZT8_12", "ZT_other"))

  all_selected <- vector("list", n_runs)
  importance_list <- vector("list", n_runs)

  for (i in seq_len(n_runs)) {
    set.seed(100 + i)
    if (method == "LASSO") {
      X <- as.matrix(x)
      y_bin <- ifelse(y == "ZT8_12", 1, 0)
      cvfit <- cv.glmnet(X, y_bin, family = "binomial", alpha = 1, nfolds = min(5, length(y_bin)))
      coef_mat <- as.matrix(coef(cvfit, s = "lambda.min"))
      selected <- rownames(coef_mat)[coef_mat[, 1] != 0]
      selected <- setdiff(selected, "(Intercept)")
    } else {
      # Use Boruta x/y interface to avoid formula problems with non-syntactic gene IDs.
      boruta_model <- Boruta(x = x, y = y, maxRuns = boruta_max_runs, doTrace = 0)
      selected <- getSelectedAttributes(boruta_model, withTentative = FALSE)
      imp <- attStats(boruta_model) %>%
        rownames_to_column("Gene") %>%
        mutate(Run = i)
      importance_list[[i]] <- imp
    }
    all_selected[[i]] <- selected
    if (verbose) cat(sprintf("[%s] run %d/%d: %d features selected\n", method, i, n_runs, length(selected)))
  }

  freq_df <- tibble(Gene = unlist(all_selected)) %>%
    count(Gene, name = "Frequency") %>%
    mutate(SelectionRate = Frequency / n_runs) %>%
    arrange(desc(Frequency), Gene)

  list(
    frequency = freq_df,
    selected_by_run = all_selected,
    importance = bind_rows(importance_list)
  )
}

plot_selection_frequency <- function(freq_df, n_runs, top_n = 60,
                                     output_file = "boruta_selection_frequency.pdf") {
  p <- freq_df %>%
    slice_max(order_by = Frequency, n = top_n, with_ties = FALSE) %>%
    mutate(Gene = factor(Gene, levels = rev(Gene))) %>%
    ggplot(aes(x = Gene, y = Frequency)) +
    geom_col() +
    coord_flip() +
    geom_hline(yintercept = n_runs, linetype = "dashed") +
    labs(x = "Gene", y = paste0("Selection frequency across ", n_runs, " runs"),
         title = "Boruta feature-selection stability") +
    theme_bw()
  ggsave(output_file, p, width = 8, height = 10)
  p
}

# ------------------------------------------------------------------------------
# Leakage-free classification using fixed expression-only features
# ------------------------------------------------------------------------------
run_repeated_classification <- function(data, feature_genes, model_type = c("RF", "SVM", "NB", "LR"),
                                        label_col = "time_class", nround = 50,
                                        train_ratio = 0.8, seed_offset = 1) {
  model_type <- match.arg(model_type)
  feature_genes <- intersect(feature_genes, get_expression_features(data, label_col = label_col))
  if (length(feature_genes) < 2) stop("Fewer than two valid expression features available.")

  y <- factor(data[[label_col]], levels = c("ZT8_12", "ZT_other"))
  x <- data[, feature_genes, drop = FALSE]

  res <- vector("list", nround)
  ctrl <- trainControl(
    method = "cv", number = 5,
    classProbs = TRUE,
    summaryFunction = twoClassSummary,
    savePredictions = "final"
  )

  for (i in seq_len(nround)) {
    set.seed(seed_offset + i)
    train_idx <- createDataPartition(y, p = train_ratio, list = FALSE)
    x_train <- x[train_idx, , drop = FALSE]
    x_test  <- x[-train_idx, , drop = FALSE]
    y_train <- droplevels(y[train_idx])
    y_test  <- factor(y[-train_idx], levels = levels(y_train))

    if (model_type == "RF") {
      fit <- train(x = x_train, y = y_train, method = "rf", metric = "ROC",
                   trControl = ctrl, tuneLength = 3, ntree = 500)
    } else if (model_type == "SVM") {
      fit <- train(x = x_train, y = y_train, method = "svmRadial", metric = "ROC",
                   trControl = ctrl, preProcess = c("center", "scale"), tuneLength = 5)
    } else if (model_type == "NB") {
      fit <- train(x = x_train, y = y_train, method = "nb", metric = "ROC",
                   trControl = ctrl, tuneLength = 3)
    } else if (model_type == "LR") {
      fit <- train(x = x_train, y = y_train, method = "glmnet", metric = "ROC",
                   trControl = ctrl, preProcess = c("center", "scale"), tuneLength = 5)
    }

    pred <- predict(fit, x_test)
    prob <- predict(fit, x_test, type = "prob")[["ZT8_12"]]
    cm <- confusionMatrix(pred, y_test, positive = "ZT8_12")
    auc_val <- tryCatch(as.numeric(auc(roc(response = y_test, predictor = prob,
                                          levels = c("ZT_other", "ZT8_12"), quiet = TRUE))),
                        error = function(e) NA_real_)

    res[[i]] <- data.frame(
      Accuracy = unname(cm$overall["Accuracy"]),
      Precision = unname(cm$byClass["Precision"]),
      Recall = unname(cm$byClass["Recall"]),
      F1 = unname(cm$byClass["F1"]),
      AUC = auc_val,
      MCC = mcc_binary(pred, y_test, positive = "ZT8_12")
    )
  }

  raw <- bind_rows(res)
  list(raw = raw, summary = metric_summary(raw), features = feature_genes)
}

run_all_models <- function(data, feature_genes, prefix, nround = params$n_outer_rounds) {
  models <- c("RF", "SVM", "NB", "LR")
  out <- setNames(vector("list", length(models)), models)
  for (m in models) {
    cat("Running", prefix, m, "\n")
    out[[m]] <- run_repeated_classification(data, feature_genes, model_type = m, nround = nround)
  }
  summary <- bind_rows(lapply(names(out), function(m) {
    out[[m]]$summary %>% mutate(Model = m, Analysis = prefix)
  })) %>% select(Analysis, Model, Metric, Mean, SD, CI95_low, CI95_high)
  list(models = out, summary = summary)
}


# ------------------------------------------------------------------------------
# Optional SHAP/model interpretation for final RF model
# ------------------------------------------------------------------------------
compute_rf_importance_and_shap <- function(data, feature_genes, output_prefix = "retina_boruta") {
  genes <- intersect(feature_genes, get_expression_features(data))
  x <- data[, genes, drop = FALSE]
  y <- factor(data$time_class, levels = c("ZT8_12", "ZT_other"))

  set.seed(999)
  rf_fit <- randomForest(x = x, y = y, ntree = 500, importance = TRUE)
  imp <- importance(rf_fit) %>%
    as.data.frame() %>%
    rownames_to_column("Gene") %>%
    arrange(desc(MeanDecreaseGini))
  write.csv(imp, paste0(output_prefix, "_RF_feature_importance.csv"), row.names = FALSE)

  # SHAP using iml can be slow; this computes SHAP for top features and a subset of samples.
  pred_fun <- function(model, newdata) {
    predict(model, newdata = newdata, type = "prob")[["ZT8_12"]]
  }
  predictor <- Predictor$new(model = rf_fit, data = x, y = y, predict.function = pred_fun, class = "ZT8_12")
  top_genes <- head(imp$Gene, min(15, nrow(imp)))
  shap_list <- lapply(seq_len(min(10, nrow(x))), function(i) {
    s <- Shapley$new(predictor, x.interest = x[i, , drop = FALSE])
    s$results %>% mutate(Sample = rownames(x)[i])
  })
  shap_df <- bind_rows(shap_list) %>% filter(feature %in% top_genes)
  write.csv(shap_df, paste0(output_prefix, "_SHAP_top_features.csv"), row.names = FALSE)
  list(rf_model = rf_fit, importance = imp, shap = shap_df)
}

# ============================================================================== 
# 3. Load and preprocess onset data
# ============================================================================== 
cat("\n=== Loading onset dataset ===\n")
onset <- preprocess_rnaseq(
  counts_file = params$onset_counts_file,
  pheno_file = params$onset_pheno_file,
  sample_id_col = params$sample_id_col,
  time_col = params$time_col,
  group_col = params$treatment_col
)

retina_onset <- make_tissue_dataset(onset$logCPM, onset$pheno, tissue = "retina",
                                    sample_id_col = params$sample_id_col,
                                    tissue_col = params$tissue_col)
choroid_onset <- make_tissue_dataset(onset$logCPM, onset$pheno, tissue = "choroid",
                                     sample_id_col = params$sample_id_col,
                                     tissue_col = params$tissue_col)

# ============================================================================== 
# 4. Rhythm-aware analysis: harmonic modelling of circadian periodicity
# ============================================================================== 
cat("\n=== Harmonic/cosinor-like gene-level modelling ===\n")
retina_expr_only <- retina_onset$data[, get_expression_features(retina_onset$data), drop = FALSE]
choroid_expr_only <- choroid_onset$data[, get_expression_features(choroid_onset$data), drop = FALSE]

harmonic_retina <- fit_harmonic_models(
  expr_df = retina_expr_only,
  pheno = retina_onset$pheno,
  covariates = c(params$treatment_col, params$sex_col),
  label = "retina_onset"
)
plot_harmonic_summary(harmonic_retina, "harmonic_periodicity_retina_onset.pdf")

harmonic_choroid <- fit_harmonic_models(
  expr_df = choroid_expr_only,
  pheno = choroid_onset$pheno,
  covariates = c(params$treatment_col, params$sex_col),
  label = "choroid_onset"
)
plot_harmonic_summary(harmonic_choroid, "harmonic_periodicity_choroid_onset.pdf")

# ============================================================================== 
# 5. Stage 1: leakage-free feature discovery in retina onset data
# ============================================================================== 
cat("\n=== Stage 1: expression-only feature selection ===\n")
set.seed(42)
train_idx <- createDataPartition(retina_onset$data$time_class, p = params$train_ratio, list = FALSE)
retina_train <- retina_onset$data[train_idx, , drop = FALSE]
retina_test  <- retina_onset$data[-train_idx, , drop = FALSE]

boruta_fs <- repeat_feature_selection(
  method = "Boruta",
  data = retina_train,
  n_runs = params$n_feature_runs,
  boruta_max_runs = params$boruta_max_runs
)

lasso_fs <- repeat_feature_selection(
  method = "LASSO",
  data = retina_train,
  n_runs = params$n_feature_runs,
  boruta_max_runs = params$boruta_max_runs
)

write.csv(boruta_fs$frequency, "boruta_selection_frequency_retina_onset_expression_only.csv", row.names = FALSE)
write.csv(lasso_fs$frequency, "lasso_selection_frequency_retina_onset_expression_only.csv", row.names = FALSE)
write.csv(boruta_fs$importance, "boruta_importance_all_runs_retina_onset.csv", row.names = FALSE)
plot_selection_frequency(boruta_fs$frequency, params$n_feature_runs, output_file = "boruta_selection_frequency_retina_onset.pdf")

boruta_selected_genes <- boruta_fs$frequency %>%
  filter(Frequency == params$n_feature_runs) %>%
  pull(Gene)
lasso_selected_genes <- lasso_fs$frequency %>%
  filter(Frequency == params$n_feature_runs) %>%
  pull(Gene)

if (length(boruta_selected_genes) < 2) {
  warning("Few/no genes selected in all Boruta runs. Using genes selected in >=80% of runs for downstream sensitivity analysis.")
  boruta_selected_genes <- boruta_fs$frequency %>%
    filter(SelectionRate >= 0.80) %>%
    pull(Gene)
}

writeLines(boruta_selected_genes, "final_boruta_signature_expression_only.txt")
writeLines(lasso_selected_genes, "final_lasso_signature_expression_only.txt")

# PCA visualization requested by reviewer: variance explained + confidence ellipse
plot_pca_signature(retina_onset$data, boruta_selected_genes,
                   output_file = "PCA_retina_onset_boruta_variance_ellipse.pdf")

# ============================================================================== 
# 6. Stage 1 classification as targeted sensitivity analysis, not primary time model
# ============================================================================== 
cat("\n=== Stage 1: targeted ZT8/12 classification sensitivity analysis ===\n")
stage1_boruta <- run_all_models(retina_onset$data, boruta_selected_genes,
                                prefix = "Stage1_retina_onset_Boruta",
                                nround = params$n_outer_rounds)
write.csv(stage1_boruta$summary, "stage1_retina_onset_boruta_expression_only_summary.csv", row.names = FALSE)

if (length(lasso_selected_genes) >= 2) {
  stage1_lasso <- run_all_models(retina_onset$data, lasso_selected_genes,
                                 prefix = "Stage1_retina_onset_LASSO",
                                 nround = params$n_outer_rounds)
  write.csv(stage1_lasso$summary, "stage1_retina_onset_lasso_expression_only_summary.csv", row.names = FALSE)
}

boruta_selected_genes <- readLines("boruta_selected_genes.txt")
boruta_selected_genes <- unique(boruta_selected_genes)

# ============================================================================== 
# 7. Reviewer 3 control analyses
# ============================================================================== 
cat("\n=== Reviewer 3 control analyses ===\n")

# 7.1 Time-label permutation control
perm_control <- permutation_label_control(
  data = retina_onset$data,
  feature_genes = boruta_selected_genes,
  model_type = "RF",
  nperm = params$n_permutations,
  nround = params$control_nround
)
write.csv(perm_control, "control_time_label_permutation_rf.csv", row.names = FALSE)


# 7.2 Random gene-set controls
random_control <- random_gene_set_control(
  data = retina_onset$data,
  signature_size = length(boruta_selected_genes),
  model_type = "RF",
  nsets = params$n_random_gene_sets,
  nround = params$control_nround
)
write.csv(random_control, "control_random_gene_sets_rf.csv", row.names = FALSE)
observed_auc <- stage1_boruta$models$RF$summary %>% filter(Metric == "AUC") %>% pull(Mean)
plot_control_auc(random_control, observed_auc, output_file = "control_random_gene_set_auc_distribution.pdf")

# 7.3 Housekeeping/reference-gene controls
# If the expression matrix uses Ensembl IDs, this converts the candidate symbols to
# Gallus gallus Ensembl gene IDs before intersecting with available features.
housekeeping_features <- housekeeping_6

hk_control <- housekeeping_gene_control(
  data = retina_onset$data,
  housekeeping_genes = housekeeping_features,
  model_type = "RF",
  nround = params$control_nround
)
if (!is.null(hk_control)) write.csv(hk_control, "control_housekeeping_genes_rf.csv", row.names = FALSE)

# 7.4 Single-time-point vs combined-window analysis
single_vs_window <- single_timepoint_window_analysis(
  data = retina_onset$data,
  feature_genes = boruta_selected_genes,
  model_type = "RF",
  nround = params$control_nround
)
write.csv(single_vs_window, "control_single_timepoint_vs_ZT8_12_window_rf.csv", row.names = FALSE)


# Add the real negative control: non-circadian genes
non_circadian_pool <- all_gene_time_tests %>%
  dplyr::filter(
    FDR_time >= 0.8
  ) %>%
  dplyr::pull(gene) %>%
  intersect(colnames(retina_onset$data))

length(non_circadian_pool)

non_circadian_control <- random_gene_set_control_from_pool(
  data = retina_onset$data,
  gene_pool = non_circadian_pool,
  signature_size = length(boruta_selected_genes),
  model_type = "RF",
  nsets = 100,
  nround = params$control_nround
)

write.csv(
  non_circadian_control,
  "control_non_circadian_53gene_sets_rf.csv",
  row.names = FALSE
)
# 7.5 non_circadian_control
random_gene_set_control_from_pool <- function(data,
                                              gene_pool,
                                              signature_size,
                                              model_type = "RF",
                                              nsets = 100,
                                              nround = 20,
                                              seed_base = 9000) {
  gene_pool <- intersect(gene_pool, colnames(data))
  
  if (length(gene_pool) < signature_size) {
    stop("Gene pool has fewer genes than signature_size.")
  }
  
  out <- vector("list", nsets)
  gene_lists <- vector("list", nsets)
  
  for (b in seq_len(nsets)) {
    set.seed(seed_base + b)
    
    genes <- sample(
      gene_pool,
      size = signature_size,
      replace = FALSE
    )
    
    fit <- run_repeated_classification(
      data = data,
      feature_genes = genes,
      model_type = model_type,
      nround = nround,
      seed_offset = seed_base + b * 10
    )
    
    out[[b]] <- data.frame(
      Set = b,
      GenePool = "non_circadian",
      N_features = length(genes),
      AUC = mean(fit$raw$AUC, na.rm = TRUE),
      Accuracy = mean(fit$raw$Accuracy, na.rm = TRUE),
      F1 = mean(fit$raw$F1, na.rm = TRUE),
      MCC = mean(fit$raw$MCC, na.rm = TRUE)
    )
    
    gene_lists[[b]] <- data.frame(
      Set = b,
      Gene_ID = genes
    )
    
    cat(
      "Non-circadian random set", b, "/", nsets,
      "AUC =", round(out[[b]]$AUC, 3), "\n"
    )
  }
  
  list(
    performance = dplyr::bind_rows(out),
    gene_lists = dplyr::bind_rows(gene_lists)
  )
}

non_circadian_control <- random_gene_set_control_from_pool(
  data = retina_onset$data,
  gene_pool = non_circadian_pool,
  signature_size = length(boruta_selected_genes),
  model_type = "RF",
  nsets = 100,
  nround = params$control_nround
)

write.csv(
  non_circadian_control$performance,
  "control_non_circadian_53gene_sets_rf.csv",
  row.names = FALSE
)

write.csv(
  non_circadian_control$gene_lists,
  "control_non_circadian_53gene_sets_gene_lists.csv",
  row.names = FALSE
)

non_circadian_summary <- non_circadian_control$performance |>
  dplyr::summarise(
    Mean_AUC = mean(AUC, na.rm = TRUE),
    SD_AUC = sd(AUC, na.rm = TRUE),
    Median_AUC = median(AUC, na.rm = TRUE),
    Min_AUC = min(AUC, na.rm = TRUE),
    Max_AUC = max(AUC, na.rm = TRUE),
    Mean_MCC = mean(MCC, na.rm = TRUE),
    SD_MCC = sd(MCC, na.rm = TRUE)
  )

non_circadian_summary


# ============================================================================== 
# 8. Stage 2: Cross-tissue validation using fixed retina-derived signature
# ============================================================================== 
cat("\n=== Stage 2: cross-tissue validation, choroid onset ===\n")
stage2_choroid <- run_all_models(choroid_onset$data, boruta_selected_genes,
                                 prefix = "Stage2_choroid_onset_Boruta",
                                 nround = params$n_outer_rounds)
write.csv(stage2_choroid$summary, "stage2_choroid_onset_boruta_expression_only_summary.csv", row.names = FALSE)

# ============================================================================== 
# 9. Stage 3: Cross-stage validation using fixed retina-derived signature
# ============================================================================== 
cat("\n=== Stage 3: cross-stage validation, progression dataset ===\n")
progression <- preprocess_rnaseq(
  counts_file = params$progression_counts_file,
  pheno_file = params$progression_pheno_file,
  sample_id_col = params$sample_id_col,
  time_col = params$time_col,
  group_col = params$treatment_col
)

progression_retina <- make_tissue_dataset(progression$logCPM, progression$pheno, tissue = "retina",
                                          sample_id_col = params$sample_id_col,
                                          tissue_col = params$tissue_col)

stage3_progression <- run_all_models(progression_retina$data, boruta_selected_genes,
                                     prefix = "Stage3_retina_progression_Boruta",
                                     nround = params$n_outer_rounds)
write.csv(stage3_progression$summary, "stage3_retina_progression_boruta_expression_only_summary.csv", row.names = FALSE)

# ============================================================================== 
# 10. Export combined summary and session information
# ============================================================================== 
combined_summary <- bind_rows(
  stage1_boruta$summary,
  if (exists("stage1_lasso")) stage1_lasso$summary else NULL,
  stage2_choroid$summary,
  stage3_progression$summary
)
write.csv(combined_summary, "combined_model_summary_expression_only_reviewer3.csv", row.names = FALSE)

sink("sessionInfo_reviewer3_revised.txt")
print(sessionInfo())
sink()

save.image("myopia_model_stages1_to_3_reviewer3_revised_v2.rdata")
cat("\nRevised Reviewer 3 workflow complete.\n")

