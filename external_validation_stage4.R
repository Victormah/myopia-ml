# ============================================================
# external validation script (nested CV, permutation, bootstrap CIs)

# NOTES:
#  - This is computationally heavier than the original LOOCV.
#  - Set nperm and nboot to smaller values for quick checks, increase for final results.
#  - Input files (counts_external.csv, external_metadata.csv, my_53_gene_signature.txt)
#    are the same as in your original script.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
  library(e1071)     # SVM + tune
  library(pROC)      # AUC
  library(glmnet)    # regularized logistic
  library(Matrix)    # sparse matrices for glmnet
})

# -----------------------------
# 0) Load data (unchanged)
# -----------------------------
counts_raw <- read.csv("GSE203604_counts_external.csv", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
gene_col <- colnames(counts_raw)[1]

counts_clean <- counts_raw %>%
  group_by(!!sym(gene_col)) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop")

counts <- as.data.frame(counts_clean)
rownames(counts) <- counts[[gene_col]]
counts[[gene_col]] <- NULL

meta <- read.csv("GSE203604_external_metadata.csv", header = TRUE, stringsAsFactors = FALSE)
stopifnot(all(c("sample_id", "group", "day") %in% colnames(meta)))

meta <- meta %>%
  mutate(
    group = factor(group, levels = c("C", "L")),
    day   = factor(day, levels = c("D1", "D3", "D6"), ordered = TRUE),
    day_num = c(D1 = 1, D3 = 3, D6 = 6)[day]
  )

stopifnot(all(meta$sample_id %in% colnames(counts)))
counts <- counts[, meta$sample_id, drop = FALSE]

sig_genes <- readLines("my_53_gene_signature.txt")
sig_genes <- unique(sig_genes)

# -----------------------------
# 1) Expression transfer: presence + variability
# -----------------------------
keep <- rowSums(counts >= 10) >= 3
counts_f <- counts[keep, , drop = FALSE]

# log2 transform for exploratory analyses (we will use z-scored values for classifiers)
expr <- log2(counts_f + 1)

sig_present <- intersect(sig_genes, rownames(expr))
cat("Signature genes detected:", length(sig_present), "of", length(sig_genes), "\n")
if (length(sig_present) < 10) warning("Very few signature genes matched. Check gene ID type (ENSGALG vs symbols).")

sig_sd <- apply(expr[sig_present, , drop = FALSE], 1, sd)
cat("Signature gene SD summary:\n"); print(summary(sig_sd))
low_var <- sort(sig_sd)[1:min(10, length(sig_sd))]
cat("Lowest-variance signature genes:\n"); print(low_var)

# Prepare z-scored expression matrix for signature genes (genes x samples)
expr_sig <- expr[sig_present, , drop = FALSE]
expr_sig_z <- t(scale(t(expr_sig)))   # gene-wise z

# transpose for modeling: samples x features
X_full <- t(expr_sig_z)
y_full <- meta$group
n_samples <- nrow(X_full)

# -----------------------------
# Helper: nested LOOCV function
# -----------------------------
# Performs outer LOOCV. For each left-out sample:
#  - Fit glmnet (inner CV) on training set to choose lambda
#  - Fit linear SVM (inner CV tune cost) on training set
# Returns predicted probabilities for left-out samples for both models.
nested_loocv <- function(X, y, svm_costs = c(0.01, 0.1, 1, 10), seed = 1) {
  set.seed(seed)
  n <- nrow(X)
  prob_logit <- rep(NA_real_, n)
  prob_svm   <- rep(NA_real_, n)
  
  for (i in seq_len(n)) {
    train_idx <- setdiff(seq_len(n), i)
    test_idx  <- i
    
    x_train <- X[train_idx, , drop = FALSE]
    y_train <- y[train_idx]
    x_test  <- X[test_idx, , drop = FALSE]
    
    # --- GLMNET (logistic) inner CV to choose lambda ---
    x_train_mat <- as.matrix(x_train)
    y_train_bin <- ifelse(y_train == "L", 1, 0)
    
    cvfit <- cv.glmnet(
      x = x_train_mat,
      y = y_train_bin,
      family = "binomial",
      alpha = 1,
      nfolds = min(5, length(train_idx))
    )
    prob_logit[test_idx] <- as.numeric(predict(cvfit, newx = as.matrix(x_test), s = "lambda.min", type = "response"))
    
    # --- Linear SVM inner tuning on training set ---
    # Ensure factor levels are consistent
    y_train_fac <- factor(y_train, levels = c("C", "L"))
    
    tune_out <- tryCatch({
      tune.svm(
        x = x_train,
        y = y_train_fac,
        kernel = "linear",
        cost = svm_costs,
        tunecontrol = tune.control(sampling = "cross", cross = min(5, length(train_idx)))
      )
    }, error = function(e) {
      message("SVM tuning failed in fold ", i, "; will fallback. Error: ", e$message)
      NULL
    })
    
    # Determine best cost
    if (!is.null(tune_out) && !is.null(tune_out$best.parameters$cost)) {
      best_cost <- tune_out$best.parameters$cost
    } else {
      # fallback default
      best_cost <- 1
    }
    
    # Retrain final SVM on full training set with probability = TRUE
    svm_fit <- tryCatch({
      svm(x = x_train, y = y_train_fac, kernel = "linear", cost = best_cost, probability = TRUE, scale = FALSE)
    }, error = function(e) {
      message("SVM training failed in fold ", i, "; fallback to default SVM without probability. Error: ", e$message)
      NULL
    })
    
    if (is.null(svm_fit)) {
      # If training failed, set NA and continue
      prob_svm[test_idx] <- NA_real_
      next
    }
    
    # Predict with probability = TRUE (guaranteed because we trained with probability = TRUE)
    pred <- predict(svm_fit, x_test, probability = TRUE)
    probs <- attr(pred, "probabilities")
    
    # Defensive: ensure column for class "L" exists
    if (!is.null(probs) && "L" %in% colnames(probs)) {
      prob_svm[test_idx] <- probs[, "L"]
    } else if (!is.null(probs) && "1" %in% colnames(probs)) {
      # sometimes probabilities are named by factor numeric; try the second column
      prob_svm[test_idx] <- probs[, ncol(probs)]
    } else {
      # As a last resort, use class prediction (0/1)
      pred_class <- as.character(pred)
      prob_svm[test_idx] <- ifelse(pred_class == "L", 1, 0)
    }
  }
  
  list(prob_logit = prob_logit, prob_svm = prob_svm)
}

# -----------------------------
# 3) Run nested LOOCV (this is the corrected, leakage-free LOOCV)
# -----------------------------
cat("Running nested LOOCV (this may take a few minutes)...\n")
nested_res <- nested_loocv(X_full, y_full, svm_costs = c(0.01, 0.1, 1, 10), seed = 123)

loocv_logit_prob <- nested_res$prob_logit
loocv_svm_prob   <- nested_res$prob_svm

roc_logit <- roc(response = y_full, predictor = loocv_logit_prob, levels = c("C", "L"), quiet = TRUE)
roc_svm   <- roc(response = y_full, predictor = loocv_svm_prob,   levels = c("C", "L"), quiet = TRUE)

cat("Nested LOOCV Logistic AUC:", round(as.numeric(auc(roc_logit)), 3), "\n")
cat("Nested LOOCV Linear SVM AUC:", round(as.numeric(auc(roc_svm)), 3), "\n")

# Plot ROC curves
plot(roc_logit, col = "#377EB8", lwd = 2, main = "Nested LOOCV ROC (53-gene signature)")
lines(roc_svm, col = "#E41A1C", lwd = 2)
legend("bottomright",
       legend = c(paste0("Logit AUC=", round(auc(roc_logit),3)),
                  paste0("SVM AUC=", round(auc(roc_svm),3))),
       col = c("#377EB8","#E41A1C"), lwd = 2, bty = "n")

# -----------------------------
# 4) Permutation test for empirical p-value (for SVM and Logit)
#    WARNING: can be slow. Set nperm smaller for quick checks.
# -----------------------------
# -----------------------------
# REPLACEMENT BLOCK: Faster permutation + bootstrap (sections 4 & 5)
# - Uses a fast gene-set permutation on the composite signature score
# - Runs a reduced, parallel nested LOOCV for classifier-level permutations
# - Runs bootstrap CIs in parallel
# - Uses a smaller SVM grid and fewer inner folds to reduce model fits
# - Includes defensive checks and progress messages
# Drop this block in place of your original sections 4 and 5.
# -----------------------------

suppressPackageStartupMessages({
  library(doParallel)
  library(foreach)
  library(glmnet)
  library(e1071)
  library(pROC)
})

# ---- Parameters you can tune ----
nperm_fast    <- 10000   # fast gene-set permutation on signature score (cheap)
nperm_nested  <- 400     # nested LOOCV permutations (parallel). Use 200-500 for diagnostics
nboot         <- 500     # bootstrap iterations (parallel). Use 200 for quick checks
cores         <- max(1, parallel::detectCores() - 1)
svm_costs_small <- c(0.1, 1)  # smaller SVM grid for inner tuning
inner_folds     <- 3          # inner CV folds for tuning (reduced from 5)
set.seed(1)

cat("Parallel workers:", cores, "\n")
cat("Fast perm:", nperm_fast, "Nested perm:", nperm_nested, "Bootstrap:", nboot, "\n\n")

# ---- Quick gene-set permutation on composite signature score (very fast) ----
# Compute composite signature score if not already present
if (!"signature_score" %in% colnames(meta)) {
  meta$signature_score <- rowMeans(X_full, na.rm = TRUE)
}

obs_diff <- abs(median(meta$signature_score[meta$group == "L"]) -
                  median(meta$signature_score[meta$group == "C"]))

cat("Running fast gene-set permutation (n =", nperm_fast, ") ...\n")
perm_diff_fast <- replicate(nperm_fast, {
  yperm <- sample(meta$group)
  abs(median(meta$signature_score[yperm == "L"]) -
        median(meta$signature_score[yperm == "C"]))
})
emp_p_fast <- mean(perm_diff_fast >= obs_diff)
cat("Fast signature-score permutation p =", emp_p_fast, "\n\n")

# If emp_p_fast is not small, you may decide not to run expensive nested permutations.
# Proceed to nested permutations only if emp_p_fast is small (or if you want classifier-level evidence).

# ---- Faster nested LOOCV function (accepts inner_folds and reduced SVM grid) ----
nested_loocv_fast <- function(X, y, svm_costs = svm_costs_small, inner_folds = inner_folds, seed = 1) {
  set.seed(seed)
  n <- nrow(X)
  prob_logit <- rep(NA_real_, n)
  prob_svm   <- rep(NA_real_, n)
  
  for (i in seq_len(n)) {
    train_idx <- setdiff(seq_len(n), i)
    test_idx  <- i
    
    x_train <- X[train_idx, , drop = FALSE]
    y_train <- y[train_idx]
    x_test  <- X[test_idx, , drop = FALSE]
    
    # Skip fold if training set has only one class (shouldn't happen in LOOCV but check)
    if (length(unique(y_train)) < 2) {
      prob_logit[test_idx] <- NA_real_
      prob_svm[test_idx]   <- NA_real_
      next
    }
    
    # --- GLMNET (logistic) inner CV to choose lambda (reduced folds) ---
    x_train_mat <- as.matrix(x_train)
    y_train_bin <- ifelse(y_train == "L", 1, 0)
    nfolds_glm <- min(inner_folds, length(train_idx))
    cvfit <- tryCatch({
      cv.glmnet(x = x_train_mat, y = y_train_bin, family = "binomial",
                alpha = 1, nfolds = nfolds_glm)
    }, error = function(e) {
      message("cv.glmnet failed in fold ", i, "; using ridge fallback. Error: ", e$message)
      NULL
    })
    if (!is.null(cvfit)) {
      prob_logit[test_idx] <- as.numeric(predict(cvfit, newx = as.matrix(x_test), s = "lambda.min", type = "response"))
    } else {
      # fallback: simple logistic on top 10 PCA components of training set
      pcs <- prcomp(x_train, center = TRUE, scale. = FALSE)
      k <- min(10, ncol(x_train), nrow(x_train) - 1)
      df_train <- data.frame(group = factor(y_train, levels = c("C","L")), pcs$x[, 1:k, drop = FALSE])
      df_test  <- data.frame(t(predict(pcs, newdata = x_test)[1:k]))
      fit_lm <- tryCatch(glm(group ~ ., data = df_train, family = binomial()), error = function(e) NULL)
      if (!is.null(fit_lm)) prob_logit[test_idx] <- predict(fit_lm, newdata = df_test, type = "response") else prob_logit[test_idx] <- NA_real_
    }
    
    # --- Linear SVM inner tuning on training set (reduced grid & folds) ---
    y_train_fac <- factor(y_train, levels = c("C", "L"))
    nfolds_svm <- min(inner_folds, length(train_idx))
    
    tune_out <- tryCatch({
      tune.svm(x = x_train, y = y_train_fac, kernel = "linear",
               cost = svm_costs, tunecontrol = tune.control(sampling = "cross", cross = nfolds_svm))
    }, error = function(e) {
      message("SVM tuning failed in fold ", i, "; fallback to cost=1. Error: ", e$message)
      NULL
    })
    
    best_cost <- if (!is.null(tune_out) && !is.null(tune_out$best.parameters$cost)) tune_out$best.parameters$cost else 1
    
    svm_fit <- tryCatch({
      svm(x = x_train, y = y_train_fac, kernel = "linear", cost = best_cost, probability = TRUE, scale = FALSE)
    }, error = function(e) {
      message("SVM training failed in fold ", i, "; setting NA. Error: ", e$message)
      NULL
    })
    
    if (is.null(svm_fit)) {
      prob_svm[test_idx] <- NA_real_
      next
    }
    
    pred <- predict(svm_fit, x_test, probability = TRUE)
    probs <- attr(pred, "probabilities")
    if (!is.null(probs) && "L" %in% colnames(probs)) {
      prob_svm[test_idx] <- probs[, "L"]
    } else if (!is.null(probs) && ncol(probs) >= 1) {
      prob_svm[test_idx] <- probs[, ncol(probs)]
    } else {
      prob_svm[test_idx] <- ifelse(as.character(pred) == "L", 1, 0)
    }
  }
  list(prob_logit = prob_logit, prob_svm = prob_svm)
}

# ---- Run nested LOOCV once (observed AUCs) ----
cat("Running nested LOOCV (observed) with reduced inner folds and grid...\n")
obs_res <- nested_loocv_fast(X_full, y_full, svm_costs = svm_costs_small, inner_folds = inner_folds, seed = 123)
loocv_logit_prob <- obs_res$prob_logit
loocv_svm_prob   <- obs_res$prob_svm

roc_logit <- roc(response = y_full, predictor = loocv_logit_prob, levels = c("C", "L"), quiet = TRUE)
roc_svm   <- roc(response = y_full, predictor = loocv_svm_prob,   levels = c("C", "L"), quiet = TRUE)

cat("Observed nested LOOCV Logistic AUC:", round(as.numeric(auc(roc_logit)), 3), "\n")
cat("Observed nested LOOCV SVM AUC:", round(as.numeric(auc(roc_svm)), 3), "\n\n")

# ---- Parallel nested permutations (classifier-level) ----
cat("Running parallel nested LOOCV permutations (n =", nperm_nested, ") ...\n")
cl <- makeCluster(cores); registerDoParallel(cl)

perm_res <- foreach(b = 1:nperm_nested, .combine = rbind, .packages = c("glmnet","e1071","pROC")) %dopar% {
  set.seed(1000 + b)
  y_perm <- sample(y_full)
  res_perm <- nested_loocv_fast(X_full, y_perm, svm_costs = svm_costs_small, inner_folds = inner_folds, seed = 1000 + b)
  auc_l <- tryCatch(as.numeric(auc(roc(response = y_perm, predictor = res_perm$prob_logit, levels = c("C","L")))), error = function(e) NA_real_)
  auc_s <- tryCatch(as.numeric(auc(roc(response = y_perm, predictor = res_perm$prob_svm,   levels = c("C","L")))), error = function(e) NA_real_)
  c(auc_l, auc_s)
}
stopCluster(cl)

perm_auc_logit <- perm_res[,1]; perm_auc_svm <- perm_res[,2]
emp_p_logit <- mean(perm_auc_logit >= as.numeric(auc(roc_logit)), na.rm = TRUE)
emp_p_svm   <- mean(perm_auc_svm   >= as.numeric(auc(roc_svm)),   na.rm = TRUE)

cat("Permutation empirical p-value (Logit):", emp_p_logit, "\n")
cat("Permutation empirical p-value (SVM):", emp_p_svm, "\n\n")

# ---- Parallel bootstrap for LOOCV AUC CIs ----
cat("Running parallel bootstrap (nboot =", nboot, ") ...\n")
cl <- makeCluster(cores); registerDoParallel(cl)

boot_res <- foreach(b = 1:nboot, .combine = rbind, .packages = c("glmnet","e1071","pROC")) %dopar% {
  set.seed(3000 + b)
  idx <- sample(seq_len(n_samples), replace = TRUE)
  Xb <- X_full[idx, , drop = FALSE]
  yb <- y_full[idx]
  if (length(unique(yb)) < 2) return(c(NA_real_, NA_real_))
  res_b <- nested_loocv_fast(Xb, yb, svm_costs = svm_costs_small, inner_folds = inner_folds, seed = 3000 + b)
  auc_b_l <- tryCatch(as.numeric(auc(roc(response = yb, predictor = res_b$prob_logit, levels = c("C","L")))), error = function(e) NA_real_)
  auc_b_s <- tryCatch(as.numeric(auc(roc(response = yb, predictor = res_b$prob_svm,   levels = c("C","L")))), error = function(e) NA_real_)
  c(auc_b_l, auc_b_s)
}
stopCluster(cl)

boot_auc_logit <- na.omit(boot_res[,1]); boot_auc_svm <- na.omit(boot_res[,2])
ci_logit <- if (length(boot_auc_logit) > 0) quantile(boot_auc_logit, probs = c(0.025, 0.975)) else c(NA, NA)
ci_svm   <- if (length(boot_auc_svm) > 0)   quantile(boot_auc_svm,   probs = c(0.025, 0.975)) else c(NA, NA)

cat("Bootstrap 95% CI (Logit AUC):", round(ci_logit[1],3), "-", round(ci_logit[2],3), "\n")
cat("Bootstrap 95% CI (SVM AUC):",   round(ci_svm[1],3),   "-", round(ci_svm[2],3),   "\n\n")

# ---- Save permutation/bootstrap summaries ----
summary_list <- list(
  n_samples = n_samples,
  n_genes_detected = length(sig_present),
  auc_logit = as.numeric(auc(roc_logit)),
  auc_svm = as.numeric(auc(roc_svm)),
  emp_p_fast = emp_p_fast,
  perm_p_logit = emp_p_logit,
  perm_p_svm = emp_p_svm,
  boot_ci_logit = ci_logit,
  boot_ci_svm = ci_svm
)


saveRDS(summary_list, file = "validation_summary_fastPerm_nestedCV.rds")
cat("Saved summary to validation_summary_fastPerm_nestedCV.rds\n")

save.image("validation_nestedCV.rdata")

# End of replacement block
# -----------------------------
