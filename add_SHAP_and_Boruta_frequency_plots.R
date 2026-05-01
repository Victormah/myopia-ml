
library(dplyr)
library(ggplot2)
library(randomForest)

if (!requireNamespace("iml", quietly = TRUE)) {
  stop("Package 'iml' is required. Install with: install.packages('iml')")
}

out_dir <- "reviewer3_feature_stability_SHAP"
dir.create(out_dir, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1) Boruta selection-frequency plot
# ------------------------------------------------------------------------------

boruta_freq_file <- "boruta_summary_retina_onset_timeclass.csv"

if (file.exists(boruta_freq_file)) {
  boruta_freq <- read.csv(boruta_freq_file, stringsAsFactors = FALSE)
} else if (file.exists("boruta_selection_frequency_retina_onset_expression_only.csv")) {
  boruta_freq <- read.csv("boruta_selection_frequency_retina_onset_expression_only.csv", stringsAsFactors = FALSE)
} else {
  stop("Boruta selection-frequency file not found.")
}

# Standardize column names if needed
if (!"Gene" %in% colnames(boruta_freq)) {
  names(boruta_freq)[1] <- "Gene"
}
if (!"Frequency" %in% colnames(boruta_freq)) {
  names(boruta_freq)[2] <- "Frequency"
}

boruta_freq <- boruta_freq %>%
  mutate(
    Consensus_50of50 = Frequency == 50,
    Selection_rate = Frequency / 50
  ) %>%
  arrange(desc(Frequency), Gene)

write.csv(
  boruta_freq,
  file.path(out_dir, "SuppData_Boruta_selection_frequency.csv"),
  row.names = FALSE
)

# Plot top genes by selection frequency. This shows both 50/50 and near-consensus genes.
plot_n <- min(80, nrow(boruta_freq))
boruta_freq_plot <- boruta_freq %>%
  slice_head(n = plot_n) %>%
  mutate(Gene = factor(Gene, levels = rev(Gene)))

p_freq <- ggplot(boruta_freq_plot, aes(x = Frequency, y = Gene, fill = Consensus_50of50)) +
  geom_col(width = 0.75) +
  geom_vline(xintercept = 50, linetype = "dashed") +
  scale_fill_manual(values = c("grey70", "black"), labels = c("<50/50", "50/50")) +
  labs(
    title = "Boruta feature-selection stability",
    subtitle = "Genes selected across 50 repeated Boruta runs",
    x = "Selection frequency out of 50 runs",
    y = "Gene ID",
    fill = "Consensus"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.y = element_text(size = 6))

print(p_freq)

ggsave(
  file.path(out_dir, "SuppFig_Boruta_selection_frequency.pdf"),
  p_freq, width = 7, height = 10
)

ggsave(
  file.path(out_dir, "SuppFig_Boruta_selection_frequency.png"),
  p_freq, width = 7, height = 10, dpi = 300
)

# ------------------------------------------------------------------------------
# 2) Train final RF model using expression-only Boruta genes
# ------------------------------------------------------------------------------

label_col <- "time_class"
if (!label_col %in% colnames(logCPM_retina)) {
  stop("Column 'time_class' was not found in logCPM_retina.")
}

boruta_features <- intersect(boruta_selected_genes, colnames(logCPM_retina))
if (length(boruta_features) < 2) {
  stop("Fewer than two Boruta-selected genes were found in logCPM_retina.")
}

# Check metadata leakage
bad_cols <- c("Time_num", "Time_rad", "sin_time", "cos_time", "time_class",
              "ZT", "ZT_class", "SampleID", "Sample_ID", "Sex", "eye_treatment")
leaked <- intersect(bad_cols, boruta_features)
if (length(leaked) > 0) {
  stop("Metadata leakage detected in feature list: ", paste(leaked, collapse = ", "))
}

X <- logCPM_retina[, boruta_features, drop = FALSE]
y <- factor(logCPM_retina[[label_col]])

# Choose positive class automatically
# Prefer ZT_812 if present, otherwise ZT8_12, otherwise the first level.
if ("ZT_812" %in% levels(y)) {
  positive_class <- "ZT_812"
} else if ("ZT8_12" %in% levels(y)) {
  positive_class <- "ZT8_12"
} else {
  positive_class <- levels(y)[1]
}

# Put positive class first; randomForest keeps this for probability column names.
y <- factor(y, levels = c(positive_class, setdiff(levels(y), positive_class)))

set.seed(123)
rf_final <- randomForest::randomForest(
  x = X,
  y = y,
  ntree = 1000,
  importance = TRUE
)

saveRDS(rf_final, file.path(out_dir, "RF_final_Boruta53_model.rds"))

# Built-in RF importance, useful as a quick companion to SHAP
rf_imp <- randomForest::importance(rf_final)
rf_imp_df <- data.frame(
  Gene = rownames(rf_imp),
  rf_imp,
  row.names = NULL,
  check.names = FALSE
) %>%
  arrange(desc(MeanDecreaseGini))

write.csv(
  rf_imp_df,
  file.path(out_dir, "SuppData_RF_feature_importance.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 3) SHAP analysis using iml::Shapley
# ------------------------------------------------------------------------------

# Prediction function returning probability for the positive class
pred_fun <- function(model, newdata) {
  as.numeric(predict(model, newdata = newdata, type = "prob")[, positive_class])
}

predictor <- iml::Predictor$new(
  model = rf_final,
  data = X,
  y = y,
  predict.function = pred_fun,
  class = positive_class
)

# To keep runtime reasonable, compute SHAP for a balanced subset of samples.
# Increase n_per_class if you want smoother plots.
n_per_class <- 10
set.seed(456)
idx_by_class <- split(seq_len(nrow(X)), y)
shap_idx <- unlist(lapply(idx_by_class, function(idx) {
  sample(idx, size = min(n_per_class, length(idx)))
}))

X_shap <- X[shap_idx, , drop = FALSE]
y_shap <- y[shap_idx]

shap_list <- vector("list", nrow(X_shap))

for (i in seq_len(nrow(X_shap))) {
  cat("Computing SHAP", i, "/", nrow(X_shap), "\n")
  sh <- iml::Shapley$new(
    predictor,
    x.interest = X_shap[i, , drop = FALSE]
  )
  tmp <- sh$results
  tmp$Sample_index <- rownames(X_shap)[i]
  tmp$Class <- as.character(y_shap[i])
  shap_list[[i]] <- tmp
}

shap_values <- bind_rows(shap_list)

# Try to extract numeric feature values from iml's feature.value column
shap_values <- shap_values %>%
  mutate(
    feature_value_numeric = suppressWarnings(as.numeric(sub("^.*=", "", feature.value)))
  )

write.csv(
  shap_values,
  file.path(out_dir, "SuppData_SHAP_values_sample_subset.csv"),
  row.names = FALSE
)

# Global mean absolute SHAP importance
shap_importance <- shap_values %>%
  group_by(feature) %>%
  summarise(
    mean_abs_SHAP = mean(abs(phi), na.rm = TRUE),
    mean_SHAP = mean(phi, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_abs_SHAP))

write.csv(
  shap_importance,
  file.path(out_dir, "SuppData_SHAP_global_importance.csv"),
  row.names = FALSE
)

# Plot top SHAP features
n_top_shap <- min(20, nrow(shap_importance))
top_shap_genes <- shap_importance$feature[seq_len(n_top_shap)]

p_shap_bar <- shap_importance %>%
  slice_head(n = n_top_shap) %>%
  mutate(feature = factor(feature, levels = rev(feature))) %>%
  ggplot(aes(x = mean_abs_SHAP, y = feature)) +
  geom_col() +
  labs(
    title = "SHAP global feature importance",
    subtitle = paste0("Mean absolute SHAP value for predicting ", positive_class),
    x = "Mean absolute SHAP value",
    y = "Gene ID"
  ) +
  theme_bw(base_size = 12)

ggsave(
  file.path(out_dir, "SuppFig_SHAP_global_importance_top20.pdf"),
  p_shap_bar, width = 7, height = 5
)

ggsave(
  file.path(out_dir, "SuppFig_SHAP_global_importance_top20.png"),
  p_shap_bar, width = 7, height = 5, dpi = 300
)

# SHAP summary dot plot for top genes
p_shap_dot <- shap_values %>%
  filter(feature %in% top_shap_genes) %>%
  mutate(feature = factor(feature, levels = rev(top_shap_genes))) %>%
  ggplot(aes(x = phi, y = feature, color = feature_value_numeric)) +
  geom_jitter(height = 0.18, width = 0, alpha = 0.8, size = 1.6) +
  labs(
    title = "SHAP summary plot",
    subtitle = paste0("Top ", n_top_shap, " Boruta genes; prediction target = ", positive_class),
    x = "SHAP value",
    y = "Gene ID",
    color = "Expression"
  ) +
  theme_bw(base_size = 12)

ggsave(
  file.path(out_dir, "SuppFig_SHAP_summary_dotplot_top20.pdf"),
  p_shap_dot, width = 7, height = 6
)

ggsave(
  file.path(out_dir, "SuppFig_SHAP_summary_dotplot_top20.png"),
  p_shap_dot, width = 7, height = 6, dpi = 300
)

cat("\nSHAP and Boruta stability outputs saved to:", out_dir, "\n")

save.image("add_SHAP_and_Boruta_frequency_plots.rdata")
