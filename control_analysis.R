# ==============================================================================
# control analyses
# ==============================================================================

cat("\n===  control analyses ===\n")

library(dplyr)
library(tidyr)
library(ggplot2)

dir.create("controls", showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Helper: random gene-set control from a defined gene pool
# ------------------------------------------------------------------------------

random_gene_set_control_from_pool <- function(data,
                                              gene_pool,
                                              signature_size,
                                              model_type = "RF",
                                              nsets = 50,
                                              nround = 10,
                                              seed_base = 9000) {
  gene_pool <- intersect(gene_pool, colnames(data))
  
  if (length(gene_pool) < signature_size) {
    stop("Gene pool has fewer genes than signature_size.")
  }
  
  out <- vector("list", nsets)
  gene_lists <- vector("list", nsets)
  
  for (b in seq_len(nsets)) {
    set.seed(seed_base + b)
    
    genes <- sample(gene_pool, size = signature_size, replace = FALSE)
    
    fit <- tryCatch(
      run_repeated_classification(
        data = data,
        feature_genes = genes,
        model_type = model_type,
        nround = nround,
        seed_offset = seed_base + b * 10
      ),
      error = function(e) {
        message("Set ", b, " failed: ", conditionMessage(e))
        return(NULL)
      }
    )
    
    if (is.null(fit)) {
      out[[b]] <- data.frame(
        Set = b,
        GenePool = "non_circadian",
        N_features = length(genes),
        AUC = NA,
        Accuracy = NA,
        F1 = NA,
        MCC = NA
      )
    } else {
      out[[b]] <- data.frame(
        Set = b,
        GenePool = "non_circadian",
        N_features = length(genes),
        AUC = mean(fit$raw$AUC, na.rm = TRUE),
        Accuracy = mean(fit$raw$Accuracy, na.rm = TRUE),
        F1 = mean(fit$raw$F1, na.rm = TRUE),
        MCC = mean(fit$raw$MCC, na.rm = TRUE)
      )
    }
    
    gene_lists[[b]] <- data.frame(Set = b, Gene_ID = genes)
    
    cat(
      "Non-circadian random set", b, "/", nsets,
      "AUC =", round(out[[b]]$AUC, 3), "\n"
    )
  }
  
  list(
    performance = bind_rows(out),
    gene_lists = bind_rows(gene_lists)
  )
}

# ------------------------------------------------------------------------------
# Observed RF AUC for the Boruta 53-gene signature
# ------------------------------------------------------------------------------

observed_auc <- stage1_boruta$models$RF$summary %>%
  filter(Metric == "AUC") %>%
  pull(Mean)

cat("Observed Boruta RF AUC =", observed_auc, "\n")


# ==============================================================================
# 7.1 Time-label permutation control
# ==============================================================================

cat("\n--- 7.1 Time-label permutation control ---\n")

perm_control <- permutation_label_control(
  data = retina_onset$data,
  feature_genes = boruta_selected_genes,
  model_type = "RF",
  nperm = params$n_permutations,
  nround = params$control_nround
)

write.csv(
  perm_control,
  "controls/7_1_time_label_permutation_rf.csv",
  row.names = FALSE
)


# ==============================================================================
# 7.2 Random 53-gene set control
# ==============================================================================

cat("\n--- 7.2 Random 53-gene set control ---\n")

random_control <- random_gene_set_control(
  data = retina_onset$data,
  signature_size = length(boruta_selected_genes),
  model_type = "RF",
  nsets = params$n_random_gene_sets,
  nround = params$control_nround
)

write.csv(
  random_control,
  "controls/7_2_random_53gene_sets_rf.csv",
  row.names = FALSE
)

plot_control_auc(
  random_control,
  observed_auc,
  output_file = "controls/7_2_random_53gene_auc_distribution.pdf"
)


# ==============================================================================
# 7.3 Housekeeping/reference-gene check
# ==============================================================================

cat("\n--- 7.3 Housekeeping/reference-gene check ---\n")

housekeeping_features <- intersect(housekeeping_6, colnames(retina_onset$data))

cat("Number of housekeeping/reference genes found:",
    length(housekeeping_features), "\n")

# 7.3A Test whether these genes are time-associated
hk_expr <- retina_onset$data[, c("time_class", "Time_num", housekeeping_features), drop = FALSE]

hk_long <- hk_expr %>%
  pivot_longer(
    cols = all_of(housekeeping_features),
    names_to = "gene",
    values_to = "expression"
  )

hk_time_tests <- hk_long %>%
  group_by(gene) %>%
  do({
    fit <- lm(expression ~ factor(Time_num), data = .)
    an <- anova(fit)
    data.frame(p_time = an$`Pr(>F)`[1])
  }) %>%
  ungroup() %>%
  mutate(FDR_time = p.adjust(p_time, method = "BH")) %>%
  arrange(FDR_time)

write.csv(
  hk_time_tests,
  "controls/7_3_housekeeping_time_association.csv",
  row.names = FALSE
)

# 7.3B Optional: model performance using housekeeping/reference genes
if (length(housekeeping_features) >= 2) {
  hk_control <- housekeeping_gene_control(
    data = retina_onset$data,
    housekeeping_genes = housekeeping_features,
    model_type = "RF",
    nround = params$control_nround
  )
  
  write.csv(
    hk_control,
    "controls/7_3_housekeeping_rf_performance.csv",
    row.names = FALSE
  )
}


# ==============================================================================
# 7.4 Single-time-point versus combined-window analysis
# ==============================================================================

cat("\n--- 7.4 Single-time-point vs combined ZT8/12 window ---\n")

single_vs_window <- single_timepoint_window_analysis(
  data = retina_onset$data,
  feature_genes = boruta_selected_genes,
  model_type = "RF",
  nround = params$control_nround
)

write.csv(
  single_vs_window,
  "controls/7_4_single_timepoint_vs_ZT8_12_window_rf.csv",
  row.names = FALSE
)


# ==============================================================================
# 7.5 Non-circadian gene control
# ==============================================================================

cat("\n--- 7.5 Non-circadian gene control ---\n")

# Create all-gene time-association table if it does not already exist
if (!exists("all_gene_time_tests")) {
  
  gene_cols <- grep("^ENSGALG", colnames(retina_onset$data), value = TRUE)
  
  all_gene_time_tests <- lapply(gene_cols, function(gene) {
    df <- data.frame(
      expression = retina_onset$data[[gene]],
      Time_num = retina_onset$data$Time_num
    )
    
    fit <- lm(expression ~ factor(Time_num), data = df)
    an <- anova(fit)
    
    data.frame(
      gene = gene,
      p_time = an$`Pr(>F)`[1],
      stringsAsFactors = FALSE
    )
  }) %>%
    bind_rows() %>%
    mutate(FDR_time = p.adjust(p_time, method = "BH")) %>%
    arrange(FDR_time)
}

write.csv(
  all_gene_time_tests,
  "controls/7_5_all_gene_time_association.csv",
  row.names = FALSE
)

# Select weakly time-associated genes as non-circadian control pool
non_circadian_pool <- all_gene_time_tests %>%
  filter(!is.na(FDR_time), FDR_time >= 0.8) %>%
  pull(gene) %>%
  intersect(colnames(retina_onset$data))

cat("Number of non-circadian control genes:", length(non_circadian_pool), "\n")

write.csv(
  data.frame(Gene_ID = non_circadian_pool),
  "controls/7_5_non_circadian_gene_pool.csv",
  row.names = FALSE
)

# Run only if enough genes are available
if (length(non_circadian_pool) >= length(boruta_selected_genes)) {
  
  non_circadian_control <- random_gene_set_control_from_pool(
    data = retina_onset$data,
    gene_pool = non_circadian_pool,
    signature_size = length(boruta_selected_genes),
    model_type = "RF",
    nsets = params$n_random_gene_sets,
    nround = params$control_nround
  )
  
  write.csv(
    non_circadian_control$performance,
    "controls/7_5_non_circadian_53gene_sets_rf.csv",
    row.names = FALSE
  )
  
  write.csv(
    non_circadian_control$gene_lists,
    "controls/7_5_non_circadian_53gene_gene_lists.csv",
    row.names = FALSE
  )
  
  # Simple figure
  pdf("controls/7_5_non_circadian_53gene_auc_distribution.pdf",
      width = 6, height = 4)
  
  print(
    ggplot(non_circadian_control$performance, aes(x = AUC)) +
      geom_histogram(bins = 20, color = "white") +
      geom_vline(xintercept = observed_auc, linetype = "dashed", linewidth = 1) +
      labs(
        title = "Non-circadian 53-gene control",
        subtitle = "Dashed line = observed Boruta 53-gene RF AUC",
        x = "Mean AUC from non-circadian 53-gene sets",
        y = "Frequency"
      ) +
      theme_bw()
  )
  
  dev.off()
  
  non_circadian_summary <- non_circadian_control$performance %>%
    summarise(
      Mean_AUC = mean(AUC, na.rm = TRUE),
      SD_AUC = sd(AUC, na.rm = TRUE),
      Median_AUC = median(AUC, na.rm = TRUE),
      Min_AUC = min(AUC, na.rm = TRUE),
      Max_AUC = max(AUC, na.rm = TRUE),
      Mean_MCC = mean(MCC, na.rm = TRUE),
      SD_MCC = sd(MCC, na.rm = TRUE)
    )
  
  write.csv(
    non_circadian_summary,
    "controls/7_5_non_circadian_summary.csv",
    row.names = FALSE
  )
  
} else {
  warning("Not enough non-circadian genes for 53-gene control.")
}

save.image("control_analysis.rdata")
cat("\n===  control analyses completed ===\n")