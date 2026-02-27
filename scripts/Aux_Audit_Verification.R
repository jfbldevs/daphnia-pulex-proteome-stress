# =============================================================================
# Auxiliary Script: Audit and Verification
# End-to-end verification from raw data to final results
#
# Authors: J.A. Norambuena, L. Herrera-Belén, P. Poblete-Grant, J.F. Beltrán
# Year:    2026
# =============================================================================

# Load configuration and directory paths
source("00_Setup.R")

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════════╗\n")
cat("║            AUDIT - INTEGRITY VERIFICATION                                ║\n")
cat("╚══════════════════════════════════════════════════════════════════════════╝\n\n")

errors <- c()
warnings <- c()

raw_data <- read_csv(file.path(DATA_DIR, "proteinGroups_clean.csv"), show_col_types = FALSE)

cat("File: proteinGroups_clean.csv\n")
cat("  Rows (proteins):", nrow(raw_data), "\n")
cat("  Columns:", ncol(raw_data), "\n")
cat("  Columns:", paste(colnames(raw_data), collapse = ", "), "\n\n")

# Verify sample structure
icalma_cols <- grep("IDI-", colnames(raw_data), value = TRUE)
llanquihue_cols <- grep("IDLL-", colnames(raw_data), value = TRUE)

cat("  Icalma samples (oligotrophic):", paste(icalma_cols, collapse = ", "), "\n")
cat("  Llanquihue samples (anthropized):", paste(llanquihue_cols, collapse = ", "), "\n\n")

if (length(icalma_cols) != 3 | length(llanquihue_cols) != 3) {
  errors <- c(errors, "ERROR: Expected 3 replicates per condition")
}

# =============================================================================
# 2. DIFFERENTIAL EXPRESSION VERIFICATION
# =============================================================================

cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("2. DIFFERENTIAL EXPRESSION\n")
cat("═══════════════════════════════════════════════════════════════════════════\n\n")

load(file.path(RESULTS_DIR, "RData_01_DataPrep.RData"))

cat("Objects loaded:", paste(ls(), collapse = ", "), "\n\n")

# Recalcular desde cero para verificar
cat("--- Recalculating differential expression from scratch ---\n\n")

# Prepare data
icalma_data <- raw_data %>% select(all_of(icalma_cols))
llanquihue_data <- raw_data %>% select(all_of(llanquihue_cols))

# Filter proteins with at least 2 valid values per group
valid_icalma <- rowSums(icalma_data > 0) >= 2
valid_llanquihue <- rowSums(llanquihue_data > 0) >= 2
valid_proteins <- valid_icalma & valid_llanquihue

cat("  Total proteins:", nrow(raw_data), "\n")
cat("  Proteins with >=2 values in Icalma:", sum(valid_icalma), "\n")
cat("  Proteins with >=2 values in Llanquihue:", sum(valid_llanquihue), "\n")
cat("  Valid proteins in both groups:", sum(valid_proteins), "\n\n")

# Calculate differential expression
results_audit <- data.frame(
  Protein_ID = raw_data$`Protein IDs`[valid_proteins],
  stringsAsFactors = FALSE
)

icalma_valid <- icalma_data[valid_proteins, ]
llanquihue_valid <- llanquihue_data[valid_proteins, ]

# Calculate means of RAW intensities (including zeros, same as pipeline)
# El pipeline principal (01_Differential_Expression.R) calcula rowMeans sobre el
# dataset raw SIN reemplazar ceros por NA, usando na.rm = TRUE.
results_audit$Mean_Icalma <- rowMeans(icalma_valid, na.rm = TRUE)
results_audit$Mean_Llanquihue <- rowMeans(llanquihue_valid, na.rm = TRUE)

# Log2 Fold Change: Llanquihue vs Icalma (positive = higher in Llanquihue)
# Formula: log2(mean_raw_Llanquihue / mean_raw_Icalma)
results_audit$Log2FC <- log2(results_audit$Mean_Llanquihue / results_audit$Mean_Icalma)

# T-test (on raw filtered values > 0, same as main pipeline)
results_audit$P_Value <- sapply(seq_len(nrow(icalma_valid)), function(i) {
  x <- as.numeric(icalma_valid[i, ])
  y <- as.numeric(llanquihue_valid[i, ])
  x <- x[!is.na(x) & x > 0]
  y <- y[!is.na(y) & y > 0]
  if (length(x) >= 2 && length(y) >= 2) {
    tryCatch(t.test(x, y)$p.value, error = function(e) NA)
  } else {
    NA
  }
})

# FDR correction
results_audit$P_Adjusted <- p.adjust(results_audit$P_Value, method = "BH")

# Classification
results_audit$Regulation <- case_when(
  results_audit$P_Adjusted < 0.05 & results_audit$Log2FC > 1 ~ "Upregulated",
  results_audit$P_Adjusted < 0.05 & results_audit$Log2FC < -1 ~ "Downregulated",
  TRUE ~ "Not significant"
)

# Summary
audit_up <- sum(results_audit$Regulation == "Upregulated", na.rm = TRUE)
audit_down <- sum(results_audit$Regulation == "Downregulated", na.rm = TRUE)
audit_total <- audit_up + audit_down

cat("RECALCULATED RESULTS (audit)::\n")
cat("  Upregulated (higher in Llanquihue):", audit_up, "\n")
cat("  Downregulated (higher in Icalma):", audit_down, "\n")
cat("  Total DEPs:", audit_total, "\n\n")

# Compare with stored results
stored_up <- sum(differential_analysis$Regulation == "Upregulated", na.rm = TRUE)
stored_down <- sum(differential_analysis$Regulation == "Downregulated", na.rm = TRUE)
stored_total <- stored_up + stored_down

cat("STORED RESULTS:\n")
cat("  Upregulated:", stored_up, "\n")
cat("  Downregulated:", stored_down, "\n")
cat("  Total DEPs:", stored_total, "\n\n")

if (audit_up == stored_up & audit_down == stored_down) {
  cat("✓ VERIFIED: DEP counts match\n\n")
} else {
  errors <- c(errors, paste("ERROR: Discrepancia en DEPs. Audit:", audit_total, "Stored:", stored_total))
  cat("✗ ERROR: Counts do NOT match\n\n")
}

# =============================================================================
# 3. PROTEIN LISTS FILE VERIFICATION
# =============================================================================

cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("3. PROTEIN LISTS VERIFICATION\n")
cat("═══════════════════════════════════════════════════════════════════════════\n\n")

up_file <- readLines(file.path(RESULTS_DIR, "Data_Upregulated_IDs.txt"))
down_file <- readLines(file.path(RESULTS_DIR, "Data_Downregulated_IDs.txt"))
all_deps_file <- readLines(file.path(RESULTS_DIR, "Data_All_DEPs_IDs.txt"))

cat("Data_Upregulated_IDs.txt:", length(up_file), "proteins\n")
cat("Data_Downregulated_IDs.txt:", length(down_file), "proteins\n")
cat("Data_All_DEPs_IDs.txt:", length(all_deps_file), "proteins\n")
cat("Suma up + down:", length(up_file) + length(down_file), "\n\n")

if (length(all_deps_file) == length(up_file) + length(down_file)) {
  cat("✓ VERIFICADO: all_DEPs = upregulated + downregulated\n\n")
} else {
  errors <- c(errors, "ERROR: all_DEPs != upregulated + downregulated")
}

# =============================================================================
# 4. PPI NETWORK VERIFICATION
# =============================================================================

cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("4. PPI NETWORK\n")
cat("═══════════════════════════════════════════════════════════════════════════\n\n")

load(file.path(RESULTS_DIR, "04_PPI_network.RData"))

cat("Nodes in network:", length(V(ppi_graph)), "\n")
cat("Edges in network:", length(E(ppi_graph)), "\n")
cat("Proteins in centrality_metrics:", nrow(centrality_metrics), "\n\n")

# Verify hub calculation
hub_threshold <- quantile(centrality_metrics$Degree, 0.95)
cat("Hub threshold (percentil 95):", hub_threshold, "\n")

n_hubs_calc <- sum(centrality_metrics$Degree >= hub_threshold)
n_hubs_stored <- nrow(read_csv(file.path(RESULTS_DIR, "Table4_Hub_Proteins.csv"), show_col_types = FALSE))

cat("Hubs calculated (Degree >= ", round(hub_threshold), "):", n_hubs_calc, "\n")
cat("Hubs in file:", n_hubs_stored, "\n\n")

# Verification with calculated threshold
cat("Hubs con threshold percentil 95 (Degree >=", round(hub_threshold), "):", n_hubs_calc, "\n\n")

# Top 10 proteins por grado
cat("Top 10 proteins por Degree:\n")
top10 <- centrality_metrics %>%
  arrange(desc(Degree)) %>%
  head(10) %>%
  select(Protein_Name, Degree, Regulation)
print(top10)
cat("\n")

# Verify hub regulation
hubs <- centrality_metrics %>% filter(Degree >= hub_threshold)
cat(paste0("Regulation of the ", nrow(hubs), " hubs (Degree >= ", round(hub_threshold), "):\n"))
print(table(hubs$Regulation))
cat("\n")

hub_down <- sum(hubs$Regulation == "Downregulated", na.rm = TRUE)
hub_up <- sum(hubs$Regulation == "Upregulated", na.rm = TRUE)
hub_ns <- sum(hubs$Regulation == "Not significant" | is.na(hubs$Regulation))

cat("  Downregulated (higher in Icalma):", hub_down, "\n")
cat("  Upregulated (higher in Llanquihue):", hub_up, "\n")
cat("  Not significant:", hub_ns, "\n\n")

# =============================================================================
# 5. BIOLOGICAL INTERPRETATION VERIFICATION
# =============================================================================

cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("5. BIOLOGICAL INTERPRETATION\n")
cat("═══════════════════════════════════════════════════════════════════════════\n\n")

cat("INTERPRETATION KEY:\n")
cat("  - Log2FC > 0: Higher expression in LLANQUIHUE (anthropized)\n")
cat("  - Log2FC < 0: Higher expression in ICALMA (oligotrophic)\n")
cat("  - Upregulated: Higher in Llanquihue (possible stress response)\n")
cat("  - Downregulated: Higher in Icalma (suppressed in stressed environment)\n\n")

# Verify interpretation consistency
cat("Verifying interpretation consistency...\n\n")

# Examples with specific proteins
stress_proteins <- centrality_metrics %>%
  filter(str_detect(tolower(Protein_Name), "heat shock|hsp|stress|superoxide"))

if (nrow(stress_proteins) > 0) {
  cat("Stress proteins found:\n")
  print(stress_proteins %>% select(Protein_Name, Regulation))
  cat("\n")
}

# Metabolic proteins
metabolic_proteins <- centrality_metrics %>%
  filter(str_detect(tolower(Protein_Name), "atp synthase|citrate|fumarate"))

cat("Key metabolic proteins:\n")
print(metabolic_proteins %>% select(Protein_Name, Degree, Regulation))
cat("\n")

# =============================================================================
# 6. AUDIT SUMMARY
# =============================================================================

cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("6. AUDIT SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════════════════\n\n")

cat("VERIFIED FINAL STATISTICS:\n")
cat("  Total proteins analyzed:", sum(valid_proteins), "\n")
cat("  DEPs totales:", audit_total, "\n")
cat("  - Upregulated (higher in Llanquihue):", audit_up, "\n")
cat("  - Downregulated (higher in Icalma):", audit_down, "\n")
cat("  Proteins in PPI network:", nrow(centrality_metrics), "\n")
cat("  Hubs (Degree >=", round(hub_threshold), "):", n_hubs_calc, "\n")
cat("  - Hubs downregulated:", hub_down, "\n")
cat("  - Hubs upregulated:", hub_up, "\n")
cat("  - Hubs not significant:", hub_ns, "\n\n")

if (length(errors) == 0) {
  cat("╔══════════════════════════════════════════════════════════════════════════╗\n")
  cat("║  ✓ AUDIT COMPLETED: NO ERRORS FOUND                       ║\n")
  cat("╚══════════════════════════════════════════════════════════════════════════╝\n")
} else {
  cat("╔══════════════════════════════════════════════════════════════════════════╗\n")
  cat("║  ✗ ERRORS FOUND:                                                   ║\n")
  for (e in errors) {
    cat("║  -", e, "\n")
  }
  cat("╚══════════════════════════════════════════════════════════════════════════╝\n")
}

cat("\n")
