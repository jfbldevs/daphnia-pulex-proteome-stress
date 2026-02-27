# =============================================================================
# Script 01: Data Preparation and Differential Expression Analysis
# Paper: Integrative bioinformatics analysis reveals disrupted metabolic
#        pathways and key hub proteins in the Daphnia pulex proteome
#        under anthropogenic stress
#
# This script loads raw proteinGroups data, performs filtering, calculates
# differential expression (t-test + BH correction), and generates
# volcano plot, MA plot, and hierarchical clustering heatmap.
#
# Authors: J.A. Norambuena, L. Herrera-Belén, P. Poblete-Grant, J.F. Beltrán
# Year:    2026
# =============================================================================

# Load configuration
source("00_setup.R")

# =============================================================================
# 1. LOAD AND CLEAN DATA
# =============================================================================

cat("=== LOADING DATA ===\n")

# Read data
protein_data <- read_csv(file.path(DATA_DIR, "proteinGroups_clean.csv"),
  show_col_types = FALSE
)

cat("Total proteins loaded:", nrow(protein_data), "\n")

# -----------------------------------------------------------------------------
# 1.1 Extraer UniProt IDs
# -----------------------------------------------------------------------------
# Extract the first ID from each protein group (Majority protein ID)
protein_data <- protein_data %>%
  mutate(
    # Extract first UniProt ID
    UniProt_ID = sapply(strsplit(`Majority protein IDs`, ";"), `[`, 1),

    # Extract protein name from FASTA header
    Protein_Name = sapply(`Fasta headers`, function(x) {
      if (is.na(x)) {
        return(NA)
      }
      # Extract part between | and OS=
      match <- regmatches(x, regexpr("\\|[^|]+\\|([^|]+)\\s+OS=", x))
      if (length(match) > 0) {
        name <- sub(".*\\|([^|]+)\\s+OS=.*", "\\1", match)
        return(trimws(name))
      }
      return(NA)
    }),

    # Clean name (remove code prefix)
    Protein_Name_Clean = sapply(Protein_Name, function(x) {
      if (is.na(x)) {
        return(NA)
      }
      # Remove pattern like "E9GMM9_DAPPU " from the start
      sub("^[A-Z0-9]+_DAPPU\\s+", "", x)
    })
  )

cat("UniProt IDs extracted:", sum(!is.na(protein_data$UniProt_ID)), "\n")

# =============================================================================
# 2. PREPARE EXPRESSION MATRIX
# =============================================================================

cat("\n=== PREPARING EXPRESSION MATRIX ===\n")

# LFQ intensity columns
icalma_cols <- c("LFQ intensity IDI-1", "LFQ intensity IDI-2", "LFQ intensity IDI-3")
llanquihue_cols <- c("LFQ intensity IDLL-1", "LFQ intensity IDLL-2", "LFQ intensity IDLL-3")

# Create expression matrix
expression_matrix <- protein_data %>%
  dplyr::select(UniProt_ID, all_of(c(icalma_cols, llanquihue_cols))) %>%
  column_to_rownames("UniProt_ID")

# Renombrar columnas para simplicidad
colnames(expression_matrix) <- c(
  "Icalma_1", "Icalma_2", "Icalma_3",
  "Llanquihue_1", "Llanquihue_2", "Llanquihue_3"
)

# Replace 0 with NA for analysis
expression_matrix[expression_matrix == 0] <- NA

cat("Expression matrix dimensions:", dim(expression_matrix), "\n")

# =============================================================================
# 3. DATA FILTERING
# =============================================================================

cat("\n=== FILTERING DATA ===\n")

# Criterion: Protein must have values in at least 2 of 3 replicates in at least one group
filter_proteins <- function(mat) {
  icalma_valid <- rowSums(!is.na(mat[, 1:3])) >= 2
  llanquihue_valid <- rowSums(!is.na(mat[, 4:6])) >= 2
  return(icalma_valid | llanquihue_valid)
}

valid_proteins <- filter_proteins(expression_matrix)
expression_filtered <- expression_matrix[valid_proteins, ]

cat("Proteins before filtering:", nrow(expression_matrix), "\n")
cat("Proteins after filtering:", nrow(expression_filtered), "\n")

# Proteins with data in BOTH groups (for direct comparison)
both_groups <- rowSums(!is.na(expression_filtered[, 1:3])) >= 2 &
  rowSums(!is.na(expression_filtered[, 4:6])) >= 2
expression_both <- expression_filtered[both_groups, ]

cat("Proteins with data in both lakes:", nrow(expression_both), "\n")

# =============================================================================
# 4. DIFFERENTIAL ANALYSIS
# =============================================================================

cat("\n=== DIFFERENTIAL ANALYSIS ===\n")

# Calculate statistics for proteins with data in both groups
differential_analysis <- protein_data %>%
  filter(UniProt_ID %in% rownames(expression_both)) %>%
  mutate(
    # Means
    Mean_Icalma = rowMeans(dplyr::select(., all_of(icalma_cols)), na.rm = TRUE),
    Mean_Llanquihue = rowMeans(dplyr::select(., all_of(llanquihue_cols)), na.rm = TRUE),

    # Standard deviations
    SD_Icalma = apply(dplyr::select(., all_of(icalma_cols)), 1, sd, na.rm = TRUE),
    SD_Llanquihue = apply(dplyr::select(., all_of(llanquihue_cols)), 1, sd, na.rm = TRUE),

    # Log2 Fold Change (Llanquihue vs Icalma)
    Log2FC = log2(Mean_Llanquihue / Mean_Icalma)
  )

# Calculate p-values with t-test
calculate_pvalue <- function(row, cols_a, cols_b) {
  a <- as.numeric(row[cols_a])
  b <- as.numeric(row[cols_b])
  a <- a[!is.na(a) & a > 0]
  b <- b[!is.na(b) & b > 0]

  if (length(a) < 2 | length(b) < 2) {
    return(NA)
  }

  tryCatch(
    {
      t.test(a, b)$p.value
    },
    error = function(e) NA
  )
}

# Apply t-test
differential_analysis$P_Value <- apply(
  differential_analysis[, c(icalma_cols, llanquihue_cols)],
  1,
  calculate_pvalue,
  cols_a = icalma_cols,
  cols_b = llanquihue_cols
)

# P-value adjustment (FDR)
differential_analysis$P_Adjusted <- p.adjust(differential_analysis$P_Value, method = "BH")

# Classify proteins with STANDARD criteria:
# - FDR-adjusted p-value < 0.05 (false positive control)
# - |Log2FC| >= 1 (at least 2-fold change)
differential_analysis <- differential_analysis %>%
  mutate(
    Regulation = case_when(
      P_Adjusted < 0.05 & Log2FC >= 1 ~ "Upregulated",
      P_Adjusted < 0.05 & Log2FC <= -1 ~ "Downregulated",
      TRUE ~ "Not significant"
    )
  )

# Summary
cat("\n--- Differential expression summary (FDR < 0.05 & |log2FC| >= 1) ---\n")
print(table(differential_analysis$Regulation))

# =============================================================================
# 5. CREATE PROTEIN LISTS FOR DOWNSTREAM ANALYSES
# =============================================================================

cat("\n=== CREATING PROTEIN LISTS ===\n")

# Upregulated proteins (p < 0.05)
upregulated <- differential_analysis %>%
  filter(Regulation == "Upregulated") %>%
  pull(UniProt_ID)

# Downregulated proteins (p < 0.05)
downregulated <- differential_analysis %>%
  filter(Regulation == "Downregulated") %>%
  pull(UniProt_ID)

# All differentially expressed proteins
DEPs <- c(upregulated, downregulated)

# Background proteins (all analyzed)
background <- differential_analysis$UniProt_ID

cat("Upregulated proteins:", length(upregulated), "\n")
cat("Downregulated proteins:", length(downregulated), "\n")
cat("Total DEPs:", length(DEPs), "\n")
cat("Background:", length(background), "\n")

# =============================================================================
# 6. SAVE RESULTS
# =============================================================================

cat("\n=== SAVING RESULTS ===\n")

# Save full differential analysis table
write_csv(differential_analysis, file.path(RESULTS_DIR, "Table1_Differential_Expression.csv"))

# Save protein ID lists for downstream analyses
write_lines(upregulated, file.path(RESULTS_DIR, "Data_Upregulated_IDs.txt"))
write_lines(downregulated, file.path(RESULTS_DIR, "Data_Downregulated_IDs.txt"))
write_lines(DEPs, file.path(RESULTS_DIR, "Data_All_DEPs_IDs.txt"))
write_lines(background, file.path(RESULTS_DIR, "Data_Background_IDs.txt"))

# Save filtered expression matrix
expression_both %>%
  rownames_to_column("UniProt_ID") %>%
  write_csv(file.path(RESULTS_DIR, "Data_Expression_Matrix.csv"))

cat("Files saved to:", RESULTS_DIR, "\n")

# =============================================================================
# 7. VISUALIZATIONS
# =============================================================================

cat("\n=== GENERATING VISUALIZATIONS ===\n")

# -----------------------------------------------------------------------------
# 7.1 Volcano Plot
# -----------------------------------------------------------------------------
volcano_plot <- ggplot(differential_analysis, aes(x = Log2FC, y = -log10(P_Value))) +
  geom_point(aes(color = Regulation), alpha = 0.6, size = 2) +
  scale_color_manual(values = regulation_colors) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray40") +
  labs(
    title = "Volcano Plot: Llanquihue vs Icalma",
    subtitle = paste("Upregulated:", length(upregulated), "| Downregulated:", length(downregulated)),
    x = expression(Log[2] ~ Fold ~ Change),
    y = expression(-Log[10] ~ P - value),
    color = "Regulation"
  ) +
  theme_paper +
  xlim(c(-8, 8))

ggsave(file.path(FIGURES_DIR, "Fig1A_Volcano_Plot.pdf"), volcano_plot,
  width = 8, height = 6, dpi = 300
)
ggsave(file.path(FIGURES_DIR, "Fig1A_Volcano_Plot.png"), volcano_plot,
  width = 8, height = 6, dpi = 300
)

# -----------------------------------------------------------------------------
# 7.2 MA Plot
# -----------------------------------------------------------------------------
differential_analysis <- differential_analysis %>%
  mutate(Average_Intensity = log10((Mean_Icalma + Mean_Llanquihue) / 2))

ma_plot <- ggplot(differential_analysis, aes(x = Average_Intensity, y = Log2FC)) +
  geom_point(aes(color = Regulation), alpha = 0.6, size = 2) +
  scale_color_manual(values = regulation_colors) +
  geom_hline(yintercept = 0, linetype = "solid", color = "gray40") +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "gray40") +
  labs(
    title = "MA Plot: Llanquihue vs Icalma",
    x = expression(Log[10] ~ Average ~ Intensity),
    y = expression(Log[2] ~ Fold ~ Change),
    color = "Regulation"
  ) +
  theme_paper

ggsave(file.path(FIGURES_DIR, "Fig1B_MA_Plot.pdf"), ma_plot,
  width = 8, height = 6, dpi = 300
)
ggsave(file.path(FIGURES_DIR, "Fig1B_MA_Plot.png"), ma_plot,
  width = 8, height = 6, dpi = 300
)

# -----------------------------------------------------------------------------
# 7.3 Heatmap of differentially expressed proteins
# -----------------------------------------------------------------------------
# Prepare matrix for heatmap
heatmap_data <- expression_both[DEPs, ]
heatmap_data <- log2(heatmap_data + 1) # Log transformation

# Impute NA with row median
for (i in 1:nrow(heatmap_data)) {
  row_median <- median(as.numeric(heatmap_data[i, ]), na.rm = TRUE)
  heatmap_data[i, is.na(heatmap_data[i, ])] <- row_median
}

# Remove rows with all NA
heatmap_data <- heatmap_data[complete.cases(heatmap_data), ]

# Column annotation
annotation_col <- data.frame(
  Lake = c(rep("Icalma", 3), rep("Llanquihue", 3))
)
rownames(annotation_col) <- colnames(heatmap_data)

# Annotation colors
ann_colors <- list(Lake = lake_colors)

# Generate heatmap
pdf(file.path(FIGURES_DIR, "Fig1C_Heatmap_DEPs.pdf"), width = 8, height = 10)
pheatmap(
  heatmap_data,
  scale = "row",
  clustering_method = "ward.D2",
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  show_rownames = FALSE,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  main = "Differentially Expressed Proteins\n(Llanquihue vs Icalma)",
  fontsize = 10
)
dev.off()

cat("Visualizations saved to:", FIGURES_DIR, "\n")

# =============================================================================
# 8. FINAL SUMMARY
# =============================================================================

cat("\n")
cat("=======================================================\n")
cat("               ANALYSIS SUMMARY\n")
cat("=======================================================\n")
cat("Total proteins in dataset:          ", nrow(protein_data), "\n")
cat("Filtered proteins (valid):          ", nrow(expression_filtered), "\n")
cat("Proteins in both lakes:             ", nrow(expression_both), "\n")
cat("-------------------------------------------------------\n")
cat("Upregulated proteins (p<0.05):      ", length(upregulated), "\n")
cat("Downregulated proteins (p<0.05):    ", length(downregulated), "\n")
cat("Total DEPs:                         ", length(DEPs), "\n")
cat("=======================================================\n")

# Save objects for downstream scripts
save(differential_analysis, upregulated, downregulated, DEPs, background,
  expression_both, expression_filtered, protein_data,
  file = file.path(RESULTS_DIR, "RData_01_DataPrep.RData")
)

cat("\nScript completed successfully!\n")
cat("Next step: Run 02_GO_Enrichment_ORA.R\n")
