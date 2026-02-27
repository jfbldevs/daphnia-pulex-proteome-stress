# =============================================================================
# Script 03: GSEA (Gene Set Enrichment Analysis)
# Enrichment analysis using ALL ranked proteins (not just DEPs)
#
# Unlike ORA (Script 02), GSEA does not require an arbitrary cutoff.
# It uses the full ranked list of proteins sorted by Log2FC, providing
# greater statistical power to detect coordinated expression changes.
#
# Authors: J.A. Norambuena, L. Herrera-Belén, P. Poblete-Grant, J.F. Beltrán
# Year:    2026
# =============================================================================

# Load configuration and directory paths
source("00_Setup.R")

cat("\n")
cat("====================================================================\n")
cat("         GSEA - Gene Set Enrichment Analysis\n")
cat("====================================================================\n\n")

# Load data
diff_results <- read_csv(file.path(RESULTS_DIR, "Table1_Differential_Expression.csv"),
  show_col_types = FALSE
)
term2gene_bp <- read_csv(file.path(RESULTS_DIR, "Data_Term2Gene_BP.csv"),
  show_col_types = FALSE
)
term2gene_mf <- read_csv(file.path(RESULTS_DIR, "Data_Term2Gene_MF.csv"),
  show_col_types = FALSE
)
term2gene_cc <- read_csv(file.path(RESULTS_DIR, "Data_Term2Gene_CC.csv"),
  show_col_types = FALSE
)

cat("Total proteins:", nrow(diff_results), "\n")
cat("GO BP terms:", n_distinct(term2gene_bp$TERM), "\n")
cat("GO MF terms:", n_distinct(term2gene_mf$TERM), "\n")
cat("GO CC terms:", n_distinct(term2gene_cc$TERM), "\n\n")

# =============================================================================
# 1. PREPARE RANKED LIST FOR GSEA
# =============================================================================

cat("=== PREPARING RANKED LIST ===\n\n")

# Create gene list ranked by Log2FC
# GSEA uses the ranking, not an arbitrary cutoff
gene_list <- diff_results %>%
  filter(!is.na(Log2FC) & is.finite(Log2FC)) %>%
  arrange(desc(Log2FC)) %>%
  pull(Log2FC, name = UniProt_ID)

cat("Proteins with valid Log2FC:", length(gene_list), "\n")
cat("Log2FC range:", round(min(gene_list), 2), "to", round(max(gene_list), 2), "\n")
cat("Median Log2FC:", round(median(gene_list), 2), "\n\n")

# =============================================================================
# 2. GSEA FUNCTION
# =============================================================================

run_gsea <- function(gene_list, term2gene, category_name) {
  cat(paste0("Running GSEA for ", category_name, "...\n"))

  # Filter term2gene for genes in our list
  term2gene_filtered <- term2gene %>%
    filter(GENE %in% names(gene_list))

  cat(paste0("  Gene-term associations: ", nrow(term2gene_filtered), "\n"))
  cat(paste0("  Unique terms: ", n_distinct(term2gene_filtered$TERM), "\n"))

  if (nrow(term2gene_filtered) < 50) {
    cat("  Insufficient associations for GSEA\n\n")
    return(NULL)
  }

  result <- tryCatch(
    {
      GSEA(
        geneList = gene_list,
        TERM2GENE = term2gene_filtered,
        pvalueCutoff = 0.25, # GSEA standard threshold
        minGSSize = 5,
        maxGSSize = 500,
        eps = 0,
        pAdjustMethod = "BH",
        verbose = FALSE
      )
    },
    error = function(e) {
      cat(paste0("  Error: ", e$message, "\n\n"))
      return(NULL)
    }
  )

  if (!is.null(result) && nrow(result@result) > 0) {
    # Add GO descriptions
    result@result$Description <- sapply(result@result$ID, function(go_id) {
      tryCatch(
        {
          term <- GO.db::GOTERM[[go_id]]
          if (!is.null(term)) Term(term) else go_id
        },
        error = function(e) go_id
      )
    })

    sig_up <- sum(result@result$p.adjust < 0.05 & result@result$NES > 0)
    sig_down <- sum(result@result$p.adjust < 0.05 & result@result$NES < 0)

    cat(paste0("  Enriched terms (p.adj < 0.05):\n"))
    cat(paste0("    - Positive NES (enriched in Llanquihue): ", sig_up, "\n"))
    cat(paste0("    - Negative NES (enriched in Icalma): ", sig_down, "\n\n"))
  } else {
    cat("  No significant results\n\n")
  }

  return(result)
}

# =============================================================================
# 3. RUN GSEA FOR EACH GO CATEGORY
# =============================================================================

cat("=== RUNNING GSEA ===\n\n")

gsea_bp <- run_gsea(gene_list, term2gene_bp, "Biological Process")
gsea_mf <- run_gsea(gene_list, term2gene_mf, "Molecular Function")
gsea_cc <- run_gsea(gene_list, term2gene_cc, "Cellular Component")

# =============================================================================
# 4. VISUALIZATIONS
# =============================================================================

cat("=== GENERATING VISUALIZATIONS ===\n\n")

# Function to create GSEA dotplot
create_gsea_dotplot <- function(gsea_result, title, top_n = 20) {
  if (is.null(gsea_result) || nrow(gsea_result@result) == 0) {
    return(NULL)
  }

  df <- gsea_result@result %>%
    filter(p.adjust < 0.25) %>%
    arrange(p.adjust) %>%
    head(top_n) %>%
    mutate(
      Description = str_trunc(Description, 45),
      Direction = ifelse(NES > 0, "Llanquihue", "Icalma")
    )

  if (nrow(df) == 0) {
    return(NULL)
  }

  ggplot(df, aes(
    x = NES, y = reorder(Description, NES),
    color = p.adjust, size = setSize
  )) +
    geom_point() +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    scale_color_gradient(
      low = "#E63946", high = "#457B9D",
      name = "Adj. P-value"
    ) +
    scale_size_continuous(range = c(3, 8), name = "Gene Set\nSize") +
    labs(
      title = title,
      subtitle = "NES > 0: Enriched in Llanquihue | NES < 0: Enriched in Icalma",
      x = "Normalized Enrichment Score (NES)",
      y = NULL
    ) +
    theme_paper +
    theme(
      axis.text.y = element_text(size = 9),
      plot.subtitle = element_text(size = 9, color = "gray40")
    )
}

# Create plots
if (!is.null(gsea_bp)) {
  p_bp <- create_gsea_dotplot(gsea_bp, "GSEA: GO Biological Process")
  if (!is.null(p_bp)) {
    ggsave(file.path(FIGURES_DIR, "FigS2_GSEA_BP.pdf"), p_bp, width = 11, height = 9)
    ggsave(file.path(FIGURES_DIR, "FigS2_GSEA_BP.png"), p_bp,
      width = 11, height = 9, dpi = 600, bg = "white"
    )
    cat("Saved: FigS2_GSEA_BP.pdf/png\n")
  }
}

if (!is.null(gsea_mf)) {
  p_mf <- create_gsea_dotplot(gsea_mf, "GSEA: GO Molecular Function")
  if (!is.null(p_mf)) {
    ggsave(file.path(FIGURES_DIR, "FigS3_GSEA_MF.pdf"), p_mf, width = 11, height = 9)
    ggsave(file.path(FIGURES_DIR, "FigS3_GSEA_MF.png"), p_mf,
      width = 11, height = 9, dpi = 600, bg = "white"
    )
    cat("Saved: FigS3_GSEA_MF.pdf/png\n")
  }
}

if (!is.null(gsea_cc)) {
  p_cc <- create_gsea_dotplot(gsea_cc, "GSEA: GO Cellular Component")
  if (!is.null(p_cc)) {
    ggsave(file.path(FIGURES_DIR, "FigS4_GSEA_CC.pdf"), p_cc, width = 11, height = 9)
    ggsave(file.path(FIGURES_DIR, "FigS4_GSEA_CC.png"), p_cc,
      width = 11, height = 9, dpi = 600, bg = "white"
    )
    cat("Saved: FigS4_GSEA_CC.pdf/png\n")
  }
}

# =============================================================================
# 5. SAVE RESULTS
# =============================================================================

cat("\n=== SAVING RESULTS ===\n\n")

if (!is.null(gsea_bp) && nrow(gsea_bp@result) > 0) {
  write_csv(gsea_bp@result, file.path(RESULTS_DIR, "TableS2_GSEA_BP.csv"))
  cat("Saved: TableS2_GSEA_BP.csv\n")
}

if (!is.null(gsea_mf) && nrow(gsea_mf@result) > 0) {
  write_csv(gsea_mf@result, file.path(RESULTS_DIR, "TableS3_GSEA_MF.csv"))
  cat("Saved: TableS3_GSEA_MF.csv\n")
}

if (!is.null(gsea_cc) && nrow(gsea_cc@result) > 0) {
  write_csv(gsea_cc@result, file.path(RESULTS_DIR, "TableS4_GSEA_CC.csv"))
  cat("Saved: TableS4_GSEA_CC.csv\n")
}

# Save R objects
save(gsea_bp, gsea_mf, gsea_cc, gene_list,
  file = file.path(RESULTS_DIR, "RData_10_GSEA.RData")
)

# =============================================================================
# 6. SUMMARY
# =============================================================================

cat("\n")
cat("====================================================================\n")
cat("                    GSEA SUMMARY\n")
cat("====================================================================\n")

summarize_gsea <- function(result, name) {
  if (is.null(result) || nrow(result@result) == 0) {
    cat(paste0(name, ": No results\n"))
    return()
  }

  total <- nrow(result@result)
  sig_005 <- sum(result@result$p.adjust < 0.05)
  sig_025 <- sum(result@result$p.adjust < 0.25)

  cat(paste0(name, ":\n"))
  cat(paste0("  Total terms: ", total, "\n"))
  cat(paste0("  Significant (p.adj < 0.05): ", sig_005, "\n"))
  cat(paste0("  Significant (p.adj < 0.25): ", sig_025, "\n"))
}

summarize_gsea(gsea_bp, "Biological Process")
summarize_gsea(gsea_mf, "Molecular Function")
summarize_gsea(gsea_cc, "Cellular Component")

cat("\n====================================================================\n")
cat("GSEA complete! Next step: Run 04_Pathway_Analysis.R\n")
cat("====================================================================\n")
