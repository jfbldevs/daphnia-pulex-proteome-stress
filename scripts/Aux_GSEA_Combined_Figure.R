# =============================================================================
# Auxiliary Script: GSEA Combined Figure
# Creates a combined barplot of all GSEA results across GO categories
#
# Authors: J.A. Norambuena, L. Herrera-Belén, P. Poblete-Grant, J.F. Beltrán
# Year:    2026
# =============================================================================

# Load configuration and directory paths
source("00_Setup.R")

# Load GSEA results
gsea_bp <- read_csv(file.path(RESULTS_DIR, "TableS2_GSEA_BP.csv"), show_col_types = FALSE)
gsea_mf <- read_csv(file.path(RESULTS_DIR, "TableS3_GSEA_MF.csv"), show_col_types = FALSE)
gsea_cc <- read_csv(file.path(RESULTS_DIR, "TableS4_GSEA_CC.csv"), show_col_types = FALSE)

# Add category labels
gsea_bp$Category <- "Biological Process"
gsea_mf$Category <- "Molecular Function"
gsea_cc$Category <- "Cellular Component"

# Combine all results
all_gsea <- bind_rows(gsea_bp, gsea_mf, gsea_cc) %>%
  filter(p.adjust < 0.25) %>%
  mutate(
    Description = str_trunc(Description, 40),
    Direction = ifelse(NES > 0, "Llanquihue\n(anthropized)", "Icalma\n(oligotrophic)"),
    Significance = case_when(
      p.adjust < 0.01 ~ "***",
      p.adjust < 0.05 ~ "**",
      p.adjust < 0.1 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  arrange(Category, NES)

cat("Total terms for visualization:", nrow(all_gsea), "\n")

# Colors
colors_direction <- c(
  "Llanquihue\n(anthropized)" = "#E63946",
  "Icalma\n(oligotrophic)" = "#457B9D"
)

# Combined figure - horizontal barplot
p_combined <- ggplot(all_gsea, aes(x = NES, y = reorder(Description, NES), fill = Direction)) +
  geom_col(alpha = 0.85, width = 0.7) +
  geom_vline(xintercept = 0, linetype = "solid", color = "gray30", linewidth = 0.5) +
  geom_text(
    aes(
      label = Significance,
      x = ifelse(NES > 0, NES + 0.1, NES - 0.1)
    ),
    hjust = ifelse(all_gsea$NES > 0, 0, 1),
    size = 4, fontface = "bold"
  ) +
  facet_grid(Category ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = colors_direction, name = "Enriched in") +
  labs(
    title = "Gene Set Enrichment Analysis (GSEA)",
    subtitle = "GO terms enriched in each lake environment",
    x = "Normalized Enrichment Score (NES)",
    y = NULL,
    caption = "*** p.adj < 0.01, ** p.adj < 0.05, * p.adj < 0.1"
  ) +
  theme_paper +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "gray95", color = NA),
    axis.text.y = element_text(size = 9),
    plot.caption = element_text(size = 8, hjust = 1, color = "gray50")
  )

ggsave(file.path(FIGURES_DIR, "Fig2_GSEA_Enrichment.pdf"), p_combined, width = 11, height = 8)
ggsave(file.path(FIGURES_DIR, "Fig2_GSEA_Enrichment.png"), p_combined,
  width = 11, height = 8, dpi = 600, bg = "white"
)

cat("Saved: Fig2_GSEA_Enrichment.pdf/png\n")

# =============================================================================
# Summary table for the paper
# =============================================================================

summary_table <- all_gsea %>%
  select(Category, Description, NES, pvalue, p.adjust, setSize, Direction) %>%
  arrange(p.adjust) %>%
  mutate(
    NES = round(NES, 2),
    pvalue = formatC(pvalue, format = "e", digits = 2),
    p.adjust = formatC(p.adjust, format = "e", digits = 2)
  )

write_csv(summary_table, file.path(RESULTS_DIR, "Table2_GSEA_Enrichment.csv"))
cat("Saved: Table2_GSEA_Enrichment.csv\n")

print(summary_table)
