# =============================================================================
# Script 04: Metabolic Pathway Analysis
# Functional classification of DEPs into metabolic pathways
# (Keyword-based classification from UniProt protein names)
#
# Authors: J.A. Norambuena, L. Herrera-Belén, P. Poblete-Grant, J.F. Beltrán
# Year:    2026
# =============================================================================

# Load configuration and directory paths
source("00_Setup.R")

cat("\n")
cat("====================================================================\n")
cat("         Metabolic Pathway Analysis\n")
cat("====================================================================\n\n")

# Load data
diff_results <- read_csv(file.path(RESULTS_DIR, "Table1_Differential_Expression.csv"), show_col_types = FALSE)
load(file.path(RESULTS_DIR, "RData_01_DataPrep.RData"))

cat("Total proteins:", nrow(diff_results), "\n")
cat("DEPs upregulated:", length(upregulated), "\n")
cat("DEPs downregulated:", length(downregulated), "\n\n")

# =============================================================================
# 1. OBTENER KEGG ANNOTATIONS DESDE UNIPROT
# =============================================================================

cat("=== RETRIEVING KEGG ANNOTATIONS (automatic attempt) ===\n\n")

# Function to retrieve KEGG from UniProt
get_uniprot_kegg <- function(uniprot_ids, batch_size = 100) {
  all_results <- data.frame()
  n_batches <- ceiling(length(uniprot_ids) / batch_size)

  for (i in 1:n_batches) {
    start_idx <- (i - 1) * batch_size + 1
    end_idx <- min(i * batch_size, length(uniprot_ids))
    batch_ids <- uniprot_ids[start_idx:end_idx]

    cat(sprintf("Procesando lote %d/%d...\n", i, n_batches))

    query_ids <- paste(batch_ids, collapse = "+OR+accession:")
    url <- paste0(
      "https://rest.uniprot.org/uniprotkb/search?",
      "query=accession:", query_ids,
      "&fields=accession,xref_kegg",
      "&format=tsv&size=", batch_size
    )

    tryCatch(
      {
        response <- readLines(url, warn = FALSE)

        if (length(response) > 1) {
          header <- strsplit(response[1], "\t")[[1]]
          data_lines <- response[-1]
          data_lines <- data_lines[nchar(data_lines) > 0]

          if (length(data_lines) > 0) {
            batch_df <- do.call(rbind, lapply(data_lines, function(line) {
              values <- strsplit(line, "\t")[[1]]
              length(values) <- length(header)
              return(values)
            }))

            batch_df <- as.data.frame(batch_df, stringsAsFactors = FALSE)
            colnames(batch_df) <- header
            all_results <- bind_rows(all_results, batch_df)
          }
        }

        Sys.sleep(0.3)
      },
      error = function(e) {
        cat(sprintf("  Error en lote %d: %s\n", i, e$message))
      }
    )
  }

  return(all_results)
}

# Check if KEGG annotations file already exists
kegg_file <- file.path(RESULTS_DIR, "Data_UniProt_KEGG_Annotations.csv")

if (file.exists(kegg_file)) {
  cat("Loading previously downloaded KEGG annotations...\n")
  kegg_annotations <- read_csv(kegg_file, show_col_types = FALSE)
} else {
  cat("Descargando anotaciones KEGG de UniProt...\n")
  kegg_annotations <- get_uniprot_kegg(background)

  if (nrow(kegg_annotations) > 0) {
    write_csv(kegg_annotations, kegg_file)
    cat("Annotations saved\n")
  }
}

cat("Proteins with KEGG data:", nrow(kegg_annotations), "\n\n")

# =============================================================================
# 2. CREAR TERM2GENE PARA KEGG
# =============================================================================

cat("=== PREPARANDO DATOS KEGG ===\n\n")

# Identificar columna KEGG
kegg_col <- colnames(kegg_annotations)[grep("kegg", colnames(kegg_annotations), ignore.case = TRUE)][1]
id_col <- colnames(kegg_annotations)[grep("entry|accession", colnames(kegg_annotations), ignore.case = TRUE)][1]

cat("Columna ID:", id_col, "\n")
cat("Columna KEGG:", kegg_col, "\n\n")

# Extraer pathways KEGG
term2gene_kegg <- data.frame(TERM = character(), GENE = character())

if (!is.na(kegg_col) && kegg_col %in% colnames(kegg_annotations)) {
  for (i in 1:nrow(kegg_annotations)) {
    protein_id <- kegg_annotations[[id_col]][i]
    kegg_text <- kegg_annotations[[kegg_col]][i]

    if (!is.na(kegg_text) && nchar(kegg_text) > 0) {
      # Extraer KEGG pathway IDs (formato: dpu:XXXXX o pathway:dpuXXXXX)
      # Also extract direct pathway IDs
      pathways <- unlist(str_extract_all(kegg_text, "dpu\\d+"))

      if (length(pathways) > 0) {
        term2gene_kegg <- rbind(term2gene_kegg, data.frame(
          TERM = pathways,
          GENE = rep(protein_id, length(pathways))
        ))
      }
    }
  }

  term2gene_kegg <- unique(term2gene_kegg)
}

cat("KEGG associations extracted:", nrow(term2gene_kegg), "\n")
cat("Unique pathways:", n_distinct(term2gene_kegg$TERM), "\n\n")

# =============================================================================
# 3. ALTERNATIVE KEGG ANALYSIS - Using manual annotation
# =============================================================================

cat("=== METABOLIC PATHWAY ANALYSIS ===\n\n")

# If insufficient direct KEGG annotations, use functional classification
# based on protein names

# Load protein names
protein_info <- diff_results %>%
  select(UniProt_ID, Protein_Name_Clean, Log2FC, P_Adjusted, Regulation) %>%
  filter(!is.na(Protein_Name_Clean))

# Define metabolic pathways manually
pathway_keywords <- list(
  "TCA Cycle" = c(
    "citrate", "isocitrate", "aconitase", "succinate", "fumarate",
    "malate", "oxoglutarate", "succinyl"
  ),
  "Glycolysis/Gluconeogenesis" = c(
    "glyceraldehyde", "phosphoglycerate", "enolase",
    "pyruvate", "hexokinase", "aldolase", "phosphofructokinase"
  ),
  "Oxidative Phosphorylation" = c(
    "ATP synthase", "cytochrome", "NADH", "ubiquinone",
    "electron transfer", "complex I", "complex II",
    "complex III", "complex IV"
  ),
  "Fatty Acid Metabolism" = c(
    "acyl-CoA", "fatty acid", "carnitine", "enoyl",
    "hydroxyacyl", "ketoacyl"
  ),
  "Amino Acid Metabolism" = c(
    "aminotransferase", "dehydrogenase", "synthetase",
    "glutamate", "glutamine", "aspartate"
  ),
  "Protein Folding/Chaperones" = c(
    "heat shock", "hsp", "chaperone", "disulfide",
    "isomerase", "folding"
  ),
  "Ribosome/Translation" = c("ribosom", "60S", "40S", "translation", "elongation factor"),
  "Cytoskeleton" = c(
    "actin", "myosin", "tubulin", "tropomyosin", "spectrin",
    "filamin", "cofilin"
  ),
  "Signal Transduction" = c(
    "kinase", "phosphatase", "GTPase", "G protein",
    "receptor", "calmodulin"
  ),
  "Antioxidant Defense" = c(
    "superoxide", "peroxidase", "catalase", "thioredoxin",
    "glutathione"
  )
)

# Classify proteins
classify_protein <- function(protein_name, keywords_list) {
  if (is.na(protein_name)) {
    return(NA)
  }
  protein_lower <- tolower(protein_name)

  for (pathway in names(keywords_list)) {
    keywords <- keywords_list[[pathway]]
    for (kw in keywords) {
      if (grepl(tolower(kw), protein_lower)) {
        return(pathway)
      }
    }
  }
  return("Other")
}

protein_info$Pathway <- sapply(protein_info$Protein_Name_Clean,
  classify_protein,
  keywords_list = pathway_keywords
)

# Summary by pathway
pathway_summary <- protein_info %>%
  filter(Pathway != "Other" & !is.na(Pathway)) %>%
  group_by(Pathway, Regulation) %>%
  summarise(
    Count = n(),
    Proteins = paste(str_trunc(Protein_Name_Clean, 30), collapse = "; "),
    Mean_Log2FC = mean(Log2FC, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Pathway, Regulation)

cat("Proteins classified by pathway:\n")
print(pathway_summary %>% select(Pathway, Regulation, Count, Mean_Log2FC))
cat("\n")

# =============================================================================
# 4. PATHWAY VISUALIZATION
# =============================================================================

cat("=== GENERATING VISUALIZATION ===\n\n")

# Prepare data for visualization
pathway_plot_data <- protein_info %>%
  filter(Pathway != "Other" & !is.na(Pathway)) %>%
  group_by(Pathway) %>%
  summarise(
    Total = n(),
    Upregulated = sum(Regulation == "Upregulated"),
    Downregulated = sum(Regulation == "Downregulated"),
    Not_Significant = sum(Regulation == "Not significant"),
    Mean_Log2FC = mean(Log2FC, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(Total >= 2) %>% # At least 2 proteins
  arrange(desc(Total))

# Reshape para stacked bar
pathway_long <- pathway_plot_data %>%
  select(Pathway, Upregulated, Downregulated) %>%
  pivot_longer(
    cols = c(Upregulated, Downregulated),
    names_to = "Regulation",
    values_to = "Count"
  ) %>%
  filter(Count > 0)

# Colores
colors_reg <- c("Upregulated" = "#E63946", "Downregulated" = "#457B9D")

# theme_paper is loaded from 00_Setup.R

# Plot
p_pathway <- ggplot(pathway_long, aes(x = reorder(Pathway, Count), y = Count, fill = Regulation)) +
  geom_col(position = "stack", width = 0.7, alpha = 0.85) +
  geom_text(aes(label = Count),
    position = position_stack(vjust = 0.5),
    size = 3.5, color = "white", fontface = "bold"
  ) +
  coord_flip() +
  scale_fill_manual(
    values = colors_reg,
    labels = c(
      "Downregulated" = "Suppressed in Llanquihue",
      "Upregulated" = "Induced in Llanquihue"
    )
  ) +
  labs(
    title = "Metabolic Pathway Analysis",
    subtitle = "Distribution of differentially expressed proteins across functional pathways",
    x = NULL,
    y = "Number of DEPs",
    fill = "Expression Pattern"
  ) +
  theme_paper +
  theme(
    axis.text.y = element_text(size = 10),
    legend.title = element_text(face = "bold")
  )

ggsave(file.path(FIGURES_DIR, "Fig3A_Metabolic_Pathways.pdf"), p_pathway, width = 10, height = 7)
ggsave(file.path(FIGURES_DIR, "Fig3A_Metabolic_Pathways.png"), p_pathway, width = 10, height = 7, dpi = 600, bg = "white")

cat("Saved: Fig3A_Metabolic_Pathways.pdf/png\n")

# =============================================================================
# 5. HEATMAP DE PATHWAYS
# =============================================================================

# Crear heatmap de Log2FC promedio por pathway
pathway_heatmap_data <- protein_info %>%
  filter(Pathway != "Other" & !is.na(Pathway) & Regulation != "Not significant") %>%
  group_by(Pathway) %>%
  summarise(
    Mean_Log2FC = mean(Log2FC, na.rm = TRUE),
    N_proteins = n(),
    .groups = "drop"
  ) %>%
  filter(N_proteins >= 2)

p_heatmap <- ggplot(
  pathway_heatmap_data,
  aes(x = 1, y = reorder(Pathway, Mean_Log2FC), fill = Mean_Log2FC)
) +
  geom_tile(width = 0.8, height = 0.8) +
  geom_text(aes(label = paste0("n=", N_proteins, "\n", round(Mean_Log2FC, 2))),
    size = 3.5, color = "white", fontface = "bold"
  ) +
  scale_fill_gradient2(
    low = "#457B9D", mid = "gray90", high = "#E63946",
    midpoint = 0, name = "Mean\nLog2FC"
  ) +
  labs(
    title = "Pathway Expression Changes",
    subtitle = "Mean Log2FC of DEPs in each pathway\n(Positive = higher in Llanquihue, Negative = higher in Icalma)",
    x = NULL,
    y = NULL
  ) +
  theme_paper +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank()
  )

ggsave(file.path(FIGURES_DIR, "Fig3B_Pathway_Heatmap.pdf"), p_heatmap, width = 6, height = 7)
ggsave(file.path(FIGURES_DIR, "Fig3B_Pathway_Heatmap.png"), p_heatmap, width = 6, height = 7, dpi = 600, bg = "white")

cat("Saved: Fig3B_Pathway_Heatmap.pdf/png\n")

# =============================================================================
# 6. GUARDAR RESULTADOS
# =============================================================================

cat("\n=== GUARDANDO RESULTADOS ===\n\n")

# Save protein classification
protein_pathway_classification <- protein_info %>%
  filter(Pathway != "Other") %>%
  select(UniProt_ID, Protein_Name_Clean, Pathway, Log2FC, P_Adjusted, Regulation) %>%
  arrange(Pathway, Log2FC)

write_csv(protein_pathway_classification, file.path(RESULTS_DIR, "Metabolic_protein_pathway_classification.csv"))
cat("Saved: Metabolic_protein_pathway_classification.csv\n")

# Save summary
write_csv(pathway_summary, file.path(RESULTS_DIR, "Metabolic_pathway_summary.csv"))
cat("Saved: Metabolic_pathway_summary.csv\n")

# =============================================================================
# 7. RESUMEN
# =============================================================================

cat("\n")
cat("====================================================================\n")
cat("            RESUMEN METABOLIC PATHWAY ANALYSIS\n")
cat("====================================================================\n\n")

cat("Pathways con DEPs:\n")
print(pathway_plot_data %>% select(Pathway, Total, Downregulated, Upregulated, Mean_Log2FC))

cat("\n")
cat("Interpretation:\n")
cat("- Mean_Log2FC negativo = suprimido en Llanquihue (antropizado)\n")
cat("- Mean_Log2FC positivo = inducido en Llanquihue\n")
cat("\n")

cat("====================================================================\n")
cat("Metabolic Pathway analysis completed!\n")
cat("====================================================================\n")
