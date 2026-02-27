# =============================================================================
# Script 02: Gene Ontology Over-Representation Analysis (ORA)
# Paper: Integrative bioinformatics analysis reveals disrupted metabolic
#        pathways and key hub proteins in the Daphnia pulex proteome
#        under anthropogenic stress
#
# This script downloads GO annotations from UniProt API, builds TERM2GENE
# mappings, and performs ORA for upregulated, downregulated, and all DEPs
# across BP, MF, and CC ontologies.
#
# Authors: J.A. Norambuena, L. Herrera-Belén, P. Poblete-Grant, J.F. Beltrán
# Year:    2026
# =============================================================================

# Load configuration
source("00_setup.R")

# Load data from previous script
load(file.path(RESULTS_DIR, "RData_01_DataPrep.RData"))

cat("=== GENE ONTOLOGY ENRICHMENT ANALYSIS ===\n\n")

# =============================================================================
# 1. RETRIEVE GO ANNOTATIONS FROM UNIPROT (robust method)
# =============================================================================

cat("=== OBTENIENDO ANOTACIONES GO DESDE UNIPROT ===\n")

# Improved function to retrieve GO terms de UniProt
get_uniprot_go_robust <- function(uniprot_ids, batch_size = 50) {
  all_results <- data.frame()
  n_batches <- ceiling(length(uniprot_ids) / batch_size)

  for (i in 1:n_batches) {
    start_idx <- (i - 1) * batch_size + 1
    end_idx <- min(i * batch_size, length(uniprot_ids))
    batch_ids <- uniprot_ids[start_idx:end_idx]

    cat(sprintf("Procesando lote %d/%d (%d IDs)...\n", i, n_batches, length(batch_ids)))

    # Construir query para UniProt API
    query_ids <- paste(batch_ids, collapse = "+OR+accession:")
    url <- paste0(
      "https://rest.uniprot.org/uniprotkb/search?",
      "query=accession:", query_ids,
      "&fields=accession,go_p,go_c,go_f,protein_name",
      "&format=tsv&size=", batch_size
    )

    tryCatch(
      {
        # Usar readLines como alternativa robusta
        response <- tryCatch(
          {
            readLines(url, warn = FALSE)
          },
          error = function(e) {
            # Fallback con httr
            resp <- GET(url)
            if (status_code(resp) == 200) {
              strsplit(rawToChar(resp$content), "\n")[[1]]
            } else {
              character(0)
            }
          }
        )

        if (length(response) > 1) {
          # Parsear TSV manualmente
          header <- strsplit(response[1], "\t")[[1]]
          data_lines <- response[-1]
          data_lines <- data_lines[nchar(data_lines) > 0]

          if (length(data_lines) > 0) {
            batch_df <- do.call(rbind, lapply(data_lines, function(line) {
              values <- strsplit(line, "\t")[[1]]
              # Ensure same number of columns
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

# Check if annotations file already exists
go_annotations_file <- file.path(RESULTS_DIR, "Data_UniProt_GO_Annotations.csv")

if (file.exists(go_annotations_file)) {
  cat("Loading previously downloaded GO annotations...\n")
  go_annotations <- read_csv(go_annotations_file, show_col_types = FALSE)
} else {
  cat("\nRetrieving GO annotations for", length(background), "proteins...\n")
  go_annotations <- get_uniprot_go_robust(background)

  if (nrow(go_annotations) > 0) {
    write_csv(go_annotations, go_annotations_file)
    cat("Annotations saved to:", go_annotations_file, "\n")
  }
}

cat("Proteins with GO annotations:", nrow(go_annotations), "\n\n")

# =============================================================================
# 2. PREPARAR TERM2GENE PARA ENRICHMENT
# =============================================================================

cat("=== PREPARANDO DATOS PARA ENRICHMENT ===\n")

# Identificar columnas de GO
go_cols <- colnames(go_annotations)
bp_col <- go_cols[grep("biological|go_p", go_cols, ignore.case = TRUE)][1]
mf_col <- go_cols[grep("molecular|go_f", go_cols, ignore.case = TRUE)][1]
cc_col <- go_cols[grep("cellular|go_c", go_cols, ignore.case = TRUE)][1]
id_col <- go_cols[grep("entry|accession", go_cols, ignore.case = TRUE)][1]

cat("Columns identified:\n")
cat("  ID:", id_col, "\n")
cat("  BP:", bp_col, "\n")
cat("  MF:", mf_col, "\n")
cat("  CC:", cc_col, "\n\n")

# Function to extract GO terms and create TERM2GENE
extract_term2gene <- function(annotations, id_column, go_column) {
  if (is.na(go_column) || !go_column %in% colnames(annotations)) {
    return(data.frame(TERM = character(), GENE = character()))
  }

  term2gene <- data.frame(TERM = character(), GENE = character())

  for (i in 1:nrow(annotations)) {
    protein_id <- annotations[[id_column]][i]
    go_text <- annotations[[go_column]][i]

    if (!is.na(go_text) && nchar(go_text) > 0) {
      # Extraer GO IDs (formato GO:0000000)
      go_ids <- unlist(regmatches(go_text, gregexpr("GO:\\d+", go_text)))

      if (length(go_ids) > 0) {
        term2gene <- rbind(term2gene, data.frame(
          TERM = go_ids,
          GENE = rep(protein_id, length(go_ids))
        ))
      }
    }
  }

  return(unique(term2gene))
}

# Create TERM2GENE for each category
if (nrow(go_annotations) > 0 && !is.na(id_col)) {
  cat("Extrayendo GO terms...\n")

  term2gene_bp <- extract_term2gene(go_annotations, id_col, bp_col)
  term2gene_mf <- extract_term2gene(go_annotations, id_col, mf_col)
  term2gene_cc <- extract_term2gene(go_annotations, id_col, cc_col)

  cat("  Biological Process:", nrow(term2gene_bp), "asociaciones\n")
  cat("  Molecular Function:", nrow(term2gene_mf), "asociaciones\n")
  cat("  Cellular Component:", nrow(term2gene_cc), "asociaciones\n\n")
}

# =============================================================================
# 3. REALIZAR GO ENRICHMENT
# =============================================================================

cat("=== REALIZANDO GO ENRICHMENT ===\n")

# Function to perform enrichment and add descriptions
run_go_enrichment <- function(genes, term2gene, background_genes, ontology_name) {
  if (nrow(term2gene) == 0) {
    cat(paste0("  ", ontology_name, ": Sin datos disponibles\n"))
    return(NULL)
  }

  # Filtrar term2gene para incluir solo genes del background
  term2gene_filtered <- term2gene %>%
    filter(GENE %in% background_genes)

  if (nrow(term2gene_filtered) < 10) {
    cat(paste0("  ", ontology_name, ": Insuficientes asociaciones\n"))
    return(NULL)
  }

  result <- tryCatch(
    {
      enricher(
        gene = genes,
        universe = background_genes,
        TERM2GENE = term2gene_filtered,
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2,
        minGSSize = 3,
        maxGSSize = 500
      )
    },
    error = function(e) {
      cat(paste0("  ", ontology_name, " error: ", e$message, "\n"))
      return(NULL)
    }
  )

  if (!is.null(result) && nrow(result@result) > 0) {
    # Add descriptions from GO.db
    result@result$Description <- sapply(result@result$ID, function(go_id) {
      tryCatch(
        {
          term <- GO.db::GOTERM[[go_id]]
          if (!is.null(term)) Term(term) else go_id
        },
        error = function(e) go_id
      )
    })

    sig_count <- sum(result@result$p.adjust < 0.05)
    cat(paste0("  ", ontology_name, ": ", sig_count, " significant terms (p.adj < 0.05)\n"))
  } else {
    cat(paste0("  ", ontology_name, ": No enriched terms\n"))
  }

  return(result)
}

# --- Enrichment para UPREGULATED ---
cat("\n--- Upregulated Proteins ---\n")

if (exists("term2gene_bp") && nrow(term2gene_bp) > 0) {
  ora_up_bp <- run_go_enrichment(upregulated, term2gene_bp, background, "BP")
  ora_up_mf <- run_go_enrichment(upregulated, term2gene_mf, background, "MF")
  ora_up_cc <- run_go_enrichment(upregulated, term2gene_cc, background, "CC")
}

# --- Enrichment para DOWNREGULATED ---
cat("\n--- Downregulated Proteins ---\n")

if (exists("term2gene_bp") && nrow(term2gene_bp) > 0) {
  ora_down_bp <- run_go_enrichment(downregulated, term2gene_bp, background, "BP")
  ora_down_mf <- run_go_enrichment(downregulated, term2gene_mf, background, "MF")
  ora_down_cc <- run_go_enrichment(downregulated, term2gene_cc, background, "CC")
}

# --- Enrichment para TODAS las DEPs ---
cat("\n--- Todas las DEPs ---\n")

if (exists("term2gene_bp") && nrow(term2gene_bp) > 0) {
  ora_all_bp <- run_go_enrichment(DEPs, term2gene_bp, background, "BP")
  ora_all_mf <- run_go_enrichment(DEPs, term2gene_mf, background, "MF")
  ora_all_cc <- run_go_enrichment(DEPs, term2gene_cc, background, "CC")
}

# =============================================================================
# 4. VISUALIZACIONES
# =============================================================================

cat("\n=== GENERANDO VISUALIZACIONES ===\n")

# Function to create barplot
create_go_barplot <- function(enrichment_result, title, fill_color, top_n = 15) {
  if (is.null(enrichment_result) || nrow(enrichment_result@result) == 0) {
    return(NULL)
  }

  df <- enrichment_result@result %>%
    filter(p.adjust < 0.05) %>%
    head(top_n) %>%
    mutate(
      Description = str_trunc(Description, 50),
      Description = fct_reorder(Description, -log10(p.adjust))
    )

  if (nrow(df) == 0) {
    return(NULL)
  }

  ggplot(df, aes(x = -log10(p.adjust), y = Description)) +
    geom_col(fill = fill_color, alpha = 0.8) +
    geom_text(aes(label = Count), hjust = -0.2, size = 3) +
    labs(title = title, x = expression(-Log[10] ~ Adjusted ~ P - value), y = NULL) +
    theme_paper +
    theme(axis.text.y = element_text(size = 9)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.15)))
}

# Generate plots for BP
if (exists("ora_up_bp") && !is.null(ora_up_bp)) {
  p <- create_go_barplot(ora_up_bp, "GO Biological Process - Upregulated", "#E63946")
  if (!is.null(p)) {
    ggsave(file.path(FIGURES_DIR, "04_GO_BP_upregulated.pdf"), p, width = 10, height = 8)
    ggsave(file.path(FIGURES_DIR, "04_GO_BP_upregulated.png"), p, width = 10, height = 8, dpi = 300)
    cat("Guardado: 04_GO_BP_upregulated.pdf\n")
  }
}

if (exists("ora_down_bp") && !is.null(ora_down_bp)) {
  p <- create_go_barplot(ora_down_bp, "GO Biological Process - Downregulated", "#457B9D")
  if (!is.null(p)) {
    ggsave(file.path(FIGURES_DIR, "05_GO_BP_downregulated.pdf"), p, width = 10, height = 8)
    ggsave(file.path(FIGURES_DIR, "05_GO_BP_downregulated.png"), p, width = 10, height = 8, dpi = 300)
    cat("Guardado: 05_GO_BP_downregulated.pdf\n")
  }
}

if (exists("ora_all_bp") && !is.null(ora_all_bp)) {
  p <- create_go_barplot(ora_all_bp, "GO Biological Process - All DEPs", "#6B4C9A")
  if (!is.null(p)) {
    ggsave(file.path(FIGURES_DIR, "06_GO_BP_all_DEPs.pdf"), p, width = 10, height = 8)
    ggsave(file.path(FIGURES_DIR, "06_GO_BP_all_DEPs.png"), p, width = 10, height = 8, dpi = 300)
    cat("Guardado: 06_GO_BP_all_DEPs.pdf\n")
  }
}

# Combined dotplot if results exist
if (exists("ora_all_bp") && !is.null(ora_all_bp) && nrow(ora_all_bp@result) > 0) {
  sig_results <- ora_all_bp@result %>% filter(p.adjust < 0.05)
  if (nrow(sig_results) > 0) {
    ora_plot <- ora_all_bp
    ora_plot@result <- head(sig_results, 20)

    p <- dotplot(ora_plot, showCategory = 20) +
      labs(title = "GO Biological Process Enrichment - All DEPs") +
      theme_paper

    ggsave(file.path(FIGURES_DIR, "07_GO_BP_dotplot.pdf"), p, width = 10, height = 10)
    cat("Guardado: 07_GO_BP_dotplot.pdf\n")
  }
}

# =============================================================================
# 5. GUARDAR RESULTADOS
# =============================================================================

cat("\n=== GUARDANDO RESULTADOS ===\n")

# Function to save results
save_if_exists <- function(result, filename) {
  if (!is.null(result) && nrow(result@result) > 0) {
    write_csv(result@result, file.path(RESULTS_DIR, filename))
    cat("Guardado:", filename, "\n")
  }
}

if (exists("ora_up_bp")) save_if_exists(ora_up_bp, "GO_BP_upregulated.csv")
if (exists("ora_down_bp")) save_if_exists(ora_down_bp, "GO_BP_downregulated.csv")
if (exists("ora_all_bp")) save_if_exists(ora_all_bp, "GO_BP_all_DEPs.csv")
if (exists("ora_up_mf")) save_if_exists(ora_up_mf, "GO_MF_upregulated.csv")
if (exists("ora_down_mf")) save_if_exists(ora_down_mf, "GO_MF_downregulated.csv")
if (exists("ora_up_cc")) save_if_exists(ora_up_cc, "GO_CC_upregulated.csv")
if (exists("ora_down_cc")) save_if_exists(ora_down_cc, "GO_CC_downregulated.csv")

# Save TERM2GENE for reference
if (exists("term2gene_bp")) write_csv(term2gene_bp, file.path(RESULTS_DIR, "Data_Term2Gene_BP.csv"))
if (exists("term2gene_mf")) write_csv(term2gene_mf, file.path(RESULTS_DIR, "Data_Term2Gene_MF.csv"))
if (exists("term2gene_cc")) write_csv(term2gene_cc, file.path(RESULTS_DIR, "Data_Term2Gene_CC.csv"))

# Save R objects
go_objects <- c(
  "ora_up_bp", "ora_down_bp", "ora_all_bp",
  "ora_up_mf", "ora_down_mf", "ora_up_cc", "ora_down_cc",
  "go_annotations", "term2gene_bp", "term2gene_mf", "term2gene_cc"
)
go_objects <- go_objects[sapply(go_objects, exists)]

save(list = go_objects, file = file.path(RESULTS_DIR, "RData_02_GO.RData"))

cat("\nScript 02 completed!\n")
cat("Next step: Run 03_GSEA_Analysis.R\n")
