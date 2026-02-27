# =============================================================================
# Script 04: Protein-Protein Interaction Network Analysis
# Paper: Systems biology analysis of Daphnia pulex under anthropogenic stress
# =============================================================================

# Load configuration
source("00_setup.R")

# Load data from previous scripts
load(file.path(RESULTS_DIR, "RData_01_DataPrep.RData"))

cat("=== PROTEIN-PROTEIN INTERACTION NETWORK ANALYSIS ===\n\n")

# =============================================================================
# 1. CONFIGURAR STRING DATABASE
# =============================================================================

cat("--- Configurando STRING Database ---\n")

# STRING tiene Daphnia pulex (taxon ID: 6669)
# Sin embargo, la cobertura puede ser limitada
# Alternativa: usar la API REST directamente

daphnia_taxon <- 6669
drosophila_taxon <- 7227 # Drosophila como backup

use_string_package <- FALSE
string_db <- NULL

# Intentar con Daphnia primero
cat("Intentando conectar con STRING para Daphnia pulex (taxon:", daphnia_taxon, ")...\n")

tryCatch(
  {
    string_db <- STRINGdb$new(
      version = "12.0",
      species = daphnia_taxon,
      score_threshold = 400, # Medium confidence
      network_type = "full",
      input_directory = RESULTS_DIR
    )
    cat("Successful connection to STRING for Daphnia pulex!\n\n")
    use_string_package <- TRUE
  },
  error = function(e) {
    cat("STRINGdb para Daphnia no disponible:", e$message, "\n")
    cat("Intentando con Drosophila como proxy...\n")

    tryCatch(
      {
        string_db <<- STRINGdb$new(
          version = "12.0",
          species = drosophila_taxon,
          score_threshold = 400,
          network_type = "full",
          input_directory = RESULTS_DIR
        )
        cat("Successful connection to STRING for Drosophila!\n\n")
        use_string_package <<- TRUE
      },
      error = function(e2) {
        cat("STRINGdb no disponible, usando API REST directamente...\n\n")
      }
    )
  }
)

# =============================================================================
# 2. ROBUST FUNCTION FOR STRING API
# =============================================================================

# Generic function to make requests to STRING API
string_api_request <- function(url, method = "GET") {
  response <- tryCatch(
    {
      if (method == "GET") {
        lines <- readLines(url, warn = FALSE)
      } else {
        resp <- POST(url)
        if (status_code(resp) == 200) {
          strsplit(rawToChar(resp$content), "\n")[[1]]
        } else {
          character(0)
        }
      }
      lines
    },
    error = function(e) {
      tryCatch(
        {
          resp <- GET(url)
          if (status_code(resp) == 200) {
            strsplit(rawToChar(resp$content), "\n")[[1]]
          } else {
            character(0)
          }
        },
        error = function(e2) {
          character(0)
        }
      )
    }
  )

  response <- response[nchar(response) > 0]
  return(response)
}

# Function to retrieve interactions from STRING API
get_string_interactions <- function(protein_ids, species = 6669, score_threshold = 400) {
  cat("Obteniendo interacciones desde STRING API...\n")

  # STRING API para network
  base_url <- "https://string-db.org/api/tsv/network"

  # Split into batches if there are many proteins
  batch_size <- 100
  all_interactions <- data.frame()

  n_batches <- ceiling(length(protein_ids) / batch_size)

  for (i in 1:n_batches) {
    start_idx <- (i - 1) * batch_size + 1
    end_idx <- min(i * batch_size, length(protein_ids))
    batch_ids <- protein_ids[start_idx:end_idx]

    cat(sprintf("  Procesando batch %d/%d...\n", i, n_batches))

    # Construir URL
    identifiers <- paste(batch_ids, collapse = "%0d")
    url <- paste0(
      base_url, "?identifiers=", identifiers,
      "&species=", species,
      "&required_score=", score_threshold
    )

    tryCatch(
      {
        lines <- string_api_request(url)

        if (length(lines) > 1) {
          # Parsear TSV
          header <- strsplit(lines[1], "\t")[[1]]
          data_lines <- lines[-1]

          if (length(data_lines) > 0) {
            batch_df <- do.call(rbind, lapply(data_lines, function(line) {
              strsplit(line, "\t")[[1]]
            }))

            if (!is.null(batch_df) && nrow(batch_df) > 0) {
              batch_df <- as.data.frame(batch_df, stringsAsFactors = FALSE)
              colnames(batch_df) <- header[1:ncol(batch_df)]
              all_interactions <- bind_rows(all_interactions, batch_df)
            }
          }
        }

        Sys.sleep(0.5) # Rate limiting
      },
      error = function(e) {
        cat(sprintf("  Error en batch %d: %s\n", i, e$message))
      }
    )
  }

  return(all_interactions)
}

# Function to map proteins to STRING IDs
get_string_mapping <- function(protein_ids, species = 6669) {
  cat("Mapping proteins to STRING IDs...\n")

  base_url <- "https://string-db.org/api/tsv/get_string_ids"

  mapping <- data.frame()
  batch_size <- 100
  n_batches <- ceiling(length(protein_ids) / batch_size)

  for (i in 1:n_batches) {
    start_idx <- (i - 1) * batch_size + 1
    end_idx <- min(i * batch_size, length(protein_ids))
    batch_ids <- protein_ids[start_idx:end_idx]

    cat(sprintf("  Mapeando batch %d/%d...\n", i, n_batches))

    identifiers <- paste(batch_ids, collapse = "%0d")
    url <- paste0(
      base_url, "?identifiers=", identifiers,
      "&species=", species, "&limit=1"
    )

    tryCatch(
      {
        lines <- string_api_request(url)

        if (length(lines) > 1) {
          header <- strsplit(lines[1], "\t")[[1]]
          data_lines <- lines[-1]

          if (length(data_lines) > 0) {
            batch_df <- do.call(rbind, lapply(data_lines, function(line) {
              strsplit(line, "\t")[[1]]
            }))

            if (!is.null(batch_df) && nrow(batch_df) > 0) {
              batch_df <- as.data.frame(batch_df, stringsAsFactors = FALSE)
              colnames(batch_df) <- header[1:ncol(batch_df)]
              mapping <- bind_rows(mapping, batch_df)
            }
          }
        }

        Sys.sleep(0.3)
      },
      error = function(e) {
        cat(sprintf("  Error en batch %d: %s\n", i, e$message))
      }
    )
  }

  return(mapping)
}

# =============================================================================
# 3. MAP PROTEINS TO STRING IDs
# =============================================================================

cat("=== MAPPING PROTEINS TO STRING ===\n")

# Prepare dataframe with protein information
proteins_for_mapping <- differential_analysis %>%
  dplyr::select(UniProt_ID, Protein_Name_Clean, Log2FC, P_Value, Regulation) %>%
  filter(!is.na(UniProt_ID))

# Mapping file saved
string_mapping_file <- file.path(RESULTS_DIR, "Data_STRING_Mapping.csv")

if (file.exists(string_mapping_file)) {
  cat("Loading previously saved STRING mapping...\n")
  string_mapping <- read_csv(string_mapping_file, show_col_types = FALSE)
} else {
  if (use_string_package && !is.null(string_db)) {
    # Usar paquete STRINGdb
    cat("Usando STRINGdb package para mapeo...\n")
    mapped_proteins <- tryCatch(
      {
        string_db$map(proteins_for_mapping, "UniProt_ID", removeUnmappedRows = TRUE)
      },
      error = function(e) {
        cat("Error con STRINGdb$map:", e$message, "\n")
        NULL
      }
    )

    if (!is.null(mapped_proteins) && nrow(mapped_proteins) > 0) {
      string_mapping <- mapped_proteins
    } else {
      string_mapping <- get_string_mapping(proteins_for_mapping$UniProt_ID, daphnia_taxon)
    }
  } else {
    # Usar API REST
    string_mapping <- get_string_mapping(proteins_for_mapping$UniProt_ID, daphnia_taxon)
  }

  if (nrow(string_mapping) > 0) {
    write_csv(string_mapping, string_mapping_file)
  }
}

# Combine with protein information
cat("Columns in string_mapping:", paste(colnames(string_mapping), collapse = ", "), "\n")
cat("Filas en string_mapping:", nrow(string_mapping), "\n")

mapped_proteins <- data.frame()

if (nrow(string_mapping) > 0) {
  if ("STRING_id" %in% colnames(string_mapping)) {
    mapped_proteins <- string_mapping
  } else if ("stringId" %in% colnames(string_mapping)) {
    # API REST usa "stringId" y "queryIndex"
    string_mapping <- string_mapping %>%
      dplyr::rename(STRING_id = stringId)

    # FIX: Extract UniProt_ID directly from STRING_id
    # El formato STRING_id es: taxon.UNIPROT_ID (ej: 6669.E9FQP0)
    # This avoids the bug where multiple STRING_ids with the same queryIndex
    # received the same UniProt_ID incorrectly

    # Extract UniProt_ID from STRING_id (part after the dot)
    string_mapping$UniProt_ID_from_STRING <- sub("^\\d+\\.", "", string_mapping$STRING_id)

    # Verify which match our proteins
    string_mapping_valid <- string_mapping %>%
      filter(UniProt_ID_from_STRING %in% proteins_for_mapping$UniProt_ID) %>%
      dplyr::rename(UniProt_ID = UniProt_ID_from_STRING) %>%
      distinct(STRING_id, UniProt_ID, .keep_all = TRUE)

    cat("  Proteins with valid STRING mapping:", nrow(string_mapping_valid), "\n")

    if (nrow(string_mapping_valid) > 0) {
      mapped_proteins <- proteins_for_mapping %>%
        inner_join(string_mapping_valid %>% dplyr::select(UniProt_ID, STRING_id),
          by = "UniProt_ID"
        )
    }
  }
}

cat("Proteins successfully mapped:", nrow(mapped_proteins), "\n")

if (nrow(mapped_proteins) > 0) {
  # Separar upregulated y downregulated
  mapped_up <- mapped_proteins %>% filter(Regulation == "Upregulated")
  mapped_down <- mapped_proteins %>% filter(Regulation == "Downregulated")

  cat("  - Upregulated mapeadas:", nrow(mapped_up), "\n")
  cat("  - Downregulated mapeadas:", nrow(mapped_down), "\n\n")
}

# =============================================================================
# 4. OBTENER RED DE INTERACCIONES
# =============================================================================

cat("=== OBTENIENDO RED DE INTERACCIONES ===\n")

interactions_file <- file.path(RESULTS_DIR, "TableS7_PPI_Interactions.csv")

if (nrow(mapped_proteins) > 0) {
  if (file.exists(interactions_file)) {
    cat("Cargando interacciones previamente descargadas...\n")
    interactions <- read_csv(interactions_file, show_col_types = FALSE)
  } else {
    # Obtener STRING IDs
    string_ids <- mapped_proteins$STRING_id

    if (use_string_package && !is.null(string_db)) {
      cat("Usando STRINGdb package para interacciones...\n")
      interactions <- tryCatch(
        {
          string_db$get_interactions(string_ids)
        },
        error = function(e) {
          cat("Error con STRINGdb:", e$message, "\n")
          get_string_interactions(mapped_proteins$UniProt_ID, daphnia_taxon)
        }
      )
    } else {
      interactions <- get_string_interactions(mapped_proteins$UniProt_ID, daphnia_taxon)
    }

    if (nrow(interactions) > 0) {
      write_csv(interactions, interactions_file)
    }
  }

  cat("Interacciones obtenidas:", nrow(interactions), "\n\n")
} else {
  cat("No mapped proteins to retrieve interactions.\n\n")
  interactions <- data.frame()
}

# =============================================================================
# 5. CONSTRUIR Y ANALIZAR RED CON IGRAPH
# =============================================================================

cat("=== CONSTRUYENDO RED DE INTERACCIONES ===\n")

if (nrow(interactions) > 0) {
  # Identificar columnas de from/to
  from_col <- intersect(c("from", "stringId_A", "preferredName_A"), colnames(interactions))[1]
  to_col <- intersect(c("to", "stringId_B", "preferredName_B"), colnames(interactions))[1]
  score_col <- intersect(c("combined_score", "score"), colnames(interactions))[1]

  if (!is.na(from_col) && !is.na(to_col)) {
    # Crear grafo
    ppi_graph <- graph_from_data_frame(
      interactions[, c(from_col, to_col)],
      directed = FALSE
    )

    # Add weights if score exists
    if (!is.na(score_col)) {
      scores <- as.numeric(interactions[[score_col]])
      # Normalizar si es necesario
      if (max(scores, na.rm = TRUE) > 1) {
        scores <- scores / 1000
      }
      E(ppi_graph)$weight <- scores
    }

    # Simplificar (remover loops y edges duplicados)
    ppi_graph <- igraph::simplify(ppi_graph)

    cat("Nodos en la red:", vcount(ppi_graph), "\n")
    cat("Edges en la red:", ecount(ppi_graph), "\n")

    # Add node attributes
    node_info <- mapped_proteins %>%
      dplyr::select(STRING_id, UniProt_ID, Protein_Name_Clean, Log2FC, Regulation) %>%
      distinct(STRING_id, .keep_all = TRUE)

    # Mapear atributos a nodos
    V(ppi_graph)$uniprot_id <- node_info$UniProt_ID[match(V(ppi_graph)$name, node_info$STRING_id)]
    V(ppi_graph)$protein_name <- node_info$Protein_Name_Clean[match(V(ppi_graph)$name, node_info$STRING_id)]
    V(ppi_graph)$log2fc <- node_info$Log2FC[match(V(ppi_graph)$name, node_info$STRING_id)]
    V(ppi_graph)$regulation <- node_info$Regulation[match(V(ppi_graph)$name, node_info$STRING_id)]

    # -----------------------------------------------------------------------------
    # 5.1 Calculate centrality metrics
    # -----------------------------------------------------------------------------
    cat("\n--- Calculating centrality metrics ---\n")

    # Degree centrality
    V(ppi_graph)$degree <- degree(ppi_graph)

    # Betweenness centrality
    V(ppi_graph)$betweenness <- betweenness(ppi_graph, normalized = TRUE)

    # Closeness centrality (only if graph is connected)
    V(ppi_graph)$closeness <- tryCatch(
      {
        closeness(ppi_graph, normalized = TRUE)
      },
      error = function(e) {
        rep(NA, vcount(ppi_graph))
      }
    )

    # Eigenvector centrality
    V(ppi_graph)$eigenvector <- tryCatch(
      {
        eigen_centrality(ppi_graph)$vector
      },
      error = function(e) {
        rep(NA, vcount(ppi_graph))
      }
    )

    # PageRank
    V(ppi_graph)$pagerank <- tryCatch(
      {
        page_rank(ppi_graph)$vector
      },
      error = function(e) {
        rep(NA, vcount(ppi_graph))
      }
    )

    # Create centrality metrics dataframe
    centrality_metrics <- data.frame(
      STRING_id = V(ppi_graph)$name,
      UniProt_ID = V(ppi_graph)$uniprot_id,
      Protein_Name = V(ppi_graph)$protein_name,
      Log2FC = V(ppi_graph)$log2fc,
      Regulation = V(ppi_graph)$regulation,
      Degree = V(ppi_graph)$degree,
      Betweenness = V(ppi_graph)$betweenness,
      Closeness = V(ppi_graph)$closeness,
      Eigenvector = V(ppi_graph)$eigenvector,
      PageRank = V(ppi_graph)$pagerank
    )

    # Ordenar por degree
    centrality_metrics <- centrality_metrics %>%
      arrange(desc(Degree))

    # Top 20 hub proteins
    cat("\n--- Top 20 Hub Proteins (por Degree) ---\n")
    print(head(centrality_metrics, 20))

    # Save metrics
    write_csv(centrality_metrics, file.path(RESULTS_DIR, "TableS6_PPI_Centrality.csv"))
    cat("\nCentrality metrics saved.\n")

    # -----------------------------------------------------------------------------
    # 5.2 Identificar hubs
    # -----------------------------------------------------------------------------
    cat("\n--- Identificando Hub Proteins ---\n")

    # Define hubs as proteins with degree >= 95th percentile
    # (Barabási & Oltvai, 2004)
    degree_threshold <- quantile(centrality_metrics$Degree, 0.95, na.rm = TRUE)

    hub_proteins <- centrality_metrics %>%
      filter(Degree >= degree_threshold)

    cat("Threshold de degree para hubs:", round(degree_threshold, 2), "\n")
    cat("Number of hub proteins:", nrow(hub_proteins), "\n")

    if (nrow(hub_proteins) > 0) {
      cat("\nHub proteins identificadas:\n")
      for (i in 1:min(10, nrow(hub_proteins))) {
        cat(sprintf(
          "  %d. %s (Degree: %d, %s)\n",
          i,
          ifelse(is.na(hub_proteins$Protein_Name[i]),
            hub_proteins$UniProt_ID[i],
            hub_proteins$Protein_Name[i]
          ),
          hub_proteins$Degree[i],
          ifelse(is.na(hub_proteins$Regulation[i]),
            "N/A",
            hub_proteins$Regulation[i]
          )
        ))
      }

      write_csv(hub_proteins, file.path(RESULTS_DIR, "Table4_Hub_Proteins.csv"))
    }

    # -----------------------------------------------------------------------------
    # 5.3 Community/module detection
    # -----------------------------------------------------------------------------
    cat("\n--- Detecting functional modules ---\n")

    # Usar algoritmo de Louvain
    communities <- tryCatch(
      {
        cluster_louvain(ppi_graph)
      },
      error = function(e) {
        cat("Error con Louvain, usando fast_greedy...\n")
        cluster_fast_greedy(ppi_graph)
      }
    )

    cat("Number of modules detected:", length(communities), "\n")
    cat("Modularidad:", round(modularity(communities), 3), "\n")

    # Add community membership
    V(ppi_graph)$community <- membership(communities)
    centrality_metrics$Community <- membership(communities)[match(
      centrality_metrics$STRING_id,
      V(ppi_graph)$name
    )]

    # Size of each community
    community_sizes <- table(membership(communities))
    cat("\nModule sizes:\n")
    print(sort(community_sizes, decreasing = TRUE))

    # Save metrics actualizadas
    write_csv(centrality_metrics, file.path(RESULTS_DIR, "TableS6_PPI_Centrality.csv")) # Update with module info
  } else {
    cat("No se pudieron identificar columnas from/to en las interacciones.\n")
    ppi_graph <- NULL
  }
} else {
  cat("No hay interacciones para construir la red.\n")
  ppi_graph <- NULL
}

# =============================================================================
# 6. VISUALIZACIONES DE LA RED
# =============================================================================

cat("\n=== GENERANDO VISUALIZACIONES DE RED ===\n")

if (exists("ppi_graph") && !is.null(ppi_graph) && vcount(ppi_graph) > 0) {
  # -----------------------------------------------------------------------------
  # 6.1 Red completa con ggraph
  # -----------------------------------------------------------------------------
  cat("Generating full network visualization...\n")

  # Convertir a tbl_graph para ggraph
  ppi_tbl <- as_tbl_graph(ppi_graph)

  # Plot de red
  set.seed(42)
  p_network <- ggraph(ppi_tbl, layout = "fr") +
    geom_edge_link(alpha = 0.2, color = "gray60") +
    geom_node_point(aes(size = degree, color = regulation), alpha = 0.7) +
    scale_color_manual(values = regulation_colors, na.value = "#CCCCCC") +
    scale_size_continuous(range = c(2, 10), name = "Degree") +
    labs(
      title = "Protein-Protein Interaction Network",
      subtitle = paste("Nodes:", vcount(ppi_graph), "| Edges:", ecount(ppi_graph)),
      color = "Regulation"
    ) +
    theme_void() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5),
      legend.position = "right"
    )

  ggsave(file.path(FIGURES_DIR, "Fig4A_PPI_Network_Metabolic.pdf"), p_network, width = 12, height = 10)
  ggsave(file.path(FIGURES_DIR, "Fig4A_PPI_Network_Metabolic.png"), p_network, width = 12, height = 10, dpi = 300)
  cat("Guardado: Fig4A_PPI_Network_Metabolic.pdf\n")

  # -----------------------------------------------------------------------------
  # 6.2 Red con labels para hubs
  # -----------------------------------------------------------------------------
  if (exists("hub_proteins") && nrow(hub_proteins) > 0 && exists("degree_threshold")) {
    # Add labels only for hubs
    V(ppi_graph)$label <- ifelse(
      V(ppi_graph)$degree >= degree_threshold,
      ifelse(is.na(V(ppi_graph)$protein_name),
        V(ppi_graph)$uniprot_id,
        V(ppi_graph)$protein_name
      ),
      NA
    )

    ppi_tbl <- as_tbl_graph(ppi_graph)

    set.seed(42)
    p_network_hubs <- ggraph(ppi_tbl, layout = "fr") +
      geom_edge_link(alpha = 0.15, color = "gray60") +
      geom_node_point(aes(size = degree, color = regulation), alpha = 0.7) +
      geom_node_text(aes(label = label), repel = TRUE, size = 3, max.overlaps = 20) +
      scale_color_manual(values = regulation_colors, na.value = "#CCCCCC") +
      scale_size_continuous(range = c(2, 12), name = "Degree") +
      labs(
        title = "PPI Network with Hub Proteins Labeled",
        subtitle = paste("Hubs (degree >=", round(degree_threshold), "):", nrow(hub_proteins)),
        color = "Regulation"
      ) +
      theme_void() +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5),
        legend.position = "right"
      )

    ggsave(file.path(FIGURES_DIR, "Fig5B_Hub_Proteins.pdf"), p_network_hubs, width = 14, height = 12)
    ggsave(file.path(FIGURES_DIR, "Fig5B_Hub_Proteins.png"), p_network_hubs, width = 14, height = 12, dpi = 300)
    cat("Guardado: Fig5B_Hub_Proteins.pdf\n")
  }

  # -----------------------------------------------------------------------------
  # 6.3 Network colored by modules
  # -----------------------------------------------------------------------------
  if (exists("communities")) {
    set.seed(42)
    p_network_modules <- ggraph(ppi_tbl, layout = "fr") +
      geom_edge_link(alpha = 0.15, color = "gray60") +
      geom_node_point(aes(size = degree, color = factor(community)), alpha = 0.7) +
      scale_color_viridis_d(option = "turbo", name = "Module") +
      scale_size_continuous(range = c(2, 10), name = "Degree") +
      labs(
        title = "PPI Network - Functional Modules",
        subtitle = paste("Modules detected:", length(communities))
      ) +
      theme_void() +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5),
        legend.position = "right"
      )

    ggsave(file.path(FIGURES_DIR, "FigS1_Functional_Categories.pdf"), p_network_modules, width = 12, height = 10)
    ggsave(file.path(FIGURES_DIR, "FigS1_Functional_Categories.png"), p_network_modules, width = 12, height = 10, dpi = 300)
    cat("Guardado: FigS1_Functional_Categories.pdf\n")
  }

  # -----------------------------------------------------------------------------
  # 6.4 Degree distribution
  # -----------------------------------------------------------------------------
  if (exists("centrality_metrics") && nrow(centrality_metrics) > 0) {
    p_degree <- ggplot(centrality_metrics, aes(x = Degree)) +
      geom_histogram(fill = "#2E86AB", color = "white", bins = 30) +
      geom_vline(xintercept = degree_threshold, linetype = "dashed", color = "#E63946", linewidth = 1) +
      annotate("text",
        x = degree_threshold + 5, y = Inf, vjust = 2,
        label = paste("Hub threshold:", round(degree_threshold)), color = "#E63946"
      ) +
      labs(
        title = "Degree Distribution of PPI Network",
        x = "Degree (number of connections)",
        y = "Frequency"
      ) +
      theme_paper

    ggsave(file.path(FIGURES_DIR, "Fig5A_Degree_Distribution.pdf"), p_degree, width = 8, height = 6)
    cat("Guardado: Fig5A_Degree_Distribution.pdf\n")

    # -----------------------------------------------------------------------------
    # 6.5 Scatter plot Degree vs Betweenness
    # -----------------------------------------------------------------------------
    p_centrality <- ggplot(centrality_metrics, aes(x = Degree, y = Betweenness)) +
      geom_point(aes(color = Regulation, size = Eigenvector), alpha = 0.6) +
      geom_text_repel(
        data = centrality_metrics %>% filter(Degree >= degree_threshold),
        aes(label = ifelse(is.na(Protein_Name), UniProt_ID, Protein_Name)),
        size = 3, max.overlaps = 15
      ) +
      scale_color_manual(values = regulation_colors, na.value = "#CCCCCC") +
      scale_size_continuous(range = c(2, 8), name = "Eigenvector\nCentrality") +
      labs(
        title = "Network Centrality Analysis",
        subtitle = "Hub proteins labeled",
        x = "Degree Centrality",
        y = "Betweenness Centrality"
      ) +
      theme_paper

    ggsave(file.path(FIGURES_DIR, "Fig5A_Centrality_Analysis.pdf"), p_centrality, width = 10, height = 8)
    cat("Guardado: Fig5A_Centrality_Analysis.pdf\n")
  }
} else {
  cat("No hay red para visualizar.\n")
}

# =============================================================================
# 7. FUNCTIONAL MODULE ANALYSIS
# =============================================================================

cat("\n=== FUNCTIONAL MODULE ANALYSIS ===\n")

if (exists("centrality_metrics") && "Community" %in% colnames(centrality_metrics)) {
  # Summary by module
  module_summary <- centrality_metrics %>%
    group_by(Community) %>%
    summarise(
      N_proteins = n(),
      N_upregulated = sum(Regulation == "Upregulated", na.rm = TRUE),
      N_downregulated = sum(Regulation == "Downregulated", na.rm = TRUE),
      Mean_Degree = mean(Degree, na.rm = TRUE),
      Mean_Log2FC = mean(Log2FC, na.rm = TRUE),
      Top_Proteins = paste(head(na.omit(Protein_Name)[order(-Degree[!is.na(Protein_Name)])], 3), collapse = ", ")
    ) %>%
    arrange(desc(N_proteins))

  cat("\nModule summary:\n")
  print(module_summary)

  write_csv(module_summary, file.path(RESULTS_DIR, "PPI_module_summary.csv"))
} else {
  cat("No module information available.\n")
}

# =============================================================================
# 8. GUARDAR RESULTADOS
# =============================================================================

cat("\n=== GUARDANDO RESULTADOS ===\n")

# Save graph
if (exists("ppi_graph") && !is.null(ppi_graph)) {
  saveRDS(ppi_graph, file.path(RESULTS_DIR, "PPI_graph.rds"))
  cat("Guardado: PPI_graph.rds\n")
}

# Save R objects
ppi_objects <- c(
  "ppi_graph", "centrality_metrics", "hub_proteins", "communities",
  "mapped_proteins", "interactions"
)
ppi_objects <- ppi_objects[sapply(ppi_objects, exists)]

if (length(ppi_objects) > 0) {
  save(list = ppi_objects, file = file.path(RESULTS_DIR, "04_PPI_network.RData"))
}

cat("\nScript 05 completed!\n")
cat("Next step: Run Aux scripts or audit.\n")
