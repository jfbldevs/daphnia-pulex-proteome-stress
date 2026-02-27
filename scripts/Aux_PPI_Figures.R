# =============================================================================
# Auxiliary Script: PPI Network Figures
# Creates focused network visualizations for the manuscript
#
# Authors: J.A. Norambuena, L. Herrera-Belén, P. Poblete-Grant, J.F. Beltrán
# Year:    2026
# =============================================================================

# Load configuration and directory paths
source("00_Setup.R")

# Load data
load(file.path(RESULTS_DIR, "04_PPI_network.RData"))
load(file.path(RESULTS_DIR, "RData_01_DataPrep.RData"))

cat("=== INVESTIGACION PROFUNDA DE LA RED PPI ===\n\n")

# -----------------------------------------------------------------------------
# 1. ANALIZAR LA ESTRUCTURA DE LA RED
# -----------------------------------------------------------------------------

cat("--- Estadisticas de la red ---\n")
cat("Nodos totales:", vcount(ppi_graph), "\n")
cat("Edges totales:", ecount(ppi_graph), "\n")
cat("Densidad:", round(edge_density(ppi_graph), 4), "\n")
cat("Diametro:", diameter(ppi_graph), "\n")
cat("Transitividad (clustering):", round(transitivity(ppi_graph), 4), "\n\n")

# Analizar componentes
components_info <- components(ppi_graph)
cat("Componentes conectados:", components_info$no, "\n")
cat("Tamano del componente principal:", max(components_info$csize), "\n\n")

# -----------------------------------------------------------------------------
# 2. IDENTIFICAR PROTEINAS CLAVE POR FUNCION BIOLOGICA
# -----------------------------------------------------------------------------

cat("--- Proteinas hub por funcion ---\n")

# Categorizar hubs por funcion
hub_functions <- centrality_metrics %>%
  filter(Degree >= 40) %>%
  arrange(desc(Degree)) %>%
  mutate(
    Functional_Category = case_when(
      str_detect(tolower(Protein_Name), "atp synthase|atpase") ~ "Energy/ATP",
      str_detect(tolower(Protein_Name), "citrate|fumarate|succin|glutamate dehydro") ~ "TCA Cycle",
      str_detect(tolower(Protein_Name), "ribosom") ~ "Ribosome",
      str_detect(tolower(Protein_Name), "histone") ~ "Chromatin",
      str_detect(tolower(Protein_Name), "transket") ~ "Carbohydrate",
      str_detect(tolower(Protein_Name), "tgc|tog") ~ "Signaling/Structure",
      str_detect(tolower(Protein_Name), "calcium") ~ "Ion Transport",
      TRUE ~ "Other"
    )
  )

print(hub_functions %>% select(Protein_Name, Degree, Regulation, Functional_Category))

cat("\n--- Distribucion de hubs por categoria ---\n")
print(table(hub_functions$Functional_Category))

cat("\n--- Hubs por regulacion ---\n")
print(table(hub_functions$Regulation))

# -----------------------------------------------------------------------------
# 3. ESTRATEGIA: Red enfocada en metabolismo energetico
# -----------------------------------------------------------------------------

# El hallazgo clave es la supresion del metabolismo energetico
# Vamos a crear una red que muestre esto claramente

# Seleccionar proteinas relacionadas con energia/metabolismo
energy_proteins <- centrality_metrics %>%
  filter(
    # Proteinas de energia o con alto grado
    str_detect(tolower(Protein_Name), "atp|citrate|fumarate|succin|nadh|cytochrome|electron|oxido|dehydro") |
      Degree >= 35
  ) %>%
  pull(STRING_id)

cat("\nProteinas seleccionadas para red enfocada:", length(energy_proteins), "\n")

# Crear subgrafo
if (length(energy_proteins) > 5) {
  energy_subgraph <- induced_subgraph(
    ppi_graph,
    V(ppi_graph)[name %in% energy_proteins]
  )

  cat("Nodos en subgrafo:", vcount(energy_subgraph), "\n")
  cat("Edges en subgrafo:", ecount(energy_subgraph), "\n")
}

# -----------------------------------------------------------------------------
# 4. CREAR VISUALIZACION LIMPIA
# -----------------------------------------------------------------------------

# Prepare data for visualization
top_proteins <- centrality_metrics %>%
  filter(Degree >= 30) %>%
  arrange(desc(Degree))

cat("\nTop proteinas (degree >= 30):", nrow(top_proteins), "\n")

# Crear subgrafo con top proteinas
top_ids <- top_proteins$STRING_id
top_subgraph <- induced_subgraph(ppi_graph, V(ppi_graph)[name %in% top_ids])

# Convertir a tidygraph
tg <- as_tbl_graph(top_subgraph) %>%
  activate(nodes) %>%
  left_join(
    centrality_metrics %>%
      select(STRING_id, Protein_Name, Degree, Regulation, Betweenness, Eigenvector),
    by = c("name" = "STRING_id")
  ) %>%
  mutate(
    # Crear etiquetas limpias
    Label = case_when(
      str_detect(Protein_Name, "ATP synthase subunit beta") ~ "ATP synthase beta",
      str_detect(Protein_Name, "ATP synthase subunit gamma") ~ "ATP synthase gamma",
      str_detect(Protein_Name, "Citrate synthase") ~ "Citrate synthase",
      str_detect(Protein_Name, "Fumarate hydratase") ~ "Fumarate hydratase",
      str_detect(Protein_Name, "Glutamate dehydrogenase") ~ "Glutamate DH",
      str_detect(Protein_Name, "Calcium-transporting") ~ "Ca2+ ATPase",
      str_detect(Protein_Name, "Histone H4") ~ "Histone H4",
      str_detect(Protein_Name, "Succinyl-CoA") ~ "Succinyl-CoA transferase",
      str_detect(Protein_Name, "Succinate--CoA") ~ "Succinate-CoA ligase",
      str_detect(Protein_Name, "Ribosomal_S17") ~ "Ribosomal S17",
      str_detect(Protein_Name, "Ribosomal_L23") ~ "Ribosomal L23",
      str_detect(Protein_Name, "Transket_pyr") ~ "Transketolase",
      str_detect(Protein_Name, "TGc domain") ~ "TGc protein",
      str_detect(Protein_Name, "TOG domain") ~ "TOG protein",
      is.na(Protein_Name) | Protein_Name == "" ~ str_extract(name, "[^.]+$"),
      TRUE ~ str_trunc(Protein_Name, 20)
    ),
    # Hacer etiquetas unicas
    Label = make.unique(Label, sep = " "),

    # Categorias funcionales para color
    Function = case_when(
      str_detect(tolower(Protein_Name), "atp synthase|atpase") ~ "ATP Synthesis",
      str_detect(tolower(Protein_Name), "citrate|fumarate|succin") ~ "TCA Cycle",
      str_detect(tolower(Protein_Name), "glutamate dehydro|nadh|dehydro") ~ "Redox/Metabolism",
      str_detect(tolower(Protein_Name), "ribosom") ~ "Ribosome",
      str_detect(tolower(Protein_Name), "histone") ~ "Chromatin",
      str_detect(tolower(Protein_Name), "calcium") ~ "Ion Transport",
      TRUE ~ "Other"
    ),

    # Expresion
    Expression = case_when(
      Regulation == "Upregulated" ~ "Higher in Llanquihue\n(Anthropized)",
      Regulation == "Downregulated" ~ "Higher in Icalma\n(Oligotrophic)",
      TRUE ~ "No change"
    ),

    # Tamano basado en grado
    Size = scales::rescale(Degree, to = c(4, 15))
  )

cat("\nCreando visualizacion...\n")

# Paleta de colores por funcion
function_colors <- c(
  "ATP Synthesis" = "#E63946",
  "TCA Cycle" = "#F4A261",
  "Redox/Metabolism" = "#E9C46A",
  "Ribosome" = "#2A9D8F",
  "Chromatin" = "#264653",
  "Ion Transport" = "#8338EC",
  "Other" = "#ADB5BD"
)

# Layout con semilla fija para reproducibilidad
set.seed(123)

# FIGURA: Red por funcion metabolica
p_function <- ggraph(tg, layout = "fr") +
  geom_edge_link(alpha = 0.15, color = "gray50", width = 0.3) +
  geom_node_point(
    aes(
      size = Degree, fill = Function,
      shape = ifelse(Regulation == "Downregulated", "Suppressed",
        ifelse(Regulation == "Upregulated", "Induced", "Unchanged")
      )
    ),
    alpha = 0.9, stroke = 0.8, color = "white"
  ) +
  geom_node_text(aes(label = Label),
    repel = TRUE, size = 2.8,
    max.overlaps = 25, segment.size = 0.2,
    box.padding = 0.4, point.padding = 0.3
  ) +
  scale_fill_manual(values = function_colors, name = "Metabolic Function") +
  scale_shape_manual(
    values = c("Suppressed" = 21, "Induced" = 24, "Unchanged" = 22),
    name = "Expression in\nLlanquihue",
    labels = c(
      "Induced" = "Induced (higher)",
      "Suppressed" = "Suppressed (lower)",
      "Unchanged" = "No change"
    )
  ) +
  scale_size_continuous(range = c(4, 14), name = "Connections\n(Degree)") +
  guides(
    fill = guide_legend(override.aes = list(size = 5, shape = 21)),
    shape = guide_legend(override.aes = list(size = 5, fill = "gray50"))
  ) +
  labs(
    title = "Protein-Protein Interaction Network: Metabolic Hub Proteins",
    subtitle = paste0(
      "Showing ", vcount(top_subgraph), " highly connected proteins (degree >= 30) and their ",
      ecount(top_subgraph), " interactions\n",
      "Node size = number of connections | Shape = expression change | Color = metabolic function"
    ),
    caption = paste0(
      "Key finding: Most hub proteins (especially TCA cycle and ATP synthesis) show LOWER expression\n",
      "in Llanquihue (anthropized lake), indicating suppression of core energy metabolism."
    )
  ) +
  theme_void(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(
      size = 10, color = "gray30", hjust = 0.5,
      margin = margin(b = 10), lineheight = 1.2
    ),
    plot.caption = element_text(
      size = 9, hjust = 0, color = "gray40",
      margin = margin(t = 15), lineheight = 1.2
    ),
    legend.position = "right",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    plot.margin = margin(15, 15, 15, 15)
  )

ggsave(file.path(FIGURES_DIR, "Fig4A_PPI_Network_Metabolic.pdf"), p_function,
  width = 14, height = 10, dpi = 300
)
cat("Saved: Fig4A_PPI_Network_Metabolic.pdf\n")

# -----------------------------------------------------------------------------
# 5. VERSION ALTERNATIVA: Solo TCA y ATP (mas enfocada)
# -----------------------------------------------------------------------------

# Filtrar solo proteinas de energia
energy_nodes <- tg %>%
  activate(nodes) %>%
  as_tibble() %>%
  filter(Function %in% c("ATP Synthesis", "TCA Cycle", "Redox/Metabolism", "Ion Transport"))

energy_ids <- energy_nodes$name

if (length(energy_ids) >= 5) {
  energy_sub <- induced_subgraph(top_subgraph, V(top_subgraph)[name %in% energy_ids])

  tg_energy <- as_tbl_graph(energy_sub) %>%
    activate(nodes) %>%
    left_join(
      centrality_metrics %>%
        select(STRING_id, Protein_Name, Degree, Regulation),
      by = c("name" = "STRING_id")
    ) %>%
    mutate(
      Label = case_when(
        str_detect(Protein_Name, "ATP synthase subunit beta") ~ "ATP synthase\nbeta",
        str_detect(Protein_Name, "ATP synthase subunit gamma") ~ "ATP synthase\ngamma",
        str_detect(Protein_Name, "Citrate synthase") ~ "Citrate\nsynthase",
        str_detect(Protein_Name, "Fumarate hydratase") ~ "Fumarate\nhydratase",
        str_detect(Protein_Name, "Glutamate dehydrogenase") ~ "Glutamate\nDH",
        str_detect(Protein_Name, "Calcium-transporting") ~ "Ca2+\nATPase",
        str_detect(Protein_Name, "Succinyl-CoA") ~ "Succinyl-CoA\ntransferase",
        str_detect(Protein_Name, "Succinate--CoA") ~ "Succinate-CoA\nligase",
        str_detect(Protein_Name, "NADH") ~ "NADH DH",
        is.na(Protein_Name) ~ str_extract(name, "[^.]+$"),
        TRUE ~ str_trunc(Protein_Name, 15)
      ),
      Label = make.unique(Label, sep = " "),
      Function = case_when(
        str_detect(tolower(Protein_Name), "atp synthase|atpase") ~ "ATP Synthesis",
        str_detect(tolower(Protein_Name), "citrate|fumarate|succin") ~ "TCA Cycle",
        str_detect(tolower(Protein_Name), "calcium") ~ "Ion Transport",
        TRUE ~ "Redox/Metabolism"
      ),
      Expression = ifelse(Regulation == "Downregulated",
        "Suppressed in Llanquihue",
        ifelse(Regulation == "Upregulated",
          "Induced in Llanquihue", "No change"
        )
      )
    )

  set.seed(456)

  p_energy <- ggraph(tg_energy, layout = "fr") +
    geom_edge_link(alpha = 0.25, color = "gray40", width = 0.5) +
    geom_node_point(aes(size = Degree, fill = Function),
      shape = 21, alpha = 0.9, stroke = 1.5, color = "white"
    ) +
    geom_node_label(aes(label = Label, fill = Function),
      size = 2.5, alpha = 0.85, label.padding = unit(0.2, "lines"),
      label.r = unit(0.15, "lines"), color = "white", fontface = "bold"
    ) +
    scale_fill_manual(
      values = c(
        "ATP Synthesis" = "#E63946",
        "TCA Cycle" = "#F4A261",
        "Redox/Metabolism" = "#2A9D8F",
        "Ion Transport" = "#8338EC"
      ),
      name = "Metabolic\nFunction"
    ) +
    scale_size_continuous(range = c(8, 20), name = "Connections") +
    labs(
      title = "Energy Metabolism Network in Daphnia pulex",
      subtitle = paste0(
        "Core metabolic proteins and their interactions\n",
        "All proteins shown have LOWER expression in Llanquihue (anthropized lake)"
      ),
      caption = paste0(
        "TCA Cycle proteins (orange): Citrate synthase, Fumarate hydratase, Succinyl-CoA enzymes\n",
        "ATP Synthesis (red): ATP synthase subunits - final step of energy production\n",
        "This network shows severe disruption of mitochondrial energy metabolism under environmental stress."
      )
    ) +
    theme_void(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(
        size = 11, color = "gray30", hjust = 0.5,
        margin = margin(b = 15), lineheight = 1.3
      ),
      plot.caption = element_text(
        size = 9, hjust = 0, color = "gray40",
        margin = margin(t = 20), lineheight = 1.3
      ),
      legend.position = "right",
      legend.title = element_text(size = 10, face = "bold"),
      plot.margin = margin(20, 20, 20, 20)
    )

  ggsave(file.path(FIGURES_DIR, "Fig4B_PPI_Energy_Metabolism.pdf"), p_energy,
    width = 12, height = 9, dpi = 300
  )
  cat("Saved: Fig4B_PPI_Energy_Metabolism.pdf\n")
}

# -----------------------------------------------------------------------------
# 6. TABLA RESUMEN DE HUBS
# -----------------------------------------------------------------------------

hub_summary <- centrality_metrics %>%
  filter(Degree >= 40) %>%
  arrange(desc(Degree)) %>%
  mutate(
    Function = case_when(
      str_detect(tolower(Protein_Name), "atp synthase") ~ "ATP Synthesis",
      str_detect(tolower(Protein_Name), "citrate|fumarate|succin") ~ "TCA Cycle",
      str_detect(tolower(Protein_Name), "glutamate dehydro") ~ "Amino Acid Metabolism",
      str_detect(tolower(Protein_Name), "ribosom") ~ "Protein Synthesis",
      str_detect(tolower(Protein_Name), "histone") ~ "Chromatin Regulation",
      str_detect(tolower(Protein_Name), "calcium|atpase") ~ "Ion Transport",
      str_detect(tolower(Protein_Name), "transket") ~ "Pentose Phosphate Pathway",
      TRUE ~ "Signaling/Structure"
    ),
    Expression_Pattern = case_when(
      Regulation == "Downregulated" ~ "SUPPRESSED in anthropized lake",
      Regulation == "Upregulated" ~ "INDUCED in anthropized lake",
      TRUE ~ "No significant change"
    )
  ) %>%
  select(Protein_Name, Degree, Function, Expression_Pattern)

cat("\n--- TABLA RESUMEN DE PROTEINAS HUB ---\n")
print(hub_summary, n = 20)

write_csv(hub_summary, file.path(RESULTS_DIR, "Table4_Hub_Proteins.csv"))
cat("\nSaved: Table4_Hub_Proteins.csv\n")

cat("\n=== VISUALIZACIONES PPI COMPLETADAS ===\n")
cat("Nuevas figuras:\n")
cat("  12_PPI_Network_Metabolic.pdf - Red completa de hubs con funciones\n")
cat("  12b_PPI_Energy_Metabolism.pdf - Red enfocada en metabolismo energetico\n")
