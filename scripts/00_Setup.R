# =============================================================================
# Script 00: Setup - Package Installation and Configuration
# Paper: Systems biology analysis of Daphnia pulex under anthropogenic stress
#
# This script installs and loads all required R packages, defines global
# plotting parameters, and sets directory paths for reproducible execution.
# All subsequent scripts source this file for consistent configuration.
#
# Authors: J.A. Norambuena, L. Herrera-Belén, P. Poblete-Grant, J.F. Beltrán
# Year:    2026
# =============================================================================

# Function to install CRAN packages if not already installed
install_if_missing <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  if (length(new_packages)) {
    install.packages(new_packages, repos = "https://cloud.r-project.org/")
  }
}

# Function to install Bioconductor packages if not already installed
install_bioc_if_missing <- function(packages) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  if (length(new_packages)) {
    BiocManager::install(new_packages, ask = FALSE, update = FALSE)
  }
}

# -----------------------------------------------------------------------------
# CRAN Packages
# -----------------------------------------------------------------------------
cran_packages <- c(
  "tidyverse", # Data manipulation and visualization
  "readxl", # Read Excel files
  "writexl", # Write Excel files
  "ggplot2", # Visualization
  "pheatmap", # Heatmaps
  "RColorBrewer", # Color palettes
  "ggrepel", # Label placement in plots
  "gridExtra", # Multiple plot arrangement
  "scales", # ggplot scales
  "igraph", # Network analysis
  "ggraph", # Network visualization with ggplot
  "tidygraph", # Tidy graph manipulation
  "httr", # HTTP requests (for APIs)
  "jsonlite", # JSON parsing
  "VennDiagram", # Venn diagrams
  "UpSetR", # UpSet plots
  "cowplot", # Plot combination
  "viridis", # Color palettes
  "here" # Portable project-relative paths
)

cat("Installing CRAN packages...\n")
install_if_missing(cran_packages)

# -----------------------------------------------------------------------------
# Bioconductor Packages
# -----------------------------------------------------------------------------
bioc_packages <- c(
  "clusterProfiler", # GO and KEGG enrichment
  "enrichplot", # Enrichment visualization
  "DOSE", # Disease Ontology
  "pathview", # KEGG pathway visualization
  "AnnotationDbi", # Annotation framework
  "GO.db", # Gene Ontology database
  "org.Dm.eg.db", # Drosophila annotation (arthropod proxy)
  "STRINGdb", # Protein-protein interaction networks
  "biomaRt", # BioMart access
  "limma" # Differential expression analysis
)

cat("Installing Bioconductor packages...\n")
install_bioc_if_missing(bioc_packages)

# -----------------------------------------------------------------------------
# Load all packages
# -----------------------------------------------------------------------------
cat("\nLoading packages...\n")

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(writexl)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(ggrepel)
  library(gridExtra)
  library(scales)
  library(igraph)
  library(ggraph)
  library(tidygraph)
  library(httr)
  library(jsonlite)
  library(VennDiagram)
  library(UpSetR)
  library(cowplot)
  library(viridis)
  library(clusterProfiler)
  library(enrichplot)
  library(DOSE)
  library(pathview)
  library(AnnotationDbi)
  library(GO.db)
  library(STRINGdb)
  library(biomaRt)
  library(limma)
})

cat("All packages loaded successfully!\n")

# -----------------------------------------------------------------------------
# Global configuration
# -----------------------------------------------------------------------------
# Custom ggplot theme for publication
theme_paper <- theme_bw() +
  theme(
    text = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    panel.grid.minor = element_blank()
  )

theme_set(theme_paper)

# Colors for lake comparison
lake_colors <- c("Icalma" = "#2E86AB", "Llanquihue" = "#A23B72")

# Colors for regulation status
regulation_colors <- c(
  "Upregulated" = "#E63946", "Downregulated" = "#457B9D",
  "Not significant" = "#CCCCCC"
)

# -----------------------------------------------------------------------------
# Project directory detection (portable across machines)
# -----------------------------------------------------------------------------
# Priority: 1) RStudio active document, 2) here::here(), 3) script location
BASE_DIR <- tryCatch(
  {
    # Try RStudio detection first
    d <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
    if (!is.null(d) && d != "" && d != ".") d else stop("fallback")
  },
  error = function(e) {
    # Fallback: use here::here() which finds the project root
    # (looks for .git, .Rproj, .here, etc.)
    tryCatch(
      {
        here::here()
      },
      error = function(e2) {
        # Last resort: infer from this script's own location
        # Works when called via source("scripts/00_Setup.R") from project root
        if (exists("ofile", envir = parent.frame(4))) {
          dirname(dirname(get("ofile", envir = parent.frame(4))))
        } else {
          stop("Cannot detect project root. Please set BASE_DIR manually or run from the project directory.")
        }
      }
    )
  }
)

DATA_DIR <- file.path(BASE_DIR, "data")
RESULTS_DIR <- file.path(BASE_DIR, "results")
FIGURES_DIR <- file.path(BASE_DIR, "figures")
SCRIPTS_DIR <- file.path(BASE_DIR, "scripts")

# Create output directories if they don't exist
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)

cat("\nConfiguration complete.\n")
cat("Project root:", BASE_DIR, "\n")
cat("Data dir:    ", DATA_DIR, "\n")
cat("Results dir: ", RESULTS_DIR, "\n")
cat("Figures dir: ", FIGURES_DIR, "\n")

# -----------------------------------------------------------------------------
# Capture session info for reproducibility
# -----------------------------------------------------------------------------
capture_session_info <- function() {
  info <- sessionInfo()
  sink(file.path(RESULTS_DIR, "session_info.txt"))
  cat("Session info captured on:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  print(info)
  sink()
  cat("Session info saved to:", file.path(RESULTS_DIR, "session_info.txt"), "\n")
}
