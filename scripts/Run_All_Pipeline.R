# =============================================================================
# Run All Pipeline — Master Script
# Paper: Integrative bioinformatics analysis reveals disrupted metabolic
#        pathways and key hub proteins in the Daphnia pulex proteome
#        under anthropogenic stress
#
# This script runs the entire analysis pipeline sequentially.
# Run from the project root directory:
#   Rscript scripts/Run_All_Pipeline.R
#
# Authors: J.A. Norambuena, L. Herrera-Belén, P. Poblete-Grant, J.F. Beltrán
# Year:    2026
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════════╗\n")
cat("║  Daphnia pulex Systems Biology Pipeline                                ║\n")
cat("╚══════════════════════════════════════════════════════════════════════════╝\n\n")

# Set working directory to scripts/ so source() calls work
script_dir <- tryCatch(
  {
    dirname(rstudioapi::getActiveDocumentContext()$path)
  },
  error = function(e) {
    # When called via Rscript, use commandArgs to find location
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
      dirname(normalizePath(sub("--file=", "", file_arg)))
    } else {
      file.path(getwd(), "scripts")
    }
  }
)

setwd(script_dir)
cat("Working directory:", getwd(), "\n\n")

# =============================================================================
# Pipeline execution
# =============================================================================

pipeline_steps <- list(
  list(
    name = "00 - Setup & Configuration",
    script = "00_Setup.R"
  ),
  list(
    name = "01 - Data Preparation & Differential Expression",
    script = "01_Differential_Expression.R"
  ),
  list(
    name = "02 - GO Enrichment (ORA)",
    script = "02_GO_Enrichment_ORA.R"
  ),
  list(
    name = "03 - GSEA Analysis",
    script = "03_GSEA_Analysis.R"
  ),
  list(
    name = "04 - Metabolic Pathway Analysis",
    script = "04_Pathway_Analysis.R"
  ),
  list(
    name = "05 - PPI Network Analysis",
    script = "05_PPI_Network.R"
  )
)

results <- list()
start_time <- Sys.time()

for (step in pipeline_steps) {
  cat("\n")
  cat("====================================================================\n")
  cat(paste0("  STEP: ", step$name, "\n"))
  cat("====================================================================\n\n")

  step_start <- Sys.time()

  result <- tryCatch(
    {
      source(step$script, local = FALSE)
      list(status = "OK", time = Sys.time() - step_start)
    },
    error = function(e) {
      list(status = "FAILED", error = e$message, time = Sys.time() - step_start)
    }
  )

  results[[step$name]] <- result

  if (result$status == "OK") {
    cat(paste0(
      "\n✓ ", step$name, " completed in ",
      round(as.numeric(result$time, units = "secs"), 1), "s\n"
    ))
  } else {
    cat(paste0("\n✗ ", step$name, " FAILED: ", result$error, "\n"))
    cat("Pipeline stopped. Fix the error and re-run.\n")
    break
  }
}

# =============================================================================
# Capture session info for reproducibility
# =============================================================================

if (exists("capture_session_info")) {
  capture_session_info()
}

# =============================================================================
# Summary
# =============================================================================

total_time <- Sys.time() - start_time

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════════╗\n")
cat("║                        PIPELINE SUMMARY                                ║\n")
cat("╚══════════════════════════════════════════════════════════════════════════╝\n\n")

for (name in names(results)) {
  r <- results[[name]]
  status <- ifelse(r$status == "OK", "✓", "✗")
  time_str <- round(as.numeric(r$time, units = "secs"), 1)
  cat(sprintf("  %s %-50s %6.1fs\n", status, name, time_str))
}

cat(sprintf("\n  Total time: %.1f seconds\n", as.numeric(total_time, units = "secs")))

n_ok <- sum(sapply(results, function(x) x$status == "OK"))
n_total <- length(results)

if (n_ok == n_total) {
  cat("\n  All steps completed successfully!\n")
  cat("\n  Outputs:\n")
  cat("    results/  — Tables, data files, R objects\n")
  cat("    figures/  — Publication-quality figures\n")
  cat("    results/session_info.txt — R session info\n")
} else {
  cat(paste0("\n  ", n_ok, "/", n_total, " steps completed. Check errors above.\n"))
}

cat("\n")
