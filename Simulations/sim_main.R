# ============================================================================
# FILE: sim_main.R
# ============================================================================
# Main simulation script for RSO comparison study
# ============================================================================

rm(list = ls())
library(RSO)      
source("sim_funcs.R")

# ----------------------------------------------------------------------------
# Simulation parameters
# ----------------------------------------------------------------------------
set.seed(2024)

# Define (n, p) pairs
np_list <- list(
  c(60, 10),
  c(60, 20),
  c(100, 20),
  c(100, 40),
  c(500, 40),
  c(500, 60)
)

# Other parameters
sparsity_types <- c("all", "p/5", "p/10")
sigma_vals <- c(1, 5)

# Fixed parameters
n_test <- 200
nsim <- 100
tau <- 1.0
corr <- 0.5

# ----------------------------------------------------------------------------
# Main execution
# ----------------------------------------------------------------------------
{
  cat("========================================\n")
  cat("RSO Simulation Study\n")
  cat("========================================\n")
  cat(sprintf("Number of (n,p) pairs: %d\n", length(np_list)))
  cat(sprintf("Sparsity types: %s\n", paste(sparsity_types, collapse=", ")))
  cat(sprintf("Sigma values: %s\n", paste(sigma_vals, collapse=", ")))
  cat(sprintf("Total scenarios: %d\n", length(np_list) * length(sparsity_types) * length(sigma_vals)))
  cat(sprintf("Replications per scenario: %d\n", nsim))
  cat("========================================\n")
  
  all_results <- run_all_scenarios(np_list, sparsity_types, sigma_vals, 
                                   n_test, tau, corr, nsim = nsim)
  
  # Convert to table
  results_table <- results_to_table(all_results)
  
  # Save results
  write.table(results_table, "RSO_simulation_results.csv", 
              row.names = FALSE,
              sep = ";")
  
  # Also save as RDS
  saveRDS(results_table, "RSO_simulation_results.rds")
  
  cat("\n========================================\n")
  cat("Simulation completed!\n")
  cat("Results saved to:\n")
  cat("  - RSO_simulation_results.csv\n")
  cat("  - RSO_simulation_results.rds\n")
  cat("========================================\n")
  
  # Print summary
  print(results_table)
}


