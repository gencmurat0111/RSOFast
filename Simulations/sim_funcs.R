# ============================================================================
# FILE: sim_funcs.R
# ============================================================================
# Helper functions for RSO simulation study
# ============================================================================

library(MASS)
library(RSO)  
source("RSOPAPER.R")

# ----------------------------------------------------------------------------
# Data generation function 
# ----------------------------------------------------------------------------
simulate_data <- function(n_train, n_test, corr, sigma, beta_true) {
  p <- length(beta_true)
  
  # correlation structure (AR(1))
  Sigma <- corr^abs(outer(1:p, 1:p, "-"))
  
  # Generate predictors
  X_all <- MASS::mvrnorm(n_train + n_test, mu = rep(0, p), Sigma = Sigma)
  
  # Split into train/test
  X_train <- X_all[1:n_train, ]
  X_test <- X_all[(n_train + 1):(n_train + n_test), ]
  
  # Generate response
  y_train <- X_train %*% beta_true + sigma * rnorm(n_train)
  y_test <- X_test %*% beta_true + sigma * rnorm(n_test)
  
  # scale variables (as detailed in paper)
  X_train <- scale(X_train) * sqrt(n_train/(n_train-1))
  X_test <- scale(X_test) * sqrt(n_test/(n_test-1))
  y_train <- y_train - mean(y_train)
  y_test <- y_test - mean(y_test)
  
  return(list(
    trainx = X_train, trainy = y_train,
    testx = X_test, testy = y_test
  ))
}

# ----------------------------------------------------------------------------
# Compute classification metrics for variable selection
# ----------------------------------------------------------------------------
compute_metrics <- function(beta_est, beta_true) {
  TP <- sum(beta_est != 0 & beta_true != 0)
  FP <- sum(beta_est != 0 & beta_true == 0)
  TN <- sum(beta_est == 0 & beta_true == 0)
  FN <- sum(beta_est == 0 & beta_true != 0)
  
  fpr <- FP / (FP + TN)
  fnr <- FN / (TP + FN)
  precision <- if (TP + FP > 0) TP / (TP + FP) else NA
  recall <- TP / (TP + FN)
  f1 <- if (!is.na(precision) && precision + recall > 0) {
    2 * precision * recall / (precision + recall)
  } else NA
  accuracy <- (TP + TN) / (TP + TN + FP + FN)
  
  return(c(TP = TP, FP = FP, TN = TN, FN = FN,
           fpr = fpr, fnr = fnr, 
           precision = precision, recall = recall, 
           f1 = f1, accuracy = accuracy))
}

# ----------------------------------------------------------------------------
# Run RSO_0 and return results (with timing)
# ----------------------------------------------------------------------------
run_rso0 <- function(data, tau, penalty_factor = NULL, beta_true = NULL) {
  X <- data$trainx
  y <- data$trainy
  X_test <- data$testx
  y_test <- data$testy
  
  p <- ncol(X)
  if (is.null(penalty_factor)) penalty_factor <- rep(1, p)
  
  # Run RSO_0 with timing
  time_start <- proc.time()
  result <- RSO_0(X, y, tau = tau, penalty.factor = penalty_factor)
  time_end <- proc.time()
  elapsed_time <- (time_end - time_start)[3]
  
  beta_est <- result$coef
  
  # Test MSE
  pred_test <- X_test %*% beta_est
  test_mse <- mean((y_test - pred_test)^2)
  
  # Number of active variables
  active <- sum(beta_est != 0)
  
  # Metrics if true beta provided
  metrics <- NULL
  if (!is.null(beta_true)) {
    metrics <- compute_metrics(beta_est, beta_true)
  }
  
  return(list(
    beta = beta_est,
    active = active,
    test_mse = test_mse,
    time = elapsed_time,
    metrics = metrics
  ))
}

# ----------------------------------------------------------------------------
# Run RSOFast and return results (with timing)
# ----------------------------------------------------------------------------
run_rsofast <- function(data, tau, penalty_factor = NULL, beta_true = NULL) {
  X <- data$trainx
  y <- data$trainy
  X_test <- data$testx
  y_test <- data$testy
  
  p <- ncol(X)
  if (is.null(penalty_factor)) penalty_factor <- rep(1, p)
  
  # Run RSOFast with timing
  time_start <- proc.time()
  result <- RSOFast(X, y, tau = tau, penalty.factor = penalty_factor)
  time_end <- proc.time()
  elapsed_time <- (time_end - time_start)[3]
  
  beta_est <- result$coef
  
  # Test MSE
  pred_test <- X_test %*% beta_est
  test_mse <- mean((y_test - pred_test)^2)
  
  # Number of active variables
  active <- sum(beta_est != 0)
  
  # Metrics if true beta provided
  metrics <- NULL
  if (!is.null(beta_true)) {
    metrics <- compute_metrics(beta_est, beta_true)
  }
  
  return(list(
    beta = beta_est,
    active = active,
    test_mse = test_mse,
    time = elapsed_time,
    metrics = metrics
  ))
}

# ----------------------------------------------------------------------------
# Compare RSO_0 and RSOFast on the same dataset
# ----------------------------------------------------------------------------
compare_methods <- function(data, tau, beta_true, penalty_factor = NULL) {
  res0 <- run_rso0(data, tau, penalty_factor, beta_true)
  res_fast <- run_rsofast(data, tau, penalty_factor, beta_true)
  
  return(list(rso0 = res0, rsofast = res_fast))
}


# ----------------------------------------------------------------------------
# Generate true beta vector based on sparsity type
# ----------------------------------------------------------------------------
generate_beta_true <- function(p, sparsity_type) {
  beta_true <- rep(0, p)
  
  if (sparsity_type == "all") {
    n_important <- p
  } else if (sparsity_type == "p/5") {
    n_important <- max(1, round(p / 5))
  } else if (sparsity_type == "p/10") {
    n_important <- max(1, round(p / 10))
  }
  
  # Random non-zero coefficients between 1 and 3, with random signs
  if (n_important > 0) {
    beta_true[1:n_important] <- runif(n_important, 1, 3) * 
      sample(c(-1, 1), n_important, replace = TRUE)
  }
  
  return(beta_true)
}

# ----------------------------------------------------------------------------
# Run single scenario
# ----------------------------------------------------------------------------
run_scenario <- function(n, p, sparsity_type, sigma, n_test, tau, corr, nsim = 100) {
  
  beta_true <- generate_beta_true(p, sparsity_type)
  n_important <- sum(beta_true != 0)
  
  cat(sprintf("\n=== n=%d, p=%d, sparsity=%s (%d important), sigma=%d ===\n", 
              n, p, sparsity_type, n_important, sigma))
  
  # Result vectors
  active_rso0 <- numeric(nsim)
  active_rsofast <- numeric(nsim)
  test_mse_rso0 <- numeric(nsim)
  test_mse_rsofast <- numeric(nsim)
  time_rso0 <- numeric(nsim)
  time_rsofast <- numeric(nsim)
  fpr_rso0 <- numeric(nsim)
  fpr_rsofast <- numeric(nsim)
  fnr_rso0 <- numeric(nsim)
  fnr_rsofast <- numeric(nsim)
  f1_rso0 <- numeric(nsim)
  f1_rsofast <- numeric(nsim)
  
  for (i in 1:nsim) {
    if (i %% 20 == 0) cat(sprintf("  Replication %d/%d\n", i, nsim))
    
    # Generate data
    data <- simulate_data(n_train = n, n_test = n_test, 
                          corr = corr, sigma = sigma, 
                          beta_true = beta_true)
    
    # Compare methods
    comp <- compare_methods(data, tau = tau, beta_true = beta_true)
    
    # Store results
    active_rso0[i] <- comp$rso0$active
    active_rsofast[i] <- comp$rsofast$active
    test_mse_rso0[i] <- comp$rso0$test_mse
    test_mse_rsofast[i] <- comp$rsofast$test_mse
    time_rso0[i] <- comp$rso0$time
    time_rsofast[i] <- comp$rsofast$time
    
    if (!is.null(comp$rso0$metrics)) {
      fpr_rso0[i] <- comp$rso0$metrics["fpr"]
      fnr_rso0[i] <- comp$rso0$metrics["fnr"]
      f1_rso0[i] <- comp$rso0$metrics["f1"]
    }
    
    if (!is.null(comp$rsofast$metrics)) {
      fpr_rsofast[i] <- comp$rsofast$metrics["fpr"]
      fnr_rsofast[i] <- comp$rsofast$metrics["fnr"]
      f1_rsofast[i] <- comp$rsofast$metrics["f1"]
    }
  }
  
  # Create summary row
  summary_row <- c(
    active_mean_rso0 = mean(active_rso0, na.rm = TRUE),
    active_sd_rso0 = sd(active_rso0, na.rm = TRUE),
    active_mean_rsofast = mean(active_rsofast, na.rm = TRUE),
    active_sd_rsofast = sd(active_rsofast, na.rm = TRUE),
    
    test_mse_mean_rso0 = mean(test_mse_rso0, na.rm = TRUE),
    test_mse_sd_rso0 = sd(test_mse_rso0, na.rm = TRUE),
    test_mse_mean_rsofast = mean(test_mse_rsofast, na.rm = TRUE),
    test_mse_sd_rsofast = sd(test_mse_rsofast, na.rm = TRUE),
    
    time_mean_rso0 = mean(time_rso0, na.rm = TRUE),
    time_sd_rso0 = sd(time_rso0, na.rm = TRUE),
    time_mean_rsofast = mean(time_rsofast, na.rm = TRUE),
    time_sd_rsofast = sd(time_rsofast, na.rm = TRUE),
    speedup = mean(time_rso0, na.rm = TRUE) / mean(time_rsofast, na.rm = TRUE),
    
    fpr_mean_rso0 = mean(fpr_rso0, na.rm = TRUE),
    fpr_sd_rso0 = sd(fpr_rso0, na.rm = TRUE),
    fpr_mean_rsofast = mean(fpr_rsofast, na.rm = TRUE),
    fpr_sd_rsofast = sd(fpr_rsofast, na.rm = TRUE),
    
    fnr_mean_rso0 = mean(fnr_rso0, na.rm = TRUE),
    fnr_sd_rso0 = sd(fnr_rso0, na.rm = TRUE),
    fnr_mean_rsofast = mean(fnr_rsofast, na.rm = TRUE),
    fnr_sd_rsofast = sd(fnr_rsofast, na.rm = TRUE),
    
    f1_mean_rso0 = mean(f1_rso0, na.rm = TRUE),
    f1_sd_rso0 = sd(f1_rso0, na.rm = TRUE),
    f1_mean_rsofast = mean(f1_rsofast, na.rm = TRUE),
    f1_sd_rsofast = sd(f1_rsofast, na.rm = TRUE)
  )
  
  return(list(
    params = c(n = n, p = p, sparsity = sparsity_type, 
               n_important = n_important, sigma = sigma),
    results = summary_row
  ))
}

# ----------------------------------------------------------------------------
# Run all scenarios from np list
# ----------------------------------------------------------------------------
run_all_scenarios <- function(np_list, sparsity_types, sigma_vals, 
                              n_test, tau, corr, nsim = 100) {
  all_results <- list()
  idx <- 1
  
  total_scenarios <- length(np_list) * length(sparsity_types) * length(sigma_vals)
  
  cat(sprintf("\nTotal scenarios: %d\n", total_scenarios))
  cat(sprintf("Replications per scenario: %d\n\n", nsim))
  
  for (np in np_list) {
    n <- np[1]
    p <- np[2]
    
    for (sp_type in sparsity_types) {
      for (sigma in sigma_vals) {
        cat("\n", rep("=", 60), "\n", sep = "")
        cat(sprintf("SCENARIO %d: n=%d, p=%d, sparsity=%s, sigma=%d\n", 
                    idx, n, p, sp_type, sigma))
        cat(rep("=", 60), "\n", sep = "")
        
        res <- run_scenario(n, p, sp_type, sigma, n_test, tau, corr, nsim = nsim)
        all_results[[idx]] <- res
        idx <- idx + 1
      }
    }
  }
  
  return(all_results)
}

# ----------------------------------------------------------------------------
# Convert results to data frame (same metrics side by side)
# ----------------------------------------------------------------------------
results_to_table <- function(all_results) {
  n_scenarios <- length(all_results)
  
  table_data <- data.frame(
    scenario_id = 1:n_scenarios,
    n = integer(n_scenarios),
    p = integer(n_scenarios),
    sparsity = character(n_scenarios),
    n_important = integer(n_scenarios),
    sigma = integer(n_scenarios),
    
    # Active
    rso0_active = numeric(n_scenarios),
    rso0_active_sd = numeric(n_scenarios),
    rsofast_active = numeric(n_scenarios),
    rsofast_active_sd = numeric(n_scenarios),
    
    # Test MSE
    rso0_test_mse = numeric(n_scenarios),
    rso0_test_mse_sd = numeric(n_scenarios),
    rsofast_test_mse = numeric(n_scenarios),
    rsofast_test_mse_sd = numeric(n_scenarios),
    
    # Time
    rso0_time = numeric(n_scenarios),
    rso0_time_sd = numeric(n_scenarios),
    rsofast_time = numeric(n_scenarios),
    rsofast_time_sd = numeric(n_scenarios),
    speedup = numeric(n_scenarios),
    
    # FPR
    rso0_fpr = numeric(n_scenarios),
    rso0_fpr_sd = numeric(n_scenarios),
    rsofast_fpr = numeric(n_scenarios),
    rsofast_fpr_sd = numeric(n_scenarios),
    
    # FNR
    rso0_fnr = numeric(n_scenarios),
    rso0_fnr_sd = numeric(n_scenarios),
    rsofast_fnr = numeric(n_scenarios),
    rsofast_fnr_sd = numeric(n_scenarios),
    
    # F1
    rso0_f1 = numeric(n_scenarios),
    rso0_f1_sd = numeric(n_scenarios),
    rsofast_f1 = numeric(n_scenarios),
    rsofast_f1_sd = numeric(n_scenarios)
  )
  
  for (i in 1:n_scenarios) {
    res <- all_results[[i]]
    table_data[i, "n"] <- res$params["n"]
    table_data[i, "p"] <- res$params["p"]
    table_data[i, "sparsity"] <- res$params["sparsity"]
    table_data[i, "n_important"] <- res$params["n_important"]
    table_data[i, "sigma"] <- res$params["sigma"]
    
    # Fill results
    table_data[i, "rso0_active"] <- res$results["active_mean_rso0"]
    table_data[i, "rso0_active_sd"] <- res$results["active_sd_rso0"]
    table_data[i, "rsofast_active"] <- res$results["active_mean_rsofast"]
    table_data[i, "rsofast_active_sd"] <- res$results["active_sd_rsofast"]
    
    table_data[i, "rso0_test_mse"] <- res$results["test_mse_mean_rso0"]
    table_data[i, "rso0_test_mse_sd"] <- res$results["test_mse_sd_rso0"]
    table_data[i, "rsofast_test_mse"] <- res$results["test_mse_mean_rsofast"]
    table_data[i, "rsofast_test_mse_sd"] <- res$results["test_mse_sd_rsofast"]
    
    table_data[i, "rso0_time"] <- res$results["time_mean_rso0"]
    table_data[i, "rso0_time_sd"] <- res$results["time_sd_rso0"]
    table_data[i, "rsofast_time"] <- res$results["time_mean_rsofast"]
    table_data[i, "rsofast_time_sd"] <- res$results["time_sd_rsofast"]
    table_data[i, "speedup"] <- res$results["speedup"]
    
    table_data[i, "rso0_fpr"] <- res$results["fpr_mean_rso0"]
    table_data[i, "rso0_fpr_sd"] <- res$results["fpr_sd_rso0"]
    table_data[i, "rsofast_fpr"] <- res$results["fpr_mean_rsofast"]
    table_data[i, "rsofast_fpr_sd"] <- res$results["fpr_sd_rsofast"]
    
    table_data[i, "rso0_fnr"] <- res$results["fnr_mean_rso0"]
    table_data[i, "rso0_fnr_sd"] <- res$results["fnr_sd_rso0"]
    table_data[i, "rsofast_fnr"] <- res$results["fnr_mean_rsofast"]
    table_data[i, "rsofast_fnr_sd"] <- res$results["fnr_sd_rsofast"]
    
    table_data[i, "rso0_f1"] <- res$results["f1_mean_rso0"]
    table_data[i, "rso0_f1_sd"] <- res$results["f1_sd_rso0"]
    table_data[i, "rsofast_f1"] <- res$results["f1_mean_rsofast"]
    table_data[i, "rsofast_f1_sd"] <- res$results["f1_sd_rsofast"]
  }
  
  return(table_data)
}