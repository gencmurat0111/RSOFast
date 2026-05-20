# ============================================================================
# Real Data Application: Superconductivity Dataset (n_train = 200)
# Solution Path: X-axis = tau, Y-axis = beta_hat (coefficient estimates)
# Following Wu (2021) methodology
# ============================================================================

library(RSO)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(xlsx)

# ----------------------------------------------------------------------------
# 1. Load and prepare the data
#
# https://archive.ics.uci.edu/ml/datasets/Superconductivty+Data
# ----------------------------------------------------------------------------
data <- read.csv("superconduct.csv", header = TRUE)

y <- data[, 1]
X <- as.matrix(data[, -1])

# Center and scale (required for RSO)
X <- scale(X, center = TRUE, scale = TRUE)
y <- y - mean(y)

n_total <- nrow(X)
p <- ncol(X)

# ----------------------------------------------------------------------------
# 2. Fixed training size (Wu, 2021)
# ----------------------------------------------------------------------------
set.seed(123)
n_train <- 200
n_test <- n_total - n_train

train_idx <- sample(1:n_total, size = n_train, replace = FALSE)
test_idx <- setdiff(1:n_total, train_idx)

X_train <- X[train_idx, ]
y_train <- y[train_idx]
X_test <- X[test_idx, ]
y_test <- y[test_idx]


# ----------------------------------------------------------------------------
# 3. Define tau grid
# ----------------------------------------------------------------------------
tau_grid <- exp(seq(log(0.01), log(20), length.out = 20))

# ----------------------------------------------------------------------------
# 4. Tau selection function for RSO_0
# ----------------------------------------------------------------------------
select_tau_rso0 <- function(X, y, penalty.factor, tau_grid, eps = 1e-6) {
  n <- length(y)
  bic_values <- numeric(length(tau_grid))
  
  for (i in seq_along(tau_grid)) {
    tau <- tau_grid[i]
    result <- RSO_0(X, y, tau = tau, penalty.factor = penalty.factor)
    df <- sum(abs(result$coef) > eps)
    
    if (df > 0) {
      rss <- sum((y - X %*% result$coef)^2)
      bic_values[i] <- n * log(rss / n) + df * log(n)
    } else {
      bic_values[i] <- Inf
    }
  }
  
  best_idx <- which.min(bic_values)
  return(list(tau = tau_grid[best_idx], bic_values = bic_values))
}

# ----------------------------------------------------------------------------
# 5. Tau selection function for RSOFast
# ----------------------------------------------------------------------------
select_tau_rsofast <- function(X, y, penalty.factor, tau_grid, eps = 1e-6) {
  n <- length(y)
  bic_values <- numeric(length(tau_grid))
  
  for (i in seq_along(tau_grid)) {
    tau <- tau_grid[i]
    result <- RSOFast(X, y, tau = tau, penalty.factor = penalty.factor)
    df <- sum(abs(result$coef) > eps)
    
    if (df > 0) {
      rss <- sum((y - X %*% result$coef)^2)
      bic_values[i] <- n * log(rss / n) + df * log(n)
    } else {
      bic_values[i] <- Inf
    }
  }
  
  best_idx <- which.min(bic_values)
  return(list(tau = tau_grid[best_idx], bic_values = bic_values))
}

# ----------------------------------------------------------------------------
# 6. Solution path function (beta_hat vs tau)
# ----------------------------------------------------------------------------
plot_solution_path <- function(X, y, penalty.factor, tau_grid, method = "RSOFast", 
                                    max_path_points = 50) {
  
  if (length(tau_grid) > max_path_points) {
    tau_path <- tau_grid[seq(1, length(tau_grid), length.out = max_path_points)]
  } else {
    tau_path <- tau_grid
  }
  
  p <- ncol(X)
  beta_path <- matrix(NA, nrow = length(tau_path), ncol = p)
  
  for (i in seq_along(tau_path)) {
    tau <- tau_path[i]
    
    if (method == "RSO_0") {
      result <- RSO_0(X, y, tau = tau, penalty.factor = penalty.factor)
      beta_vec <- as.numeric(result$coef)
    } else {
      result <- RSOFast(X, y, tau = tau, penalty.factor = penalty.factor)
      beta_vec <- as.numeric(result$coef)
    }
    
    if (length(beta_vec) != p) {
      cat(sprintf("Warning: tau = %.4f için beta vektörü alınamadı\n", tau))
      beta_vec <- rep(NA, p)
    }
    
    beta_path[i, ] <- beta_vec
  }
  
  # Convert to data frame for ggplot
  df_path <- data.frame()
  for (j in 1:p) {
    df_path <- rbind(df_path, data.frame(
      tau = tau_path,
      beta = beta_path[, j],
      variable = as.character(j)
    ))
  }
  
  return(df_path)
}

# ----------------------------------------------------------------------------
# 7. Comparison
# ----------------------------------------------------------------------------
penalty.factor <- rep(1, p)

# RSO_0
cat("--- RSO_0 ---\n")
cat("Performing BIC-based tau selection...\n")
time0 <- system.time({
  tau_selection0 <- select_tau_rso0(X_train, y_train, penalty.factor, tau_grid)
  tau0 <- tau_selection0$tau
  cat(sprintf("Selected tau: %.4f\n", tau0))
  
  cat("Fitting final model...\n")
  result0 <- RSO_0(X_train, y_train, tau = tau0, penalty.factor = penalty.factor)
})
cat(sprintf("Total time: %.3f seconds\n\n", time0[3]))

# RSOFast
cat("--- RSOFast ---\n")
cat("Performing BIC-based tau selection...\n")
time_fast <- system.time({
  tau_selection_fast <- select_tau_rsofast(X_train, y_train, penalty.factor, tau_grid)
  tau_fast <- tau_selection_fast$tau
  cat(sprintf("Selected tau: %.4f\n", tau_fast))
  
  cat("Fitting final model...\n")
  result_fast <- RSOFast(X_train, y_train, tau = tau_fast, penalty.factor = penalty.factor)
})
cat(sprintf("Total time: %.3f seconds\n\n", time_fast[3]))

# ----------------------------------------------------------------------------
# 8. Generate solution paths (beta_hat vs tau)
# ----------------------------------------------------------------------------
cat("========================================\n")
cat("Generating Solution Paths...\n")
cat("========================================\n")

tau_path_grid <- exp(seq(log(0.01), log(20), length.out = 20))

cat("Computing RSO_0 solution path...\n")
df_path0 <- plot_solution_path(X_train, y_train, penalty.factor, 
                                    tau_path_grid, method = "RSO_0")

cat("Computing RSOFast solution path...\n")
df_path_fast <- plot_solution_path(X_train, y_train, penalty.factor, 
                                        tau_path_grid, method = "RSOFast")

# ----------------------------------------------------------------------------
# 9. Plot RSO_0 solution path
# ----------------------------------------------------------------------------
p_path0 <- ggplot(df_path0, aes(x = tau, y = beta, group = variable, color = variable)) +
  geom_line(alpha = 0.7) +
  scale_x_log10() +
  labs(x = expression(tau), y = expression(hat(beta)[j]),
       title = "RSO_0 Solution Path") +
  theme_bw() +
  theme(legend.position = "none") +
  geom_vline(xintercept = tau0, linetype = "dashed", color = "red", alpha = 0.7)

# ----------------------------------------------------------------------------
# 10. Plot RSOFast solution path
# ----------------------------------------------------------------------------
p_path_fast <- ggplot(df_path_fast, aes(x = tau, y = beta, group = variable, color = variable)) +
  geom_line(alpha = 0.7) +
  scale_x_log10() +
  labs(x = expression(tau), y = expression(hat(beta)[j]),
       title = "RSOFast Solution Path") +
  theme_bw() +
  theme(legend.position = "none") +
  geom_vline(xintercept = tau_fast, linetype = "dashed", color = "red", alpha = 0.7)

# ----------------------------------------------------------------------------
# 11. Evaluate on test set
# ----------------------------------------------------------------------------
pred0 <- X_test %*% result0$coef
mse0 <- mean((y_test - pred0)^2)
active0 <- sum(abs(result0$coef) != 0)

pred_fast <- X_test %*% result_fast$coef
mse_fast <- mean((y_test - pred_fast)^2)
active_fast <- sum(abs(result_fast$coef) != 0)

# ----------------------------------------------------------------------------
# 12. Plot BIC curves
# ----------------------------------------------------------------------------
bic_df <- data.frame(
  tau = tau_grid,
  RSO_0 = tau_selection0$bic_values,
  RSOFast = tau_selection_fast$bic_values
)

bic_melted <- melt(bic_df, id.vars = "tau", 
                   variable.name = "Method", 
                   value.name = "BIC")

p_bic <- ggplot(bic_melted, aes(x = tau, y = BIC, color = Method)) +
  geom_line() +
  geom_point() +
  geom_vline(data = data.frame(Method = "RSO_0", tau = tau0), 
             aes(xintercept = tau, color = Method), linetype = "dashed") +
  geom_vline(data = data.frame(Method = "RSOFast", tau = tau_fast), 
             aes(xintercept = tau, color = Method), linetype = "dashed") +
  scale_x_log10() +
  labs(x = expression(tau), y = "BIC", 
       title = "BIC Curves for Tau Selection") +
  theme_bw() +
  theme(legend.position = "bottom")

# ----------------------------------------------------------------------------
# 13. Results summary
# ----------------------------------------------------------------------------
cat("\n========================================\n")
cat("RESULTS SUMMARY\n")
cat("========================================\n")

cat("\n--- Computation Time (including tau selection) ---\n")
cat(sprintf("RSO_0:     %.3f seconds\n", time0[3]))
cat(sprintf("RSOFast:   %.3f seconds\n", time_fast[3]))
cat(sprintf("Speedup:   %.1fx\n", time0[3] / time_fast[3]))

cat("\n--- Selected Tau Values ---\n")
cat(sprintf("RSO_0:     %.4f\n", tau0))
cat(sprintf("RSOFast:   %.4f\n", tau_fast))

cat("\n--- Prediction Performance (Test Set, n_test = %d) ---\n", n_test)
cat(sprintf("RSO_0:     MSE = %.6f\n", mse0))
cat(sprintf("RSOFast:   MSE = %.6f\n", mse_fast))
cat(sprintf("Difference: %.2e\n", abs(mse0 - mse_fast)))

cat("\n--- Variable Selection ---\n")
cat(sprintf("RSO_0:     %d active variables selected\n", active0))
cat(sprintf("RSOFast:   %d active variables selected\n", active_fast))

# ----------------------------------------------------------------------------
# 14. Coefficient comparison
# ----------------------------------------------------------------------------
nonzero_idx <- which(abs(result0$coef) != 0 | abs(result_fast$coef) != 0)

if (length(nonzero_idx) > 0) {
  coef_comparison <- data.frame(
    Variable = nonzero_idx,
    RSO_0 = as.numeric(result0$coef[nonzero_idx]),
    RSOFast = as.numeric(result_fast$coef[nonzero_idx])
  )
  
  cat("\n--- Coefficient Comparison (Non-zero Variables) ---\n")
  print(coef_comparison)
  
  cor_coef <- cor(result0$coef[nonzero_idx], result_fast$coef[nonzero_idx])
  cat(sprintf("\nCorrelation between coefficient estimates: %.6f\n", cor_coef))
}

# ----------------------------------------------------------------------------
# 15. Save plots
# ---------------------------------------------------------------------------
ggsave("bic_curves.pdf", p_bic, width = 5, height = 3)
ggsave("solution_path_rso0.pdf", p_path0, width = 5, height = 3)
ggsave("solution_path_rsofast.pdf", p_path_fast, width = 5, height = 3)


# ----------------------------------------------------------------------------
# 16. Display plots
# ----------------------------------------------------------------------------
print(p_bic)
print(p_path0)
print(p_path_fast)

# ----------------------------------------------------------------------------
# 17. Final summary table
# ----------------------------------------------------------------------------
summary_table <- data.frame(
  Metric = c("Selected Tau", "Total Time (s)", "Test MSE", "Active Variables", "Speedup"),
  RSO_0 = c(round(tau0, 4), round(time0[3], 3), round(mse0, 6), active0, NA),
  RSOFast = c(round(tau_fast, 4), round(time_fast[3], 3), round(mse_fast, 6), active_fast, 
              round(time0[3] / time_fast[3], 1))
)

print(summary_table)

write.xlsx(summary_table, "summary_table.xlsx")

