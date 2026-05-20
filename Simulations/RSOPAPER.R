#' Ridge Regression for RSO
#'
#' Computes the sum of squared residuals (SSR) for ridge regression with individual penalty factors.
#' Used internally by the RSO algorithm.
#'
#' @param x Centered predictor matrix (n x p)
#' @param y Centered response vector (n x 1)
#' @param gam Ridge penalty parameters vector (p x 1) where gam_j = 1/nu_j
#' @param penalty.factor Adaptive penalty factors (p x 1), default = 1 for RSO
#'
#' @return Sum of squared residuals (SSR)
#'
#' @export
ridgereg_0 <- function(x, y, gam, penalty.factor) {
  H <- x %*%
    diag(sqrt(gam)) %*%
    solve(diag(sqrt(gam)) %*%
            t(x) %*%
            x %*%
            diag(sqrt(gam)) + diag(penalty.factor)) %*%
    diag(sqrt(gam)) %*%
    t(x)
  return(sum((y - H %*% y)^2))
}

#' Ridge Regression with Diagnostics for RSO
#'
#' Computes the sum of squared residuals (SSR), degrees of freedom (trace of hat matrix),
#' and coefficient estimates for ridge regression with individual penalty factors.
#'
#' @param x Centered predictor matrix (n x p)
#' @param y Centered response vector (n x 1)
#' @param gam Ridge penalty parameters vector (p x 1) where gam_j = 1/nu_j
#' @param penalty.factor Adaptive penalty factors (p x 1), default = 1 for RSO
#'
#' @return A list containing:
#' \item{ssr}{Sum of squared residuals}
#' \item{df}{Degrees of freedom (trace of hat matrix)}
#' \item{coef}{Coefficient estimates}
#'
#' @export
ridgeregWdf_0 <- function(x, y, gam, penalty.factor) {
  H <- x %*%
    diag(sqrt(gam)) %*%
    solve(diag(sqrt(gam)) %*%
            t(x) %*%
            x %*%
            diag(sqrt(gam)) + diag(penalty.factor)) %*%
    diag(sqrt(gam)) %*%
    t(x)
  coef <- diag(sqrt(gam)) %*%
    solve(diag(sqrt(gam)) %*%
            t(x) %*%
            x %*%
            diag(sqrt(gam)) + diag(penalty.factor)) %*%
    diag(sqrt(gam)) %*%
    t(x) %*%
    y
  return(list(ssr = sum((y - H %*% y)^2), df = sum(diag(H)), coef = coef))
}

#' BIC-based Tau Selection for RSO
#'
#' Selects the optimal regularization parameter tau for RSO (Ridge Selection Operator)
#' using Bayesian Information Criterion (BIC). Tests a grid of tau values and returns
#' the one that minimizes BIC.
#'
#' @param x Centered predictor matrix (n x p)
#' @param y Centered response vector (n x 1)
#' @param penalty.factor Adaptive penalty factors (p x 1), default = 1 for RSO
#' @param tau_grid Numeric vector of tau values to test. If NULL, automatically generates a grid.
#' @param gaminit Initial lambda values (optional)
#' @param eps Threshold for considering coefficients as zero (default = 1e-6)
#' @param max_iter Maximum iterations for RSO (default = 20, not used directly)
#'
#' @return Optimal tau value that minimizes BIC
#'
#' @export
#' @examples
#' \dontrun{
#' n <- 100; p <- 10
#' x <- matrix(rnorm(n*p), n, p)
#' x <- scale(x, scale = FALSE)
#' y <- rnorm(n)
#' y <- y - mean(y)
#' best_tau <- RSO_tau_bic(x, y, rep(1, p))
#' }
RSO_tau_bic_0 <- function(x, y, penalty.factor = rep(1, ncol(x)),
                        tau_grid = NULL, gaminit = NULL,
                        eps = 1e-6, max_iter = 20) {

  p <- ncol(x)
  n <- length(y)

  # Create default tau grid if not provided
  if (is.null(tau_grid)) {
    tau_max <- sum(penalty.factor) * 0.5
    tau_grid <- seq(0.001, tau_max, length = 50)
  }

  bic_values <- numeric(length(tau_grid))

  for (i in seq_along(tau_grid)) {
    tau <- tau_grid[i]

    # Run RSO
    temp <- RSO(x, y, tau = tau, penalty.factor = penalty.factor, gaminit = gaminit)

    # Degrees of freedom (number of non-zero coefficients)
    df <- sum(abs(temp$coef) > eps)

    # Residual sum of squares
    rss <- sum((y - x %*% temp$coef)^2)

    # Calculate BIC
    bic_values[i] <- n * log(rss / n) + df * log(n)
  }

  # Find tau with minimum BIC
  best_tau <- tau_grid[which.min(bic_values)]

  return(best_tau)
}

#' Ridge Selection Operator (RSO)
#'
#' Implements the modified coordinate descent algorithm for RSO as described in
#' Wu (2021). Solves for optimal lambda parameters that minimize SSR subject to
#' sum(lambda) = tau and lambda >= 0.
#'
#' @param x Centered predictor matrix (n x p)
#' @param y Centered response vector (n x 1)
#' @param tau Regularization parameter (sum of lambda)
#' @param penalty.factor Adaptive penalty factors (p x 1), default = 1 for RSO
#' @param gaminit Initial lambda values (optional). If NULL, uses tau/p for all.
#'
#' @return A list containing:
#' \item{lamhat}{Optimal lambda vector (p x 1)}
#' \item{coef}{Coefficient estimates (p x 1)}
#'
#' @export
#' @examples
#' \dontrun{
#' n <- 100; p <- 10
#' x <- matrix(rnorm(n*p), n, p)
#' x <- scale(x, scale = FALSE)
#' y <- rnorm(n)
#' y <- y - mean(y)
#' result <- RSO(x, y, tau = 1.0, penalty.factor = rep(1, p))
#' }
RSO_0 <- function(x, y, tau, penalty.factor, gaminit = NULL) {

  ### predictors x: n x p, column centered
  ### response y: n x 1, centered

  p <- ncol(x)

  # Initialize gamma (lambda) values
  if (is.null(gaminit)) {
    gamcur <- rep(tau / p, p)
  } else {
    gamcur <- tau * gaminit / sum(gaminit)
  }

  kkk <- 20
  mygrid <- (0:kkk) / kkk
  test <- 1
  count <- 1

  while (test) {

    oldgam <- gamcur

    for (j in 1:p) {
      gamcandj <- rep(0, p)
      gamcandj[j] <- 1

      gamcandMj <- gamcur
      gamcandMj[j] <- 0

      if (sum(gamcandMj) > 0.01 * tau) {
        gamcandMj <- gamcandMj / sum(gamcandMj)

        fffun <- function(ttt) {
          return(ridgereg(x, y, (gamcandj * ttt + gamcandMj * (1 - ttt)) * tau, penalty.factor))
        }
        tttmin <- optimize(fffun, c(0, 1), tol = 0.0001)
        ttt <- tttmin$minimum
        gamcur <- (gamcandj * ttt + gamcandMj * (1 - ttt)) * tau
      } # end if
    }  # end for over j

    count <- count + 1
    if (max(abs(gamcur - oldgam)) < 0.001 * tau) {
      test <- 0
    }
    if (count > 20) {
      test <- 0
      print("Takes more than 20 loops to converge in modified coordinate descent!!!!")
    }
  } # end while over test

  # Set very small values to zero (variable selection)
  gamcur[which(gamcur <= 2 * 0.001 * tau)] <- 0

  lamhat <- gamcur

  # Compute final results
  temp <- ridgeregWdf_0(x, y, lamhat, penalty.factor = penalty.factor)

  return(list(lamhat = lamhat, coef = temp$coef))
}
