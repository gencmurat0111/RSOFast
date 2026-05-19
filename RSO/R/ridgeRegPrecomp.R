#' Individually Penalized Ridge Regression with Precomputed Matrices
#'
#' Computes the sum of squared residuals (SSR) for individually penalized
#' ridge regression using precomputed X'X and X'y matrices.
#' This is optimized for repeated calls
#' with the same x matrix but different gam parameters (e.g., during grid search).
#'
#' @param x Centered predictor matrix (n x p)
#' @param y Centered response vector (n x 1)
#' @param gam Ridge penalty parameters vector (p x 1) where gam_j = 1/nu_j
#' @param penalty.factor Adaptive penalty factors (p x 1), default = 1 for RSO
#' @param precomp List containing precomputed values:
#'   \itemize{
#'     \item xtx: X'X matrix (p x p)
#'     \item xty: X'y vector (p x 1)
#'   }
#'
#' @return Sum of squared residuals (SSR)
#'
#' @examples
#' \dontrun{
#' n <- 100; p <- 10
#' x <- matrix(rnorm(n*p), n, p)
#' x <- scale(x, scale = FALSE)
#' y <- rnorm(n)
#' y <- y - mean(y)
#' precomp <- list(xtx = crossprod(x), xty = crossprod(x,y))
#' gam <- runif(p)
#' result <- ridgeRegPrecomp(x, y, gam, rep(1, p), precomp)
#' }
#'
#' @export
ridgeRegPrecomp <- function(x, y, gam, penalty.factor = rep(1, length(gam)), precomp) {
  # Input validation
  if (!is.matrix(x)) {
    stop("x must be a matrix")
  }
  if (!is.numeric(y)) {
    stop("y must be numeric")
  }
  if (length(y) != nrow(x)) {
    stop("Number of rows in x must match length of y")
  }
  if (length(gam) != ncol(x)) {
    stop("Length of gam must match number of columns in x")
  }
  if (length(penalty.factor) != ncol(x)) {
    stop("Length of penalty.factor must match number of columns in x")
  }
  if (!is.list(precomp) || !all(c("xtx", "xty") %in% names(precomp))) {
    stop("precomp must be a list with 'xtx' and 'xty' elements")
  }

  ridgeRegPrecomp_arma(x, y, gam, penalty.factor, precomp)
}
