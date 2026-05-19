#' Individually Penalized Ridge Regression
#'
#' Computes the sum of squared residuals (SSR) for individually penalized
#' ridge regression without explicitly forming the hat matrix.
#' Uses Cholesky decomposition for speed and numerical stability.
#'
#' @param x Centered predictor matrix (n x p)
#' @param y Centered response vector (n x 1)
#' @param gam Ridge penalty parameters vector (p x 1) where gam_j = 1/nu_j
#' @param penalty.factor Adaptive penalty factors (p x 1), default = 1 for RSO
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
#' gam <- runif(p)
#' result <- ridgereg(x, y, gam, rep(1, p))
#' }
#'
#' @export
ridgereg <- function(x, y, gam, penalty.factor = rep(1, length(gam))) {
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

  # Call C++ function
  ridgereg_fast_arma(x, y, gam, penalty.factor)
}



