#' Ridge Regression with Degrees of Freedom
#'
#' Computes the sum of squared residuals (SSR), degrees of freedom (trace of hat matrix),
#' and coefficient estimates for ridge regression. Uses Cholesky decomposition for speed
#' and numerical stability.
#'
#' @param x Centered predictor matrix (n x p)
#' @param y Centered response vector (n x 1)
#' @param gam Ridge penalty parameters vector (p x 1) where gam_j = 1/nu_j
#' @param penalty.factor Adaptive penalty factors (p x 1), default = 1 for RSO
#'
#' @return A list containing:
#' \item{ssr}{Sum of squared residuals}
#' \item{df}{Degrees of freedom (trace of hat matrix)}
#' \item{coef}{Coefficient estimates (p x 1)}
#'
#' @examples
#' \dontrun{
#' n <- 100; p <- 10
#' x <- matrix(rnorm(n*p), n, p)
#' x <- scale(x, scale = FALSE)
#' y <- rnorm(n)
#' y <- y - mean(y)
#' gam <- runif(p)
#' result <- ridgereg_df(x, y, gam, rep(1, p))
#' }
#'
#' @export
ridgereg_df <- function(x, y, gam, penalty.factor = rep(1, length(gam))) {
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

  # C++ fonksiyonunu çağır
  ridgeregWdf_fast_arma(x, y, gam, penalty.factor)
}
