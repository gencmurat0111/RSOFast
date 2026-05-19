#' Fast RSO (Ridge Selection Operator)
#'
#' Implements the modified coordinate descent algorithm for RSO as described in
#' Wu (2021). Solves for optimal lambda parameters that minimize SSR subject to
#' sum(lambda) = tau and lambda >= 0.
#'
#' @param x Centered predictor matrix (n x p)
#' @param y Centered response vector (n x 1)
#' @param tau Regularization parameter (sum of lambda)
#' @param penalty.factor Adaptive penalty factors (p x 1), default = 1 for RSO
#' @param gaminit Initial lambda values. If NULL, uses tau/p for all.
#'
#' @return A list containing:
#' \item{gamma}{Optimal lambda vector (p x 1)}
#' \item{lambda}{Alias for gamma}
#' \item{coefficients}{Coefficient estimates (p x 1)}
#' \item{coef}{Alias for coefficients}
#' \item{df}{Degrees of freedom}
#' \item{ssr}{Sum of squared residuals}
#' \item{iterations}{Number of iterations}
#' \item{converged}{Convergence status}
#' \item{n_selected}{Number of selected variables}
#' \item{tau}{Original tau parameter}
#'
#' @references
#' Wu, Y. (2021). Can't ridge regression perform variable selection?.
#' Technometrics, 63(2), 263-271. \doi{10.1080/00401706.2020.1791254}
#'
#' @examples
#' \dontrun{
#' n <- 100; p <- 10
#' x <- matrix(rnorm(n*p), n, p)
#' x <- scale(x, scale = FALSE)
#' y <- rnorm(n)
#' y <- y - mean(y)
#' result <- RSOFast(x, y, tau = 1.0, rep(1, p))
#' }
#'
#' @export
RSOFast <- function(x, y, tau, penalty.factor = rep(1, ncol(x)), gaminit = NULL) {

  #============================================================================
  # INPUT VALIDATION
  #============================================================================
  if (!is.matrix(x)) {
    stop("x must be a matrix")
  }
  if (!is.numeric(y)) {
    stop("y must be numeric")
  }
  if (length(y) != nrow(x)) {
    stop("Number of rows in x must match length of y")
  }
  if (length(penalty.factor) != ncol(x)) {
    stop("Length of penalty.factor must match number of columns in x")
  }
  if (!is.numeric(tau) || length(tau) != 1 || tau <= 0) {
    stop("tau must be a single positive number")
  }

  #============================================================================
  # GAMINIT - NULL ise tau/p ile doldur
  #============================================================================
  if (is.null(gaminit)) {
    gaminit <- rep(tau / ncol(x), ncol(x))
  } else {
    if (length(gaminit) != ncol(x)) {
      stop("Length of gaminit must match number of columns in x")
    }
    if (any(gaminit < 0)) {
      stop("gaminit must contain non-negative values")
    }
  }

  #============================================================================
  # CALL C++ FUNCTION
  #============================================================================
  RSOFast_arma(
    x = x,
    y = y,
    tau = tau,
    penalty_factor = penalty.factor,
    gaminit = gaminit  # artık NULL değil
  )
}
