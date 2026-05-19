#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
double ridgereg_fast_arma(const arma::mat& x,
                          const arma::vec& y,
                          const arma::vec& gam,
                          const arma::vec& penalty_factor) {

  int n = x.n_rows;
  int p = x.n_cols;

  arma::vec sqrt_gam = arma::sqrt(gam);

  // X scaled olustur (bir kere)
  arma::mat X_scaled(n, p);
  for (int j = 0; j < p; ++j) {
    X_scaled.col(j) = x.col(j) * sqrt_gam(j);
  }

  arma::mat M = X_scaled.t() * X_scaled;
  M.diag() += penalty_factor;

  arma::vec Xty = X_scaled.t() * y;

  // Cholesky
  arma::mat R;
  try {
    R = arma::chol(M);
  } catch (const std::runtime_error& e) {
    Rcpp::stop("Cholesky decomposition failed. Matrix might not be positive definite.");
  }

  arma::vec beta_tilde = arma::solve(arma::trimatu(R),
                                     arma::solve(arma::trimatl(R.t()), Xty));

  // H'y = X_scaled * beta_tilde
  arma::vec Hy = X_scaled * beta_tilde;

  // SSR
  arma::vec resid = y - Hy;
  double ssr = arma::sum(resid % resid);

  // Numerik stabilite
  if (ssr < 0 && ssr > -1e-10) ssr = 0;

  return ssr;
}
