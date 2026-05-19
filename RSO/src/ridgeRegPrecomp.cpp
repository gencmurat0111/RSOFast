#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
double ridgeRegPrecomp_arma(const arma::mat& x,
                             const arma::vec& y,
                             const arma::vec& gam,
                             const arma::vec& penalty_factor,
                             const Rcpp::List& precomp) {

  int n = x.n_rows;
  int p = x.n_cols;

  arma::vec sqrt_gam = arma::sqrt(gam);

  // =========================
  // CASE 1: PRIMAL (n >= p)
  // =========================
  if (n >= p) {

    arma::mat xtx = Rcpp::as<arma::mat>(precomp["xtx"]);
    arma::vec xty = Rcpp::as<arma::vec>(precomp["xty"]);

    arma::mat M = xtx;
    M.each_col() %= sqrt_gam;
    M.each_row() %= sqrt_gam.t();
    M.diag() += penalty_factor;

    arma::vec B = sqrt_gam % xty;

    arma::mat R = arma::chol(M);

    arma::vec v = arma::solve(arma::trimatl(R.t()), B);
    arma::vec beta_tilde = arma::solve(arma::trimatu(R), v);

    arma::vec beta = sqrt_gam % beta_tilde;
    arma::vec Hy = x * beta;

    arma::vec resid = y - Hy;
    double ssr = arma::dot(resid, resid);

    if (ssr < 0 && ssr > -1e-10) ssr = 0;
    return ssr;
  }

  // =========================
  // CASE 2: DUAL (n < p)
  // =========================
  else {

    // w = sqrt(gam) / sqrt(penalty)
    arma::vec w = sqrt_gam / arma::sqrt(penalty_factor);

    // Z = X * diag(w)
    arma::mat Z = x;
    Z.each_row() %= w.t();

    // K = Z Z^T + I
    arma::mat K = Z * Z.t();
    K.diag() += 1.0;

    arma::mat R = arma::chol(K);

    // solve: (K) alpha = y
    arma::vec v = arma::solve(arma::trimatl(R.t()), y);
    arma::vec alpha = arma::solve(arma::trimatu(R), v);

    // beta_tilde = (sqrt(gam)/penalty) % (X^T alpha)
    arma::vec beta_tilde = (sqrt_gam / penalty_factor) % (x.t() * alpha);

    arma::vec beta = sqrt_gam % beta_tilde;

    arma::vec Hy = x * beta;

    arma::vec resid = y - Hy;
    double ssr = arma::dot(resid, resid);

    if (ssr < 0 && ssr > -1e-10) ssr = 0;
    return ssr;
  }
}
