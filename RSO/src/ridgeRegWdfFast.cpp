#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]

Rcpp::List ridgeregWdf_fast_arma(const arma::mat& x,
                                         const arma::vec& y,
                                         const arma::vec& gam,
                                         const arma::vec& penalty_factor) {

  const arma::uword n = x.n_rows;
  const arma::uword p = x.n_cols;

  arma::vec sqrt_gam = arma::sqrt(gam);

  // =========================================================
  // CASE 1: PRIMAL (p <= n)
  // =========================================================
  if (p <= n) {

    arma::mat Xs = x;
    Xs.each_row() %= sqrt_gam.t();

    arma::mat M = Xs.t() * Xs;
    M.diag() += penalty_factor;

    arma::vec Xty = Xs.t() * y;

    arma::mat R = arma::chol(M);

    arma::vec tmp = arma::solve(arma::trimatl(R.t()), Xty);
    arma::vec beta_tilde = arma::solve(arma::trimatu(R), tmp);

    arma::vec coef = sqrt_gam % beta_tilde;

    arma::vec Hy = Xs * beta_tilde;
    arma::vec resid = y - Hy;

    double ssr = arma::dot(resid, resid);
    if (ssr < 0 && ssr > -1e-10) ssr = 0.0;

    // 🔥 FAST + STABLE df
    arma::vec diag_Minv = arma::diagvec(arma::inv_sympd(M));

    double df = (double)p - arma::accu(diag_Minv % penalty_factor);
    if (df < 0 && df > -1e-10) df = 0.0;

    return Rcpp::List::create(
      Rcpp::Named("ssr")  = ssr,
      Rcpp::Named("df")   = df,
      Rcpp::Named("coef") = coef
    );
  }

  // =========================================================
  // CASE 2: DUAL (p > n)
  // =========================================================
  else {

    arma::vec w = sqrt_gam / arma::sqrt(penalty_factor);

    arma::mat Z = x;
    Z.each_row() %= w.t();

    arma::mat K = Z * Z.t();
    K.diag() += 1.0;

    arma::mat R = arma::chol(K);

    arma::vec tmp = arma::solve(arma::trimatl(R.t()), y);
    arma::vec alpha = arma::solve(arma::trimatu(R), tmp);

    arma::vec beta_tilde =
      (sqrt_gam / penalty_factor) % (x.t() * alpha);

    arma::vec coef = sqrt_gam % beta_tilde;

    arma::vec Hy = x * coef;
    arma::vec resid = y - Hy;

    double ssr = arma::dot(resid, resid);
    if (ssr < 0 && ssr > -1e-10) ssr = 0.0;

    // 🔥 EXACT DF (dual, stabil)
    arma::mat V =
      arma::solve(arma::trimatu(R),
                  arma::solve(arma::trimatl(R.t()), Z));

    arma::vec diag_Minv =
      (1.0 / penalty_factor)
      - (arma::sum(Z % V, 0).t() / penalty_factor);

    double df = (double)p - arma::accu(diag_Minv % penalty_factor);
    if (df < 0 && df > -1e-10) df = 0.0;

    return Rcpp::List::create(
      Rcpp::Named("ssr")  = ssr,
      Rcpp::Named("df")   = df,
      Rcpp::Named("coef") = coef
    );
  }
}
