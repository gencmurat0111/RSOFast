#include <RcppArmadillo.h>
#include "RSO_types.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
Rcpp::List RSOFast_arma(const arma::mat& x,
                        const arma::vec& y,
                        double tau,
                        const arma::vec& penalty_factor,
                        Rcpp::Nullable<Rcpp::NumericVector> gaminitNV = R_NilValue) {

  //============================================================================
  // INPUT VALIDATION
  //============================================================================
  arma::uword n = x.n_rows;
  arma::uword p = x.n_cols;

  if (n != y.n_elem) {
    Rcpp::stop("Number of rows in x must match length of y");
  }
  if (p != penalty_factor.n_elem) {
    Rcpp::stop("Number of columns in x must match length of penalty_factor");
  }
  if (tau <= 0) {
    Rcpp::stop("tau must be positive");
  }

  //============================================================================
  // INITIALIZATION
  //============================================================================
  // Precomputed values (for speed)
  arma::mat xtx = x.t() * x;
  arma::vec xty = x.t() * y;

  Rcpp::List precomp = Rcpp::List::create(
    Rcpp::Named("xtx") = xtx,
    Rcpp::Named("xty") = xty
  );

  // Initialize gamma (lambda) values
  arma::vec gam(p);

  if (gaminitNV.isNotNull()) {
    // Scale initial values to sum to tau
    arma::vec gaminit = Rcpp::as<arma::vec>(gaminitNV);
    double sum_init = arma::sum(gaminit);

    if (sum_init <= 0) {
      Rcpp::warning("Sum of initial gam values is <= 0. Using equal distribution.");
      gam.fill(tau / p);
    } else {
      gam = tau * gaminit / sum_init;
    }
  } else {
    // Default: equal distribution
    gam.fill(tau / p);
  }

  //============================================================================
  // ALGORITHM PARAMETERS
  //============================================================================
  const double tol = 0.001 * tau;           // Convergence tolerance
  const double zero_threshold = 2.0 * tol;  // Threshold for variable selection
  const int max_iter = 50;                   // Maximum iterations
  const int golden_iter = 15;                // Golden section iterations
  const double gr = (std::sqrt(5.0) - 1.0) / 2.0;  // Golden ratio

  //============================================================================
  // COORDINATE DESCENT
  //============================================================================
  int iter = 0;
  bool converged = false;
  double current_ssr = ridgeRegPrecomp_arma(x, y, gam, penalty_factor, precomp);

  while (!converged && iter < max_iter) {
    iter++;
    arma::vec old_gam = gam;
    double old_ssr = current_ssr;

    // Loop over each coordinate
    for (arma::uword j = 0; j < p; j++) {

      // Skip if penalty_factor[j] == 0 (unpenalized variable)
      if (penalty_factor(j) == 0.0) continue;

      // Current value and sum of others
      double gam_j = gam(j);
      double other_sum = arma::sum(gam) - gam_j;

      // Skip if other_sum is too small
      if (other_sum <= 1e-10 * tau) continue;

      //----------------------------------------------------------------------
      // Golden section search for optimal t in [0, 1]
      // t = 0: gam(j) = 0, others scaled to sum to tau
      // t = 1: gam(j) = tau, others = 0
      //----------------------------------------------------------------------
      double a = 0.0;
      double b = 1.0;
      double c = b - gr * (b - a);
      double d = a + gr * (b - a);

      arma::vec gam_test(p);

      // Function to evaluate SSR for a given t
      auto eval_ssr = [&](double t) -> double {
        gam_test = gam;

        // Set j-th coordinate to t * tau
        gam_test(j) = t * tau;

        // Scale others proportionally to maintain sum = tau
        if (other_sum > 0) {
          double scale = (tau - gam_test(j)) / other_sum;
          for (arma::uword i = 0; i < p; i++) {
            if (i != j) {
              gam_test(i) = gam(i) * scale;
            }
          }
        }

        return ridgeRegPrecomp_arma(x, y, gam_test, penalty_factor, precomp);
      };

      double fc = eval_ssr(c);
      double fd = eval_ssr(d);

      // Golden section loop
      for (int gi = 0; gi < golden_iter; gi++) {
        if (fc < fd) {
          b = d;
          d = c;
          fd = fc;
          c = b - gr * (b - a);
          fc = eval_ssr(c);
        } else {
          a = c;
          c = d;
          fc = fd;
          d = a + gr * (b - a);
          fd = eval_ssr(d);
        }

        if (std::abs(b - a) < 1e-8) break;
      }

      double best_t = (a + b) / 2.0;
      double best_ssr = std::min(fc, fd);

      // Update if improvement found
      if (best_ssr < current_ssr - 1e-10) {
        // Update gam with optimal value
        gam(j) = best_t * tau;

        // Scale others proportionally
        if (other_sum > 0) {
          double scale = (tau - gam(j)) / other_sum;
          for (arma::uword i = 0; i < p; i++) {
            if (i != j) {
              gam(i) *= scale;
            }
          }
        }

        current_ssr = best_ssr;
      }
    }

    // Check convergence
    double max_change = arma::max(arma::abs(gam - old_gam));
    double ssr_change = std::abs(current_ssr - old_ssr);

    if (max_change < tol || ssr_change < tol * 1e-3) {
      converged = true;
    }
  }

  //============================================================================
  // POST-PROCESSING
  //============================================================================
  // Warning if not converged
  if (!converged) {
    Rcpp::warning("RSOFast did not converge in %d iterations", max_iter);
  }

  // Zero out very small values (variable selection)
  gam.elem(arma::find(gam <= zero_threshold)).zeros();

  // Renormalize to sum to tau (optional, but safe)
  double current_sum = arma::sum(gam);
  if (current_sum > 0) {
    gam *= tau / current_sum;
  }

  //============================================================================
  // FINAL CALCULATIONS
  //============================================================================
  // Compute final coefficients and degrees of freedom
  Rcpp::List temp = ridgeregWdf_fast_arma(x, y, gam, penalty_factor);
  arma::vec coef = Rcpp::as<arma::vec>(temp["coef"]);
  double df = Rcpp::as<double>(temp["df"]);

  // Count selected variables (non-zero gamma)
  arma::uword n_selected = arma::accu(gam > zero_threshold * 0.1);

  //============================================================================
  // RETURN RESULTS
  //============================================================================
  return Rcpp::List::create(
    Rcpp::Named("gamma") = gam,           // Regularization parameters
    Rcpp::Named("lambda") = gam,          // Alias for gamma
    Rcpp::Named("coefficients") = coef,   // Regression coefficients
    Rcpp::Named("coef") = coef,           // Alias for coefficients
    Rcpp::Named("df") = df,               // Degrees of freedom
    Rcpp::Named("ssr") = current_ssr,     // Sum of squared residuals
    Rcpp::Named("iterations") = iter,     // Number of iterations
    Rcpp::Named("converged") = converged, // Convergence status
    Rcpp::Named("n_selected") = n_selected, // Number of selected variables
    Rcpp::Named("tau") = tau               // Original tau parameter
  );
}
