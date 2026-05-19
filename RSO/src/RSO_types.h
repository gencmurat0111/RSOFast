#ifndef RSO_TYPES_H
#define RSO_TYPES_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Fonksiyon prototipleri - varsayılan değer YOK
double ridgeRegPrecomp_arma(const arma::mat& x,
                       const arma::vec& y,
                       const arma::vec& gam,
                       const arma::vec& penalty_factor,
                       const Rcpp::List& precomp);

// İSİM DÜZELTİLDİ: ridgeregWdf_fast_arma ile aynı olmalı
Rcpp::List ridgeregWdf_fast_arma(const arma::mat& x,
                                 const arma::vec& y,
                                 const arma::vec& gam,
                                 const arma::vec& penalty_factor);

Rcpp::List RSO_2(const arma::mat& x,
                   const arma::vec& y,
                   double tau,
                   const arma::vec& penalty_factor,
                   Rcpp::Nullable<Rcpp::NumericVector> gaminitNV);  // varsayılan değer YOK

#endif
