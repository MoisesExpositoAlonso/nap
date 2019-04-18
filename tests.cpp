#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <random>
#include <chrono>
#include <ctime>

#include <cstdio>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <list>
#include <iostream>
#include <string>


// when armadillo is loaded, remove this below
//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace std;
////// using namespace arma; // I think this might be causing problems


#include <bigmemory/MatrixAccessor.hpp>
#include <bigmemory/isna.hpp>

// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(bigmemory)]]
// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

////////////////////////////////////////////////////////////////////////////////
/// Profiling utilities
////////////////////////////////////////////////////////////////////////////////

// RcppExport SEXP start_profiler(SEXP str) {
//   ProfilerStart(as<const char*>(str));
//   return R_NilValue;
// }

// RcppExport SEXP stop_profiler() {
//   ProfilerStop();
//   return R_NilValue;
// }


////////////////////////////////////////////////////////////////////////////////
/// Utilities
////////////////////////////////////////////////////////////////////////////////
//
//
// // [[Rcpp::export]]
// arma::mat test(){
//   arma::mat s_chain(100,1000);
//
//   arma::vec news,s;
//   news=Rcpp::rnorm(100,0,0.1);
//   s=news;
//   s_chain.col(0)=s;
//   return(s_chain);
// }

// // [[Rcpp::export]]
// bool worksrunif(double Paccept=0.5){
//   return(Rcpp::runif(1)(0)<Paccept);
// }


// // [[Rcpp::export]]
// arma::vec repC(arma::ivec myrows, int rep ){
//   return Rcpp::rep(myrows,rep);
// }

// // [[Rcpp::export]]
// arma::mat writereadtable(){
//   arma::mat A;
//   A << 0.165300 << 0.454037 << 0.995795 << 0.124098 << 0.047084 << arma::endr
//     << 0.688782 << 0.036549 << 0.552848 << 0.937664 << 0.866401 << arma::endr
//     << 0.348740 << 0.479388 << 0.506228 << 0.145673 << 0.491547 << arma::endr
//     << 0.148678 << 0.682258 << 0.571154 << 0.874724 << 0.444632 << arma::endr
//     << 0.245726 << 0.595218 << 0.409327 << 0.367827 << 0.385736 << arma::endr;
//   A.save("A.txt", arma::raw_ascii);
//   arma::mat B;
//   B.load("A.txt");
//   return(B);
// }
// // [[Rcpp::export]]
// arma::mat loadvec(){
//   arma::mat A;
//   A.load("databig/s_svar05.txt");
//   return(A);
// }


// // [[Rcpp::export]]
// bool accessmatrix(const SEXP A,
//                  const arma::uvec & mycols,
//                  const arma::uvec & myrows) {
//
//   Rcpp::XPtr<BigMatrix> bigMat(A);
//   MatrixAccessor<double> macc(*bigMat);
//
//   for (int j = 0; j <mycols.n_elem; j++) {
//       for (int i = 0; i < myrows.n_elem; i++) {
//         cout << macc[mycols(j)-1][myrows(i)-1] << endl;
//   }}
//   return wrap(true);
// }


// [[Rcpp::export]]
double medianC(const arma::vec& ar, double burnin=0.1){
  arma::vec tosub;
  tosub = arma::linspace((int)round(ar.n_elem * burnin), ar.n_elem-1);
  arma::uvec indices = arma::conv_to<arma::uvec>::from(tosub);
  return(arma::median(ar.elem( indices)));
}

// [[Rcpp::export]]
double tral(const arma::mat& ar, double burnin=0.1){
  return(medianC(ar.row(0), burnin) );
}

// [[Rcpp::export]]
arma::rowvec t2(arma::mat ar, double burnin=0.1){
  return(ar.row(0));
}


// [[Rcpp::export]]
arma::vec medianCmat(const arma::mat& ar,double burnin=0.1){
  if(ar.n_cols<10) burnin=0;
  arma::vec res(ar.n_rows);
  res.fill(0.0);
  for(int i=0; i<ar.n_rows; i++){
    double tmp=medianC(arma::trans(ar.row(i)),burnin);
    res(i) = tmp;
  }
  return(res);
}
