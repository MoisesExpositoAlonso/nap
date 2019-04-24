#include <stdio.h>
#include <math.h>
#include <vector>
#include <list>
#include <iostream>
#include <string>

#define ARMA_64BIT_WORD 1
//// https://stackoverflow.com/questions/40592054/large-matrices-in-rcpparmadillo-via-the-arma-64bit-word-define

// when armadillo is loaded, remove this below
//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace std;
using namespace arma;

// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(BH, bigmemory)]]
// [[Rcpp::plugins(cpp11)]]

#include <bigmemory/MatrixAccessor.hpp>
#include <bigmemory/isna.hpp>


// // [[Rcpp::export]]
// arma::Mat<double> BMsubset(SEXP A, const arma::uvec & myrows, const arma::uvec & mycols ){
//       Rcpp::XPtr<BigMatrix> bigMat(A);
//       arma::Mat<double> X0((double*) bigMat->matrix(), bigMat->nrow(), bigMat->ncol(), false, false);
//                                       // consider saying true, perhaps is faster
//       // Subset matrix
//     	if(myrows.n_elem == X0.n_rows){
//     		X0=X0.cols(mycols);
//     	}else if(mycols.n_elem == X0.n_rows){
//     	  X0=X0.rows(myrows);
//     	}else{
//     		X0=X0.submat(myrows,mycols);
//     	}
//       return(X0);
// }


// [[Rcpp::export]]
arma::vec upperTmat(const arma::mat mymat){
  arma::vec output((mymat.n_cols*(mymat.n_cols-1))/2);
  arma::mat::const_iterator it = mymat.begin() + mymat.n_rows; //iterator at rows to skip in every step, starts at second column
  long toSkipInVec = 0;
  for(int i = 1; i < mymat.n_cols; i++) //Starts with 1 to skip the diagonal
  {
    std::copy(it, it + i, output.begin() + toSkipInVec);
    toSkipInVec += i;
    it += mymat.n_rows;
  }
  return output;
}

// [[Rcpp::export]]
arma::mat Xmcenter(arma::mat X){
  arma::mat newX(X.n_rows,X.n_cols);
  for(int j=0; j<X.n_cols; j++){
   newX.col(j) = (X.col(j) - arma::mean( X.col(j))) ;
  }
  return(newX);
}

// [[Rcpp::export]]
arma::mat Xmvcenter(arma::mat X){
  arma::mat newX(X.n_rows,X.n_cols);
  for(int j=0; j<X.n_cols; j++){
   newX.col(j) = (X.col(j) - arma::mean( X.col(j))) /arma::stddev( X.col(j));
  }
  return(newX);
}

// [[Rcpp::export]]
arma::mat LDrelative(SEXP A, arma::uvec  m, bool debug = false){

  Rcpp::XPtr<BigMatrix> bigMat(A);
  if(bigMat->matrix_type() !=8) stop("Big matrix is not of type double");

  // Read the genome matrix from address
  arma::Mat<double> X((double*) bigMat->matrix(), bigMat->nrow(), bigMat->ncol(), false, false);
  X=X.cols(m);

  // mean and var center for LD calculation
  X=Xmvcenter(X);
  if(debug) cout << X << endl;

  // Get the relative LD for the proposals
  arma::mat R2 =  arma::trans(X)*X ;
  if(debug)  cout << R2 << endl;
  R2 = R2/ arma::sum(upperTmat(R2));
  if(debug)  cout << arma::sum(upperTmat(R2)) << endl;
  return(R2);
}

arma::mat LDrelative(arma::mat X, bool debug = false){

  // mean and var center for LD calculation
  X=Xmvcenter(X);
  if(debug) cout << X << endl;

  // Get the relative LD for the proposals
  arma::mat R2 =  arma::trans(X)*X ;
  if(debug)  cout << R2 << endl;
  R2 = R2/ arma::sum(upperTmat(R2));
  if(debug)  cout << arma::sum(upperTmat(R2)) << endl;
  return(R2);
}
