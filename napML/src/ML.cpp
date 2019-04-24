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

#include <fstream>
#include <sstream>



// when armadillo is loaded, remove this below
//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace std;
// using namespace arma;

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
/// Constants
////////////////////////////////////////////////////////////////////////////////

// #define MIN_NUM = std::numeric_limits<float>::min(); // problem is that it does not know the type
const double MIN_NUM = std::numeric_limits<float>::min();
// #define PIraw = 3.14159265358979323846;
// const double PI= 3.14159265358979323846;
// # define PI 3.14159265358979323846  /* pi */
//const double MAX_NUM = std::numeric_limits<float>::max();


////////////////////////////////////////////////////////////////////////////////
/// Utilities with matrices
////////////////////////////////////////////////////////////////////////////////


// [[Rcpp::export]]
bool BMrecode(SEXP A, double & a, const double & b ){
    Rcpp::XPtr<BigMatrix> bigMat(A);
    MatrixAccessor<double> macc(*bigMat);
    int i, j;
    double val=0;
    for (j = 0; j <bigMat->ncol(); j++) {
      for (i = 0; i < bigMat->nrow(); i++) {
        val= (macc[j][i]);
        if(val == a){
        macc[j][i] = b;
       }
      }
    }
    return Rcpp::wrap(true);
}

// [[Rcpp::export]]
bool BMsimulate(SEXP A){
      Rcpp::XPtr<BigMatrix> bigMat(A);
      MatrixAccessor<double> macc(*bigMat);

      NumericVector maf = Rcpp::runif( bigMat->ncol(), 0,0.49);

      int i, j;
      for (j = 0; j <bigMat->ncol(); j++) {
        for (i = 0; i < bigMat->nrow(); i++) {
          if(Rcpp::runif(1)(0) < maf(j)) macc[j][i] = 2;
        }
      }
      return Rcpp::wrap(true);
}

// [[Rcpp::export]]
bool BMwrite012(SEXP A, std::string outfile){
     // Accessor
      Rcpp::XPtr<BigMatrix> bigMat(A);
      MatrixAccessor<double> macc(*bigMat);
      // Open file
      ofstream myfile;
      myfile.open(outfile);

      // header looks like ._A
      std::string head;
      std::string spacer=" ";
      for(int l=0;l< bigMat->ncol(); l++) head.append("._A");
      myfile<<head<<"\n";

      // the matrix
      for (int i = 0; i < bigMat->nrow(); i++) {
        for (int j = 0; j <bigMat->ncol(); j++) {
          myfile<<macc[j][i]<<" ";
        }
        myfile<<"\n";
      }
      return Rcpp::wrap(true);
}


// [[Rcpp::export]]
bool BMwritePED(SEXP A, std::string outfile){
     // Accessor
      Rcpp::XPtr<BigMatrix> bigMat(A);
      MatrixAccessor<double> macc(*bigMat);
      // Open file
      ofstream myfile;
      myfile.open(outfile);

      // the matrix
      double val=0;
      for (int i = 0; i < bigMat->nrow(); i++) {
        for (int j = 0; j <bigMat->ncol(); j++) {
          val= (macc[j][i]);
          if (val==0) myfile<< "C C" <<" ";
          else myfile<< "A A" <<" ";
        }
        myfile<<"\n";
      }
      return Rcpp::wrap(true);
}

// [[Rcpp::export]]
arma::Mat<double> BMsubset(SEXP A,
                           const arma::uvec & myrows,
                           const arma::uvec & mycols ){
      Rcpp::XPtr<BigMatrix> bigMat(A);
      arma::Mat<double> X0((double*) bigMat->matrix(), bigMat->nrow(), bigMat->ncol(), false, false);
                                      // consider saying true, perhaps is faster
      // Subset matrix
    	if(myrows.n_elem == X0.n_rows){
    		X0=X0.cols(mycols-1); // BUG CHECK - position 0 vs 1 in R
    	}else if(mycols.n_elem == X0.n_rows){
    	  X0=X0.rows(myrows-1);// BUG CHECK - position 0 vs 1 in R
    	}else{
    		X0=X0.submat(myrows-1,mycols-1);// BUG CHECK - position 0 vs 1 in R
    	}
      return(X0);
}

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
arma::mat LDrelative(SEXP A, arma::uvec  mycols, bool debug = false){

  Rcpp::XPtr<BigMatrix> bigMat(A);
  if(bigMat->matrix_type() !=8) stop("Big matrix is not of type double");

  // Read the genome matrix from address
  arma::Mat<double> X((double*) bigMat->matrix(), bigMat->nrow(), bigMat->ncol(), false, false);
  X=X.cols(mycols);

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

////////////////////////////////////////////////////////////////////////////////
/// GWA
////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export()]]
arma::vec rankC(arma::vec x) {
  // not a typo, rankC -> order R equivalent rankC(rankC) -> rank R equivalent
    arma::uvec rankarma  = arma::stable_sort_index(x);
    arma::vec z = arma::conv_to<arma::vec>::from(rankarma);
    return z;
}


// [[Rcpp::export]]
List accuracies (const arma::vec & y,
                         const arma::vec & what){
  // linear model
  arma::mat xs(what.n_elem, 2);
  xs.fill(1);
  xs.col(1)=y;
  arma::colvec coef = arma::solve(xs, what);
  double R2= 1- (sum(pow(y-xs*coef,2))/
                 sum(pow(y-arma::mean(y),2))
                 );
  // Pearson r
  arma::mat m1;
  m1.insert_cols(m1.n_cols, y);
  m1.insert_cols(m1.n_cols, what);
  arma::mat r_=arma::cor(m1);
  double r = r_(0,1);
  // Spearman rho
  arma::mat m2;
  m2.insert_cols(m2.n_cols, rankC(rankC(y))); // not a typo, rankC -> order R equivalent rankC(rankC) -> rank R equivalent
  m2.insert_cols(m2.n_cols, rankC(rankC(what)));
  arma::mat rho_=arma::cor(m2);
  double rho = rho_(0,1);

  return List::create(Named("a") = coef(0),
                      Named("b")=coef(1),
                      Named("R2")=R2,
                      Named("rho")=rho,
                      Named("r")=r
                      );
}
// [[Rcpp::export]]
arma::colvec My(const arma::colvec & y, const arma::colvec & h){
  /*
  * Mean trait per genotype
  */
  // Declarations
  arma::colvec hunique = unique(h);
  arma::colvec m(hunique.n_elem);
  // cout << hunique << endl;

  for(int i=0; i< hunique.n_elem; i++ ){
    // Create temporal vector
    arma::colvec ytmp;
    // Fill with all values corresponding to the same genotype
    for(int j=0; j<y.n_elem;j++){
      if(h(j) == hunique(i)) {
        ytmp.resize(ytmp.size()+1);
        ytmp(ytmp.size()-1) = y(j);
      }
    }
    // Compute variance
      if(ytmp.n_elem ==1){
       // v(i)=0;
       	m(i)=ytmp(0);
      }else{
      	m(i)=arma::mean( ytmp );
      }
  }
  return(m);
}


// [[Rcpp::export]]
arma::colvec BMs(const SEXP A,
                 const arma::colvec & y,
                 const arma::uvec & mycols,
                 const arma::uvec & myrows,
                 bool debug=false) {
  Rcpp::XPtr<BigMatrix> bigMat(A);
  MatrixAccessor<double> macc(*bigMat);
  arma::vec coef (mycols.n_elem);

  for (int j = 0; j <mycols.n_elem; j++) {
    int n1=0, n0=0;
    double m1=0, m0=0;
      for (int i = 0; i < myrows.n_elem; i++) {
        if(macc[mycols(j)-1][myrows(i)-1] == 0){
          n0++;
          m0 += y(myrows(i)-1);
        }else{
          n1++;
          m1 += y(myrows(i)-1);
        }
      }
      coef(j) = (m1/n1) - (m0/n0);
      if(isna(coef(j))) coef(j)=0;
  }
  double sds = sqrt(arma::var(coef));
  for (int j = 0; j <mycols.n_elem; j++) if(abs(coef(j))>sds*3) coef(j)=0;
  for (int j = 0; j <mycols.n_elem; j++) if(coef(j) < -1) coef(j)=-0.999;
  return coef;
}


// [[Rcpp::export]]
arma::colvec BMmgwa( arma::mat X0,
                     arma::colvec & yraw,
                    bool debug=false) {
      /////////////////////////
 	// Centering
  arma::mat X=Xmcenter(X0);
  arma::vec y = (yraw/arma::mean(yraw) ) - 1;
  // Ordineary Least Squares
  arma::colvec coef(X.n_cols);
	arma::colvec val;
	// cout << "calculate effects" << endl;
  for(int i=0; i< X.n_cols; i++){
		arma::mat Xsub= X.col(i);
		arma::vec val = solve(Xsub,arma::mat(y));
		if(debug) cout << val << endl;
		coef(i) = val(0);
	}
 return(coef);
}


////////////////////////////////////////////////////////////////////////////////
/// Likelihood, Probabilities, Proposals
////////////////////////////////////////////////////////////////////////////////

////// PRIOR
// double uniform(const arma::colvec & s, const double min, const double max){ // uniform probability
//   double L= 0;
//   int N=s.n_elem;
//   for(int i=0;i<N ;i++){
//     L+= R::dunif(s(i),min,max,true);
//   }
//   return L;
// }
// double loggaussian(const arma::colvec & s, double svar){ //
//   arma::vec x = log(1+s);
//   // double L = -.5*n*log(2*PI) -.5*n*log(svar) -(1/(2*svar))*sum(arma::pow((x-mean),2));
//   double L=0;
//   for(int i = 0; i<s.n_elem; i++){
//     L+= R::dnorm(x(i),0,svar,true);
//   }
//   return L;
// }
// double logmixgaussian(const arma::colvec & s, double svar, double ss){
//     arma::vec x = log(1+s);
//     double L=0;
//     for(int i = 0; i<s.n_elem; i++){
//       if(x(i)==0){
//         L += log(ss  + (1-ss) * R::dnorm(x(i),0,svar,false)) ;
//       }else{
//         L += (1-ss) * R::dnorm(x(i),0,svar,true);  // L+= R::dnorm(x(i),0,svar,true);
//       }
//     }
//     return L;
//   }
///// prior

// //[[Rcpp::export]]
// double logmixgaussian(const double & s, double svar, double ss){
//     double x = log(1+s);
//     double L;
//       if(s==0){
//         L= (ss  + (1-ss) * R::dnorm(s,0,svar,false)) ;
//       }else{
//         L= ((1-ss) * R::dnorm(s,0,svar,false));  // L+= R::dnorm(x(i),0,svar,true);
//       }
//     return log(L);
// }

/////  LIKELIHOOD
// [[Rcpp::export]]
double LLGaussMix(const double & y, const double & w, const double & v, const double & p){
  double LL;
  if(y==0){
    LL = p  + (1-p) *  R::pnorm(0,w,v,true,false) ;
  }else{
    LL = (1-p) * R::dnorm(y,w,v,false);
  }
return log(LL);
}

// [[Rcpp::export]]
double LIKELIHOOD(const arma::vec & y, // worked in previous iteration
                  const arma::vec & w,
                  const double & b, const double & a,
                  const double & p, const double & mu,
                  const double & epi,
                  bool verbose=false,
                  bool printall=false){
    double L=0;
    double LL;
      for(int i=0; i< y.n_elem ; i ++){ // check for infinity
        LL= LLGaussMix(y(i)/mu,w(i),w(i)*b+a,p);
        if(!std::isnan(LL)) L += LL; // if there is one nan, all sum is nan
      }
  if(verbose) cout<< L << endl;
  if(std::isinf(L)) L=MIN_NUM;
  return(L);
}


////////////////////////////////////////////////////////////////////////////////
/// Fitness
////////////////////////////////////////////////////////////////////////////////

   // fitness mode
  double wmode(const double &s , const double &x, const int & FITmode){
    double w_;
    switch(FITmode){
      case 1:  // additive
        w_=  (s * x) ;
        break;
      case 2:  // multiplicative
        w_= 1 + (s * x) ;
        break;
      case 3:  // inverse multiplicative
        w_= pow((1 + s),x);
        break;
    }
    return(w_ );
  }
  // add fitness contribution of one SNP to overall fitness
  void wupdate(double &prevw, const double &s,const double &x, const int & FITmode) {     // void operator +*=(double w, double s,int x,int mode)  // operator not needed unless class
    switch(FITmode){
      case 1:  // additive
        prevw= prevw + wmode(s,x, FITmode);
        break;
      default:  // non-additive
        prevw= prevw * wmode(s,x, FITmode);
        break;
    }
  }
  // remove fitness contribution of one SNP to overall fitness (necessary to avoid repeating computations)
  void unwupdate(double &prevw, const double &s,const double &x, const int & FITmode) {
    switch(FITmode){
      case 1 :  // additive
        prevw= prevw - wmode(s,x, FITmode);
        break;
      default:  // non-additive
        prevw= prevw / wmode(s,x, FITmode);
        break;
    }
  }

// [[Rcpp::export]]
arma::vec wC(
               const arma::Mat<double> X,
               const arma::vec & s,
               const int & mode,
               double epi=1,
               bool verbose=false){
    arma::vec w(X.n_rows); // w= new arma::vec(X.n_rows)
    w.fill(1); // need to fill otherwise unstable floats
    for (int i = 0; i < X.n_cols; i ++) {
        for(int j=0; j < X.n_rows ; j++){
          wupdate(w(j),s(i),X(j,i), mode );
        }
    }
  return(pow(w,epi));
}

// [[Rcpp::export]]
arma::vec wCBM(SEXP A, // worked in previous iteration
                   const arma::vec & s,
                   const arma::uvec & mycols,
                   const arma::uvec & myrows,
                   const int & mode,
                   double epi=1,
                   bool verbose=false){

  Rcpp::XPtr<BigMatrix> bigMat(A);
  MatrixAccessor<double> macc(*bigMat);
  arma::vec w(myrows.n_elem); // w= new arma::vec(X.n_rows)
    w.fill(1); // need to fill otherwise unstable floats
  int i, j;
  double val=0;
  for (j = 0; j <mycols.n_elem; j++) {
    for (i = 0; i < myrows.n_elem; i++) {
      val= (macc[mycols(j)-1][myrows(i)-1] );
      // val= (macc[j][i] );
      wupdate(w(i),s(j), val, mode );
    }
  }
   for (i = 0; i < myrows.n_elem; i++)  if(w(i)<0) w(i)=0; // CHECK BUGS
  return(pow(w,epi));
}



// [[Rcpp::export]]
void wCBMupdate(arma::vec wnew,
                SEXP A,
                const arma::uvec & mycols,
                const arma::uvec & myrows,
                const int & indecito,
                const double & s0,
                const double & s1,
                const int & mode){
        Rcpp::XPtr<BigMatrix> bigMat(A);
        MatrixAccessor<double> macc(*bigMat);
        double val=0;
        for(int j=0; j < myrows.n_elem ; j++){
          val= (macc[indecito][myrows(j)-1]) ;
          unwupdate(wnew(j),s0, val, mode);
          wupdate(wnew(j),s1,val, mode);
         }
  }
// [[Rcpp::export]]
arma::vec wCBMupdateforR(arma::vec wnew,
                SEXP A,
                const arma::uvec & mycols,
                const arma::uvec & myrows,
                const int & indecito,
                const double & s0,
                const double & s1,
                const int & mode){
  Rcpp::XPtr<BigMatrix> bigMat(A);
        MatrixAccessor<double> macc(*bigMat);
        double val=0;
        for(int j=0; j < myrows.n_elem ; j++){
          val= (macc[indecito][myrows(j)-1]) ;
          unwupdate(wnew(j),s0, val, mode);
          wupdate(wnew(j),s1,val, mode);
         }
  return(wnew);
}
// [[Rcpp::export]]
arma::vec sampleWC(
                    const arma::vec & w,
                    const double & a,
                    const double & b,
                    const double & p,
                    const int & rep
                    ){
  arma::vec y(w.n_elem * rep);
  y.fill(0);
  double val=0;
  int count=0;
  for(int i=0; i<w.n_elem ; i++){
    for(int j=0; j<rep;j++){
      val = Rcpp::rnorm(1,w(i), abs(a+(w(i)*b)))(0);
      if(val<0) val=0;
      if(Rcpp::runif(1)(0) < p ) val=0;
      y[i] =val;
      count++;
    }
  }
  return(y);
}

// [[Rcpp::export]]
bool ssaveC(arma::vec s, std::string path){
  s.save(path,arma::raw_ascii);
  return(wrap(true));
}


// [[Rcpp::export]]
arma::vec ssimC(arma::uvec snps,
                double svar){
  arma::vec news;
  news.load("databig/s_svar01.txt");
  news= exp( log(1+news) * (svar/0.01) )-1 ;

  arma::vec subnews(snps.n_elem,arma::fill::zeros);
  for( int i=0; i<snps.n_elem;i++){
    subnews(i) = news(snps(i)-1);
  }
  // if(svar==0.05){
  //   news.load("databig/s_svar05.txt");
  // }else if(svar==0.25){
  //   news.load("databig/s_svar25.txt");
  // }else if(svar==0.5){
  //   news.load("databig/s_svar50.txt");
  // } else{
  //   cout << "S values not found for accuracy calculations!" << endl;
  //   news=exp(Rcpp::rnorm(nsnp,0,svar))-1;
  // }
  return subnews;
}
