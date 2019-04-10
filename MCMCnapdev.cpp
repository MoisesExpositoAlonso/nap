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



#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

void printProgress (double percentage)
{
  int val = (int) (percentage * 100);
  int lpad = (int) (percentage * PBWIDTH);
  int rpad = PBWIDTH - lpad;
  printf ("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
  fflush (stdout);
}


template <typename T>
Rcpp::NumericVector arma2vec(const T& x) {
  return Rcpp::NumericVector(x.begin(), x.end());
}

template <typename T>
arma::vec vec2arma(const T& x) {
  return  Rcpp::as<arma::vec>(x);
}

// [[Rcpp::export]]
NumericVector sample_num( NumericVector x,
                          int size,
                          bool replace,
                          NumericVector prob = NumericVector::create()
)
{
  NumericVector ret = RcppArmadillo::sample(x, size, replace, prob) ;
  return ret ;
}


// #define MIN_NUM = std::numeric_limits<float>::min(); // problem is that it does not know the type
// #define PIraw = 3.14159265358979323846;
// const double PI= 3.14159265358979323846;
// # define PI 3.14159265358979323846  /* pi */
//const double MAX_NUM = std::numeric_limits<float>::max();

const double MIN_NUM = std::numeric_limits<float>::min();

////////////////////////////////////////////////////////////////////////////////
/// Matrix and LD utilities
////////////////////////////////////////////////////////////////////////////////


// [[Rcpp::export]]
arma::Mat<double> BMsubset(SEXP A, const arma::uvec & myrows, const arma::uvec & mycols ){
      Rcpp::XPtr<BigMatrix> bigMat(A);
      arma::Mat<double> X0((double*) bigMat->matrix(), bigMat->nrow(), bigMat->ncol(), false, false);
                                      // consider saying true, perhaps is faster
      // Subset matrix
      if(myrows.n_elem == X0.n_rows){
        X0=X0.cols(mycols);
      }else if(mycols.n_elem == X0.n_rows){
        X0=X0.rows(myrows);
      }else{
        X0=X0.submat(myrows,mycols);
      }
      return(X0);
}


// [[Rcpp::export]]
arma::vec upperTmat(const arma::mat mat){
  arma::vec output((mat.n_cols*(mat.n_cols-1))/2);
  arma::mat::const_iterator it = mat.begin() + mat.n_rows; //iterator at rows to skip in every step, starts at second column
  long toSkipInVec = 0;
  for(int i = 1; i < mat.n_cols; i++) //Starts with 1 to skip the diagonal
  {
    std::copy(it, it + i, output.begin() + toSkipInVec);
    toSkipInVec += i;
    it += mat.n_rows;
  }
  return output;
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



////////////////////////////////////////////////////////////////////////////////
/// Likelihood, Probabilities, Proposals
////////////////////////////////////////////////////////////////////////////////


////// Utils proposals //////

// [[Rcpp::export]]
double runif_reflect(double minhere,double maxhere,double min,double max){
  // int counter=1;
  double newval;
  if( min == max){
    newval= min; // if values want to be kept constant
  }else{
    newval =Rcpp::runif(1,minhere,maxhere)(0);
    if(newval<min){newval = (min-newval) + min;}
    else if(newval>max){newval = max- (newval-max);}
  }
  if(newval < min || newval>max){   // Check if it is out of bounds
      newval=(max-min)/2;
  }
  return(newval);
}

////////////////////////////////////////////////////////////////////////////////
////// Hyperparameter Proposals //////

class GPROPOSAL{
  private:
    double b; double bmin; double bmax;
    double a; double amin; double amax;
    double p; double pmin=0; double pmax=1;
    double mu; double mumin; double mumax;
    double epi; double epimin; double epimax;
    double svar; double svarmin; double svarmax;
    double ss; double ssmin; double ssmax;
    int nupdates=1;
    bool verbose=false;
    double bw=0.1;
  public:
    GPROPOSAL(
            double b_=0.5,double bmin_=0,double bmax_=1,
            double a_=0.1,double amin_=0,double amax_=1,
            double p_=0.5,
            double mu_=1,double mumin_=0, double mumax_=50,
            double epi_=1,double epimin_=0, double epimax_=5,
            double svar_=0.5,double svarmin_=0, double svarmax_=5,
            double ss_=0.1,double ssmin_=0, double ssmax_=1
            ){
      b=b_;bmin=bmin_;bmax=bmax_;
      a=a_;amin=amin_;amax=amax_;
      p=p_;
      mu=mu_;mumin=mumin_;mumax=mumax_;
      epi=epi_;epimin=epimin_;epimax=epimax_;
      svar=svar_;svarmin=svarmin_;svarmax=svarmax_;
      ss=ss_;ssmin=ssmin_;ssmax=ssmax_;
    }
    void setbandwidth(double bandwidth){bw=bandwidth;}
    void setupdatesnum(int ups){nupdates=ups;}
    void setverbose(bool verbose_){verbose=verbose_;}
    void printatributes(){
      cout <<"bw = " << bw << endl;
      cout <<"b = " << b << " [" << bmin << " " << bmax << "]" << endl;
      cout <<"a = " << a << " [" << amin << " " << amax << "]" << endl;
      cout <<"p = " << p << " [" << pmin << " " << pmax << "]" << endl;
      cout <<"mu = " << mu << " [" << mumin << " " << mumax << "]" << endl;
      cout <<"epi = " << epi << " [" << epimin << " " << epimax << "]" << endl;
      cout <<"svar = " << svar << " [" << svarmin << " " << svarmax << "]" << endl;
      cout <<"ss = " << ss << " [" << ssmin << " " << ssmax << "]" << endl;
    }
    arma::vec fn(arma::vec g){
      // New proposal
      arma::vec news=g;
      // Update one position
      double minhere,maxhere;
      double newval;
      if(verbose) cout << "Loop to substitute position" << endl;
      for(int j=0; j< nupdates; j++){
        int randomIndex = rand() % g.size();
        switch(randomIndex){
            if(verbose) cout << g(randomIndex) << endl;
            case 0:
              minhere=g(randomIndex)- (bw *(bmax-bmin)) ;
              maxhere=g(randomIndex)+ (bw *(bmax-bmin));
              newval= runif_reflect(minhere,maxhere,bmin,bmax); break;
            case 1:
              minhere=g(randomIndex)- (bw *(amax-amin)) ;
              maxhere=g(randomIndex)+ (bw *(amax-amin));
              newval= runif_reflect(minhere,maxhere,amin,amax); break;
            case 2:
              minhere=g(randomIndex)- (bw *(pmax-pmin)) ;
              maxhere=g(randomIndex)+ (bw *(pmax-pmin));
              newval= runif_reflect(minhere,maxhere,pmin,pmax); break;
            case 3:
              minhere=g(randomIndex)- (bw *(mumax-mumin)) ;
              maxhere=g(randomIndex)+ (bw *(mumax-mumin));
              newval= runif_reflect(minhere,maxhere,mumin,mumax); break;
            case 4:
              minhere=g(randomIndex)- (bw *(epimax-epimin)) ;
              maxhere=g(randomIndex)+ (bw *(epimax-epimin));
              newval= runif_reflect(minhere,maxhere,epimin,epimax); break;
            case 5:
              minhere=g(randomIndex)- (bw *(svarmax-svarmin)) ;
              maxhere=g(randomIndex)+ (bw *(svarmax-svarmin));
              newval= runif_reflect(minhere,maxhere,svarmin,svarmax); break;
            case 6:
              minhere=g(randomIndex)- (bw *(ssmax-ssmin)) ;
              maxhere=g(randomIndex)+ (bw *(ssmax-ssmin));
              newval= runif_reflect(minhere,maxhere,ssmin,ssmax); break;
        }
        if(verbose) cout << newval << endl;
        news(randomIndex) = newval;
      }
      if(verbose) cout << "End loop" << endl;
      return(news);
  }
};


// [[Rcpp::export]]
arma::mat test_GPROPOSAL(double b=1,
                          double a=1,
                          double p=1,
                          double mu=1,
                          double epi=1,
                          double svar=1,
                          double ss=0.1,
                          double bw=0.5,
                          int iter=1000,
                          bool verbose = false,
                          int nupdates=7
                          ){

  arma::vec g(7);
  g(0)= b;
  g(1)= a;
  g(2)= p;
  g(3)= mu;
  g(4)= epi;
  g(5)= svar;
  g(6)= ss;
  GPROPOSAL GProp;

  GProp.setupdatesnum(nupdates);
  GProp.setbandwidth(bw);
  GProp.setverbose(verbose);
  GProp.printatributes();

  arma::mat res(7,iter);
  res.col(0)=g;

  cout << "Runing several proposals "<< endl;
  for(int i=1; i<iter; i++){
    g=GProp.fn(g);
    res.col(i)=g;
  }
  return(res);
}


////////////////////////////////////////////////////////////////////////////////
////// Hyperparameter Proposals //////

class SPROPOSAL {
  private:
    arma::vec s;
    int nsnp;
    double min;
    double max;
    int mode;
    bool verbose;
    double svar;
    double bw=0.5;
  public:
  SPROPOSAL(
            int mode_=1,
            int nsnp_=100,
            double svar_= 0.1,
            double min_=-0.98039215686274505668,
            double max_=50,
            bool verbose_= false
            ){
   nsnp=nsnp_;min=min_;max=max_;mode=mode_;verbose=verbose_;svar=svar_;
   s=start(nsnp_, mode_,svar_,min_,max_);
  }
  arma::vec start(int nsnp,int mode,double svar,double min, double max){
    // arma::vec news(nsnp);
    arma::vec news(nsnp);
    switch(mode){
    case 1:
      news=Rcpp::runif(nsnp,min,max);
      break;
    case 2:
      news=Rcpp::rnorm(nsnp,0,svar);
      break;
    }
    return(news);
  }
  void printatributes(){
    cout <<"bw = " << bw << endl;
    cout <<"min = " << min << endl;
    cout <<"max = " <<max << endl;
    cout <<"mode = " <<mode << endl;
    cout <<"verbose = " <<verbose << endl;
  }
  void setS(arma::vec s_){s=s_;}
  void setbandwidth(double bandwidth){bw=bandwidth;}
  arma::vec getS(){return(s);}
  void update(int mode=2, double svar=0.1){
    switch(mode){
      case 1:
        return flat();
        break;
      case 2:
        return logn(svar);
        break;
      default:
        return logn(svar);
        break;
    }
  }
   void flat(){
      arma::vec news=s;
      double minhere,maxhere,newval;
      int randomIndex = rand() % s.size();
      minhere=s(randomIndex)-bw;
      maxhere=s(randomIndex)+bw;
      newval = runif_reflect(minhere,maxhere,min,max);
      news(randomIndex) = newval;
      // return(news);
      s=news;
   }
   void logn(double svar){
      arma::vec news= s;
      double meanhere,newval;
      int randomIndex = rand() % s.size();
      meanhere=log(1+s(randomIndex)); //*** transform mean to log dimension
      if(std::isinf(meanhere)) meanhere= 0; // check for infinity
      newval = Rcpp::rnorm(1,meanhere,svar)(0); // sample with the mean being the last value.
      news(randomIndex) = exp(newval)-1; //*** back transform to natural dimension
      // return(news);
      s=news;
  }
};

  //  arma::vec flat(int nproposals){
  //     arma::vec news=s;
  //     double minhere,maxhere,newval;
  //     if(verbose) cout << "Loop to substitute position" << endl;
  //     for(int j=0; j< nproposals; j++){
  //       int randomIndex = rand() % s.size();
  //       minhere=s(randomIndex)-bw;
  //       maxhere=s(randomIndex)+bw;
  //       newval = runif_reflect(minhere,maxhere,min,max);
  //       news(randomIndex) = newval;
  //     }
  //     if(verbose) cout << "End loop" << endl;
  //     return(news);
  // }
  //   arma::vec logn(int nproposals,double svar){
  //     arma::vec news= s;
  //     double meanhere,newval;
  //     if(verbose) cout << "Loop to substitute position" << endl;
  //     for(int j=0; j< nproposals; j++){
  //       int randomIndex = rand() % s.size();
  //       meanhere=log(1+s(randomIndex)); //*** transform mean to log dimension
  //       // newval = Rcpp::rnorm(1,1,svar)(0); // this will always draw from the same distribution. The chain is not reversible anymore and the aceptance ratio has to be updated
  //       //*// if(std::isinf(meanhere)) meanhere= (max-min)/2; // check for infinity
  //       if(std::isinf(meanhere)) meanhere= 0; // check for infinity
  //       newval = Rcpp::rnorm(1,meanhere,svar)(0); // sample with the mean being the last value.
  //       //*// if(newval<log(1+min)) newval= log(1+min)+ (newval-log(1+min)); // refract in limits
  //       //*// if(newval>log(1+max)) newval= log(1+max)- (newval-log(1+max)); // refract in limits
  //       // NOTE: if any of the lines //*// the chain gets biased towards lower values
  //       news(randomIndex) = exp(newval)-1; //*** back transform to natural dimension
  //       if(verbose) cout << newval << endl;
  //       if(verbose) cout << news(randomIndex) << endl;
  //     }
  //     if(verbose) cout << "End loop" << endl;
  //     return(news);
  // }

// [[Rcpp::export]]
arma::vec PropoS(int nsnps=100,
                 int mode=2,
                 double svar = 0.1){
  SPROPOSAL ps(mode,nsnps,svar);
  return(ps.getS());
}

// [[Rcpp::export]]
arma::vec test_UPDATES(
                 int nsnps=100,
                 int nupdates=1,
                 int mode=2,
                 double svar = 0.1){
  SPROPOSAL ps(mode,nsnps,svar);
  arma::rowvec s=ps.getS();
  arma::mat res(s.n_rows,nupdates, arma::fill::zeros);
  cout << res.col(0).n_elem << endl;
  cout << res.row(0).n_elem << endl;
  for(int j=0; j<nupdates; j++){
    ps.update();
    arma::rowvec snew=ps.getS();
    res.col(0)=snew;
  }
  return(res);
}

//   NumericMatrix res(nsnps,nupdates+1);
//   res(_,0) = arma2vec(ps.getS());
//   for (int i=0; i<nupdates+1;i++){
//     ps.update();
//     res(_,i+1) = arma2vec(ps.getS());
//   }
//   // return List::create(Named("s_0") = s_0,
//   //                      Named("s_1") = s_1
//   //                        );
//   return(res);
// }

// arma::vec test_UPDATES(
//                  int nsnps=100,
//                  int nupdates=1,
//                  int mode=2,
//                  double svar = 0.1){
//   SPROPOSAL ps(mode,nsnps,svar);
//   arma::vec s_0;
//   arma::vec s_1;
//   s_0 = ps.getS();
//   for (int i=0; i<nupdates;i++){
//     ps.update();
//     cout << ps.getS() << endl;
//     s_1 = ps.getS();
//   }
//   return List::create(Named("s_0") = s_0,
//                        Named("s_1") = s_1
//                          );
// }





// // [[Rcpp::export]]
// arma::vec PropoS( // I think this function might be unnecessary
//                  int nupdates=1,
//                  double svar = 0.5){
//   arma::vec s(nupdates);
//   // s.fill(0);
//   SPROPOSAL ps(nupdates,svar);
//   return(ps.fn(s,svar));
// }
//
// // [[Rcpp::export]]
// arma::mat test_SPROPOSAL(int nupdates=1,
//                          double svar =0.5,
//                          int m=3,
//                          int iter=3000,
//                          int mode=3
//                          ){
//
//   // initialize class
//   arma::vec s(m);
//   s.fill(0);
//   SPROPOSAL PS(nupdates,mode);
//
//   // create result matrix
//   arma::mat res(m,iter);
//   res.col(0)=s;
//
//   cout << "Runing several proposals "<< endl;
//   for(int i=1; i<iter; i++){
//     s=PS.fn(s,svar);
//     res.col(i)=s;
//   }
//
//   return(res);
// }
