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
using namespace arma;

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
/// Utilities with matrices
////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
arma::Mat<double> BMsubset(SEXP A, const arma::uvec & myrows, const arma::uvec & mycols ){
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

// [[Rcpp::export]]
double modeC(const arma::vec& ar){
  return(arma::median(ar));
}

// [[Rcpp::export]]
arma::vec modeCmat(const arma::mat& ar){
  arma::vec res(ar.n_rows);
  res.fill(0.0);
  for(int i=0; i<ar.n_rows; i++){
    double tmp=arma::median(ar.row(i));
    res(i) = tmp;
  }
  return(res);
}

////////////////////////////////////////////////////////////////////////////////
/// GWA
////////////////////////////////////////////////////////////////////////////////


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

// [[Rcpp::export]]
arma::colvec BMcgwa(arma::mat X0,const arma::vec & yraw) {
   	// Centering
    arma::mat X=Xmcenter(X0);
    arma::vec y = (yraw/arma::mean(yraw) ) - 1;
    // Ordineary Least Squares
	arma::colvec coef = solve(X,y);
 return(coef);
}

// [[Rcpp::export]]
arma::vec BMridge(arma::mat X0,const arma::colvec & yraw,
             const double & lambda=1) {
 	// Centering
  arma::mat X=Xmcenter(X0);
  arma::vec y = (yraw/arma::mean(yraw) ) - 1;
  // Precompute some values
	int n = X.n_rows;
	int p = X.n_cols;

	arma::vec coef(p, fill::zeros);
	colvec y_tilde = join_cols(y, arma::zeros<colvec>(p));

		mat X_tilde = join_cols(X, sqrt(lambda) * eye(p, p));
		mat Q, R;
		qr_econ(Q, R, X_tilde);
		coef = solve(R, Q.t() * y_tilde);

  return coef;
	// return List::create(Named("coef") = coef,
	//           					Named("lambda") = lambda);
}
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
  NumericVector ret = sample(x, size, replace, prob) ;
  return ret ;
}


// #define MIN_NUM = std::numeric_limits<float>::min(); // problem is that it does not know the type
const double MIN_NUM = std::numeric_limits<float>::min();
// #define PIraw = 3.14159265358979323846;
// const double PI= 3.14159265358979323846;
// # define PI 3.14159265358979323846  /* pi */
//const double MAX_NUM = std::numeric_limits<float>::max();


////////////////////////////////////////////////////////////////////////////////
/// Matrix and LD utilities
////////////////////////////////////////////////////////////////////////////////


///// Utilities get locations of genotypes in genome matrix from IDs attached to fitness
// [[Rcpp::export]]
arma::vec hsub(const arma::vec & h){
  arma::vec hunique = unique(h);
  arma::vec hpos(h.n_elem);
  for(int i=0; i<h.n_elem;i++){
    for(int j=0; j< hunique.n_elem;j++){
      if(h(i) == hunique(j)) hpos(i) = j;
    }
  }
  return(hpos);
}


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


///// Find which selection coefficient is different in vector
// [[Rcpp::export]]
int whichchanged(const arma::vec & s1, const arma::vec & s2){
  int index= -9;
  int count=0;
  for(int i=0;i<s1.size();i++){
    if(s1(i) != s2(i)){ index=i; count++;}
  }
  if(count>1){
    cout<<"Counter of changes = "<<count<<endl;
    stop("More than one parameter updated at a time!"); // also using this for genotypes
  }
  return(index);
}

////////////////////////////////////////////////////////////////////////////////
/// Likelihood, Probabilities, Proposals
////////////////////////////////////////////////////////////////////////////////


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

//[[Rcpp::export]]
double logmixgaussian(const double & s, double svar, double ss){
    double x = log(1+s);
    double L;
      if(s==0){
        L= (ss  + (1-ss) * R::dnorm(s,0,svar,false)) ;
      }else{
        L= ((1-ss) * R::dnorm(s,0,svar,false));  // L+= R::dnorm(x(i),0,svar,true);
      }
    return log(L);
}

//[[Rcpp::export]]
double iterprior(const arma::colvec & s, const double & svar, const double & ss){
    arma::vec x = log(1+s);
    double L=0;
    for(int i = 0; i<s.n_elem; i++){
        L += logmixgaussian(x(i),svar,ss);
    }
    return L;
}

//[[Rcpp::export]]
double PRIOR (const arma::colvec & s,
              double const & par1=0.1,
              double const & par2=0.1,
              const int & type=2){
  switch(type){
    case 1: // moc mode, return 1
      return 0.0;
      break;
    case 2: // true
      return iterprior(s,par1,par2);
      break;
  }
}

//[[Rcpp::export]]
double ePRIOR(double e, double const & m=1, double const & v=0.2){
  // return(R::dlnorm(e,m,v,true));
  return(R::dnorm(e,m,v,true));
}

// update probabilities
//[[Rcpp::export]]
double uPRIOR ( double pri,
                const arma::colvec & s0,
                const arma::colvec & s1,
                const double & par10=0.1,
                const double & par11=0.1,
                const double & par20=0.1,
                const double & par21=0.1,
                const int & type=2){
  switch(type){
  case 1: // moc mode, return 1
    return 0.0;
    break;
  case 2: // true
    int position = whichchanged(s0,s1);
    if(position != -9){
      double L0, L1;
      L0 = logmixgaussian(s0(position),par10,par20);
      L1 = logmixgaussian(s1(position),par11,par21);
      pri-=L0;
      pri+=L1;
    }
    return(pri);
    break;
  }
}

/////  likelihood
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
double LIKELIHOOD(const arma::vec & y,
                  const arma::vec & hs,
                  const arma::vec & w,
                  const double & b, const double & a,
                  const double & p, const double & mu,
                  const double & epi,
                  bool verbose=false,int LIKmode=2){
  double L=0;
  switch(LIKmode){
  case 1: // moc likelihood
    break;
  case 2:
    arma::vec w_;
    if(epi==1){w_=pow(w,epi);}
    else{w_=w;}

    double LL;
      for(int i=0; i< y.n_elem ; i ++){ // check for infinity
        LL= LLGaussMix(y(i)/mu,w_(hs(i)),w_(hs(i))*b+a,p);

        // if(verbose and std::isinf(LL)){
        if(std::isinf(LL)){
          cout << "---" << endl;
          cout << i << endl;
          cout << y(i) << " "<< w_(hs(i)) << " "<< w_(hs(i))*b+a <<" "<< p << endl;
          cout << LL << endl;
        }
        L += LL;
      }
    break;
  };
  // cout<< L << endl; /////// CHECKING LIKELIHOOD!
  return(L);
}


// double uLIKELIHOOD(double LL, // probably no work around, the majority of genotypes will be affected by any change in parameters
//                   const arma::vec & y,const arma::vec & hs,
//                   const arma::vec & w0,const arma::vec & w1,
//                   const double & b0, const double & b1,
//                   const double & a0, const double & a1,
//                   const double & p0, const double & p1,
//                   const double & mu0, const double & mu1,
//                   const double & epi0, const double & epi1,
//                   bool verbose=false){
//   LLL=LIKELIHOOD(y,hs,w1,b1,a1,p1,mu1,epi1,verbose);
//   return (LLL);
// }

// double uLIKELIHOOD(double LL, // probably no work around, the majority of genotypes will be affected by any change in parameters
//                   const arma::vec & y,const arma::vec & hs,
//                   const arma::vec & w0,const arma::vec & w1,
//                   const double & b0, const double & b1,
//                   const double & a0, const double & a1,
//                   const double & p0, const double & p1,
//                   const double & mu0, const double & mu1,
//                   const double & epi0, const double & epi1,
//                   bool verbose=false){
//   // ind out whether it was s or hyperparameters that changed
//   if(epi0 != epi1){
//     LL=LIKELIHOOD(y,hs,w1,b1,a1,p1,mu1,epi1)
//   }else if(b0 != b1 ||
//            a0 != a1 ||
//            p0 != p1 ||
//            mu0 != mu1){
//
//   }
//   int i = whichchanged(w0,w1);
//   double L0, L1;
//   L0 = LLGaussMix(y(i)/mu0,w0(hs(i)),w0(hs(i))*b0+a0,p0);
//   L1 = LLGaussMix(y(i)/mu1,w1(hs(i)),w1(hs(i))*b1+a1,p1);
//   LL-=L0;
//   LL+=L1;
//       if(verbose and std::isinf(LL)){ // check for infinity
//         cout << "---" << endl;
//         cout << i << endl;
//         cout << y(i) << " "<< w1(hs(i)) << " "<< w1(hs(i))*b1+a1 <<" "<< p1 << endl;
//         cout << LL << endl;
//       }
//   return(LL);
//   }

////////////////////////////////////////////////////////////////////////////////
/// Fitness
////////////////////////////////////////////////////////////////////////////////

   // fitness mode
  double wmode(const double &s , const int &x, const int & FITmode){
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
  void wupdate(double &prevw, const double &s,const int &x, const int & FITmode) {     // void operator +*=(double w, double s,int x,int mode)  // operator not needed unless class
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
  void unwupdate(double &prevw, const double &s,const int &x, const int & FITmode) {
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
               double mu=1,
               bool verbose=false){
    arma::vec w(X.n_rows); // w= new arma::vec(X.n_rows)
    w.fill(1); // need to fill otherwise unstable floats
    for (int i = 0; i < X.n_cols; i ++) {
        for(int j=0; j < X.n_rows ; j++){
          wupdate(w(j),s(i),X(j,i), mode );
        }
    }
  return(pow(w,epi)*mu);
}
////////////////////////////////////////////////////////////////////////////////
/// MCMC
////////////////////////////////////////////////////////////////////////////////


class NAP{
  private:
    // hyperparameters
    double b; double bmin; double bmax;
    double a; double amin; double amax;
    double p; double pmin; double pmax;
    double mu; double mumin; double mumax;
    double epi; double epimin; double epimax;
    double svar; double svarmin; double svarmax;
    double ss; double ssmin; double ssmax;
    int npar=7; // update if new parameters are added

    // s parameters and fitness
    arma::vec s; int nsnp=100; double smin; double smax;
    arma::vec y;
    arma::vec hs;
    arma::vec w;
    int nind=515;
    arma::mat X;

    // misc
    int nupdates=1; // 0 for all
    bool verbose=false;
    bool test=false;
    double updateratio= -9.0;
    double bw=0.001; /// BANDWIDT OF CHANGES IS REALLY IMPORTANT HARDCODED PARAM.
    int iterations=1000;
    double burnin=0.1;
    int iter=0;
    int Smode; int PRImode; int LIKmode ;int FITmode;

    // chains
    arma::mat s_chain;
    arma::mat par_chain;
    arma::mat w_chain;
    arma::vec prob_chain;arma::vec pri_chain;arma::vec lik_chain;
    arma::vec ok_chain;

    SEXP A;
  public:
    NAP(
            const arma::vec & y_,
            const arma::vec & h,
            const SEXP & A_, // instead of arma::mat X,
            const arma::uvec & m, // the positions of SNPs // check start 0 vs 1
            const arma::uvec & n , // the positions of individuals
            arma::vec & s_,

            double b_=0.5,double bmin_=0,double bmax_=1,
            double a_=0.1,double amin_=0,double amax_=1,
            double p_=0.5,double pmin_=0,double pmax_=1,
            double mu_=1,double mumin_=0, double mumax_=50,
            double epi_=1,double epimin_=0.5, double epimax_=2,
            double svar_=0.1,double svarmin_=0, double svarmax_=5,
            double ss_=0.1,double ssmin_=0, double ssmax_=1,
            double smin_=-0.98039215686274505668, double smax_=50,

            int iterations_=1000, double burnin_=0.1,
            int Smode_=2,int LIKmode_=3, int PRImode_=2, int FITmode_=2,
            bool verbose_=false,
            bool test_=false,
            double updateratio_ = -9.0,
            double bw_=0.01
            ){
      ////////////////
      // Initialize //
      b=b_;bmin=bmin_;bmax=bmax_;
      a=a_;amin=amin_;amax=amax_;
      p=p_;pmin=pmin_;pmax=pmax_;
      mu=mu_;mumin=mumin_;mumax=mumax_;
      epi=epi_;epimin=epimin_;epimax=epimax_;
      svar=svar_;svarmin=svarmin_;svarmax=svarmax_;
      ss=ss_;ssmin=ssmin_;ssmax=ssmax_;
      smin=smin_;smax=smax_;

      A=A_;
      nsnp=m.n_elem;
      nind=n.n_elem;
      s=s_;

      y=y_;
      // hs=hsub(h); // BUG CHECK
      hs=h-1; // BUG CHECK

      iterations=iterations_;burnin=burnin_;
      Smode=Smode_;LIKmode=LIKmode_;PRImode=PRImode_;FITmode=FITmode_;
      verbose=verbose_; test=test_; updateratio=updateratio_; bw=bw_;

      /////////////
      // Chains //
      prob_chain.zeros(iterations_+1);
      pri_chain.zeros(iterations_+1);
      lik_chain.zeros(iterations_+1);
      ok_chain.zeros(iterations_+1);
      w_chain.zeros(nind,iterations_+1);
      s_chain.zeros(nsnp,iterations_+1);
      par_chain.zeros(npar, iterations_+1);
      s_chain.col(0)=s;
      par_chain(0,0)=b;// actually this does not matter, pas are selected randomly
      par_chain(1,0)=a;
      par_chain(2,0)=p;
      par_chain(3,0)=mu;
      par_chain(4,0)=epi;
      par_chain(5,0)=svar;
      par_chain(6,0)=ss; //

      /////////////////////////
     // GENOME MATRIX //
      X.zeros(n.n_elem,m.n_elem);
      if(TYPEOF(A) == EXTPTRSXP){
      cout<< "Reading external pointer and subsetting genome matrix ... "<<endl;
         X=BMsubset(A,n,m);
      }else if(TYPEOF(A) == REALSXP){
      cout<< "Matrix provided already subsetted from R"<<endl;
        // NumericMatrix Xr(A);
        // cout << " nrow= " << Xr.nrow() << " ncol= " << Xr.ncol() << endl;
        // arma::mat X(Xr.begin(), Xr.nrow(), Xr.ncol(), false); // this was the original
         // X(Xr.begin(), Xr.nrow(), Xr.ncol(), false); // without the arma::mat at the beginning does not work
         // arma::mat X=A; // no viable converstion
        NumericMatrix Xr(A);
        X = as<arma::mat>( Xr ) ;
      }
    };

    ////////////////////////////////////////////////////////////////////////////
    // functions
    ////////////////////////////////////////////////////////////////////////////
    void setverbose(bool verb){verbose=verb;}; // set verbose
    ////////////////////
    // update selection
     void Sstart(){
       if(verbose) cout << "starting S" <<endl;
        arma::vec news;
        switch(Smode){
        case 1:
          news=Rcpp::runif(nsnp,smin,smax);
          break;
        case 2:
          news=exp(Rcpp::rnorm(nsnp,0,svar))-1;
          break;
        }
       s=news;
       s_chain.col(0)=s;
      }
     void Sstartgwa(){
        arma::vec news;
        arma::vec ymeans=My(y,hs);
        news=BMridge(X,ymeans,1); ///////////// BUG CAREFULL GWA
        for(int j;j<news.n_elem;j++){
          if(news(j)< -1) news(j) = -0.99;
        }
        s=news;
        s_chain.col(0)=s;
     }
      void Scopy(){s_chain.col(iter)=s_chain.col(iter-1);}
      void Supdate(){
      s=s_chain.col(iter-1);
      if(verbose) cout << "updating S" <<endl;
        switch(Smode){
          case 1:
            Sflat();
            break;
          case 2:
            Slogn(svar);
            break;
        }
        s_chain.col(iter)=s; // update S
      }
     void Sflat(){
        double minhere,maxhere,newval;
        int randomIndex = rand() % s.size();
        minhere=s(randomIndex)-bw;
        maxhere=s(randomIndex)+bw;
        newval = runif_reflect(minhere,maxhere,smin,smax);
        if(newval< -1) newval=-0.999;
        s(randomIndex) = newval;
       if(verbose) cout << "New S:" << newval << endl;
     }
     void Slogn(double svar){
        double meanhere,newval;
        int randomIndex = rand() % s.size();
        meanhere=log(1+s(randomIndex)); //*** transform mean to log dimension
        if(std::isinf(meanhere)) meanhere= 0; // check for infinity
        newval = Rcpp::rnorm(1,meanhere,svar)(0); // sample with the mean being the last value.
        if(newval< -1) newval=-0.999;
        s(randomIndex) = exp(newval)-1; //*** back transform to natural dimension
       if(verbose) cout << "New S:" << newval << endl;
    }
    ////////////////////
    // update parameters
    void Gstart(){
      if(verbose) cout << "starting hyperparameters" <<endl;
      par_chain(0,0)=Rcpp::runif(1,bmin,bmax)(0);
      par_chain(1,0)=Rcpp::runif(1,amin,amax)(0);
      par_chain(2,0)=Rcpp::runif(1,pmin,pmax)(0);
      par_chain(3,0)=Rcpp::runif(1,mumin,mumax)(0);
      par_chain(4,0)=Rcpp::runif(1,epimin,epimax)(0);
      par_chain(5,0)=Rcpp::runif(1,svarmin,svarmax)(0);
      par_chain(6,0)=Rcpp::runif(1,ssmin,ssmax)(0);
    };
    void Gcopy(){par_chain.col(iter)=par_chain.col(iter-1);}

    void Gupdate(){
      if(verbose) cout << "updating hyperparameters" <<endl;
      double minhere,maxhere, newval;
      // next s
      par_chain.col(iter)=par_chain.col(iter-1);
      for(int j=0; j< 1; j++){ // only 1 updatepossible
        int randomIndex = rand() % npar;
        switch(randomIndex){
            if(verbose) cout << par_chain(randomIndex,iter) << endl;
            case 0:
              minhere=par_chain(randomIndex,iter)- (bw *(bmax-bmin)) ;
              maxhere=par_chain(randomIndex,iter)+ (bw *(bmax-bmin));
              newval= runif_reflect(minhere,maxhere,bmin,bmax); break;
              if(newval<0) newval=0;
            case 1:
              minhere=par_chain(randomIndex,iter)- (bw *(amax-amin)) ;
              maxhere=par_chain(randomIndex,iter)+ (bw *(amax-amin));
              newval= runif_reflect(minhere,maxhere,amin,amax); break;
              if(newval<0) newval=0;
            case 2:
              minhere=par_chain(randomIndex,iter)- (bw *(pmax-pmin)) ;
              maxhere=par_chain(randomIndex,iter)+ (bw *(pmax-pmin));
              newval= runif_reflect(minhere,maxhere,pmin,pmax); break;
              if(newval<0) newval=0;
              if(newval>1) newval=1;
            case 3:
              minhere=par_chain(randomIndex,iter)- (bw *(mumax-mumin)) ;
              maxhere=par_chain(randomIndex,iter)+ (bw *(mumax-mumin));
              newval= runif_reflect(minhere,maxhere,mumin,mumax); break;
              if(newval<0) newval=0;
            case 4:
              minhere=par_chain(randomIndex,iter)- (bw *(epimax-epimin)) ;
              maxhere=par_chain(randomIndex,iter)+ (bw *(epimax-epimin));
              newval= runif_reflect(minhere,maxhere,epimin,epimax); break;
              if(newval<0) newval=0;
            case 5:
              minhere=par_chain(randomIndex,iter)- (bw *(svarmax-svarmin)) ;
              maxhere=par_chain(randomIndex,iter)+ (bw *(svarmax-svarmin));
              newval= runif_reflect(minhere,maxhere,svarmin,svarmax); break;
              if(newval<0) newval=0;
            case 6:
              minhere=par_chain(randomIndex,iter)- (bw *(ssmax-ssmin)) ;
              maxhere=par_chain(randomIndex,iter)+ (bw *(ssmax-ssmin));
              newval= runif_reflect(minhere,maxhere,ssmin,ssmax); break;
              if(newval<0) newval=0;
        }
        if(verbose) cout << "New hyperpar:" << newval << endl;
        par_chain(randomIndex,iter) = newval;
      }
      // if(verbose) cout << "End loop" << endl;
  }
  ////////////////////
  // calculate fitness
  void FITstart(){
    if(verbose) cout << "start computing fitness from scratch" <<endl;
    arma::vec w=wC(X,s_chain.col(0),FITmode,par_chain(4,0),par_chain(3,0)); // matrix, s, mode, epi, mu
    // arma::vec w(nind,1.0); // w= new arma::vec(X.n_rows)
    // for(int j=0; j < X.n_rows ; j++){
    // for (int i = 0; i < X.n_cols; i ++) {
    //       wupdate(w(j),s_chain(i,iter),X(j,i), par_chain(4,iter) ); // epi in position 4
    //     }
    // } // THIS WAS CAUSING A BUG
    w_chain.col(iter)=w;
  }
  void FITcopy(){w_chain.col(iter)=w_chain.col(iter-1);}
  void FITupdate(){
    if(verbose) cout << "updating fitness" <<endl;
    // copy to only modify
    arma::vec wnew = w_chain.col(iter-1);
    // find which changed
    int indecito=whichchanged(s_chain.col(iter),s_chain.col(iter-1));
    // update selection
    if(indecito != -9){ // only update if s and s2 are not identical
      for(int j=0; j < X.n_rows ; j++){
        unwupdate(wnew(j),s_chain(indecito,iter-1),X(j,indecito), FITmode);
       }
      for(int j=0; j < X.n_rows ; j++){
        wupdate(wnew(j),s_chain(indecito,iter),X(j,indecito), FITmode);
       }
    }
    // update epistasis
    if(par_chain(4,iter-1) != par_chain(4,iter)){
      wnew = pow(wnew,1/par_chain(4,iter-1));
      wnew = pow(wnew,par_chain(4,iter));
    }
      w_chain.col(iter)=wnew;
  }
  ////////////////////
  // calculate probabilities
  void PROBstart(){
    if(verbose) cout << "calculating probability from scratch" <<endl;
        double newpri, newlik;
        newpri=PRIOR(s_chain.col(iter),
                      par_chain(5,iter), // svar in position 5
                      par_chain(6,iter), // ss in position 6
                      PRImode); //+ ePRIOR(par_chain(4,iter)); // epi in position 4
        newlik=LIKELIHOOD(y,hs,
                                     w_chain.col(iter),
                                     par_chain(0,iter), // b in position 0
                                     par_chain(1,iter), // a in position 1
                                     par_chain(2,iter), // p in position 2
                                     par_chain(3,iter), // mu in position 3
                                     par_chain(4,iter), // epi in position 3
                                     verbose,LIKmode);
        pri_chain(iter)=newpri;
        lik_chain(iter)=newlik;
        prob_chain(iter)=newpri + newlik;
        if(verbose) cout << "Start Prior: " << newpri << endl;
        if(verbose) cout << "Start Likelihood: " << newlik << endl;
        if(verbose) cout << "Start Posterior: " << newlik + newpri << endl;
  }
  void PROBupdate(){
    if(verbose) cout << "updating probability" <<endl;
        double newpri, newlik;
        newpri=uPRIOR(pri_chain(iter-1),
                       s_chain.col(iter-1),s_chain.col(iter),
                       par_chain(5,iter-1),par_chain(5,iter), // svar in position 5
                       par_chain(6,iter-1),par_chain(6,iter), // ss in position 6
                       PRImode) ;//+ ePRIOR(par_chain(4,iter)); // epi in position 4
        newlik=LIKELIHOOD(y,hs,
                                     w_chain.col(iter),
                                     par_chain(0,iter), // b in position 0
                                     par_chain(1,iter), // a in position 1
                                     par_chain(2,iter), // p in position 2
                                     par_chain(3,iter), // mu in position 3
                                     par_chain(4,iter), // epi in position 3
                                     verbose,LIKmode);
        pri_chain(iter)=newpri;
        lik_chain(iter)=newlik;
        prob_chain(iter)=newpri + newlik;
}

  ////////////////////
  void additerations() {iter ++;}
  // decide in chain
  void decide(){
    // Ratio of probs
    double Paccept= exp( prob_chain(iter) - prob_chain(iter-1));
    bool accept = Rcpp::runif(1)(0)<Paccept;
    if(verbose) cout<< "Accept proposal " << accept <<endl;
    if(accept){
      ok_chain(iter)=1;
    }else{
      ok_chain(iter)=0;
      s_chain.col(iter) = s_chain.col(iter-1);
      par_chain.col(iter) = par_chain.col(iter-1);
      prob_chain(iter) = prob_chain(iter-1);
      pri_chain(iter) = pri_chain(iter-1);
      lik_chain(iter) = lik_chain(iter-1);
    }
  }
    ////////////////////
  List report(){
    arma::vec shat=modeCmat(s_chain);
    arma::vec par=modeCmat(par_chain);
    arma::vec what=wC(X,shat,FITmode,par(4),par(3)); // matrix, s, mode, epi, mu
    return List::create(Named("chain") = s_chain.t(),
                        Named("shat")=shat,
                        Named("parchain") = par_chain.t(),
                        Named("par")=par,
                        Named("w")=what,
                        Named("posterior") = prob_chain,
                        Named("prior") = pri_chain,
                        Named("likelihood") = lik_chain,
                        Named("accept") = ok_chain,
                        Named("FITmode") = FITmode);
  }

  ////////////////////
  // MAIN
  ////////////////////
    void RUNTEST(bool verbose=true){
      int iter=0;
      // Sstart();
      // Gstart(); // BUG CHECK!
      FITstart();
      PROBstart();


      if(verbose) cout << "w: " << w_chain.col(0) << endl;
      if(verbose) cout << "s: " << s_chain.col(0) << endl;
      if(verbose) cout << "b: " << par_chain(0,0) << endl;
      if(verbose) cout << "a: " << par_chain(1,0) << endl;
      if(verbose) cout << "p: " << par_chain(2,0) << endl;
      if(verbose) cout << "mu: " << par_chain(3,0) << endl;
      if(verbose) cout << "epi: " << par_chain(4,0) << endl;
      if(verbose) cout << "svar: " << par_chain(5,0) << endl;
      if(verbose) cout << "ss: " << par_chain(6,0) << endl;
      if(verbose) cout << "Prior: " << pri_chain(0) << endl;
      if(verbose) cout << "Likelihood: " << lik_chain(0) << endl;
      if(verbose) cout << "Posterior: " << prob_chain(0) << endl;
    }
    ///////////////////////////////////////////////////////////////////////////
    void MAIN(){
    ///////////////////////////////////////////////////////////////////////////
      int iter=0;
      cout<< endl;
      cout<< "Arguments:"<<endl;
      cout<< "# iterations = "<< iterations <<endl;
      cout<< "----------"<<endl;
      cout<< "Total number of individual's observations = "<<  y.n_elem <<endl;
      cout<< "Total number of SNPs = "<<  s.n_elem <<endl;
      cout<< "----------"<<endl;
      cout<< "TEST run = "<< test <<endl;
      cout<< "Verbose = "<< verbose <<endl;
      cout<< "----------"<<endl;
      cout<< "Initializing & start chrono ... "<<endl;
      std::chrono::time_point<std::chrono::system_clock> start, end; // start chronometer values
      start = std::chrono::system_clock::now();
      cout<< endl;

    /// Start NAP & handle -INF probability
    ///////////////////////////////////////////////////////////////////////////
    // Sstart();
    // Sstartgwa(); // initialize s based on gwa
    // Gstart(); // BUG careful
    FITstart();
    PROBstart();
    int attemptcounter=0;
    int maxattempts=1000;
    while(std::isinf(prob_chain(0)) ||std::isnan(prob_chain(0)) ||attemptcounter==0){
      if(verbose && attemptcounter==0) cout << "Running proposals and first probability" << endl;
      if(attemptcounter==maxattempts) stop("Attempted initializing 1000 times!!! ... tune starting parameters!");
      if(attemptcounter==1) cout << "Posterior is infinite!!!. Attempting new start..." << endl<< endl;

      // Sstart();
    // Sstartgwa(); // initialize s based on gwa
    // Gstart(); // BUG careful
      FITstart();
      PROBstart();
      attemptcounter++;
    }
    if(verbose) cout << "Successful MCMC start after " << attemptcounter << " attempts " << endl<< endl;


    /// Run MCMC
    ///////////////////////////////////////////////////////////////////////////
    if(verbose) cout << "----------"<< endl;
    if(verbose) cout << "*** Starting iterations ***" << endl;

    if(updateratio== -9.0) updateratio = (double)s.n_elem/((double)s.n_elem+(double)par_chain.n_rows);
    if(verbose) cout << "updates of s/hyperp follow the ratio ratio: " <<  updateratio << endl;

    additerations(); // chain starts in 1
    for(int l=1; l<iterations+1;l++){
      if(verbose) cout << "Step: " << l << endl;
      // update
        if(Rcpp::runif(1)(0) < updateratio ){ // update selection coefficients
          Supdate();
          Gcopy();
        }else{ // update hyperparameters
          Gupdate();
          Scopy();
        }
        // if(verbose) cout << s_chain.col(iter) << endl;
        FITcopy();
        FITupdate();
        PROBupdate();
        if(verbose) cout << "w_1: " << w_chain(1,l) << endl;
        if(verbose) cout << "s_1: " << s_chain(1,l) << endl;
        if(verbose) cout << "b: " << par_chain(0,l) << endl;
        if(verbose) cout << "a: " << par_chain(1,l) << endl;
        if(verbose) cout << "p: " << par_chain(2,l) << endl;
        if(verbose) cout << "mu: " << par_chain(3,l) << endl;
        if(verbose) cout << "epi: " << par_chain(4,l) << endl;
        if(verbose) cout << "svar: " << par_chain(5,l) << endl;
        if(verbose) cout << "ss: " << par_chain(6,l) << endl;
        if(verbose) cout << "Prior: " << pri_chain(l) << endl;
        if(verbose) cout << "Likelihood: " << lik_chain(l) << endl;
        if(verbose) cout << "Posterior: " << prob_chain(l) << endl;
      // decide
        decide();
      // print
        additerations();// printProgress(iter/iterations);
        if(verbose) cout << endl;
    }

    /// END
    ///////////////////////////////////////////////////////////////////////////
    cout<< endl;
    cout<< "Summary:"<<endl;
    cout<< "Acceptance ratio = "<< sum(ok_chain) / iterations << endl;
    cout<< "Final posterior = "<< prob_chain(iterations) << endl;
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds1 = end-start;
    std::cout << "elapsed time: " << elapsed_seconds1.count() << " seconds" <<endl;
    cout<< "----------"<<endl;
    cout<< ends;
    }
};


//main function
// [[Rcpp::export]]
List napMCMC(
            arma::vec & y,
            arma::vec & h,
            const SEXP & A, // instead of arma::mat X,
            const arma::uvec m, // the positions of SNPs
            const arma::uvec n , // the positions of individuals
            arma::vec & s,
            double b=0.5,double bmin=0,double bmax=1,
            double a=0.1,double amin=0,double amax=1,
            double p=0.5, double pmin=0,double pmax=1,
            double mu=1,double mumin=0, double mumax=50,
            double epi=1,double epimin=1, double epimax=1,
            double svar=0.1,double svarmin=0.001, double svarmax=1,
            double ss=0.1,double ssmin=0, double ssmax=1,
            double smin=-0.98039215686274505668, double smax=10,
            int iterations=1e4,
            double burnin=0.1,
            int Smode=2,int LIKmode=2, int PRImode=1, int FITmode=2,
            bool verbose=false,
            bool test=false,
            double updateratio= -9.0,
            double bw= 0.001,
            std::string file2sink= "output.log",
            bool sink2file = false){

  // direct outputs
  if(sink2file) std::freopen(file2sink.c_str(), "w", stdout);

  // start class and analyses
  NAP nap(y,h,A,m,n,s,b,bmin,bmax,a,amin,amax,p,pmin,pmax,mu,mumin,mumax,epi,epimin, epimax,
          svar,svarmin, svarmax,ss,ssmin, ssmax,smin, smax,
          iterations,burnin,
          Smode,LIKmode,PRImode, FITmode,
          verbose, test,
          updateratio, bw
          );

  // run
  nap.MAIN();

  return(nap.report());
}

// [[Rcpp::export]]
List test_napMCMC(
            arma::vec & y,
            arma::vec & h,
            const SEXP & A, // instead of arma::mat X,
            const arma::uvec m, // the positions of SNPs
            const arma::uvec n , // the positions of individuals
            arma::vec & s,
            double b=0.5,double bmin=0,double bmax=1,
            double a=0.1,double amin=0,double amax=1,
            double p=0.5, double pmin=0,double pmax=1,
            double mu=1,double mumin=0, double mumax=50,
            double epi=1,double epimin=1, double epimax=1,
            double svar=0.1,double svarmin=0.001, double svarmax=1,
            double ss=0.1,double ssmin=0, double ssmax=1,
            double smin=-0.98039215686274505668, double smax=10,
            int iterations=1e4,
            double burnin=0.1,
            int Smode=2,int LIKmode=2, int PRImode=1, int FITmode=2,
            bool verbose=false,
            bool test=false,
            double updateratio= -9.0,
            double bw= 0.001,
            std::string file2sink= "output.log",
            bool sink2file = false){

  // direct outputs
  if(sink2file) std::freopen(file2sink.c_str(), "w", stdout);

  // start class and analyses
  NAP nap(y,h,A,m,n,s,b,bmin,bmax,a,amin,amax,p,pmin,pmax,mu,mumin,mumax,epi,epimin, epimax,
          svar,svarmin, svarmax,ss,ssmin, ssmax,smin, smax,
          iterations,burnin,
          Smode,LIKmode,PRImode, FITmode,
          verbose, test,
          updateratio, bw
          );

  // run
  nap.RUNTEST();

  return(nap.report());
}
