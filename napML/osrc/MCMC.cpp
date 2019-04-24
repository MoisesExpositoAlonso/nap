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

// To transform from two matrices
// mat  A = randu<mat>(5,5);
// fmat B = conv_to<fmat>::from(A);


#include <Rcpp.h>
#include <fstream>
#include <sstream>
#include <string>
using namespace Rcpp;


//Function is taking a path to a numeric file and return the same datain a NumericMatrix object
// [[Rcpp::export]]
NumericMatrix readfilecpp(std::string path, int rows, int cols){
    NumericMatrix output(cols,rows);// output matrix (specifying the size is critical otherwise R crashes)
    std::ifstream myfile(path.c_str()); //Opens the file. c_str is mandatory here so that ifstream accepts the string path
    std::string line;
    std::getline(myfile,line,'\n'); //skip the first line (col names in our case). Remove those lines if note necessary

       for (int row=0; row<rows; ++row){ // basic idea: getline() will read lines row=0:19 and for each line will put the value separated by ',' into 46749 columns
           std::string line;
           std::getline(myfile,line,'\n'); //Starts at the second line because the first one was ditched previously
           if(!myfile.good() ) //If end of rows then break
               break;
           std::stringstream iss(line); // take the line into a stringstream
           std::string val;
           std::getline(iss,val,','); ///skips the first column (row names)
           for (int col=0; col<cols; ++col ){
             std::string val;
             std::getline(iss,val,','); //reads the stringstream line and separate it into 49749 values (that were delimited by a ',' in the stringstream)
             std::stringstream convertor(val); //get the results into another stringstream 'convertor'
             convertor >> output(col,row); //put the result into our output matrix at for the actual row and col
           }
       }
    return(output);
}


//Function is taking a path to a numeric file and return the same datain a NumericMatrix object
// [[Rcpp::export]]
NumericVector readfilevector(std::string path, int rows){
    NumericVector output(rows);// output matrix (specifying the size is critical otherwise R crashes)
    std::ifstream myfile(path.c_str()); //Opens the file. c_str is mandatory here so that ifstream accepts the string path
    std::string line;
    std::getline(myfile,line,'\n'); //skip the first line (col names in our case). Remove those lines if note necessary

       for (int row=0; row<rows; ++row){ // basic idea: getline() will read lines row=0:19 and for each line will put the value separated by ',' into 46749 columns
           std::string line;
           std::getline(myfile,line,'\n'); //Starts at the second line because the first one was ditched previously
           if(!myfile.good() ) //If end of rows then break
               break;
           std::stringstream iss(line); // take the line into a stringstream
           iss >> output(row); //put the result into our output matrix at for the actual row and col
         }
    return(output);
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

// [[Rcpp::export]]
double medianC(const arma::vec& ar, double burnin=0.1){
  arma::vec tosub;
  tosub = arma::linspace((int)round(ar.n_elem * burnin), ar.n_elem-1);
  arma::uvec indices = arma::conv_to<arma::uvec>::from(tosub);
  return(arma::median(ar.elem( indices)));
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

	arma::vec coef(p, arma::fill::zeros);
	arma::colvec y_tilde = join_cols(y, arma::zeros<arma::colvec>(p));

		arma::mat X_tilde = join_cols(X, sqrt(lambda) * arma::eye(p, p));
		arma::mat Q, R;
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
  arma::vec hunique = unique(h-1);
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
double ePRIOR(double e, double const & m=1, double const & v=0.1){
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
                  bool verbose=false,
                  int LIKmode=2,
                  bool printall=false){
  double L=0;
  switch(LIKmode){
  case 1: // moc likelihood
    break;
  case 2:
    // arma::vec w_; // I think this was causing problems
    // if(epi!=1){w_=pow(w,epi);}
    // else{w_=w;}

    double LL;
      for(int i=0; i< y.n_elem ; i ++){ // check for infinity
        LL= LLGaussMix(y(i)/mu,w(hs(i)),w(hs(i))*b+a,p);

        if((verbose and (std::isinf(LL) or std::isnan(LL))) or printall ){
        // if(std::isinf(LL)){
          cout << "---" << endl;
          cout << i << endl;
          cout << y(i) << " "<< w(hs(i)) << " "<< w(hs(i))*b+a <<" "<< p << endl;
          cout << LL << endl;
        }
        if(!std::isnan(LL)) L += LL; // if there is one nan, all sum is nan
      }
    break;
  };
  if(verbose) cout<< L << endl; /////// CHECKING LIKELIHOOD!
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
arma::vec wCBM(SEXP A,
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
// // [[Rcpp::export]]
// arma::vec ssimC(int nsnp,
//                 double svar){
//   arma::vec news(nsnp);
//   if(svar==0.05){
//     news.load("databig/s_svar05.txt");
//   }else if(svar==0.25){
//     news.load("databig/s_svar25.txt");
//   }else if(svar==0.5){
//     news.load("databig/s_svar50.txt");
//   } else{
//     cout << "S values not found for accuracy calculations!" << endl;
//     news=exp(Rcpp::rnorm(nsnp,0,svar))-1;
//   }
//   return news;
// }

// [[Rcpp::export]]
arma::vec ssimC(int nsnp,
                double svar){
  arma::vec news(nsnp);
  news.load("databig/s_svar01.txt");
  news= exp( log(1+news) * (svar/0.01) )-1 ;
  //
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
  return news;
}

// // [[Rcpp::export]]
// List wsimC(SEXP A,
//                    const arma::uvec & myrows,
//                    const arma::uvec & mycols,
//                    const double & b,
//                    const double & a,
//                    const double & p,
//                    const double & mu,
//                    const double & epi,
//                    const double & svar,
//                    const int & mode,
//                    int rep=1
//                   ){
//
//   arma::vec s, Ey, y, h;
//   s = ssimC(mycols.n_elem,svar);
//   Ey=wCBM(A,s,myrows,mycols,mode,epi,mu);
//   y=sampleWC(Ey,a,b,p,rep = rep);
//   // h=Rcpp::sort(Rcpp::rep(myrows,rep)); // CHECK BUG - this needed if want to simulate replicates per genotype
//   return(List::create(Named("w")=Ey ,
//                       Named("y")=y )
//         );
// }


///////////////////////////////////////////////////////////////////////////////
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
    arma::vec s; int nsnp; double smin; double smax;
    arma::vec y;
    arma::vec hs;
    arma::vec w;
    int nind;
    arma::mat X;
    arma::uvec  mycols;
    arma::uvec  myrows;


    // misc
    int nupdates=1; // 0 for all
    bool verbose=false;
    bool test=false;
    double updateratio;
    double bw; /// BANDWIDT OF CHANGES IS REALLY IMPORTANT HARDCODED PARAM.
    int iterations;
    double burnin;
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
            const arma::uvec & mycols_, // the positions of SNPs // check start 0 vs 1
            const arma::uvec & myrows_ , // the positions of individuals
            arma::vec & s_,

            double b_,double bmin_,double bmax_,
            double a_,double amin_,double amax_,
            double p_,double pmin_,double pmax_,
            double mu_,double mumin_, double mumax_,
            double epi_,double epimin_, double epimax_,
            double svar_,double svarmin_, double svarmax_,
            double ss_,double ssmin_, double ssmax_,
            double smin_, double smax_,

            int iterations_, double burnin_,
            int Smode_,int LIKmode_, int PRImode_, int FITmode_,
            bool verbose_,
            bool test_,
            double updateratio_,
            double bw_
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
      mycols=mycols_;
      myrows=myrows_;
      nsnp=mycols.n_elem;
      nind=myrows.n_elem;
      s=s_;

      y=y_;
      // hs=hsub(h); // BUG CHECK
      hs=h-1; // BUG CHECK // old

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
     // GENOME MATRIX // Trying to get rid of the transformation
    //   X.zeros(n.n_elem,m.n_elem);
    //   if(TYPEOF(A) == EXTPTRSXP){
    //   cout<< "Reading external pointer and subsetting genome matrix ... "<<endl;
    //      X=BMsubset(A,n,m);
    //   }else if(TYPEOF(A) == REALSXP){
    //   cout<< "Matrix provided already subsetted from R"<<endl;
    //     // NumericMatrix Xr(A);
    //     // cout << " nrow= " << Xr.nrow() << " ncol= " << Xr.ncol() << endl;
    //     // arma::mat X(Xr.begin(), Xr.nrow(), Xr.ncol(), false); // this was the original
    //      // X(Xr.begin(), Xr.nrow(), Xr.ncol(), false); // without the arma::mat at the beginning does not work
    //      // arma::mat X=A; // no viable converstion
    //     NumericMatrix Xr(A);
    //     X = as<arma::mat>( Xr ) ;
    //     X=Xmvcenter(X); // CHECK FOR BUGS!
    //
    // }

    }; // end class

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
          // news=exp(Rcpp::rnorm(nsnp,0,0.01))-1;
          news=Rcpp::rnorm(nsnp,0,0.01);
          break;
        }
       s=news;
       s_chain.col(0)=s;
      }
     // void Sstartgwa(){
     //    arma::vec news;
     //    arma::vec ymeans=My(y,hs);
     //    news=BMridge(X,ymeans,1); ///////////// BUG CAREFULL GWA
     //    for(int j;j<news.n_elem;j++){
     //      if(news(j)< -1) news(j) = -0.99;
     //    }
     //    s=news;
     //    s_chain.col(0)=s;
     // }
     void Sstartdif(){
       if(verbose) cout << "starting S using w1-w0" <<endl;
        arma::vec news;
        // if(verbose) cout << "calculate average per genotype" <<endl;
        // arma::vec ymeans=My(y,hs);
        if(verbose) cout << "calculate w1-w0" <<endl;
        news=BMs(A,y,mycols,myrows); ///////////// BUG CAREFULL GWA
        if(verbose) cout << "correct s < -1 " <<endl;
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
     void Slogn(double svar){ // CHECK BUG, fixed svar below!
        double meanhere,newval;
        int randomIndex = rand() % s.size();
        meanhere=log(1+s(randomIndex)); //*** transform mean to log dimension
        if(std::isinf(meanhere)) meanhere= 0; // check for infinity
        newval = Rcpp::rnorm(1,meanhere,0.01)(0);
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
   void Ggeneralist(){
      if(verbose) cout << "starting hyperparameters" <<endl;
      par_chain(0,0)=0.1;
      par_chain(1,0)=0.1;
      par_chain(2,0)=0.1;
      par_chain(3,0)=1;
      par_chain(4,0)=1;
      par_chain(5,0)=0.1;
      par_chain(6,0)=0;
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
              if(newval<0) newval=1; // CHECK BUG
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
  }
  ////////////////////
  // calculate fitness
  void FITstart(){
    if(verbose) cout << "start computing fitness from scratch" <<endl;
    // arma::vec w=wC(X,s_chain.col(0),FITmode,par_chain(4,0)); // matrix, s, mode, epi, mu
    arma::vec w=wCBM(A,s_chain.col(0),
                     mycols,myrows, // CHECK BUG positions. Inside NAP already subtracted
                     FITmode,par_chain(4,0)); // CHECK BUG
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
      if(verbose) cout << "selection coef changed for SNP " <<indecito <<endl;
        Rcpp::XPtr<BigMatrix> bigMat(A);
        MatrixAccessor<double> macc(*bigMat);
        double val=0;
        for(int j=0; j < myrows.n_elem ; j++){
          val= (macc[indecito][myrows(j)-1]) ;
          if(val!=0){
            unwupdate(wnew(j),s_chain(indecito,iter-1), val, FITmode);
            wupdate(wnew(j),s_chain(indecito,iter),val, FITmode);
          }
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
                      PRImode); // + ePRIOR(par_chain(4,iter)); // epi in position 4
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
                       PRImode);// + ePRIOR(par_chain(4,iter)); // epi in position 4
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
    if(Paccept>1) Paccept=1;
    bool accept = Rcpp::runif(1)(0)<=Paccept;
    if(verbose) cout<< "Pt-1 " << prob_chain(iter-1) <<endl;
    if(verbose) cout<< "Pt " << prob_chain(iter) <<endl;
    if(verbose) cout<< "Accept proposal " << accept <<endl;
    if(accept){
      ok_chain(iter)=1;
      // bw=0.001; // BUG CHECK -- interesting application, if accepted you want to climb slowly
    }else{
      bw+=0.001; // BUG CHECK -- interesting application, although it is incompatible with the MH algorithm
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
    if(verbose) cout << "Summarizing the chain for report" << endl;
    // Summarize chain
    arma::vec shat=medianCmat(s_chain);
    arma::vec par=medianCmat(par_chain);
    if(verbose) cout << "Calculating final individualfitness" << endl;
    arma::vec what=wCBM(A,shat,mycols,myrows,FITmode);
    // Calculate pseudo-p
    if(verbose) cout << "Calculating pseudo-P-value" << endl;
    // arma::vec pp=;
    // Calculate ranges parameters
    // arma::vec par_up =
    // arma::vec par_low =
    // Calculate accuracies
    if(verbose) cout << "Calculating prediction accuracy" << endl;
    List res=accuracies(y,what);

    return List::create(Named("chain") = s_chain.t(),
                        Named("shat")=shat,
                        // Named("pp")=pp,
                        Named("parchain") = par_chain.t(),
                        Named("par")=par,
                        Named("par_up")=par,
                        Named("par_low")=par,
                        Named("w")=what,
                        Named("accuracy")=res,
                        Named("posterior") = prob_chain,
                        Named("prior") = pri_chain,
                        Named("likelihood") = lik_chain,
                        Named("final_likelihood") =
                          LIKELIHOOD(y,hs,
                                     what,
                                     par(0), // b in position 0
                                     par(1), // a in position 1
                                     par(2), // p in position 2
                                     par(3), // mu in position 3
                                     par(4), // epi in position 3
                                     verbose,LIKmode),
                        Named("accept") = ok_chain,
                        Named("FITmode") = FITmode);
  }


  ////////////////////
  // MAIN
  ////////////////////
    void RUNTEST(bool verbose=false){
      int iter=0;
      // S start();
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
    cout << "----------------------------------------------------------"<< endl;
    cout << "----------------------------------------------------------"<< endl;

      int iter=0;
      cout<< endl;
      cout<< "# individual's observations = "<<  y.n_elem <<endl;
      cout<< "# of SNPs = "<<  s.n_elem <<endl;
      cout<< "# iterations = "<< iterations <<endl;
      cout<< "Prior mode  = "<< PRImode <<endl;
      cout<< "Likelihood mode  = "<<  LIKmode <<endl;
      cout<< "S sampling mode  = "<< Smode <<endl;
      cout<< "Fitness mode  = "<< FITmode <<endl;
      cout<< "Verbose = "<< (bool)verbose <<endl;
      cout<< "Initializing & running MCMC ... "<<endl;
      std::chrono::time_point<std::chrono::system_clock> start, end; // start chronometer values
      start = std::chrono::system_clock::now();
      cout<< endl;

    /// Start NAP & handle -INF probability
    ///////////////////////////////////////////////////////////////////////////
    // Gstart();
    // Sstartdif();
    FITstart();
    PROBstart();

    //// Uncomment to stop if not right start
    // if(std::isinf(prob_chain(0)) ||std::isnan(prob_chain(0))){
    //   stop("Attempted initializing but Probabiligy is -Inf !!! ... tune starting parameters!");
    // }

    //// Uncomment if trying to find a good start
   int attemptcounter=0;
   int maxattempts=1000;
   while(std::isinf(prob_chain(0)) ||std::isnan(prob_chain(0))){
     if(verbose && attemptcounter==0) cout << "Running proposals and first probability" << endl;
     if(attemptcounter==maxattempts) stop("Attempted initializing 1000 times!!! ... tune starting parameters!");
     if(attemptcounter==1) cout << "Posterior is infinite!!!. Attempting new start..." << endl<< endl;
     // Sstart();
     Gstart(); // BUG careful
     // Ggeneralist(); // BUG careful
     FITstart();
     PROBstart();
     attemptcounter++;
   }
   if(verbose) cout << "Successful MCMC start after " << attemptcounter << " attempts " << endl<< endl;


    /// Run MCMC
    ///////////////////////////////////////////////////////////////////////////
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
        // if(verbose) cout << "w_1: " << w_chain.col(l) << endl;
        // if(verbose) cout << "s_1: " << s_chain.col(l) << endl;
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
    cout<< "Final prior = "<< pri_chain(iterations) << endl;
    cout<< "Final likelihood = "<< lik_chain(iterations) << endl;
    cout<< "Starting posterior = "<< prob_chain(0) << endl;
    cout<< "Final posterior = "<< prob_chain(iterations) << endl;
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds1 = end-start;
    std::cout << "elapsed time: " << elapsed_seconds1.count() << " seconds"<<endl;
    cout << "----------------------------------------------------------"<< endl;
    cout<< ends;
    }

}; // end class NAP


//main function
// [[Rcpp::export]]
List napMCMC(
            arma::vec & y,
            arma::vec & h,
            const SEXP & A, // instead of arma::mat X,
            const arma::uvec mycols, // the positions of SNPs
            const arma::uvec myrows , // the positions of individuals
            arma::vec & s,
            double b=0.1,double bmin=0,double bmax=1,
            double a=0.1,double amin=0,double amax=1,
            double p=0.1, double pmin=0,double pmax=1,
            double mu=1,double mumin=1, double mumax=1,
            double epi=1,double epimin=0.9, double epimax=1.1,
            double svar=0.1,double svarmin=0.001, double svarmax=0.5,
            double ss=0,double ssmin=0, double ssmax=0.0,
            double smin=-0.98039215686274505668, double smax=10,
            int iterations=1e4,
            double burnin=0.1,
            int Smode=1,int LIKmode=2, int PRImode=1, int FITmode=2,
            bool verbose=false,
            bool test=false,
            double updateratio= -9.0,
            double bw= 0.001,
            std::string file2sink= "output.log",
            bool sink2file = false
        ){

  // direct outputs
  if(sink2file) std::freopen(file2sink.c_str(), "w", stdout);

  // start class and analyses
  cout << "COLD MCMC CHAIN" << endl;
  NAP nap1(y,h,A,mycols,myrows,s,b,bmin,bmax,a,amin,amax,p,pmin,pmax,mu,mumin,mumax,epi,epimin, epimax,
          svar,svarmin, svarmax,ss,ssmin, ssmax,smin, smax,
          iterations,burnin,
          Smode,LIKmode,PRImode, FITmode,
          verbose, test,
          updateratio, bw
          );
  nap1.MAIN();
  List coldreport = nap1.report();

  cout << "HOT MCMC CHAIN" << endl;
  NAP nap2(y,h,A,mycols,myrows, s,b,bmin,bmax,a,amin,amax,p,pmin,pmax,mu,mumin,mumax,epi,epimin, epimax,
          svar,svarmin, svarmax,ss,ssmin, ssmax,smin, smax,
          iterations,burnin,
          Smode,LIKmode,PRImode, FITmode,
          verbose, test,
          updateratio, bw*10
          );
  nap2.MAIN();
  List hotreport = nap2.report();

  cout << "CHAINS END" << endl;

  if( (double)coldreport["final_likelihood"] > (double)hotreport["final_likelihood"] ){
    cout << "KEEPING COLD CHAIN" << endl;
    return(coldreport);
  }else{
    cout << "KEEPING HOT CHAIN" << endl;
    return(hotreport);
  }

  // return(coldreport);

}

