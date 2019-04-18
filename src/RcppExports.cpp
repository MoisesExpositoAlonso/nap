// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// readfilecpp
NumericMatrix readfilecpp(std::string path, int rows, int cols);
RcppExport SEXP _nap_readfilecpp(SEXP pathSEXP, SEXP rowsSEXP, SEXP colsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type path(pathSEXP);
    Rcpp::traits::input_parameter< int >::type rows(rowsSEXP);
    Rcpp::traits::input_parameter< int >::type cols(colsSEXP);
    rcpp_result_gen = Rcpp::wrap(readfilecpp(path, rows, cols));
    return rcpp_result_gen;
END_RCPP
}
// readfilevector
NumericVector readfilevector(std::string path, int rows);
RcppExport SEXP _nap_readfilevector(SEXP pathSEXP, SEXP rowsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type path(pathSEXP);
    Rcpp::traits::input_parameter< int >::type rows(rowsSEXP);
    rcpp_result_gen = Rcpp::wrap(readfilevector(path, rows));
    return rcpp_result_gen;
END_RCPP
}
// BMsimulate
bool BMsimulate(SEXP A);
RcppExport SEXP _nap_BMsimulate(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(BMsimulate(A));
    return rcpp_result_gen;
END_RCPP
}
// BMwrite012
bool BMwrite012(SEXP A, std::string outfile);
RcppExport SEXP _nap_BMwrite012(SEXP ASEXP, SEXP outfileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type A(ASEXP);
    Rcpp::traits::input_parameter< std::string >::type outfile(outfileSEXP);
    rcpp_result_gen = Rcpp::wrap(BMwrite012(A, outfile));
    return rcpp_result_gen;
END_RCPP
}
// BMwritePED
bool BMwritePED(SEXP A, std::string outfile);
RcppExport SEXP _nap_BMwritePED(SEXP ASEXP, SEXP outfileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type A(ASEXP);
    Rcpp::traits::input_parameter< std::string >::type outfile(outfileSEXP);
    rcpp_result_gen = Rcpp::wrap(BMwritePED(A, outfile));
    return rcpp_result_gen;
END_RCPP
}
// BMsubset
arma::Mat<double> BMsubset(SEXP A, const arma::uvec& myrows, const arma::uvec& mycols);
RcppExport SEXP _nap_BMsubset(SEXP ASEXP, SEXP myrowsSEXP, SEXP mycolsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type myrows(myrowsSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type mycols(mycolsSEXP);
    rcpp_result_gen = Rcpp::wrap(BMsubset(A, myrows, mycols));
    return rcpp_result_gen;
END_RCPP
}
// upperTmat
arma::vec upperTmat(const arma::mat mymat);
RcppExport SEXP _nap_upperTmat(SEXP mymatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type mymat(mymatSEXP);
    rcpp_result_gen = Rcpp::wrap(upperTmat(mymat));
    return rcpp_result_gen;
END_RCPP
}
// Xmcenter
arma::mat Xmcenter(arma::mat X);
RcppExport SEXP _nap_Xmcenter(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(Xmcenter(X));
    return rcpp_result_gen;
END_RCPP
}
// Xmvcenter
arma::mat Xmvcenter(arma::mat X);
RcppExport SEXP _nap_Xmvcenter(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(Xmvcenter(X));
    return rcpp_result_gen;
END_RCPP
}
// LDrelative
arma::mat LDrelative(SEXP A, arma::uvec mycols, bool debug);
RcppExport SEXP _nap_LDrelative(SEXP ASEXP, SEXP mycolsSEXP, SEXP debugSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type mycols(mycolsSEXP);
    Rcpp::traits::input_parameter< bool >::type debug(debugSEXP);
    rcpp_result_gen = Rcpp::wrap(LDrelative(A, mycols, debug));
    return rcpp_result_gen;
END_RCPP
}
// medianC
double medianC(const arma::vec& ar, double burnin);
RcppExport SEXP _nap_medianC(SEXP arSEXP, SEXP burninSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type ar(arSEXP);
    Rcpp::traits::input_parameter< double >::type burnin(burninSEXP);
    rcpp_result_gen = Rcpp::wrap(medianC(ar, burnin));
    return rcpp_result_gen;
END_RCPP
}
// medianCmat
arma::vec medianCmat(const arma::mat& ar, double burnin);
RcppExport SEXP _nap_medianCmat(SEXP arSEXP, SEXP burninSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type ar(arSEXP);
    Rcpp::traits::input_parameter< double >::type burnin(burninSEXP);
    rcpp_result_gen = Rcpp::wrap(medianCmat(ar, burnin));
    return rcpp_result_gen;
END_RCPP
}
// accuracies
List accuracies(const arma::vec& y, const arma::vec& what);
RcppExport SEXP _nap_accuracies(SEXP ySEXP, SEXP whatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type what(whatSEXP);
    rcpp_result_gen = Rcpp::wrap(accuracies(y, what));
    return rcpp_result_gen;
END_RCPP
}
// My
arma::colvec My(const arma::colvec& y, const arma::colvec& h);
RcppExport SEXP _nap_My(SEXP ySEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(My(y, h));
    return rcpp_result_gen;
END_RCPP
}
// BMs
arma::colvec BMs(const SEXP A, const arma::colvec& y, const arma::uvec& mycols, const arma::uvec& myrows, bool debug);
RcppExport SEXP _nap_BMs(SEXP ASEXP, SEXP ySEXP, SEXP mycolsSEXP, SEXP myrowsSEXP, SEXP debugSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const SEXP >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type mycols(mycolsSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type myrows(myrowsSEXP);
    Rcpp::traits::input_parameter< bool >::type debug(debugSEXP);
    rcpp_result_gen = Rcpp::wrap(BMs(A, y, mycols, myrows, debug));
    return rcpp_result_gen;
END_RCPP
}
// BMmgwa
arma::colvec BMmgwa(arma::mat X0, arma::colvec& yraw, bool debug);
RcppExport SEXP _nap_BMmgwa(SEXP X0SEXP, SEXP yrawSEXP, SEXP debugSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X0(X0SEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type yraw(yrawSEXP);
    Rcpp::traits::input_parameter< bool >::type debug(debugSEXP);
    rcpp_result_gen = Rcpp::wrap(BMmgwa(X0, yraw, debug));
    return rcpp_result_gen;
END_RCPP
}
// BMcgwa
arma::colvec BMcgwa(arma::mat X0, const arma::vec& yraw);
RcppExport SEXP _nap_BMcgwa(SEXP X0SEXP, SEXP yrawSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X0(X0SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type yraw(yrawSEXP);
    rcpp_result_gen = Rcpp::wrap(BMcgwa(X0, yraw));
    return rcpp_result_gen;
END_RCPP
}
// BMridge
arma::vec BMridge(arma::mat X0, const arma::colvec& yraw, const double& lambda);
RcppExport SEXP _nap_BMridge(SEXP X0SEXP, SEXP yrawSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X0(X0SEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type yraw(yrawSEXP);
    Rcpp::traits::input_parameter< const double& >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(BMridge(X0, yraw, lambda));
    return rcpp_result_gen;
END_RCPP
}
// sample_num
NumericVector sample_num(NumericVector x, int size, bool replace, NumericVector prob);
RcppExport SEXP _nap_sample_num(SEXP xSEXP, SEXP sizeSEXP, SEXP replaceSEXP, SEXP probSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< bool >::type replace(replaceSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type prob(probSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_num(x, size, replace, prob));
    return rcpp_result_gen;
END_RCPP
}
// hsub
arma::vec hsub(const arma::vec& h);
RcppExport SEXP _nap_hsub(SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(hsub(h));
    return rcpp_result_gen;
END_RCPP
}
// runif_reflect
double runif_reflect(double minhere, double maxhere, double min, double max);
RcppExport SEXP _nap_runif_reflect(SEXP minhereSEXP, SEXP maxhereSEXP, SEXP minSEXP, SEXP maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type minhere(minhereSEXP);
    Rcpp::traits::input_parameter< double >::type maxhere(maxhereSEXP);
    Rcpp::traits::input_parameter< double >::type min(minSEXP);
    Rcpp::traits::input_parameter< double >::type max(maxSEXP);
    rcpp_result_gen = Rcpp::wrap(runif_reflect(minhere, maxhere, min, max));
    return rcpp_result_gen;
END_RCPP
}
// whichchanged
int whichchanged(const arma::vec& s1, const arma::vec& s2);
RcppExport SEXP _nap_whichchanged(SEXP s1SEXP, SEXP s2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type s1(s1SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type s2(s2SEXP);
    rcpp_result_gen = Rcpp::wrap(whichchanged(s1, s2));
    return rcpp_result_gen;
END_RCPP
}
// logmixgaussian
double logmixgaussian(const double& s, double svar, double ss);
RcppExport SEXP _nap_logmixgaussian(SEXP sSEXP, SEXP svarSEXP, SEXP ssSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type svar(svarSEXP);
    Rcpp::traits::input_parameter< double >::type ss(ssSEXP);
    rcpp_result_gen = Rcpp::wrap(logmixgaussian(s, svar, ss));
    return rcpp_result_gen;
END_RCPP
}
// iterprior
double iterprior(const arma::colvec& s, const double& svar, const double& ss);
RcppExport SEXP _nap_iterprior(SEXP sSEXP, SEXP svarSEXP, SEXP ssSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type s(sSEXP);
    Rcpp::traits::input_parameter< const double& >::type svar(svarSEXP);
    Rcpp::traits::input_parameter< const double& >::type ss(ssSEXP);
    rcpp_result_gen = Rcpp::wrap(iterprior(s, svar, ss));
    return rcpp_result_gen;
END_RCPP
}
// PRIOR
double PRIOR(const arma::colvec& s, double const& par1, double const& par2, const int& type);
RcppExport SEXP _nap_PRIOR(SEXP sSEXP, SEXP par1SEXP, SEXP par2SEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type s(sSEXP);
    Rcpp::traits::input_parameter< double const& >::type par1(par1SEXP);
    Rcpp::traits::input_parameter< double const& >::type par2(par2SEXP);
    Rcpp::traits::input_parameter< const int& >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(PRIOR(s, par1, par2, type));
    return rcpp_result_gen;
END_RCPP
}
// ePRIOR
double ePRIOR(double e, double const& m, double const& v);
RcppExport SEXP _nap_ePRIOR(SEXP eSEXP, SEXP mSEXP, SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type e(eSEXP);
    Rcpp::traits::input_parameter< double const& >::type m(mSEXP);
    Rcpp::traits::input_parameter< double const& >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(ePRIOR(e, m, v));
    return rcpp_result_gen;
END_RCPP
}
// uPRIOR
double uPRIOR(double pri, const arma::colvec& s0, const arma::colvec& s1, const double& par10, const double& par11, const double& par20, const double& par21, const int& type);
RcppExport SEXP _nap_uPRIOR(SEXP priSEXP, SEXP s0SEXP, SEXP s1SEXP, SEXP par10SEXP, SEXP par11SEXP, SEXP par20SEXP, SEXP par21SEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type pri(priSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type s0(s0SEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type s1(s1SEXP);
    Rcpp::traits::input_parameter< const double& >::type par10(par10SEXP);
    Rcpp::traits::input_parameter< const double& >::type par11(par11SEXP);
    Rcpp::traits::input_parameter< const double& >::type par20(par20SEXP);
    Rcpp::traits::input_parameter< const double& >::type par21(par21SEXP);
    Rcpp::traits::input_parameter< const int& >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(uPRIOR(pri, s0, s1, par10, par11, par20, par21, type));
    return rcpp_result_gen;
END_RCPP
}
// LLGaussMix
double LLGaussMix(const double& y, const double& w, const double& v, const double& p);
RcppExport SEXP _nap_LLGaussMix(SEXP ySEXP, SEXP wSEXP, SEXP vSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const double& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const double& >::type v(vSEXP);
    Rcpp::traits::input_parameter< const double& >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(LLGaussMix(y, w, v, p));
    return rcpp_result_gen;
END_RCPP
}
// LIKELIHOOD
double LIKELIHOOD(const arma::vec& y, const arma::vec& hs, const arma::vec& w, const double& b, const double& a, const double& p, const double& mu, const double& epi, bool verbose, int LIKmode, bool printall);
RcppExport SEXP _nap_LIKELIHOOD(SEXP ySEXP, SEXP hsSEXP, SEXP wSEXP, SEXP bSEXP, SEXP aSEXP, SEXP pSEXP, SEXP muSEXP, SEXP epiSEXP, SEXP verboseSEXP, SEXP LIKmodeSEXP, SEXP printallSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type hs(hsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const double& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double& >::type p(pSEXP);
    Rcpp::traits::input_parameter< const double& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const double& >::type epi(epiSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< int >::type LIKmode(LIKmodeSEXP);
    Rcpp::traits::input_parameter< bool >::type printall(printallSEXP);
    rcpp_result_gen = Rcpp::wrap(LIKELIHOOD(y, hs, w, b, a, p, mu, epi, verbose, LIKmode, printall));
    return rcpp_result_gen;
END_RCPP
}
// wC
arma::vec wC(const arma::Mat<double> X, const arma::vec& s, const int& mode, double epi, bool verbose);
RcppExport SEXP _nap_wC(SEXP XSEXP, SEXP sSEXP, SEXP modeSEXP, SEXP epiSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::Mat<double> >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type s(sSEXP);
    Rcpp::traits::input_parameter< const int& >::type mode(modeSEXP);
    Rcpp::traits::input_parameter< double >::type epi(epiSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(wC(X, s, mode, epi, verbose));
    return rcpp_result_gen;
END_RCPP
}
// wCBM
arma::vec wCBM(SEXP A, const arma::vec& s, const arma::uvec& mycols, const arma::uvec& myrows, const int& mode, double epi, bool verbose);
RcppExport SEXP _nap_wCBM(SEXP ASEXP, SEXP sSEXP, SEXP mycolsSEXP, SEXP myrowsSEXP, SEXP modeSEXP, SEXP epiSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type s(sSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type mycols(mycolsSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type myrows(myrowsSEXP);
    Rcpp::traits::input_parameter< const int& >::type mode(modeSEXP);
    Rcpp::traits::input_parameter< double >::type epi(epiSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(wCBM(A, s, mycols, myrows, mode, epi, verbose));
    return rcpp_result_gen;
END_RCPP
}
// wCBMupdate
void wCBMupdate(arma::vec wnew, SEXP A, const arma::uvec& mycols, const arma::uvec& myrows, const int& indecito, const double& s0, const double& s1, const int& mode);
RcppExport SEXP _nap_wCBMupdate(SEXP wnewSEXP, SEXP ASEXP, SEXP mycolsSEXP, SEXP myrowsSEXP, SEXP indecitoSEXP, SEXP s0SEXP, SEXP s1SEXP, SEXP modeSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type wnew(wnewSEXP);
    Rcpp::traits::input_parameter< SEXP >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type mycols(mycolsSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type myrows(myrowsSEXP);
    Rcpp::traits::input_parameter< const int& >::type indecito(indecitoSEXP);
    Rcpp::traits::input_parameter< const double& >::type s0(s0SEXP);
    Rcpp::traits::input_parameter< const double& >::type s1(s1SEXP);
    Rcpp::traits::input_parameter< const int& >::type mode(modeSEXP);
    wCBMupdate(wnew, A, mycols, myrows, indecito, s0, s1, mode);
    return R_NilValue;
END_RCPP
}
// wCBMupdateforR
arma::vec wCBMupdateforR(arma::vec wnew, SEXP A, const arma::uvec& mycols, const arma::uvec& myrows, const int& indecito, const double& s0, const double& s1, const int& mode);
RcppExport SEXP _nap_wCBMupdateforR(SEXP wnewSEXP, SEXP ASEXP, SEXP mycolsSEXP, SEXP myrowsSEXP, SEXP indecitoSEXP, SEXP s0SEXP, SEXP s1SEXP, SEXP modeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type wnew(wnewSEXP);
    Rcpp::traits::input_parameter< SEXP >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type mycols(mycolsSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type myrows(myrowsSEXP);
    Rcpp::traits::input_parameter< const int& >::type indecito(indecitoSEXP);
    Rcpp::traits::input_parameter< const double& >::type s0(s0SEXP);
    Rcpp::traits::input_parameter< const double& >::type s1(s1SEXP);
    Rcpp::traits::input_parameter< const int& >::type mode(modeSEXP);
    rcpp_result_gen = Rcpp::wrap(wCBMupdateforR(wnew, A, mycols, myrows, indecito, s0, s1, mode));
    return rcpp_result_gen;
END_RCPP
}
// sampleWC
arma::vec sampleWC(const arma::vec& w, const double& a, const double& b, const double& p, const int& rep);
RcppExport SEXP _nap_sampleWC(SEXP wSEXP, SEXP aSEXP, SEXP bSEXP, SEXP pSEXP, SEXP repSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const double& >::type p(pSEXP);
    Rcpp::traits::input_parameter< const int& >::type rep(repSEXP);
    rcpp_result_gen = Rcpp::wrap(sampleWC(w, a, b, p, rep));
    return rcpp_result_gen;
END_RCPP
}
// ssaveC
bool ssaveC(arma::vec s, std::string path);
RcppExport SEXP _nap_ssaveC(SEXP sSEXP, SEXP pathSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type s(sSEXP);
    Rcpp::traits::input_parameter< std::string >::type path(pathSEXP);
    rcpp_result_gen = Rcpp::wrap(ssaveC(s, path));
    return rcpp_result_gen;
END_RCPP
}
// ssimC
arma::vec ssimC(int nsnp, double svar);
RcppExport SEXP _nap_ssimC(SEXP nsnpSEXP, SEXP svarSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nsnp(nsnpSEXP);
    Rcpp::traits::input_parameter< double >::type svar(svarSEXP);
    rcpp_result_gen = Rcpp::wrap(ssimC(nsnp, svar));
    return rcpp_result_gen;
END_RCPP
}
// napMCMC
List napMCMC(arma::vec& y, arma::vec& h, const SEXP& A, const arma::uvec mycols, const arma::uvec myrows, arma::vec& s, double b, double bmin, double bmax, double a, double amin, double amax, double p, double pmin, double pmax, double mu, double mumin, double mumax, double epi, double epimin, double epimax, double svar, double svarmin, double svarmax, double ss, double ssmin, double ssmax, double smin, double smax, int iterations, double burnin, int Smode, int LIKmode, int PRImode, int FITmode, bool verbose, bool test, double updateratio, double bw, std::string file2sink, bool sink2file);
RcppExport SEXP _nap_napMCMC(SEXP ySEXP, SEXP hSEXP, SEXP ASEXP, SEXP mycolsSEXP, SEXP myrowsSEXP, SEXP sSEXP, SEXP bSEXP, SEXP bminSEXP, SEXP bmaxSEXP, SEXP aSEXP, SEXP aminSEXP, SEXP amaxSEXP, SEXP pSEXP, SEXP pminSEXP, SEXP pmaxSEXP, SEXP muSEXP, SEXP muminSEXP, SEXP mumaxSEXP, SEXP epiSEXP, SEXP epiminSEXP, SEXP epimaxSEXP, SEXP svarSEXP, SEXP svarminSEXP, SEXP svarmaxSEXP, SEXP ssSEXP, SEXP ssminSEXP, SEXP ssmaxSEXP, SEXP sminSEXP, SEXP smaxSEXP, SEXP iterationsSEXP, SEXP burninSEXP, SEXP SmodeSEXP, SEXP LIKmodeSEXP, SEXP PRImodeSEXP, SEXP FITmodeSEXP, SEXP verboseSEXP, SEXP testSEXP, SEXP updateratioSEXP, SEXP bwSEXP, SEXP file2sinkSEXP, SEXP sink2fileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type h(hSEXP);
    Rcpp::traits::input_parameter< const SEXP& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::uvec >::type mycols(mycolsSEXP);
    Rcpp::traits::input_parameter< const arma::uvec >::type myrows(myrowsSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type bmin(bminSEXP);
    Rcpp::traits::input_parameter< double >::type bmax(bmaxSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type amin(aminSEXP);
    Rcpp::traits::input_parameter< double >::type amax(amaxSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type pmin(pminSEXP);
    Rcpp::traits::input_parameter< double >::type pmax(pmaxSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type mumin(muminSEXP);
    Rcpp::traits::input_parameter< double >::type mumax(mumaxSEXP);
    Rcpp::traits::input_parameter< double >::type epi(epiSEXP);
    Rcpp::traits::input_parameter< double >::type epimin(epiminSEXP);
    Rcpp::traits::input_parameter< double >::type epimax(epimaxSEXP);
    Rcpp::traits::input_parameter< double >::type svar(svarSEXP);
    Rcpp::traits::input_parameter< double >::type svarmin(svarminSEXP);
    Rcpp::traits::input_parameter< double >::type svarmax(svarmaxSEXP);
    Rcpp::traits::input_parameter< double >::type ss(ssSEXP);
    Rcpp::traits::input_parameter< double >::type ssmin(ssminSEXP);
    Rcpp::traits::input_parameter< double >::type ssmax(ssmaxSEXP);
    Rcpp::traits::input_parameter< double >::type smin(sminSEXP);
    Rcpp::traits::input_parameter< double >::type smax(smaxSEXP);
    Rcpp::traits::input_parameter< int >::type iterations(iterationsSEXP);
    Rcpp::traits::input_parameter< double >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< int >::type Smode(SmodeSEXP);
    Rcpp::traits::input_parameter< int >::type LIKmode(LIKmodeSEXP);
    Rcpp::traits::input_parameter< int >::type PRImode(PRImodeSEXP);
    Rcpp::traits::input_parameter< int >::type FITmode(FITmodeSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< bool >::type test(testSEXP);
    Rcpp::traits::input_parameter< double >::type updateratio(updateratioSEXP);
    Rcpp::traits::input_parameter< double >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< std::string >::type file2sink(file2sinkSEXP);
    Rcpp::traits::input_parameter< bool >::type sink2file(sink2fileSEXP);
    rcpp_result_gen = Rcpp::wrap(napMCMC(y, h, A, mycols, myrows, s, b, bmin, bmax, a, amin, amax, p, pmin, pmax, mu, mumin, mumax, epi, epimin, epimax, svar, svarmin, svarmax, ss, ssmin, ssmax, smin, smax, iterations, burnin, Smode, LIKmode, PRImode, FITmode, verbose, test, updateratio, bw, file2sink, sink2file));
    return rcpp_result_gen;
END_RCPP
}
// rowSumsC
NumericVector rowSumsC(NumericMatrix x);
RcppExport SEXP _nap_rowSumsC(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rowSumsC(x));
    return rcpp_result_gen;
END_RCPP
}
// allelesimCvec
NumericVector allelesimCvec(double mu, double nu, double m, double wAA, double wAa, double waa, double p0, double psource, double Fi, double d, int N, int tmax);
RcppExport SEXP _nap_allelesimCvec(SEXP muSEXP, SEXP nuSEXP, SEXP mSEXP, SEXP wAASEXP, SEXP wAaSEXP, SEXP waaSEXP, SEXP p0SEXP, SEXP psourceSEXP, SEXP FiSEXP, SEXP dSEXP, SEXP NSEXP, SEXP tmaxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< double >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type wAA(wAASEXP);
    Rcpp::traits::input_parameter< double >::type wAa(wAaSEXP);
    Rcpp::traits::input_parameter< double >::type waa(waaSEXP);
    Rcpp::traits::input_parameter< double >::type p0(p0SEXP);
    Rcpp::traits::input_parameter< double >::type psource(psourceSEXP);
    Rcpp::traits::input_parameter< double >::type Fi(FiSEXP);
    Rcpp::traits::input_parameter< double >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type tmax(tmaxSEXP);
    rcpp_result_gen = Rcpp::wrap(allelesimCvec(mu, nu, m, wAA, wAa, waa, p0, psource, Fi, d, N, tmax));
    return rcpp_result_gen;
END_RCPP
}
// allelesimCmat
NumericMatrix allelesimCmat(double mu, double nu, double m, double wAA, double wAa, double waa, double p0, double psource, double Fi, double d, int N, int tmax, int rep);
RcppExport SEXP _nap_allelesimCmat(SEXP muSEXP, SEXP nuSEXP, SEXP mSEXP, SEXP wAASEXP, SEXP wAaSEXP, SEXP waaSEXP, SEXP p0SEXP, SEXP psourceSEXP, SEXP FiSEXP, SEXP dSEXP, SEXP NSEXP, SEXP tmaxSEXP, SEXP repSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< double >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type wAA(wAASEXP);
    Rcpp::traits::input_parameter< double >::type wAa(wAaSEXP);
    Rcpp::traits::input_parameter< double >::type waa(waaSEXP);
    Rcpp::traits::input_parameter< double >::type p0(p0SEXP);
    Rcpp::traits::input_parameter< double >::type psource(psourceSEXP);
    Rcpp::traits::input_parameter< double >::type Fi(FiSEXP);
    Rcpp::traits::input_parameter< double >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type tmax(tmaxSEXP);
    Rcpp::traits::input_parameter< int >::type rep(repSEXP);
    rcpp_result_gen = Rcpp::wrap(allelesimCmat(mu, nu, m, wAA, wAa, waa, p0, psource, Fi, d, N, tmax, rep));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_nap_readfilecpp", (DL_FUNC) &_nap_readfilecpp, 3},
    {"_nap_readfilevector", (DL_FUNC) &_nap_readfilevector, 2},
    {"_nap_BMsimulate", (DL_FUNC) &_nap_BMsimulate, 1},
    {"_nap_BMwrite012", (DL_FUNC) &_nap_BMwrite012, 2},
    {"_nap_BMwritePED", (DL_FUNC) &_nap_BMwritePED, 2},
    {"_nap_BMsubset", (DL_FUNC) &_nap_BMsubset, 3},
    {"_nap_upperTmat", (DL_FUNC) &_nap_upperTmat, 1},
    {"_nap_Xmcenter", (DL_FUNC) &_nap_Xmcenter, 1},
    {"_nap_Xmvcenter", (DL_FUNC) &_nap_Xmvcenter, 1},
    {"_nap_LDrelative", (DL_FUNC) &_nap_LDrelative, 3},
    {"_nap_medianC", (DL_FUNC) &_nap_medianC, 2},
    {"_nap_medianCmat", (DL_FUNC) &_nap_medianCmat, 2},
    {"_nap_accuracies", (DL_FUNC) &_nap_accuracies, 2},
    {"_nap_My", (DL_FUNC) &_nap_My, 2},
    {"_nap_BMs", (DL_FUNC) &_nap_BMs, 5},
    {"_nap_BMmgwa", (DL_FUNC) &_nap_BMmgwa, 3},
    {"_nap_BMcgwa", (DL_FUNC) &_nap_BMcgwa, 2},
    {"_nap_BMridge", (DL_FUNC) &_nap_BMridge, 3},
    {"_nap_sample_num", (DL_FUNC) &_nap_sample_num, 4},
    {"_nap_hsub", (DL_FUNC) &_nap_hsub, 1},
    {"_nap_runif_reflect", (DL_FUNC) &_nap_runif_reflect, 4},
    {"_nap_whichchanged", (DL_FUNC) &_nap_whichchanged, 2},
    {"_nap_logmixgaussian", (DL_FUNC) &_nap_logmixgaussian, 3},
    {"_nap_iterprior", (DL_FUNC) &_nap_iterprior, 3},
    {"_nap_PRIOR", (DL_FUNC) &_nap_PRIOR, 4},
    {"_nap_ePRIOR", (DL_FUNC) &_nap_ePRIOR, 3},
    {"_nap_uPRIOR", (DL_FUNC) &_nap_uPRIOR, 8},
    {"_nap_LLGaussMix", (DL_FUNC) &_nap_LLGaussMix, 4},
    {"_nap_LIKELIHOOD", (DL_FUNC) &_nap_LIKELIHOOD, 11},
    {"_nap_wC", (DL_FUNC) &_nap_wC, 5},
    {"_nap_wCBM", (DL_FUNC) &_nap_wCBM, 7},
    {"_nap_wCBMupdate", (DL_FUNC) &_nap_wCBMupdate, 8},
    {"_nap_wCBMupdateforR", (DL_FUNC) &_nap_wCBMupdateforR, 8},
    {"_nap_sampleWC", (DL_FUNC) &_nap_sampleWC, 5},
    {"_nap_ssaveC", (DL_FUNC) &_nap_ssaveC, 2},
    {"_nap_ssimC", (DL_FUNC) &_nap_ssimC, 2},
    {"_nap_napMCMC", (DL_FUNC) &_nap_napMCMC, 41},
    {"_nap_rowSumsC", (DL_FUNC) &_nap_rowSumsC, 1},
    {"_nap_allelesimCvec", (DL_FUNC) &_nap_allelesimCvec, 12},
    {"_nap_allelesimCmat", (DL_FUNC) &_nap_allelesimCmat, 13},
    {NULL, NULL, 0}
};

RcppExport void R_init_nap(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
