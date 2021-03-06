// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// BMrecode
bool BMrecode(SEXP A, double& a, const double& b);
RcppExport SEXP _nap_BMrecode(SEXP ASEXP, SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type A(ASEXP);
    Rcpp::traits::input_parameter< double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double& >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(BMrecode(A, a, b));
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
// rankC
arma::vec rankC(arma::vec x);
RcppExport SEXP _nap_rankC(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rankC(x));
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
double LIKELIHOOD(const arma::vec& y, const arma::vec& w, const double& b, const double& a, const double& p, const double& mu, const double& epi, bool verbose, bool printall);
RcppExport SEXP _nap_LIKELIHOOD(SEXP ySEXP, SEXP wSEXP, SEXP bSEXP, SEXP aSEXP, SEXP pSEXP, SEXP muSEXP, SEXP epiSEXP, SEXP verboseSEXP, SEXP printallSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const double& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double& >::type p(pSEXP);
    Rcpp::traits::input_parameter< const double& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const double& >::type epi(epiSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< bool >::type printall(printallSEXP);
    rcpp_result_gen = Rcpp::wrap(LIKELIHOOD(y, w, b, a, p, mu, epi, verbose, printall));
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
// ssimC
arma::vec ssimC(arma::uvec snps, double svar);
RcppExport SEXP _nap_ssimC(SEXP snpsSEXP, SEXP svarSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uvec >::type snps(snpsSEXP);
    Rcpp::traits::input_parameter< double >::type svar(svarSEXP);
    rcpp_result_gen = Rcpp::wrap(ssimC(snps, svar));
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
    {"_nap_BMrecode", (DL_FUNC) &_nap_BMrecode, 3},
    {"_nap_BMsimulate", (DL_FUNC) &_nap_BMsimulate, 1},
    {"_nap_BMwrite012", (DL_FUNC) &_nap_BMwrite012, 2},
    {"_nap_BMwritePED", (DL_FUNC) &_nap_BMwritePED, 2},
    {"_nap_BMsubset", (DL_FUNC) &_nap_BMsubset, 3},
    {"_nap_upperTmat", (DL_FUNC) &_nap_upperTmat, 1},
    {"_nap_Xmcenter", (DL_FUNC) &_nap_Xmcenter, 1},
    {"_nap_Xmvcenter", (DL_FUNC) &_nap_Xmvcenter, 1},
    {"_nap_LDrelative", (DL_FUNC) &_nap_LDrelative, 3},
    {"_nap_rankC", (DL_FUNC) &_nap_rankC, 1},
    {"_nap_accuracies", (DL_FUNC) &_nap_accuracies, 2},
    {"_nap_My", (DL_FUNC) &_nap_My, 2},
    {"_nap_BMs", (DL_FUNC) &_nap_BMs, 5},
    {"_nap_BMmgwa", (DL_FUNC) &_nap_BMmgwa, 3},
    {"_nap_LLGaussMix", (DL_FUNC) &_nap_LLGaussMix, 4},
    {"_nap_LIKELIHOOD", (DL_FUNC) &_nap_LIKELIHOOD, 9},
    {"_nap_wC", (DL_FUNC) &_nap_wC, 5},
    {"_nap_wCBM", (DL_FUNC) &_nap_wCBM, 7},
    {"_nap_ssimC", (DL_FUNC) &_nap_ssimC, 2},
    {"_nap_rowSumsC", (DL_FUNC) &_nap_rowSumsC, 1},
    {"_nap_allelesimCvec", (DL_FUNC) &_nap_allelesimCvec, 12},
    {"_nap_allelesimCmat", (DL_FUNC) &_nap_allelesimCmat, 13},
    {NULL, NULL, 0}
};

RcppExport void R_init_nap(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
