// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// gethcpp
List gethcpp(int nc1, int cc, int nmix, int nk, int jj, int mm, const IntegerVector PIA, const IntegerVector cumss, const NumericVector Tsk, const NumericVector hk);
RcppExport SEXP _openCR_gethcpp(SEXP nc1SEXP, SEXP ccSEXP, SEXP nmixSEXP, SEXP nkSEXP, SEXP jjSEXP, SEXP mmSEXP, SEXP PIASEXP, SEXP cumssSEXP, SEXP TskSEXP, SEXP hkSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nc1(nc1SEXP);
    Rcpp::traits::input_parameter< int >::type cc(ccSEXP);
    Rcpp::traits::input_parameter< int >::type nmix(nmixSEXP);
    Rcpp::traits::input_parameter< int >::type nk(nkSEXP);
    Rcpp::traits::input_parameter< int >::type jj(jjSEXP);
    Rcpp::traits::input_parameter< int >::type mm(mmSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type PIA(PIASEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type cumss(cumssSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type Tsk(TskSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type hk(hkSEXP);
    rcpp_result_gen = Rcpp::wrap(gethcpp(nc1, cc, nmix, nk, jj, mm, PIA, cumss, Tsk, hk));
    return rcpp_result_gen;
END_RCPP
}
// makegkcpp
List makegkcpp(int cc, int kk, int mm, int detectfn, int sigmai, const NumericMatrix openval, const NumericMatrix traps, const NumericMatrix mask);
RcppExport SEXP _openCR_makegkcpp(SEXP ccSEXP, SEXP kkSEXP, SEXP mmSEXP, SEXP detectfnSEXP, SEXP sigmaiSEXP, SEXP openvalSEXP, SEXP trapsSEXP, SEXP maskSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type cc(ccSEXP);
    Rcpp::traits::input_parameter< int >::type kk(kkSEXP);
    Rcpp::traits::input_parameter< int >::type mm(mmSEXP);
    Rcpp::traits::input_parameter< int >::type detectfn(detectfnSEXP);
    Rcpp::traits::input_parameter< int >::type sigmai(sigmaiSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type openval(openvalSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type traps(trapsSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type mask(maskSEXP);
    rcpp_result_gen = Rcpp::wrap(makegkcpp(cc, kk, mm, detectfn, sigmai, openval, traps, mask));
    return rcpp_result_gen;
END_RCPP
}
// makegkParallelcpp
List makegkParallelcpp(int detectfn, int sigmai, int grain, const NumericMatrix& openval, const NumericMatrix& traps, const NumericMatrix& mask);
RcppExport SEXP _openCR_makegkParallelcpp(SEXP detectfnSEXP, SEXP sigmaiSEXP, SEXP grainSEXP, SEXP openvalSEXP, SEXP trapsSEXP, SEXP maskSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type detectfn(detectfnSEXP);
    Rcpp::traits::input_parameter< int >::type sigmai(sigmaiSEXP);
    Rcpp::traits::input_parameter< int >::type grain(grainSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type openval(openvalSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type traps(trapsSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type mask(maskSEXP);
    rcpp_result_gen = Rcpp::wrap(makegkParallelcpp(detectfn, sigmai, grain, openval, traps, mask));
    return rcpp_result_gen;
END_RCPP
}
// PCH1cpp
NumericVector PCH1cpp(int type, int x, int nc, int jj, const IntegerVector cumss, int nmix, const NumericMatrix openval0, const IntegerVector PIA0, const IntegerVector PIAJ, const NumericVector intervals);
RcppExport SEXP _openCR_PCH1cpp(SEXP typeSEXP, SEXP xSEXP, SEXP ncSEXP, SEXP jjSEXP, SEXP cumssSEXP, SEXP nmixSEXP, SEXP openval0SEXP, SEXP PIA0SEXP, SEXP PIAJSEXP, SEXP intervalsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< int >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type nc(ncSEXP);
    Rcpp::traits::input_parameter< int >::type jj(jjSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type cumss(cumssSEXP);
    Rcpp::traits::input_parameter< int >::type nmix(nmixSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type openval0(openval0SEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type PIA0(PIA0SEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type PIAJ(PIAJSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type intervals(intervalsSEXP);
    rcpp_result_gen = Rcpp::wrap(PCH1cpp(type, x, nc, jj, cumss, nmix, openval0, PIA0, PIAJ, intervals));
    return rcpp_result_gen;
END_RCPP
}
// PCH1secrcpp
NumericVector PCH1secrcpp(int type, bool individual, int x, int nc, int jj, const IntegerVector cumss, int kk, int mm, const NumericMatrix openval0, const IntegerVector PIA0, const IntegerVector PIAJ, const NumericVector gk0, int binomN, const NumericMatrix Tsk, const NumericVector intervals, const IntegerVector moveargsi, int movemodel, const CharacterVector usermodel, const IntegerMatrix kernel, const IntegerMatrix mqarray, double cellsize);
RcppExport SEXP _openCR_PCH1secrcpp(SEXP typeSEXP, SEXP individualSEXP, SEXP xSEXP, SEXP ncSEXP, SEXP jjSEXP, SEXP cumssSEXP, SEXP kkSEXP, SEXP mmSEXP, SEXP openval0SEXP, SEXP PIA0SEXP, SEXP PIAJSEXP, SEXP gk0SEXP, SEXP binomNSEXP, SEXP TskSEXP, SEXP intervalsSEXP, SEXP moveargsiSEXP, SEXP movemodelSEXP, SEXP usermodelSEXP, SEXP kernelSEXP, SEXP mqarraySEXP, SEXP cellsizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< bool >::type individual(individualSEXP);
    Rcpp::traits::input_parameter< int >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type nc(ncSEXP);
    Rcpp::traits::input_parameter< int >::type jj(jjSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type cumss(cumssSEXP);
    Rcpp::traits::input_parameter< int >::type kk(kkSEXP);
    Rcpp::traits::input_parameter< int >::type mm(mmSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type openval0(openval0SEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type PIA0(PIA0SEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type PIAJ(PIAJSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type gk0(gk0SEXP);
    Rcpp::traits::input_parameter< int >::type binomN(binomNSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type Tsk(TskSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type intervals(intervalsSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type moveargsi(moveargsiSEXP);
    Rcpp::traits::input_parameter< int >::type movemodel(movemodelSEXP);
    Rcpp::traits::input_parameter< const CharacterVector >::type usermodel(usermodelSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix >::type kernel(kernelSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix >::type mqarray(mqarraySEXP);
    Rcpp::traits::input_parameter< double >::type cellsize(cellsizeSEXP);
    rcpp_result_gen = Rcpp::wrap(PCH1secrcpp(type, individual, x, nc, jj, cumss, kk, mm, openval0, PIA0, PIAJ, gk0, binomN, Tsk, intervals, moveargsi, movemodel, usermodel, kernel, mqarray, cellsize));
    return rcpp_result_gen;
END_RCPP
}
// PCH0secrjcpp
NumericVector PCH0secrjcpp(int type, int x, int nc, int jj, const IntegerVector cumss, int kk, int mm, int cc0, const IntegerVector PIA0, const NumericVector gk0, int binomN, const NumericMatrix Tsk);
RcppExport SEXP _openCR_PCH0secrjcpp(SEXP typeSEXP, SEXP xSEXP, SEXP ncSEXP, SEXP jjSEXP, SEXP cumssSEXP, SEXP kkSEXP, SEXP mmSEXP, SEXP cc0SEXP, SEXP PIA0SEXP, SEXP gk0SEXP, SEXP binomNSEXP, SEXP TskSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< int >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type nc(ncSEXP);
    Rcpp::traits::input_parameter< int >::type jj(jjSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type cumss(cumssSEXP);
    Rcpp::traits::input_parameter< int >::type kk(kkSEXP);
    Rcpp::traits::input_parameter< int >::type mm(mmSEXP);
    Rcpp::traits::input_parameter< int >::type cc0(cc0SEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type PIA0(PIA0SEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type gk0(gk0SEXP);
    Rcpp::traits::input_parameter< int >::type binomN(binomNSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type Tsk(TskSEXP);
    rcpp_result_gen = Rcpp::wrap(PCH0secrjcpp(type, x, nc, jj, cumss, kk, mm, cc0, PIA0, gk0, binomN, Tsk));
    return rcpp_result_gen;
END_RCPP
}
// PCH1secrparallelcpp
NumericVector PCH1secrparallelcpp(int x, int type, int grain, bool individual, int jj, int mm, int nc, const IntegerVector cumss, const NumericMatrix openval0, const IntegerVector PIA0, const IntegerVector PIAJ, const NumericVector gk0, int binomN, const NumericMatrix Tsk, const NumericVector intervals, const IntegerVector moveargsi, int movemodel, const IntegerMatrix kernel, const IntegerMatrix mqarray, double cellsize);
RcppExport SEXP _openCR_PCH1secrparallelcpp(SEXP xSEXP, SEXP typeSEXP, SEXP grainSEXP, SEXP individualSEXP, SEXP jjSEXP, SEXP mmSEXP, SEXP ncSEXP, SEXP cumssSEXP, SEXP openval0SEXP, SEXP PIA0SEXP, SEXP PIAJSEXP, SEXP gk0SEXP, SEXP binomNSEXP, SEXP TskSEXP, SEXP intervalsSEXP, SEXP moveargsiSEXP, SEXP movemodelSEXP, SEXP kernelSEXP, SEXP mqarraySEXP, SEXP cellsizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< int >::type grain(grainSEXP);
    Rcpp::traits::input_parameter< bool >::type individual(individualSEXP);
    Rcpp::traits::input_parameter< int >::type jj(jjSEXP);
    Rcpp::traits::input_parameter< int >::type mm(mmSEXP);
    Rcpp::traits::input_parameter< int >::type nc(ncSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type cumss(cumssSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type openval0(openval0SEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type PIA0(PIA0SEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type PIAJ(PIAJSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type gk0(gk0SEXP);
    Rcpp::traits::input_parameter< int >::type binomN(binomNSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type Tsk(TskSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type intervals(intervalsSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type moveargsi(moveargsiSEXP);
    Rcpp::traits::input_parameter< int >::type movemodel(movemodelSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix >::type kernel(kernelSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix >::type mqarray(mqarraySEXP);
    Rcpp::traits::input_parameter< double >::type cellsize(cellsizeSEXP);
    rcpp_result_gen = Rcpp::wrap(PCH1secrparallelcpp(x, type, grain, individual, jj, mm, nc, cumss, openval0, PIA0, PIAJ, gk0, binomN, Tsk, intervals, moveargsi, movemodel, kernel, mqarray, cellsize));
    return rcpp_result_gen;
END_RCPP
}
// pradelloglikcpp
NumericVector pradelloglikcpp(int type, const IntegerVector w, int nc, int jj, int nmix, const NumericMatrix openval, const IntegerVector PIAJ, const NumericVector intervals);
RcppExport SEXP _openCR_pradelloglikcpp(SEXP typeSEXP, SEXP wSEXP, SEXP ncSEXP, SEXP jjSEXP, SEXP nmixSEXP, SEXP openvalSEXP, SEXP PIAJSEXP, SEXP intervalsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type w(wSEXP);
    Rcpp::traits::input_parameter< int >::type nc(ncSEXP);
    Rcpp::traits::input_parameter< int >::type jj(jjSEXP);
    Rcpp::traits::input_parameter< int >::type nmix(nmixSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type openval(openvalSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type PIAJ(PIAJSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type intervals(intervalsSEXP);
    rcpp_result_gen = Rcpp::wrap(pradelloglikcpp(type, w, nc, jj, nmix, openval, PIAJ, intervals));
    return rcpp_result_gen;
END_RCPP
}
// prwicpp
double prwicpp(int type, int n, int x, int nc, int jj, const IntegerVector cumss, int nmix, const IntegerVector w, const IntegerVector fi, const IntegerVector li, const NumericMatrix openval, const IntegerVector PIA, const IntegerVector PIAJ, const NumericVector intervals, int CJSp1);
RcppExport SEXP _openCR_prwicpp(SEXP typeSEXP, SEXP nSEXP, SEXP xSEXP, SEXP ncSEXP, SEXP jjSEXP, SEXP cumssSEXP, SEXP nmixSEXP, SEXP wSEXP, SEXP fiSEXP, SEXP liSEXP, SEXP openvalSEXP, SEXP PIASEXP, SEXP PIAJSEXP, SEXP intervalsSEXP, SEXP CJSp1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type nc(ncSEXP);
    Rcpp::traits::input_parameter< int >::type jj(jjSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type cumss(cumssSEXP);
    Rcpp::traits::input_parameter< int >::type nmix(nmixSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type w(wSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type fi(fiSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type li(liSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type openval(openvalSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type PIA(PIASEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type PIAJ(PIAJSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type intervals(intervalsSEXP);
    Rcpp::traits::input_parameter< int >::type CJSp1(CJSp1SEXP);
    rcpp_result_gen = Rcpp::wrap(prwicpp(type, n, x, nc, jj, cumss, nmix, w, fi, li, openval, PIA, PIAJ, intervals, CJSp1));
    return rcpp_result_gen;
END_RCPP
}
// allhistparallelcpp
NumericVector allhistparallelcpp(int x, int type, int nc, int CJSp1, int grain, const NumericMatrix pmix, const NumericVector intervals, const IntegerVector cumss, const IntegerVector w, const IntegerVector fi, const IntegerVector li, const NumericMatrix openval, const IntegerVector PIA, const IntegerVector PIAJ);
RcppExport SEXP _openCR_allhistparallelcpp(SEXP xSEXP, SEXP typeSEXP, SEXP ncSEXP, SEXP CJSp1SEXP, SEXP grainSEXP, SEXP pmixSEXP, SEXP intervalsSEXP, SEXP cumssSEXP, SEXP wSEXP, SEXP fiSEXP, SEXP liSEXP, SEXP openvalSEXP, SEXP PIASEXP, SEXP PIAJSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< int >::type nc(ncSEXP);
    Rcpp::traits::input_parameter< int >::type CJSp1(CJSp1SEXP);
    Rcpp::traits::input_parameter< int >::type grain(grainSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type pmix(pmixSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type intervals(intervalsSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type cumss(cumssSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type w(wSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type fi(fiSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type li(liSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type openval(openvalSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type PIA(PIASEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type PIAJ(PIAJSEXP);
    rcpp_result_gen = Rcpp::wrap(allhistparallelcpp(x, type, nc, CJSp1, grain, pmix, intervals, cumss, w, fi, li, openval, PIA, PIAJ));
    return rcpp_result_gen;
END_RCPP
}
// prwisecrcpp
double prwisecrcpp(int type, int n, int x, int nc, int jj, int kk, int mm, int nmix, const IntegerVector cumss, const IntegerVector w, const IntegerVector fi, const IntegerVector li, const NumericVector gk, const NumericMatrix openval, const IntegerVector PIA, const IntegerVector PIAJ, int binomN, const NumericMatrix Tsk, const NumericVector intervals, const IntegerVector moveargsi, const NumericMatrix h, const IntegerMatrix hindex, int CJSp1, int movemodel, const CharacterVector usermodel, const IntegerMatrix kernel, const IntegerMatrix mqarray, double cellsize);
RcppExport SEXP _openCR_prwisecrcpp(SEXP typeSEXP, SEXP nSEXP, SEXP xSEXP, SEXP ncSEXP, SEXP jjSEXP, SEXP kkSEXP, SEXP mmSEXP, SEXP nmixSEXP, SEXP cumssSEXP, SEXP wSEXP, SEXP fiSEXP, SEXP liSEXP, SEXP gkSEXP, SEXP openvalSEXP, SEXP PIASEXP, SEXP PIAJSEXP, SEXP binomNSEXP, SEXP TskSEXP, SEXP intervalsSEXP, SEXP moveargsiSEXP, SEXP hSEXP, SEXP hindexSEXP, SEXP CJSp1SEXP, SEXP movemodelSEXP, SEXP usermodelSEXP, SEXP kernelSEXP, SEXP mqarraySEXP, SEXP cellsizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type nc(ncSEXP);
    Rcpp::traits::input_parameter< int >::type jj(jjSEXP);
    Rcpp::traits::input_parameter< int >::type kk(kkSEXP);
    Rcpp::traits::input_parameter< int >::type mm(mmSEXP);
    Rcpp::traits::input_parameter< int >::type nmix(nmixSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type cumss(cumssSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type w(wSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type fi(fiSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type li(liSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type gk(gkSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type openval(openvalSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type PIA(PIASEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type PIAJ(PIAJSEXP);
    Rcpp::traits::input_parameter< int >::type binomN(binomNSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type Tsk(TskSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type intervals(intervalsSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type moveargsi(moveargsiSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type h(hSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix >::type hindex(hindexSEXP);
    Rcpp::traits::input_parameter< int >::type CJSp1(CJSp1SEXP);
    Rcpp::traits::input_parameter< int >::type movemodel(movemodelSEXP);
    Rcpp::traits::input_parameter< const CharacterVector >::type usermodel(usermodelSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix >::type kernel(kernelSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix >::type mqarray(mqarraySEXP);
    Rcpp::traits::input_parameter< double >::type cellsize(cellsizeSEXP);
    rcpp_result_gen = Rcpp::wrap(prwisecrcpp(type, n, x, nc, jj, kk, mm, nmix, cumss, w, fi, li, gk, openval, PIA, PIAJ, binomN, Tsk, intervals, moveargsi, h, hindex, CJSp1, movemodel, usermodel, kernel, mqarray, cellsize));
    return rcpp_result_gen;
END_RCPP
}
// allhistsecrparallelcpp
NumericVector allhistsecrparallelcpp(int x, int type, int mm, int nc, int binomN, int CJSp1, int grain, const NumericMatrix pmix, const NumericVector intervals, const IntegerVector cumss, const IntegerVector w, const IntegerVector fi, const IntegerVector li, const NumericVector gk, const NumericMatrix openval, const IntegerVector PIA, const IntegerVector PIAJ, const NumericMatrix Tsk, const NumericMatrix h, const IntegerMatrix hindex, int movemodel, const IntegerVector moveargsi, const IntegerMatrix kernel, const IntegerMatrix mqarray, double cellsize);
RcppExport SEXP _openCR_allhistsecrparallelcpp(SEXP xSEXP, SEXP typeSEXP, SEXP mmSEXP, SEXP ncSEXP, SEXP binomNSEXP, SEXP CJSp1SEXP, SEXP grainSEXP, SEXP pmixSEXP, SEXP intervalsSEXP, SEXP cumssSEXP, SEXP wSEXP, SEXP fiSEXP, SEXP liSEXP, SEXP gkSEXP, SEXP openvalSEXP, SEXP PIASEXP, SEXP PIAJSEXP, SEXP TskSEXP, SEXP hSEXP, SEXP hindexSEXP, SEXP movemodelSEXP, SEXP moveargsiSEXP, SEXP kernelSEXP, SEXP mqarraySEXP, SEXP cellsizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< int >::type mm(mmSEXP);
    Rcpp::traits::input_parameter< int >::type nc(ncSEXP);
    Rcpp::traits::input_parameter< int >::type binomN(binomNSEXP);
    Rcpp::traits::input_parameter< int >::type CJSp1(CJSp1SEXP);
    Rcpp::traits::input_parameter< int >::type grain(grainSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type pmix(pmixSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type intervals(intervalsSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type cumss(cumssSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type w(wSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type fi(fiSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type li(liSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type gk(gkSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type openval(openvalSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type PIA(PIASEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type PIAJ(PIAJSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type Tsk(TskSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type h(hSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix >::type hindex(hindexSEXP);
    Rcpp::traits::input_parameter< int >::type movemodel(movemodelSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type moveargsi(moveargsiSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix >::type kernel(kernelSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix >::type mqarray(mqarraySEXP);
    Rcpp::traits::input_parameter< double >::type cellsize(cellsizeSEXP);
    rcpp_result_gen = Rcpp::wrap(allhistsecrparallelcpp(x, type, mm, nc, binomN, CJSp1, grain, pmix, intervals, cumss, w, fi, li, gk, openval, PIA, PIAJ, Tsk, h, hindex, movemodel, moveargsi, kernel, mqarray, cellsize));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_openCR_gethcpp", (DL_FUNC) &_openCR_gethcpp, 10},
    {"_openCR_makegkcpp", (DL_FUNC) &_openCR_makegkcpp, 8},
    {"_openCR_makegkParallelcpp", (DL_FUNC) &_openCR_makegkParallelcpp, 6},
    {"_openCR_PCH1cpp", (DL_FUNC) &_openCR_PCH1cpp, 10},
    {"_openCR_PCH1secrcpp", (DL_FUNC) &_openCR_PCH1secrcpp, 21},
    {"_openCR_PCH0secrjcpp", (DL_FUNC) &_openCR_PCH0secrjcpp, 12},
    {"_openCR_PCH1secrparallelcpp", (DL_FUNC) &_openCR_PCH1secrparallelcpp, 20},
    {"_openCR_pradelloglikcpp", (DL_FUNC) &_openCR_pradelloglikcpp, 8},
    {"_openCR_prwicpp", (DL_FUNC) &_openCR_prwicpp, 15},
    {"_openCR_allhistparallelcpp", (DL_FUNC) &_openCR_allhistparallelcpp, 14},
    {"_openCR_prwisecrcpp", (DL_FUNC) &_openCR_prwisecrcpp, 28},
    {"_openCR_allhistsecrparallelcpp", (DL_FUNC) &_openCR_allhistsecrparallelcpp, 25},
    {NULL, NULL, 0}
};

RcppExport void R_init_openCR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
