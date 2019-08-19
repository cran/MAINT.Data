#ifndef EMGAUSSH
#define EMGAUSSH

#include <vector>
#include <Rcpp.h>
#include <numeric>

using namespace std;
using namespace Rcpp;

const double MINZ = 1e-300;

RcppExport
SEXP CEMGauss(SEXP X_s, SEXP k_s, SEXP Config_s, SEXP Homoc_s, SEXP maxiter_s, SEXP tautol_s, SEXP convtol_s,
  SEXP InitSolz_s, SEXP InitSoltau_s, SEXP InitSolmuk_s, SEXP InitSolSigma_s, SEXP InitSolSigmak_s, SEXP InitSolLnLik_s);       

#endif
