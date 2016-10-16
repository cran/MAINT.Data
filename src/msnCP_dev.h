#ifndef _RCppMDt_msnCP_dev_H
#define _RCppMDt_msnCP_dev_H

#include <Rcpp.h>

using namespace Rcpp ;

const double b = sqrt(2./PI);

void Rprintivctzzz(const int p,const IntegerVector& v);

template<class VT>
void Rprintv(const int p,const VT& v)
{
	for (int i=0;i<p;i++)  Rprintf("%f ",v(i));
	Rprintf("\n");
} 

template<class MT>
void RprintM(const int m,const int n,const MT& M)
{
	for (int r=0;r<m;r++) {
		for (int c=0;c<n;c++) Rprintf("%f ",M(r,c));
		Rprintf("\n");
	}
	Rprintf("\n");
}

inline int ncovp(int Config,int q,int p)
{
	switch (Config)  {
		case 1: return p*(p+1)/2;
		case 2: return p+q+q*(q-1);
		case 3: return p+q;
		case 4: return p+q*(q-1);
		case 5: return p;
	}
  return 0;
}

RcppExport SEXP msnCP_dev(SEXP param_s, SEXP y_s, SEXP grpind_s, SEXP Config_s, SEXP n_s, SEXP p_s, SEXP k_s, 
		SEXP trace_s, SEXP c2tol_s, SEXP ldRtol_s, SEXP PenF_s, SEXP PenC_s, 
		SEXP nopenalty_s, SEXP MachineEPS_s) ;

RcppExport SEXP msnCP_dev_grad(SEXP param_s, SEXP y_s, SEXP grpind_s, SEXP Config_s, SEXP n_s, SEXP p_s, SEXP k_s, 
		SEXP trace_s, SEXP c2tol_s, SEXP ldRtol_s, SEXP beta0tol_s, SEXP PenF_s, SEXP MachineEPS_s) ;

#endif
