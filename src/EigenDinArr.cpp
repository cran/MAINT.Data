#include <RcppEigen.h>

#include "msnCP_dev.cpp"
#include "msnCP_dev_grad.cpp"

using namespace Eigen ;

SEXP msnCP_dev(SEXP param_s, SEXP y_s, SEXP grpind_s, SEXP Config_s, SEXP n_s, SEXP p_s, SEXP k_s, 
  SEXP trace_s, SEXP c2tol_s, SEXP ldRtol_s, SEXP PenF_s, SEXP PenC_s,
	SEXP nopenalty_s, SEXP MachineEPS_s)
{
	int Config(as<int>(Config_s)), n(as<int>(n_s)), p(as<int>(p_s)), k(as<int>(k_s));
	double c2tol(as<double>(c2tol_s)), ldRtol(as<double>(ldRtol_s));
	double PenF(as<double>(PenF_s)), PenC(as<double>(PenC_s)), MachineEPS(as<double>(MachineEPS_s));
	bool trace(as<bool>(trace_s)), nopenalty(as<bool>(nopenalty_s));
	NumericVector param(param_s);
	NumericMatrix y(y_s);
	IntegerVector grpind(grpind_s);

	return wrap( msnCP_dev1<VectorXd,RowVectorXd,MatrixXd,MatrixXd>(
			param,y,grpind,Config,n,p,k,trace,c2tol,ldRtol,
			PenF,PenC,nopenalty,MachineEPS,false) );
}

SEXP msnCP_dev_grad(SEXP param_s, SEXP y_s, SEXP grpind_s, SEXP Config_s, SEXP n_s, SEXP p_s, SEXP k_s, 
	SEXP trace_s, SEXP c2tol_s, SEXP ldRtol_s, SEXP beta0tol_s, SEXP PenF_s,
	SEXP MachineEPS_s)
{
	int Config(as<int>(Config_s)), n(as<int>(n_s)), p(as<int>(p_s)), k(as<int>(k_s));
	double c2tol(as<double>(c2tol_s)), ldRtol(as<double>(ldRtol_s)), beta0tol(as<double>(beta0tol_s));
	double PenF(as<double>(PenF_s)), MachineEPS(as<double>(MachineEPS_s));
	bool trace(as<bool>(trace_s));
	NumericVector param(param_s);
	NumericMatrix y(y_s);
	IntegerVector grpind(grpind_s);

	return wrap( msnCP_dev_grad1<VectorXd,RowVectorXd,RowVectorXd,RowVectorXd,MatrixXd,MatrixXd>(
			param,y,grpind,Config,n,p,k,trace,c2tol,ldRtol,beta0tol,PenF,MachineEPS,false) 
		);
}



