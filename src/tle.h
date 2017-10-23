#ifndef TLEH
#define TLEH

#include <vector>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <numeric>

using namespace std;
using namespace Rcpp;
using namespace Eigen;

const double Inf = DBL_MAX;
void parcolmeans(const NumericMatrix& X,const vector<int>& Set,VectorXd& res);

class Estimate {
	public: 
		Estimate(int p) : p_(p) { muE_.resize(p); SigmaE_.resize(p,p); };
		int p(void) { return p_; }
		VectorXd& muE(void) { return muE_; }
		MatrixXd& SigmaE(void) { return SigmaE_; }
		double logLik(void) { return logLik_; }
		void setlogLik(double LogL) { logLik_ = LogL; } 
	private:
		int p_;
		VectorXd muE_;
		MatrixXd SigmaE_;
		double logLik_;		
};

RcppExport SEXP Cfasttle(SEXP X_s, SEXP n_s, SEXP p_s, SEXP Poolm_s, SEXP m_s, SEXP kdblstar_s, SEXP k_s, SEXP nrep_s,
	SEXP Cnf_s, SEXP c0_s, SEXP maxrefstps_s, SEXP ClctSt_s);

RcppExport SEXP Cfulltle(SEXP X_s, SEXP n_s, SEXP p_s, SEXP k_s, SEXP Cnf_s, SEXP c0_s);



#endif
