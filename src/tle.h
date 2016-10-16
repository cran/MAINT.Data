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
	private:
		int p_;
		VectorXd muE_;
		MatrixXd SigmaE_;
		double logLik_;		
};

#endif
