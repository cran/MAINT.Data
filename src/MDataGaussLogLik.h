#ifndef GAUSSLOGLIKH
#define GAUSSLOGLIKH

#include<vector>
#include <RcppEigen.h>

using namespace Eigen;

const double LN2PI = log(2*PI);

// void MDataGaussLogLik(const int n, const int p, const MatrixXd& X, const VectorXd& u, 
//                      const MatrixXd& Sigma, std::vector<double>& res);

void MDataGaussLogLik(const int n, const int p, const int Config, const MatrixXd& X, 
                      const VectorXd& u, const MatrixXd& Sigma, std::vector<double>& res);


#endif
