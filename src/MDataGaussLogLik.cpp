#include<cmath>
#include <RcppEigen.h>
#include "MDataGaussLogLik.h"
#include "EigenTmplFoo.cpp"

using namespace Eigen ;

void MDataGaussLogLik(const int n, const int p, const int Config, const MatrixXd& X, 
                      const VectorXd& u, const MatrixXd& Sigma, std::vector<double>& res)
{
  static const double PenF = 1e6;	           // penalty factor for numerically singular covariance matrices	
  double lndetSig, Singviol,penalty;
 
  if (Config==1) {

    MatrixXd dev,tmp;
    if (dev.rows()!=p || dev.cols()!=n)  dev.resize(p,n);
    if (tmp.rows()!=p || tmp.cols()!=n)  tmp.resize(p,n);
    for (int obs=0;obs<n;++obs) for (int j=0;j<p;++j) dev(j,obs) = X(obs,j) - u(j);

    if (!safepdsolve<MatrixXd,MatrixXd>(Sigma, dev, tmp, lndetSig, Singviol, MINLNDET, MINLNK2)) {
      penalty = -PenF*Singviol;
      for (int obs=0;obs<n;++obs) res[obs] = penalty;
      return;
    }  

    double c0 = -0.5 * ( p*LN2PI + lndetSig );
//    for (int obs=0;obs<n;++obs) {
//      res[obs] = 0.;
//      for (int j=0;j<p;++j) res[obs] -= dev(j,obs)*tmp(j,obs);
//      res[obs] /= 2.;
//      res[obs] += c0;      
//    }
     for (int obs=0;obs<n;++obs) res[obs] = c0 - static_cast<double>(dev.col(obs).transpose()*tmp.col(obs))/2;

    return;
  }
 
  else if (Config==3) {

    int q = p/2;
    double a,b,c,det,lndet,dev1,dev2;
 
    double c0 = -0.5 * p*LN2PI;
    for (int obs=0;obs<n;++obs) res[obs] = c0;

    for (int j=0;j<q;++j) {

      a = Sigma(j,j);  
      b = Sigma(q+j,q+j);  
      c = Sigma(j,q+j);  
      det = a*b - c*c;

      if (det < std::numeric_limits<double>::min()) {  
        penalty = -std::numeric_limits<double>::infinity();
        for (int obs=0;obs<n;++obs) res[obs] = penalty;
        return;
      }
      lndet = log(det); 
      if (lndet < MINLNDET) {
        penalty = -PenF*(MINLNDET-lndet);
        for (int obs=0;obs<n;++obs) res[obs] = penalty;
        return;
      }

      for (int obs=0;obs<n;++obs) {
        dev1 = X(obs,j) - u(j);  
        dev2 = X(obs,q+j) - u(q+j);  
        res[obs] -= ( lndet + (a*dev2*dev2 + b*dev1*dev1 -2*c*dev1*dev2)/det ) / 2;
      }
    }

    return;    
  }

  else if (Config==4) {

    MatrixXd dev,tmp;
    int q = p/2;
    if (dev.rows()!=q || dev.cols()!=n)  dev.resize(q,n);
    if (tmp.rows()!=q || tmp.cols()!=n)  tmp.resize(q,n);
    bool singularMat;
//    double resinc;
    double c0 = -0.5 * p*LN2PI;
    for (int obs=0;obs<n;++obs) res[obs] = c0;

    for (int block=0;block<2;++block) {   // block==0: MidPoints ; block==1: LogRanges  

      if (block==0) {
        for (int obs=0;obs<n;++obs) for (int j=0;j<q;++j) dev(j,obs) = X(obs,j) - u(j);
        singularMat = !safepdsolve<MatrixXd,MatrixXd>(Sigma.block(0,0,q,q), dev, tmp, lndetSig, Singviol, MINLNDET, MINLNK2);
      } else {
        for (int obs=0;obs<n;++obs) for (int j=0;j<q;++j) dev(j,obs) = X(obs,q+j) - u(q+j);
        singularMat = !safepdsolve<MatrixXd,MatrixXd>(Sigma.block(q,q,q,q), dev, tmp, lndetSig, Singviol, MINLNDET, MINLNK2);
      }
 
      if (singularMat) {
        penalty = -PenF*Singviol;
        for (int obs=0;obs<n;++obs) res[obs] = penalty;
        return;
      }  

//      for (int obs=0;obs<n;++obs) {
//        resinc = -lndetSig;
//        for (int j=0;j<q;++j) resinc -= dev(j,obs)*tmp(j,obs);
//        resinc /= 2.;
//        res[obs] += resinc;      
//      }
       for (int obs=0;obs<n;++obs) res[obs] -= ( lndetSig + static_cast<double>(dev.col(obs).transpose()*tmp.col(obs)) ) / 2;
    }

    return;
  }

  else if (Config==5) {

    double det,lndet,dev;
    double c0 = -0.5 * p*LN2PI;
    for (int obs=0;obs<n;++obs) res[obs] = c0;

    for (int j=0;j<p;++j) {

      det = Sigma(j,j);

      if (det < std::numeric_limits<double>::min()) {  
        penalty = -std::numeric_limits<double>::infinity();
        for (int obs=0;obs<n;++obs) res[obs] = penalty;
        return;
      }
      lndet = log(det); 
      if (lndet < MINLNDET) {
        penalty = -PenF*(MINLNDET-lndet);
        for (int obs=0;obs<n;++obs) res[obs] = penalty;
        return;
      }

      for (int obs=0;obs<n;++obs) {
        dev = X(obs,j) - u(j);  
        res[obs] -= (lndet + dev*dev/det) / 2;
      }
    }

    return;    
  }


}  


