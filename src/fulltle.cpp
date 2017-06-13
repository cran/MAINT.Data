#include "tle.h"

double loglik(const NumericMatrix& X,const int n,const int p,const int Cnf,
	double c0,const int k,const vector<int>& Set,Estimate& tmpsol);

RcppExport
SEXP Cfulltle(SEXP X_s, SEXP n_s, SEXP p_s, SEXP k_s, SEXP Cnf_s, SEXP c0_s)
{
  int n(as<int>(n_s)), p(as<int>(p_s)), k(as<int>(k_s)), Cnf(as<int>(Cnf_s));
  NumericMatrix X(X_s);
  Estimate tmpsol(p);
  vector<int> tmpSet(k),bestSet(k);
  double c0(as<double>(c0_s)),tmploglik,bestloglik(-Inf);
	
  bool start(true);
    while (TRUE)  {
    if (start) {
      for (int i=0;i<k;++i) tmpSet[i] = i;
      start = false;
    }  else {
      unsigned j(k);
      bool cont(TRUE);
      if (tmpSet[k-1]<n-1) cont=FALSE;
      while (cont) {
        if (j==0) cont=FALSE;
        else if  (tmpSet[j-1]<n-(k-(j-1))) cont=FALSE;
        if (j>0 && cont) --j;
      }
      if (j==0) break;
      ++tmpSet[j-1];
      if (j-1<k-1) for (int j1=j;j1<k;j1++)
       tmpSet[j1] = tmpSet[j1-1]+1;
    }

    tmploglik = loglik(X,n,p,Cnf,c0,k,tmpSet,tmpsol);
    if (tmploglik > bestloglik)  {
      bestSet = tmpSet;
      bestloglik = tmploglik;
    }
  }
  return List::create(Named("LogLik")=bestloglik,Named("Set")=bestSet);
}

double parcovloglik(const MatrixXd& Xdev,const vector<int>& Set,MatrixXd& Sigma,const double c0);
double parcovloglik3(const MatrixXd& Xdev,const vector<int>& Set,MatrixXd& Sigma,const double c0);
double parcovloglik4(const MatrixXd& Xdev,const vector<int>& Set,MatrixXd& Sigma,const double c0);
double parcovloglik5(const MatrixXd& Xdev,const vector<int>& Set,MatrixXd& Sigma,const double c0);

double loglik(const NumericMatrix& X,const int n,const int p,const int Cnf,
	double c0,const int k,const vector<int>& Set,Estimate& tmpsol)
{
	static MatrixXd Xdev;
	if (Xdev.rows()!=n || Xdev.cols()!=p) Xdev.resize(n,p);

	parcolmeans(X,Set,tmpsol.muE());
	double muEj;
	for (int j=0;j<p;j++) {
		muEj = tmpsol.muE()(j);
		for (int i=0;i<n;i++) Xdev(i,j) = X(i,j) - muEj;
	}
	switch (Cnf)  {
		case 1:
			return parcovloglik(Xdev,Set,tmpsol.SigmaE(),c0);
			break;    
		case 3:
			return parcovloglik3(Xdev,Set,tmpsol.SigmaE(),c0);
			break;    
		case 4:
			return parcovloglik4(Xdev,Set,tmpsol.SigmaE(),c0);
			break;    
		case 5:
			return parcovloglik5(Xdev,Set,tmpsol.SigmaE(),c0);
			break;    
	}
  return 0.;
}

double parcovloglik(const MatrixXd& Xdev,const vector<int>& Set,MatrixXd& Sigma,const double c0)
{
	int n(Xdev.rows()), p(Xdev.cols()), k(Set.size());
	double psum;

	for (int j1=0;j1<p;j1++) for (int j2=j1;j2<p;j2++) {
		psum = 0.;
		for (int i=0;i<k;i++) psum += Xdev(Set[i],j1)*Xdev(Set[i],j2);
		Sigma(j1,j2) = psum/k;
		if (j2>j1) Sigma(j2,j1) = Sigma(j1,j2);
	}

	static LLT<MatrixXd> cholSigma;
	static MatrixXd SigmaSr;
	if (SigmaSr.rows()!=p || SigmaSr.cols()!=p) SigmaSr.resize(p,p);

	cholSigma.compute(Sigma);
	SigmaSr = cholSigma.matrixL();
	double lnSigSrdet(0.);
	for (int j=0;j<p;++j) lnSigSrdet += log(SigmaSr(j,j));
	return c0-n*lnSigSrdet; 
}

double parcovloglik3(const MatrixXd& Xdev,const vector<int>& Set,MatrixXd& Sigma,const double c0)
{
	int n(Xdev.rows()), p(Xdev.cols()), q(p/2), k(Set.size());

	Sigma.setZero(p,p);
	double psumMP,psumLR,psumMPLG;
	for (int j=0;j<q;j++) {
		psumMP = psumLR = psumMPLG = 0.;
		for (int i=0;i<k;i++) {
			int i1(Set[i]);
			psumMP += Xdev(i1,j)*Xdev(i1,j);
			psumLR += Xdev(i1,q+j)*Xdev(i1,q+j);
			psumMPLG += Xdev(i1,j)*Xdev(i1,q+j);
		}
		Sigma(j,j) = psumMP/k;
		Sigma(q+j,q+j) = psumLR/k;
		Sigma(j,q+j) = Sigma(q+j,j) = psumMPLG/k;
	}

	double b,lnSigdet(0.);
	for (int j=0;j<q;j++) {
		b = Sigma(j,q+j);
		lnSigdet += log(Sigma(j,j)*Sigma(q+j,q+j)-b*b)/2; 
	}
	return c0-n*lnSigdet/2; 
}

double parcovloglik4(const MatrixXd& Xdev,const vector<int>& Set,MatrixXd& Sigma,const double c0)
{
	int n(Xdev.rows()), p(Xdev.cols()), q(p/2);
	static MatrixXd pSigma;
	if (pSigma.rows()!=q || pSigma.cols()!=q) pSigma.resize(q,q);
	double loglik(c0);

	Sigma.setZero(p,p);
	loglik += parcovloglik(Xdev.block(0,0,n,q),Set,pSigma,c0/2);
	Sigma.block(0,0,q,q) = pSigma;
	loglik += parcovloglik(Xdev.block(0,q,n,q),Set,pSigma,c0/2);
	Sigma.block(q,q,q,q) = pSigma;

	return loglik;
}

double parcovloglik5(const MatrixXd& Xdev,const vector<int>& Set,MatrixXd& Sigma,const double c0)
{
	int n(Xdev.rows()), p(Xdev.cols()), k(Set.size());

	Sigma.setZero(p,p);
	double psum,Xdevij,lnSigdet(0.);
	for (int j=0;j<p;++j) {
 		psum = 0.;
		for (int i=0;i<k;++i) {
			Xdevij = Xdev(Set[i],j);
			psum += Xdevij*Xdevij;
		}
		lnSigdet += log(Sigma(j,j)=psum/k);
	}
	return c0-n*lnSigdet/2; 
}

