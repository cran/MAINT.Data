#ifndef _EigenTmplFoo_cpp
#define _EigenTmplFoo_cpp

#include <RcppEigen.h>

using namespace Eigen ;

template<class SQMATTP> void SetIdentity(SQMATTP& Ip,int p,int* Ipdim)
{
	if (Ip.rows()!=p || Ip.cols()!=p) Ip.resize(p,p);
	for (int r=0;r<p;r++)  {
		Ip(r,r) = 1.;
		for (int c=0;c<r;c++) Ip(r,c) = Ip(c,r) = 0.;
	}
	*Ipdim = p;
}

template<class SQMATTP>
bool pdsolve(const SQMATTP& M,SQMATTP& MInv,double* logDet)
{
	/* static */ LDLT<SQMATTP> MSr;
	/* static */ SQMATTP Ip;
	/* static */ int Ipdim(0);
	int p(M.rows());

	if (Ipdim!=p) SetIdentity(Ip,p,&Ipdim);

	MSr = M.ldlt();      // Maybe replace this by a better rank-reaveling criteron
	if (MSr.info()!=Success) return false;
	MInv = MSr.solve(Ip);
	if (logDet) {
		Diagonal<const SQMATTP> MSrdiag(MSr.vectorD());      
		*logDet = log(MSrdiag(0));
		for (int i=1;i<M.rows();i++) *logDet += log(MSrdiag(i));
	}
	return true;
}

template<class SQMATTP,class VCTTP>
bool pdsolve(const SQMATTP& M, const VCTTP& b, VCTTP& res, double* logDet)
{
	/* static */ LDLT<SQMATTP> MSr;

	MSr = M.ldlt();      // Maybe replace this by a better rank-reaveling criteron
	if (MSr.info()!=Success) return false;
	res = MSr.solve(b);
	if (logDet) {
		Diagonal<const SQMATTP> MSrdiag(MSr.vectorD());      
		*logDet = log(MSrdiag(0));
		for (int i=1;i<M.rows();i++) *logDet += log(MSrdiag(i));
	}
	return true;
}

template<class VCTTP> void SetZero(VCTTP& v,const int n,bool cheksize)
{ 
	if (cheksize && v.size()!=n) v.resize(n); 
	for (int i=0;i<n;i++) v(i) = 0.;
}

template<class MATTP> void SetZero(MATTP& M,const int m,const int n,bool cheksize)
{ 
	if ( cheksize && (M.rows()!=m || M.cols()!=n) ) M.resize(m,n); 
	for (int r=0;r<m;r++) for (int c=0;c<n;c++) M(r,c) = 0.;
}

template<class SQMATTP>
bool chcksing(const SQMATTP& M, double& logDet, double& viol, const double minlndet, const double minlnk2)
{
  double det = M.determinant();
// Rprintf("chcksing 1 -- det = %g std::numeric_limits<double>::min() = %g\n",det,std::numeric_limits<double>::min());
  if (det < std::numeric_limits<double>::min()) {
    viol = std::numeric_limits<double>::infinity();
    return false;
  }  
  logDet = log(det);
// Rprintf("chcksing 2 -- logDet = %g minlndet = %g\n",logDet,minlndet);
  if (logDet < minlndet) {
    viol = minlndet-logDet;
    return false;
  }  
  double trace = M.trace();
// Rprintf("chcksing 3 -- trace = %g std::numeric_limits<double>::min() = %g\n",trace,std::numeric_limits<double>::min());
  if (trace < std::numeric_limits<double>::min()) {
    viol = std::numeric_limits<double>::infinity();
    return false;
  }  

  double logtrace = log(trace);
  int p(M.rows());
// Rprintf("chcksing 4 -- p*logtrace-logDet = %g minlnk2 = %g\n",p*logtrace-logDet,minlnk2);
  if (p*logtrace-logDet > minlnk2) {
    SelfAdjointEigenSolver<SQMATTP> VLV(M,EigenvaluesOnly);
    VectorXd egval(VLV.eigenvalues()); 
    if (egval(0) < std::numeric_limits<double>::epsilon()) {
      viol = std::numeric_limits<double>::infinity();
      return false;
    }  
    double logk2 = log(egval(p-1)/egval(0)); 
// Rprintf("chcksing 5 -- logk2 = %g minlnk2 = %g\n",logk2,minlnk2);
    if (logk2 > minlnk2) {
      viol = logk2 - minlnk2;
      return false;
    }  
  } 
  return true;
}

const double MINLNK2 = log(1e15);     //  maximum (L2 norm) condition number for a matrix to be considered numerically non-singular 
const double MINLNDET = -500;         //  minimum log-determinant for a matrix to be considered numerically non-singular 

template<class SQMATTP>
bool safepdsolve(const SQMATTP& M, SQMATTP& MInv, double& logDet, double& viol,  
  const double minlndet=MINLNDET, const double minlnk2=MINLNK2)
{
  if ( !chcksing(M, logDet, viol, minlndet, minlnk2) )  
     return false;
  
  LLT<SQMATTP> MSr = M.llt(); 
   if (MSr.info()!=Success) {
     viol =  std::numeric_limits<double>::infinity();
     return false;
   }
   SQMATTP Ip;
   static int Ipdim(0);
   int p(M.rows());
   if (Ipdim!=p) SetIdentity(Ip,p,&Ipdim);
   MInv = MSr.solve(Ip);
   return true; 
}

template<class SQMATTP, class RCTMATTP>
bool safepdsolve(const SQMATTP& M, const RCTMATTP& rhs, RCTMATTP& res, double& logDet, double& viol,  
  const double minlndet=MINLNDET, const double minlnk2=MINLNK2)
{
   if ( !chcksing(M, logDet, viol, minlndet, minlnk2) )
     return false;

// Rprintf("safepdsolve 1 -- Got here\n");
   
   LLT<SQMATTP> MSr = M.llt();  
   if (MSr.info()!=Success) {
     viol =  std::numeric_limits<double>::infinity();
     return false;
   }

// Rprintf("safepdsolve 2 -- Got here\n");

   res = MSr.solve(rhs);
   return true; 
}

#endif




