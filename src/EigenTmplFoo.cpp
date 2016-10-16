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
bool pdsolve(const SQMATTP& M,const VCTTP& b,VCTTP& res,double* logDet)
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

#endif




