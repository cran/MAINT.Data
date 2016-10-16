#ifndef _RCppMDt_AuxTmplFoo_cpp
#define _RCppMDt_AuxTmplFoo_cpp

#include <Rcpp.h>
#include "EigenTmplFoo.cpp"

using namespace Rcpp ;

template<class VCTTP,class MATTP>
void outerprod(const int p,const VCTTP& v1,const VCTTP& v2,MATTP& res)
{
	for (int r=0;r<p;r++) for(int c=0;c<p;c++) res(r,c) = v1(r)*v2(c);
	return;
}

template<class VCTTP,class MATTP>
void outerprod(const int p,const VCTTP& v,MATTP& res)
{
	for (int r=0;r<p;r++) for(int c=0;c<=r;c++) {
		res(r,c) = v(r)*v(c);
		if (c<r) res(c,r) = res(r,c);
	}
	return;
}

template<class SQMATTP>
SQMATTP RestCov(const int q,NumericVector::iterator xpos,const int Config,const bool FixedArrays)
{
	int p(2*q);
	/* static */	SQMATTP A,Sigma;
/*	int NullMatdim(0);
	SQMATTP A,Sigma,NullppMat;
	// if (NullMatdim.size()!=p) NullppMat.resize(p,p);
	if (NullMatdim!=p) { 
		SetZero<SQMATTP>(NullppMat,p,p,!FixedArrays);
 		NullMatdim = p;
	}
*/
 	switch (Config)  {
		case 1: 
//			if (A.size()!=p) A.resize(p,p);
//			A = NullppMat;  			
			SetZero<SQMATTP>(A,p,p,!FixedArrays);
			for (int i=0;i<p;i++) A(i,i) = *xpos++;
			for (int c=1;c<p;c++) for (int r=0;r<c;r++) A(r,c) = *xpos++;
			return A.transpose() * A;
		case 2: 
//			if (A.size()!=p) A.resize(p,p);
//			A = NullppMat;  			
			SetZero<SQMATTP>(A,p,p,!FixedArrays);
 			for (int i=0;i<p;i++) A(i,i) = *xpos++;
  			for (int c=1;c<p;c++) for (int r=0;r<c;r++)  {
    				if ( (r<q && c<q)  || (r>=q && c>=q) || c==r+q )  A(r,c) = *xpos++;
    				else if (r>0)  {
					double dbltmp = 0.;
					for (int i=0;i<r;i++) dbltmp -= A(i,r)*A(i,c);   
      					A(r,c) = dbltmp/A(r,r);
				}
  			}  
			return A.transpose() * A;
		case 3: 
			int qplsi;
			double A11,A12,A22;
//			if (Sigma.size()!=p) Sigma.resize(p,p);
//			Sigma = NullppMat;
			SetZero<SQMATTP>(Sigma,p,p,!FixedArrays);
			for (int i=0;i<q;i++)  {
				A11 = *(xpos+i);
				A12 = *(xpos+p+i);
				A22 = *(xpos+(qplsi=q+i));
				Sigma(i,i) = A11*A11;
				Sigma(i,qplsi) = Sigma(qplsi,i) = A11*A12;  
				Sigma(qplsi,qplsi) = A12*A12 + A22*A22;
			}
			return Sigma;  
		case 4: 
			if (A.size()!=q) {
				A.resize(q,q);
  				for (int c=0;c<q-1;c++) for (int r=c+1;r<q;r++)  A(r,c) = 0.;
			}
//			if (Sigma.size()!=p) Sigma.resize(p,p);
//			Sigma = NullppMat;
			SetZero<SQMATTP>(Sigma,p,p,!FixedArrays);
 			for (int i=0;i<q;i++) A(i,i) = *(xpos+i);
  			for (int c=1,xind=p;c<q;c++) for (int r=0;r<c;r++,xind++)  A(r,c) = *(xpos+xind);
			Sigma.block(0,0,q,q) = A.transpose() * A; 
 			for (int i=0;i<q;i++) A(i,i) = *(xpos+q+i);
  			for (int c=1,xind=p+q*(q-1)/2;c<q;c++) for (int r=0;r<c;r++,xind++)  A(r,c) = *(xpos+xind);
			Sigma.block(q,q,q,q) = A.transpose() * A;
			return Sigma;
		case 5: 
			double x;
//			if (Sigma.size()!=p) Sigma.resize(p,p);
//			Sigma = NullppMat;
			SetZero<SQMATTP>(Sigma,p,p,!FixedArrays);
			for (int i=0;i<p;i++) {
        x = *xpos++;
        Sigma(i,i) = x*x;
			}  
			return Sigma;
	}
  return Sigma;
}

#endif
