#include "msnCP_dev.h"
#include "msnCP_Aux.h"
#include <limits>

#include "EigenTmplFoo.cpp"
#include "AuxTmplFoo.cpp"

const double ln2pi = log(2.*PI);

template<class SQMATTP>
void cov2cor(const int p,const SQMATTP& S,SQMATTP& R)
{
	for (int r=0;r<p;r++)  {
		R(r,r) = 1.;
		for (int c=0;c<r;c++) {
			R(r,c) = S(r,c)/sqrt(S(r,r)*S(c,c));
			R(c,r) = R(r,c);
		}
	}
	return;
} 

template<class VCTTP,class SQMATTP>
void cnvCPtoDP(const int p,const NumericVector mu,const SQMATTP& Sigma,const NumericVector gamma1,
		VCTTP& ksi, SQMATTP& Omega, VCTTP& alpha, SQMATTP& Omegabar, VCTTP& delta,
		double* c2, bool* admissible, double* viol, const double tol,const bool FixedArrays)
{
	static VCTTP c,muz,omega,mu0,sigmaz,tmp;
	static SQMATTP mu0OtP;

	if (!FixedArrays)  {
		if (tmp.size()!=p) tmp.resize(p);
		if (c.size()!=p) c.resize(p);
		if (muz.size()!=p) muz.resize(p);
		if (omega.size()!=p) omega.resize(p);
		if (mu0.size()!=p) mu0.resize(p);
		if (sigmaz.size()!=p) sigmaz.resize(p);
		if (mu0OtP.rows()!=p || mu0OtP.cols()!=p)  mu0OtP.resize(p,p);
	}
	for (int i=0;i<p;i++) { 
		c(i) = pow(2*fabs(gamma1(i))/(4.-PI),1./3);
		if (gamma1(i)<0) c(i) *= -1;
		muz(i) = c(i)/sqrt(1.+c(i)*c(i));
		sigmaz(i) = sqrt(1.-muz(i)*muz(i));
		delta(i) = muz(i)/b;
		omega(i) = sqrt(Sigma(i,i))/sigmaz(i);
		mu0(i) = omega(i)*muz(i);
		ksi(i) = mu(i)-mu0(i);
	}
	outerprod<VCTTP,SQMATTP>(p,mu0,mu0OtP);
	Omega = Sigma+mu0OtP;
	cov2cor<SQMATTP>(p,Omega,Omegabar);

	if (!pdsolve<SQMATTP,VCTTP>(Omegabar,delta,tmp,NULL))  {
		for (int i=0;i<p;i++) alpha(i) = (double)NAN;
	  	*c2 = (double)NAN;
	  	*admissible = false;
	  	*viol = -Omegabar.determinant();
	  	return;
	}
	*c2 = 1.- delta.dot(tmp);
	if (*c2 < tol)  {
  	alpha = tmp/tol;
		*admissible = false;
		*viol = -*c2;
	}
	else  {
  	alpha = tmp/sqrt(*c2);
		*admissible = true;
		*viol = (double)NAN;
	}
	return;
} 

template<class VCTTP,class ROWVCTTP,class SQMATTP,class RCTMATTP>
double msnCP_dev1(NumericVector& param, const NumericMatrix& y, const IntegerVector& grpind, 
		const int Config, const int n, const int p, const int k, 
		const bool trace, const double c2tol, const double ldRtol, 
		const double PenF, const double PenC, const bool nopenalty,
		const double MachineEPS, const bool FixedArrays)
{
	double dbltmp,penalty,DPc2,DPviol;
	bool DPadmissible;
	int q(p/2), nSigmapar(ncovp(Config,q,p));
	NumericVector::iterator mu1ptr(param.begin());
	NumericVector::iterator beta2kptr(mu1ptr+p);
	NumericVector::iterator Sigmaptr(mu1ptr+k*p);
	NumericVector::iterator gamma1ptr(Sigmaptr+nSigmapar);

	static SQMATTP Sigma,OmegaInv,DPOmega, DPOmegabar;
	static RCTMATTP y0;
	static NumericMatrix beta2k;
	static VCTTP omega,alphoveromg,DPksi1,DPalpha,DPdelta;     
	static ROWVCTTP y0i;

	if (y0.rows()!=n || y0.cols()!=p) y0.resize(n,p);
	if (!FixedArrays)  {
		if (Sigma.rows()!=p || Sigma.cols()!=p)  Sigma.resize(p,p);
		if (OmegaInv.rows()!=p || OmegaInv.cols()!=p)  OmegaInv.resize(p,p);
		if (omega.size()!=p) omega.resize(p);
		if (alphoveromg.size()!=p) alphoveromg.resize(p);
		if (y0i.size()!=p) y0i.resize(p);
		if (DPksi1.size()!=p) DPksi1.resize(p);
		if (DPOmega.rows()!=p || DPOmega.cols()!=p)  DPOmega.resize(p,p);
		if (DPalpha.size()!=p) DPalpha.resize(p);
		if (DPOmegabar.rows()!=p || DPOmegabar.cols()!=p)  DPOmegabar.resize(p,p);
		if (DPdelta.size()!=p) DPdelta.resize(p);
	}
	Sigma = RestCov<SQMATTP>(q,Sigmaptr,Config,FixedArrays);

	double ldR(log(Sigma.determinant()));
	if ( !(nopenalty) && ldR < ldRtol ) {
    dbltmp = ldRtol-ldR;
    penalty = PenF * (PenC+dbltmp*dbltmp);
	}
  else penalty = 0.;
	NumericVector mu1(mu1ptr,mu1ptr+p);
  if (k>1)  beta2k = NumericMatrix(k-1,p,beta2kptr);
	NumericVector gamma1(gamma1ptr,gamma1ptr+p);

	cnvCPtoDP<VCTTP,SQMATTP>(p,mu1,Sigma,gamma1,
		DPksi1,DPOmega,DPalpha,DPOmegabar,DPdelta,&DPc2,&DPadmissible,&DPviol,
		MachineEPS,FixedArrays);

	if (!nopenalty && DPc2 < c2tol) {
    dbltmp = c2tol-DPc2; 
    penalty += PenF * (PenC + dbltmp*dbltmp); 
	}
  if (!DPadmissible)  {
     		if (nopenalty) return INFINITY;  
     		else return penalty;
	}       

	double logDet;
	if (!pdsolve<SQMATTP>(DPOmega,OmegaInv,&logDet)) return INFINITY;
	for (int i=0;i<p;i++) {
		omega(i) = sqrt(DPOmega(i,i));
		alphoveromg(i) = DPalpha(i)/omega(i);
	}
	for (int r=0;r<n;r++) for(int c=0;c<p;c++) 
		if (grpind(r)<0) y0(r,c) = y(r,c) - DPksi1(c); 
		else y0(r,c) = y(r,c) - DPksi1(c) - beta2k(grpind(r),c); 

	double dev(n*(p*ln2pi+logDet));
	for (int obs=0;obs<n;obs++) {
		y0i = y0.row(obs);
		dev += y0i.dot(OmegaInv*y0i.transpose()) - 2*zeta(0,y0i.dot(alphoveromg));
	}

	if(trace) { 
    		Rprintf("msnCP.dev %f\n",dev);
    		Rprintf("Centred parameters:\n");
      		Rprintf("mu1 = ") ; Rprintv<NumericVector>(p,mu1);
                if (k>1) { Rprintf("beta2k =") ; RprintM<NumericMatrix>(k-1,p,beta2k); }
    		Rprintf("gamma1 = ") ; Rprintv<NumericVector>(p,gamma1);
    		Rprintf("Sigma =\n") ; RprintM<SQMATTP>(p,p,Sigma);
    		Rprintf("Direct parameters:\n");
      		Rprintf("ksi1 = ") ; Rprintv<VCTTP>(p,DPksi1);
    		Rprintf("alpha = ") ; Rprintv<VCTTP>(p,DPalpha);
    		Rprintf("Omega =\n") ; RprintM<SQMATTP>(p,p,DPOmega);
	}

  	return dev+penalty;
}

