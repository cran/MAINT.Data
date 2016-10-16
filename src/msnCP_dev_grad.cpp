#include <limits>
#include "msnCP_dev.h"

#include "AuxTmplFoo.cpp"
#include "EigenRestCovGrad.cpp"

const double b0 = 2./(4.-PI);
const double cubrootb0 = pow(b0,(1./3));

template<class VCTTP,class ROWVCTTP,class EXTVCTTP,class SQMATTP,class RCTMATTP>
void msnCP_ll_grad(const NumericVector& mu1, const NumericMatrix& beta2k, const SQMATTP& Sigma, 
		const NumericVector& gamma1, const NumericMatrix& y, const IntegerVector& grpind, 
		const int n, const int p, const int k, const int nvcovpar, const bool trace,
		const double PenF, const double ldRtol, const double c2tol, const double beta0tol,
		EXTVCTTP& CPgrad, const double Machinetol, const bool FixedArrays)
{
	int p2(p*p);

	static VCTTP sigma,mu0,SigmaImu0,mu0bar,ksi1;
	static ROWVCTTP y0sum,sumz1eta,y0i,z1eta;
	static EXTVCTTP penaltygrad;
	static SQMATTP SigmaI,Ip,D33,Dtld33,tmpM,tmpM2; // See if there is a diagonal matrix -- for Dtld33!!!
	static MatrixXd y0,yg0sum,D23,Dtld32,D32;     	// Try sparse matrix representations for Dtld32  !!!
	static int Ipdim(0);
	static VectorXi dvecind;
	static int dvecinddim(0);
	static SQMATTP Omega,OmegaInv,mu0OtP,sclMatI,OmgbI; 
	static VCTTP omega,bomega,delta,OmgbIdelta,tmpv,eta;
	static VCTTP ksi1grad,etagrad,gamma1grad;
	static RCTMATTP beta2kgrad;
	static EXTVCTTP Omegagrad,OmegaInvgrad1,Sigmagrad;
	static MatrixXd sumgz1eta,OmegaInvgrad2;
	static SQMATTP nOmgminusS0;

	if (dvecinddim!=p2)  {
		dvecind.resize(p2);
		for (int r=0,i1=0,i2=0;r<p;r++) for (int c=0;c<p;c++,i1++)
			if (c<r) dvecind(i1) = (int)NAN;
			else dvecind(i1) = i2++;
		dvecinddim = p2;
	}

	if (y0.rows()!=n || y0.cols()!=p) y0.resize(n,p);
	SetZero<MatrixXd>(D23,nvcovpar,p,true);
	SetZero<MatrixXd>(D32,p,nvcovpar,true);
	SetZero<MatrixXd>(Dtld32,p,nvcovpar,true);
	SetZero<SQMATTP>(Dtld33,p,p,!FixedArrays);
	if (Omegagrad.size()!=nvcovpar)  Omegagrad.resize(nvcovpar);
	if (OmegaInvgrad1.size()!=nvcovpar)  OmegaInvgrad1.resize(nvcovpar);
	if (OmegaInvgrad2.rows()!=nvcovpar || OmegaInvgrad2.cols()!=nvcovpar) 
		OmegaInvgrad2.resize(nvcovpar,nvcovpar);
	if (Sigmagrad.size()!=nvcovpar)  Sigmagrad.resize(nvcovpar);
	if (y0sum.size()!=p) y0sum.resize(p);
	SetZero<MatrixXd>(yg0sum,k,p,true);
	SetZero<MatrixXd>(sumgz1eta,k,p,true);

	if (!FixedArrays)  {
		if (sigma.size()!=p) sigma.resize(p);
		if (mu0.size()!=p) mu0.resize(p);
		if (SigmaI.rows()!=p || SigmaI.cols()!=p)  SigmaI.resize(p,p);
		if (SigmaImu0.size()!=p)  SigmaImu0.resize(p);
		if (mu0bar.size()!=p)  mu0bar.resize(p);
		if (Dtld33.rows()!=p || Dtld33.cols()!=p)  Dtld33.resize(p,p);
		if (ksi1.size()!=p)  ksi1.resize(p);
		if (Omega.rows()!=p || Omega.cols()!=p)  Omega.resize(p,p);
		if (OmegaInv.rows()!=p || OmegaInv.cols()!=p)  OmegaInv.resize(p,p);
		if (mu0OtP.rows()!=p || mu0OtP.cols()!=p)  mu0OtP.resize(p,p);
		if (omega.size()!=p)  omega.resize(p);
		if (bomega.size()!=p)  bomega.resize(p);
		if (delta.size()!=p)  delta.resize(p);
		if (OmgbIdelta.size()!=p)  OmgbIdelta.resize(p);
		if (sclMatI.rows()!=p || sclMatI.cols()!=p)  sclMatI.resize(p,p);
		if (OmgbI.rows()!=p || OmgbI.cols()!=p)  OmgbI.resize(p,p);
		if (tmpM.rows()!=p || tmpM.cols()!=p)  tmpM.resize(p,p);
		if (tmpM2.rows()!=p || tmpM2.cols()!=p)  tmpM2.resize(p,p);
		if (tmpv.size()!=p)  tmpv.resize(p);
		if (eta.size()!=p)  eta.resize(p);
		if (D33.rows()!=p || D33.cols()!=p)  D33.resize(p,p);
		if (ksi1grad.size()!=p)  ksi1grad.resize(p);
		if (etagrad.size()!=p)  etagrad.resize(p);
		if (beta2kgrad.rows()!=k-1 || beta2kgrad.cols()!=p)  beta2kgrad.resize(k-1,p);
		if (nOmgminusS0.rows()!=p || nOmgminusS0.cols()!=p)  nOmgminusS0.resize(p,p);
		if (y0i.size()!=p) y0i.resize(p);
		if (sumz1eta.size()!=p) sumz1eta.resize(p);
		if (gamma1grad.size()!=p) gamma1grad.resize(p);
		if (z1eta.size()!=p) z1eta.resize(p);
	}

	int ngradpar((k+1)*p+nvcovpar);
	SetZero<EXTVCTTP>(penaltygrad,ngradpar,!FixedArrays);
	if (Ipdim!=p) SetIdentity(Ip,p,&Ipdim);

	for (int i=0;i<p;i++) { 
		sigma(i) = sqrt(Sigma(i,i));
		mu0(i) = sigma(i)*pow(b0*fabs(gamma1(i)),1./3);
		if (gamma1(i)<0) mu0(i) *= -1;
	}

	double lRdet;
	if ( !pdsolve<SQMATTP>(Sigma,SigmaI,&lRdet) || lRdet < -std::numeric_limits<double>::max() )
	{
		warning("Non-positive definite covariance matrix found during gradient computations (which returned 0).\n");
		for (int i=0;i<ngradpar;i++) CPgrad(i) = 0.;
		return;
	}
	if (lRdet < ldRtol) {
		double cnst(PenF*(ldRtol-lRdet)); 
		for (int c=0,i=k*p;c<p;c++) {
			penaltygrad(i) = SigmaI(c,c) - 1./Sigma(c,c); 
			for (int r=c;r<p;r++,i++) {
				if (r!=c) penaltygrad(i) = 2*SigmaI(r,c); 
			 	penaltygrad(i) *= cnst;
			}
		}
	}
	SigmaImu0 = SigmaI * mu0;
	double beta02(mu0.dot(SigmaImu0));
	double beta0(sqrt(beta02));
	if (beta0>beta0tol) { 
		for (int i=0;i<p;i++) mu0bar(i) = mu0(i)/(sigma(i)*beta0);
		for (int c=0;c<p;c++) for(int j=0;j<p;j++)
			if (j==c) D23(dvecind(c*p+j),c) = 2*mu0(j);
			else if (j<c) D23(dvecind(j*p+c),c) = mu0(j);
			else D23(dvecind(c*p+j),c) = mu0(j);
		for (int i=0;i<p;i++)  {
			double mu0bi(mu0bar(i));
			Dtld32(i,dvecind(i*(p+1))) = beta0*mu0bi/(2*sigma(i));
			Dtld33(i,i) = (b0/(3*beta02))*sigma(i)/(mu0bi*mu0bi);
		}
	}
	else {
		SetZero<MatrixXd>(Dtld32,p,nvcovpar,true);
		SetZero<MatrixXd>(Dtld33,p,p,true);
	}
	for (int i=0;i<p;i++) ksi1(i) = mu1(i)-mu0(i);
	outerprod<VCTTP,SQMATTP>(p,mu0,mu0OtP);
	Omega = Sigma+mu0OtP;
	double logDet;
	if (!pdsolve<SQMATTP>(Omega,OmegaInv,&logDet)) {
		warning("Non-positive definite Omega matrix found during gradient computations (which returned 0).\n");
		for (int i=0;i<p;i++) CPgrad(i) = 0.;
		return;
	}
	for (int i=0;i<p;i++)  {
		omega(i) = sqrt(Omega(i,i));
		bomega(i) = b*omega(i);
		delta(i) = mu0(i)/(bomega(i));
	}
	outerprod<VCTTP,SQMATTP>(p,omega,sclMatI);
	OmgbI = sclMatI.cwiseProduct(OmegaInv); 
	OmgbIdelta = OmgbI * delta;
	double c2(1.-delta.dot(OmgbIdelta));

	if (c2 < c2tol) {    

		static MatrixXd dOmgbdOmg,vtmpM,vdOmgb,tmpextv,dOmgdgam1;     // Try sparse matrix representations for dOmgbdOmg !!!
		static SQMATTP Omgb; 
		static VectorXi diagind;
		static int diaginddim(0);
		static VCTTP dAomgm1,Adomgm1,gamma1grp1,gamma1grp2;
		double tmpc,penc(PenF*(c2-c2tol));

		if (diaginddim!=p)  {
			diagind.resize(p);
			diagind(0) = 0;
			int cnst(p*(p+1));
			for (int i0=1,i1=p-1;i0<p;i0++,i1--) diagind(i0) = (cnst-i1*(i1+1))/2;
			diaginddim = p;
		}
		if (vtmpM.rows()!=1 || vtmpM.cols()!=nvcovpar) vtmpM.resize(1,nvcovpar);  
		if (vdOmgb.rows()!=1 || vdOmgb.cols()!=nvcovpar) vdOmgb.resize(1,nvcovpar);  
		if (tmpextv.rows()!=1 || tmpextv.cols()!=nvcovpar) tmpextv.resize(1,nvcovpar);  
		SetZero<MatrixXd>(dOmgbdOmg,nvcovpar,nvcovpar,true); 
		if (dOmgdgam1.rows()!=nvcovpar || dOmgdgam1.cols()!=p) dOmgdgam1.resize(nvcovpar,p);  
		if (!FixedArrays)  {
			if (Omgb.rows()!=p || Omgb.cols()!=p)  Omgb.resize(p,p);
			if (dAomgm1.size()!=p)  dAomgm1.resize(p);
			if (Adomgm1.size()!=p)  Adomgm1.resize(p);
			if (gamma1grp1.size()!=p)  gamma1grp1.resize(p);
  			if (gamma1grp2.size()!=p)  gamma1grp2.resize(p);
		}
		Omgb = Omega.array() / sclMatI.array(); 
		for (int c=0,i=1;c<p-1;c++,i++) for (int r=c+1;r<p;r++,i++) {
			dOmgbdOmg(i,i) = 1./sclMatI(r,c);
			dOmgbdOmg(i,diagind(r)) = -Omgb(r,c)/(2*Omega(r,r));
			dOmgbdOmg(i,diagind(c)) = -Omgb(r,c)/(2*Omega(c,c));
		}
		outerprod<VCTTP,SQMATTP>(p,OmgbIdelta,tmpM);
		for (int i=0,c=0;c<p;c++) for (int r=c;r<p;r++,i++)
			if (r==c) vtmpM(i) = tmpM(r,c);
			else vtmpM(i) = 2*tmpM(r,c);
		vdOmgb = vtmpM * dOmgbdOmg;
		tmpextv = vdOmgb * D23 * Dtld32;
		for (int i=0;i<nvcovpar;i++) 
			penaltygrad(k*p+i) -= penc*(vdOmgb(i)+tmpextv(i));
 
		dOmgdgam1 = D23 * Dtld33;		
		for (int i=0;i<p;i++) {
			tmpc = cubrootb0*sigma(i);
			dAomgm1(i) = tmpc / (3*bomega(i)*pow(gamma1(i)*gamma1(i),1./3));
			Adomgm1(i) = -( tmpc*pow(fabs(gamma1(i)),1./3) / (bomega(i)*omega(i)) ) *
				dOmgdgam1(diagind(i),i)/(2*omega(i)) ;
			if (gamma1(i) < 0.) Adomgm1(i) *= -1; 
		}
		gamma1grp1 = -2*OmgbIdelta.cwiseProduct(dAomgm1+Adomgm1);    
		gamma1grp2 = vdOmgb*dOmgdgam1; 
		for (int i=0;i<p;i++)
			penaltygrad(k*p+nvcovpar+i) -= penc*(gamma1grp1(i)+gamma1grp2(i)); 
		if (c2 < Machinetol) {
			CPgrad = penaltygrad;
			return;
		}
	}

	double b2(b*b);
	double c1(sqrt((b2-(1.-b2)*beta02)/(1.+beta02)));
	double q1(1./(c1*(1.+beta02)));    
	double q2(q1*(2*c1-q1)/2);
	double q1q2(q1*q2);
	tmpv = sqrt(fabs(q1q2)) * SigmaImu0;
	outerprod<VCTTP,SQMATTP>(p,tmpv,tmpM);
	if (q1q2 < 0.) tmpM *= -1;
	D33 = q1*SigmaI-tmpM;    
	for (int i=0;i<p;i++) for (int c=0,j=0;c<p;c++) for (int r=c;r<p;r++,j++) 
		if (r==c) D32(i,j) = -SigmaImu0(c) * D33(i,c) ;
		else D32(i,j) = -SigmaImu0(c)*D33(i,r) - SigmaImu0(r)*D33(i,c) ;
	D33 -= tmpM;    
	eta = q1*SigmaImu0;

	for (int r=0;r<n;r++) for(int c=0;c<p;c++) 
		if (grpind(r)<0) y0(r,c) = y(r,c) - ksi1(c); 
		else y0(r,c) = y(r,c) - ksi1(c) - beta2k(grpind(r),c); 

	nOmgminusS0 = n*Omega - y0.transpose()*y0;
	for (int c=0,j=0;c<p;c++) for (int r=c;r<p;r++,j++)  { 
		if (r==c) {
			OmegaInvgrad1(j) = nOmgminusS0(c,c)/2 ;
			outerprod<VCTTP,SQMATTP>(p,OmegaInv.col(c),tmpM);
			for (int c1=0,j1=0;c1<p;c1++) for (int r1=c1;r1<p;r1++,j1++)   
				OmegaInvgrad2(j1,j) = -tmpM(r1,c1);
		}
		else {
			OmegaInvgrad1(j) = nOmgminusS0(r,c) ;
			outerprod<VCTTP,SQMATTP>(p,OmegaInv.col(c),OmegaInv.col(r),tmpM);
			outerprod<VCTTP,SQMATTP>(p,OmegaInv.col(r),OmegaInv.col(c),tmpM2);
			for (int c1=0,j1=0;c1<p;c1++) for (int r1=c1;r1<p;r1++,j1++)   
				OmegaInvgrad2(j1,j) = -tmpM(r1,c1)-tmpM2(r1,c1) ;
		}
	}
	Omegagrad = OmegaInvgrad1 * OmegaInvgrad2;

	yg0sum.row(grpind(0)+1) = y0i = y0.row(0);   
	double z1(zeta(1,y0i.dot(eta)));
	z1eta = z1*eta;
	sumgz1eta.row(grpind(0)+1) = z1eta;
	etagrad = z1*y0i;
	for (int obs=1;obs<n;obs++) {
		y0i = y0.row(obs);
  		yg0sum.row(grpind(obs)+1) += y0i;
		z1 = zeta(1,y0i.dot(eta));
		z1eta = z1*eta;
		sumgz1eta.row(grpind(obs)+1) += z1eta;
		etagrad += z1*y0i;
	}
	y0sum = yg0sum.row(0);
	sumz1eta = sumgz1eta.row(0);
	for (int g=1;g<k;g++) {
		y0sum += yg0sum.row(g);
		sumz1eta += sumgz1eta.row(g);
	}

        ksi1grad = OmegaInv*VectorXd(y0sum) - VectorXd(sumz1eta);
	for (int g=1;g<k;g++) beta2kgrad.row(g-1) = OmegaInv*VectorXd(yg0sum.row(g)) - VectorXd(sumgz1eta.row(g));
	for (int i=0;i<p;i++) CPgrad(i) = penaltygrad(i) + ksi1grad(i);
	for (int i=p,j=0;j<p;j++) for (int g=1;g<k;g++,i++) CPgrad(i) = penaltygrad(i) + beta2kgrad(g-1,j);
	Sigmagrad = -ksi1grad.transpose()*Dtld32 + Omegagrad + Omegagrad*D23*Dtld32 + 
		etagrad.transpose()*(D32+D33*Dtld32) ;
	for (int i0=0,i=k*p;i0<nvcovpar;i0++,i++) CPgrad(i) = penaltygrad(i) + Sigmagrad(i0);
  	gamma1grad = -ksi1grad.transpose()*Dtld33 + Omegagrad*D23*Dtld33 + 
		etagrad.transpose()*D33*Dtld33 ;
	for (int i0=0,i=k*p+nvcovpar;i0<p;i0++,i++) CPgrad(i) = penaltygrad(i)+gamma1grad(i0);

	if(trace) {         
		Rprintf("msnCP.ll.grad -- penalty ="); Rprintv<VectorXd>(2*p+nvcovpar,penaltygrad);
		Rprintf("\nCentred parameters:");
		Rprintf("mu1 = ") ; Rprintv<NumericVector>(p,mu1);
                if (k>1) { Rprintf("beta2k = ") ; RprintM<NumericMatrix>(k-1,p,beta2k); }
		Rprintf("gamma1 = ") ; Rprintv<NumericVector>(p,gamma1);
		Rprintf("Sigma =\n") ; RprintM<SQMATTP>(p,p,Sigma);
		Rprintf("Direct parameters:\n");
		Rprintf("ksi1 = ") ; Rprintv<VCTTP>(p,ksi1);
		Rprintf("alpha = ") ; Rprintv<VCTTP>(p,omega.cwiseProduct(eta));
		Rprintf("Omega =\n") ; RprintM<SQMATTP>(p,p,Omega);
	}

	return;
}

template<class VCTTP,class ROWVCTTP,class EXTVCTTP,class DMVCTTP,class SQMATTP,class RCTMATTP>
NumericVector msnCP_dev_grad1(NumericVector& param, const NumericMatrix& y, const IntegerVector& grpind,
		const int Config, const int n, const int p, const int k,
		const bool trace, const double c2tol, const double ldRtol, const double beta0tol, 
		const double PenF, const double MachineEPS, const bool FixedArrays)
{
	int q(p/2), nvcovpar(p*(p+1)/2), nvcovsrpar(ncovp(Config,q,p));
	int CPgradl((k+1)*p+nvcovpar), gradl((k+1)*p+nvcovsrpar);
	NumericVector::iterator mu1ptr(param.begin());
	NumericVector::iterator beta2kptr(mu1ptr+p);
	NumericVector::iterator Sigmaptr(mu1ptr+k*p);
	NumericVector::iterator gamma1ptr(Sigmaptr+nvcovsrpar);
	double Machinetol(sqrt(MachineEPS));

	static SQMATTP Sigma;
	static NumericMatrix beta2k;
	static EXTVCTTP CPgrad,tmpv;
	static DMVCTTP SigmaSrgrad;
	static RCTMATTP DSigSigSr,tmpM;   // Try sparse matrix representation for DSigSigSr !!! 

	static NumericVector paramgrad;
	if (paramgrad.size()!=gradl) paramgrad = clone(param); 
	SetZero<RCTMATTP>(DSigSigSr,nvcovpar,nvcovsrpar,true); 
	if (SigmaSrgrad.size()!=nvcovsrpar)  SigmaSrgrad.resize(nvcovsrpar);
	if (tmpM.rows()!=nvcovpar || tmpM.cols()!=nvcovsrpar) tmpM.resize(nvcovpar,nvcovsrpar);
	if ( !FixedArrays)  {
		if (CPgrad.size()!=CPgradl)  CPgrad.resize(CPgradl);
		if (Sigma.rows()!=p || Sigma.cols()!=p)  Sigma.resize(p,p);
		if (tmpv.size()!=nvcovpar)  tmpv.resize(nvcovpar);
	}
  
	NumericVector mu1(mu1ptr,mu1ptr+p);
        if (k>1)  beta2k = NumericMatrix(k-1,p,beta2kptr);
	Sigma = RestCov<SQMATTP>(q,Sigmaptr,Config,FixedArrays);
	NumericVector gamma1(gamma1ptr,gamma1ptr+p);

	msnCP_ll_grad<VCTTP,ROWVCTTP,EXTVCTTP,SQMATTP,RCTMATTP>(
		mu1,beta2k,Sigma,gamma1,y,grpind,n,p,k,nvcovpar,trace,
		PenF,ldRtol,c2tol,beta0tol,CPgrad,MachineEPS,FixedArrays );

	bool zerograd(true);
	for (int i=0;zerograd&&i<CPgradl;i++)
		if (fabs(CPgrad(i)) > Machinetol) zerograd = false;
	if (zerograd) {
		for (int i=0;i<gradl;i++) paramgrad(i) = 0.;
		return paramgrad;
	}
	RestCov_grad<RCTMATTP,SQMATTP,EXTVCTTP>(
		p,q,nvcovpar,Config,param.begin()+k*p,FixedArrays,DSigSigSr);

	for (int c=0,ind=0;c<p;c++) for (int r=c;r<p;r++,ind++) {
		tmpv(ind) = CPgrad(p*k+ind);
		tmpM.row(ind) = DSigSigSr.row(utind1(c,r));
	}

	SigmaSrgrad = tmpv * tmpM; 
	for (int i=0;i<k*p;i++) paramgrad(i) = -2*CPgrad(i);
	for (int i0=0,i1=k*p;i0<nvcovsrpar;i0++,i1++) paramgrad(i1) = -2*SigmaSrgrad(i0);
	for (int i0=k*p+nvcovpar,i1=k*p+nvcovsrpar;i0<CPgradl;i0++,i1++) paramgrad(i1) = -2*CPgrad(i0);

	return paramgrad;
}

