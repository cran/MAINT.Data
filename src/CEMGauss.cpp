#include<vector>
#include <RcppEigen.h>
#include "EigenTmplFoo.cpp"
#include "CEMGauss.h"
#include "MDataGaussLogLik.h"

const double MINLIKEXP = -200.;

using namespace Eigen ;

typedef Map<VectorXd> Mvct;
typedef Map<MatrixXd> Mmat;

RcppExport
SEXP CEMGauss(SEXP X_s, SEXP k_s, SEXP Config_s, SEXP Homoc_s, SEXP maxiter_s, SEXP tautol_s, SEXP convtol_s,
    SEXP InitSolz_s, SEXP InitSoltau_s, SEXP InitSolmuk_s, SEXP InitSolSigma_s, SEXP InitSolSigmak_s, SEXP InitSolLnLik_s)        
{
   Mmat X(as<Mmat>(X_s));
   int k(as<int>(k_s)), Config(as<int>(Config_s)), maxiter(as<int>(maxiter_s));
   bool Homoc(static_cast<bool>(as<int>(Homoc_s)));
   double tautol(as<double>(tautol_s)), convtol(as<double>(convtol_s));
   int n(X.rows()),p(X.cols());
   NumericMatrix muk(p,k);
   MatrixXd z(n,k),Likk(n,k),LikExp(n,k),Sigma(p,p);
   NumericVector tau(k),nrmfct(n),Sigmak(Dimension(p,p,k));
   double Likall,LnLik1;

   bool Initz;
   NumericMatrix InitSolz,InitSolmuk,InitSolSigma;
   NumericVector InitSoltau,InitSolSigmak;

   if (!Rf_isNull(InitSolz_s))  {
     InitSolz = InitSolz_s;
     Initz = true;
   } else {
     Initz = false;
   }

   if (!Rf_isNull(InitSolmuk_s)) InitSolmuk = InitSolmuk_s;
   if (!Rf_isNull(InitSolSigma_s))  InitSolSigma = InitSolSigma_s;
   if (!Rf_isNull(InitSoltau_s))  InitSoltau = InitSoltau_s;
   if (!Rf_isNull(InitSolSigmak_s))  {
     InitSolSigmak = InitSolSigmak_s;
     InitSolSigmak.attr("dim") = IntegerVector::create(p,p,k);
   }  

   double LnLik(as<double>(InitSolLnLik_s));
   NumericMatrix Wk(p,p),wdev(n,p); 
   static vector<double> tmpvct;
   if (tmpvct.size()!=n) tmpvct.resize(n);
   if (!(Initz))  {
     tau = clone(InitSoltau);
     muk = clone(InitSolmuk);
     if (Homoc) Sigma = as<MatrixXd>(clone(InitSolSigma));
     for (int g=0;g<k;++g) {	
        Mvct mukj(muk.begin()+g*p,p); 
        if (Homoc) MDataGaussLogLik(n,p,Config,X,mukj,Sigma,tmpvct);
        else {
          Mmat Sigmakj(InitSolSigmak.begin()+g*p*p,p,p);           
          MDataGaussLogLik(n,p,Config,X,mukj,Sigmakj,tmpvct);
        }
        for (int obs=0;obs<n;++obs) LikExp(obs,g) = tmpvct[obs];
     }

     for (int obs=0;obs<n;++obs) {
       nrmfct(obs) = 0.;
       double maxLikExp = LikExp(obs,0);
       for (int g=1;g<k;++g) maxLikExp = fmax(LikExp(obs,g),maxLikExp);
       if (maxLikExp < MINLIKEXP) { 
         for (int g=0;g<k;++g) LikExp(obs,g) -= maxLikExp;
         nrmfct(obs) = maxLikExp;
       }
       for (int g=0;g<k;++g) Likk(obs,g) = tau[g] * exp(LikExp(obs,g));
     }     

     for (int obs=0;obs<n;++obs) {
       Likall = 0.;
       for (int g=0;g<k;++g) Likall += Likk(obs,g);
       if (Likall>MINZ) for (int g=0;g<k;++g) z(obs,g) = fmax(Likk(obs,g)/Likall,MINZ);
       else for (int g=0;g<k;++g) z(obs,g) = 1./k;
      }
   } else { 
     z = as<MatrixXd>(clone(InitSolz));
  }  
   
  bool converg(false),stopcycle(false);
  int iter =0;
  std::vector<double> nk(k); 

  while (!converg && iter<maxiter)  {
    for (int g=0; g<k; ++g) {
      nk[g] = 0;
      for (int obs=0;obs<n;++obs) nk[g] += z(obs,g);
      if (nk[g] < tautol) stopcycle = true;
      tau[g] = fmax(nk[g]/n,0.);
    }
    if (stopcycle) break;
      if (Homoc==TRUE) Sigma.setZero();
      for (int g=0,Sigmakp=0; g<k; ++g,  Sigmakp+=p*p)  {
        for (int j=0; j<p; ++j) {
          double tmpsum = 0;
          for (int obs=0; obs<n; ++obs) tmpsum += z(obs,g) * X(obs,j);
          muk(j,g) = tmpsum/nk[g];
        }    
        for (int obs=0;obs<n;++obs) {
          double zweight = sqrt(z(obs,g)); 
          for (int j=0; j<p; ++j) wdev(obs,j) = zweight * (X(obs,j)-muk(j,g));
        }   

        for (int j1=0;j1<p;++j1) {
          double psum = 0.;
          for (int obs=0;obs<n;obs++) psum += wdev(obs,j1)*wdev(obs,j1);
          Wk(j1,j1) = psum;
          if (Config==1) {
            for (int j2=0;j2<j1;++j2) {
              double psum = 0.;
              for (int obs=0;obs<n;obs++) psum += wdev(obs,j1)*wdev(obs,j2);
              Wk(j1,j2) = Wk(j2,j1) = psum;
            }
          } 
        }
        if (Config!=1 && Config!=5) {
          int q = p/2;
          if (Config==3) {
            for (int j1=0;j1<q;++j1) {
              double psum = 0.;
              for (int obs=0;obs<n;obs++) psum += wdev(obs,j1)*wdev(obs,q+j1);
              Wk(j1,q+j1) = Wk(q+j1,j1) =psum;
            }
          } else if (Config==4) {
            for (int j1=0;j1<q;++j1) for (int j2=0;j2<j1;++j2) {
              double psum = 0.;
              for (int obs=0;obs<n;obs++) psum += wdev(obs,j1)*wdev(obs,j2);
              Wk(j1,j2) = Wk(j2,j1) = psum;
            }
            for (int j1=q;j1<p;++j1) for (int j2=q;j2<j1;++j2) {
              double psum = 0.;
              for (int obs=0;obs<n;obs++) psum += wdev(obs,j1)*wdev(obs,j2);
              Wk(j1,j2) = Wk(j2,j1) = psum;
            }
          } 
        }

        if (Homoc==FALSE) {
          for (int j1=0;j1<p;++j1) {
            Sigmak(Sigmakp+j1*p+j1) = Wk(j1,j1)/nk[g];
            for (int j2=0;j2<j1;++j2) Sigmak(Sigmakp+j1*p+j2) = Sigmak(Sigmakp+j2*p+j1) = Wk(j1,j2)/nk[g];
          } 
          Mmat Sigmakj(Sigmak.begin()+g*p*p,p,p);           
          Mvct mukj(muk.begin()+g*p,p); 
          MDataGaussLogLik(n,p,Config,X,mukj,Sigmakj,tmpvct);
        }          
        for (int obs=0;obs<n;++obs) LikExp(obs,g) = tmpvct[obs];
        if (Homoc==TRUE) {
          for (int j1=0;j1<p;++j1) {
            Sigma(j1,j1) += Wk(j1,j1)/n;
              for (int j2=0;j2<j1;++j2) {
                Sigma(j1,j2) += Wk(j1,j2)/n;
                Sigma(j2,j1) = Sigma(j1,j2);
              }
          }
        }
      } 

      if (Homoc==TRUE) {
        for (int g=0;g<k;++g) {
          Mvct mukj(muk.begin()+g*p,p); 
          MDataGaussLogLik(n,p,Config,X,mukj,Sigma,tmpvct);
          for (int obs=0;obs<n;++obs) LikExp(obs,g) = tmpvct[obs];
        }
      }
      
      for (int obs=0;obs<n;++obs) {
        nrmfct(obs) = 0.;
        double maxLikExp = LikExp(obs,0);
        for (int g=1;g<k;++g) maxLikExp = fmax(LikExp(obs,g),maxLikExp);
        if (maxLikExp < MINLIKEXP) { 
          for (int g=0;g<k;++g) LikExp(obs,g) -= maxLikExp;
          nrmfct(obs) = maxLikExp;
        }
        for (int g=0;g<k;++g) Likk(obs,g) = tau[g] * exp(LikExp(obs,g));
      }     
      
      LnLik1 = 0.;
      for (int obs=0;obs<n;++obs) {
        Likall = 0.;
        for (int g=0;g<k;++g) Likall += Likk(obs,g);
        if (Likall>MINZ) for (int g=0;g<k;++g) z(obs,g) = fmax(Likk(obs,g)/Likall,MINZ);
        else for (int g=0;g<k;++g) z(obs,g) = 1./k;
        LnLik1 += nrmfct(obs) + log(Likall);
      }
      if (!isfinite(LnLik1)) { 
        LnLik = -std::numeric_limits<double>::infinity();
        converg = true;
      } else {
        if (LnLik1 < LnLik+convtol) converg = true; 
        LnLik = LnLik1;
        iter += 1;
      } 
    }

    IntegerVector grp(n);
    for (int obs=0;obs<n;++obs) {
      double maxz = 0.;
      for (int g=0;g<k;++g) if (z(obs,g) > maxz) {
        grp[obs] = g+1;
        maxz = z(obs,g);
      }
    }

    int npar,Sigmapar;
    if (Config==1) Sigmapar = p*(p+1)/2;	
    else if (Config==3) Sigmapar = 3*p/2;	
    else if (Config==4) Sigmapar = p*(p/2+1)/2;	
    else if (Config==5) Sigmapar = p;
    if (Homoc==TRUE) npar = p*k + Sigmapar + k-1;		
    else npar = (p+Sigmapar)*k + k-1;
    double BIC = -2*LnLik + log(n)*npar;	 		
    double AIC = -2*LnLik + 2*npar;	 		    

    if (Homoc==TRUE) return List::create(
      Named("tau")=tau,
      Named("muk")=muk,
      Named("Sigma")=Sigma,
      Named("z")=z,
      Named("clusters")=grp,
      Named("LnLik")=LnLik,
      Named("npar")=npar,
      Named("BIC")=BIC,
      Named("AIC")=AIC
    );
    else return List::create(
      Named("tau")=tau,
      Named("muk")=muk,
      Named("Sigmak")=Sigmak,
      Named("z")=z,
      Named("clusters")=grp,
      Named("LnLik")=LnLik,
      Named("npar")=npar,
      Named("BIC")=BIC,
      Named("AIC")=AIC
    );
}
