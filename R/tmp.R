CovtoParC2 <- function(Sigma,q,Sqrt)	
{					 
# 	Gets the vectorized form of covariance matrix according to Configuration 2

	intcomb <- q*(q+1)/2	   	    # Number of possible combinations between pairs of interval variables
	par <- array(dim=2*intcomb+q)
	if (Sqrt==TRUE) mat <- t(chol(Sigma))
	else mat <- Sigma
	cnt <- 1	
	for (i in 1:q) {
		for (j in 1:i)  {
			par[cnt] <- mat[i,j] 
			par[intcomb+cnt] <- mat[q+i,q+j]
			cnt <- cnt + 1
		}
		par[2*intcomb+i] <- mat[q+i,i]
	}
	par
} 

PartoCovC2 <- function(par,q)

#  Gets the covariance matrix from its vectorized form according to Configuration 2

#  Arguments

#  par    - The free Configuration 2 covariances in vector form
#  q      - Number of integer variables
  
{
	Sigma <- matrix(0.,2*q,2*q)
	intcomb <- q*(q+1)/2		# Number of possible combinations between pairs of interval variables
	cnt <- 1	
	for (i in 1:q) {
		for (j in 1:i)  {
			Sigma[i,j] <- Sigma[j,i] <- par[cnt]  
			Sigma[q+i,q+j] <- Sigma[q+j,q+i] <- par[intcomb+cnt]  
			cnt <- cnt + 1
		}
		Sigma[q+i,i] <- Sigma[i,q+i] <- par[2*intcomb+i]  
	}
	Sigma
}

logLik <- function(Data,par,q)
{
	Sigma <- PartoCovC2(par,q)
	spar <- CovtoParC2(Sigma,q,Sqrt=TRUE)
	-GC2mLogLik(spar,Data)
}	

