CEMGaus <- function(X,k,InitSol,Config,Homoc,maxiter,tautol,convtol)
{
#   cat("CEMGaus -- Got here\n") 

  .Call( "CEMGaus", X, k, InitSol, Config, Homoc, maxiter, tautol, convtol,
          PACKAGE = "MAINT.Data" )
}

