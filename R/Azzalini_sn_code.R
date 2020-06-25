# Functions copied from A.Azzalini sn package (version 1.3-0) with minor adaptations 
# 
# Original author:  A.Azzalini
# Home-page: http://azzalini.stat.unipd.it/SN

# Adaptations: Pedro Duarte Silva (PDS)

mymsn.mle <-function( param, y, x=rep(1,nrow(y)), w=rep(1,nrow(y)), trace=FALSE, 
                      algorithm=c("nlminb", "Nelder-Mead", "L-BFGS-B", "CG", "SANN"), 
			p=ifelse(is.matrix(y),ncol(y),1), 
			lbound=c(rep(-Inf,length(param))),
			ubound=c(rep(Inf,length(param))), 
#			control=list(), ...
			limlnk2, control=list(), ...
                    )
{
#  if (is.null(param)) return(NULL)
  
  y <- data.matrix(y)
  if (!is.numeric(x)) { stop("x must be numeric") }
  opt.method <- match.arg(algorithm)
  x <- data.matrix(x) 
  d <- ncol(y)  
  n <- sum(w)
  p <- ncol(x)
  if (is.null(param)) {
    fit0  <- lm.wfit(x, y, w, method="qr")
    beta  <- as.matrix(coef(fit0))
    res   <- resid(fit0)
#    a     <- msn.moment.fit(res)
    a     <- msn.moment.fit(res,limlnk2=limlnk2)
    Omega <- a$Omega
    omega <- a$omega
    alpha <- a$alpha
    if (!a$admissible) { alpha<-alpha/(1+max(abs(alpha))) }
    beta[1,] <- beta[1,]-omega*a$delta*sqrt(2/pi)  
    eta <-alpha/omega
    param <- c(beta,eta)
    if(trace)
    { 
      cat("Initial parameters:\n")
      print(cbind(t(beta),eta,Omega))
    }
  }
  dev <- mymsn.dev(param, x, y, w)    
  if(opt.method == "nlminb")
  {
#    opt <- nlminb(param, mymsn.dev, mymsn.dev.grad, control=control, 
    opt <- nlminb(param, mymsn.dev, mymsn.dev.grad, control=control, limlnk2=limlnk2, 
      x=x, y=y, w=w, trace=trace, lower=lbound, upper=ubound)

    opt$value <- opt$objective 
    opt$objective <- NULL

  }  else  {
#    opt <- optim(param, fn=mymsn.dev, gr=mymsn.dev.grad, method=opt.method, control=control,
    opt <- optim(param, fn=mymsn.dev, gr=mymsn.dev.grad, method=opt.method, control=control, limlnk2=limlnk2,
      x=x, y=y, w=w, trace=trace, lower=lbound, upper=ubound)
 } 
  if(trace)
  {
    cat("Message from function", opt.method, ":", opt$message,"\n")
    cat("Output parameters " , format(opt$par), "\n")
  }
  opt$iterations <- NULL

  opt
}

#mymsn.dev <- function(param, x, y, w=rep(1,nrow(y)), trace=FALSE)
mymsn.dev <- function(param, x, y, w=rep(1,nrow(y)), limlnk2, trace=FALSE)
{
  d <- ncol(y)
  n <- sum(w)
  p <- ncol(x)
  beta <- matrix(param[1:(p*d)],p,d)
  eta <- param[(p*d+1):(p*d+d)]
  y0 <- y-x %*% beta
  Omega <- (t(y0) %*% (y0*w))/n  
  if (all(is.finite(Omega))) {          # my (PDS) adaptation  - original Azzalini code does not handle
    D <- diag(qr(2*pi*Omega)[[1]])      # the possibility of an infinite Omega within the optimization
  } else {
    return(.Machine$double.xmax)
  }
  D <- diag(qr(2*pi*Omega)[[1]])
  logDet <- sum(log(abs(D)))
  dev <- n*logDet - 2*sum(zeta(0, y0 %*% eta) * w) + n*d
  if(trace)
  { 
    cat("\nmsn.dev:",dev,"\n","parameters:"); 
    print(rbind(beta,eta))
  }

  dev
}

#mymsn.dev.grad <- function(param, x, y, w, trace=FALSE)
mymsn.dev.grad <- function(param, x, y, w, limlnk2, trace=FALSE)
{
  d <- ncol(y)
  if(missing(w)) w <- rep(1,nrow(y))
  n <- sum(w)
  p <- ncol(x)
  beta <- matrix(param[1:(p*d)],p,d)
  eta <- param[(p*d+1):(p*d+d)]
  y0 <- y-x %*% beta
  Omega <- (t(y0) %*% (w*y0))/n
  if (!all(is.finite(Omega)))        # my (PDS) addition  - original Azzalini code does not handle
    { return(rep(0., p*d+d)) }       # the possibility that an infinite Omega is found along the optimization search
  p1 <- zeta(1,as.vector(y0 %*% eta)) * w
#  Omega.inv <- pdwt.solve(Omega, silent=TRUE)
  Omega.inv <- Safepdsolve(Omega,maxlnk2=limlnk2,scale=TRUE)
  if(is.null(Omega.inv)) return(rep(0., p*d+d))  # my (PDS) addition  - original Azzalini returned NA instead of 0. 
  Dbeta <- (t(x) %*% (y0*w) %*% Omega.inv - outer(as.vector(t(x) %*% p1), eta))
  Deta <- as.vector(t(y0) %*% p1)
  if(trace){
    cat("gradient:\n")
    print(rbind(Dbeta,Deta))}
  -2*c(Dbeta,Deta)
}

#msn.moment.fit <- function(y)
#msn.moment.fit <- function(y,singtodiag=TRUE)
msn.moment.fit <- function(y, limlnk2, singtodiag=TRUE)
{# 31-12-1997: simple fit of MSN distribution usign moments
  y     <- as.matrix(y)
  k     <- ncol(y)
  if (k<= nrow(y)) {
    if (singtodiag) var.y <- diag(apply(y,2,function(x) sd(x)^2))
    else return(NULL) 
  } else {
    var.y <- var(y)
  }
  m.y   <- apply(y, 2, mean)  
  y0    <- (t(y) - m.y)/sqrt(diag(var.y))
  gamma1<- apply(y0^3, 1, mean)
  out   <- (abs(gamma1) > 0.99527)
  gamma1[out] <- sign(gamma1[out])*0.995
  a     <- sign(gamma1)*(2*abs(gamma1)/(4-pi))^0.33333
  delta <- sqrt(pi/2)*a/sqrt(1+a^2)
  m.z   <- delta * sqrt(2/pi) 
  omega <- sqrt(diag(var.y)/(1-m.z^2))
  Omega <- var.y + outer(omega*m.z, omega*m.z) 
  xi    <- m.y-omega*m.z
  O.cor <- cov2cor(Omega)
#Oegval <- eigen(O.cor,only.values=TRUE)$values
#  O.inv <- pdwt.solve(O.cor)
#  O.inv <- pdwt.solve(O.cor,silent=TRUE)
  O.inv <- Safepdsolve(O.cor,maxlnk2=limlnk2,scale=FALSE)
  if (is.null(O.inv)) {
    if (singtodiag) O.inv <- diag(rep(1.,k))
    else return(NULL) 
  }
  tmp   <- as.vector(1 - t(delta) %*% O.inv %*% delta)
  if(tmp<=0) {
    tmp <- 0.0001; admissible <- FALSE
  } else  {
    admissible <- TRUE
  }
  alpha <- as.vector(O.inv %*% delta)/sqrt(tmp)
  list(xi=xi, Omega=Omega, alpha=alpha, Omega.cor=O.cor, omega=omega, 
       delta=delta, skewness=gamma1, admissible=admissible) 
}

#mydp2cpMv <- function(dp, family, cp.type="proper", fixed.nu=NULL, aux=FALSE, upto=NULL) 
mydp2cpMv <- function(dp, family, limlnk2, cp.type="proper", fixed.nu=NULL, aux=FALSE, upto=NULL) 
{# internal. NB: name of cp[1] must change according to dp[1]
  cp.type <- match.arg(cp.type, c("proper", "pseudo", "auto"))
  family <- toupper(family)
  if(!(family %in% c("SN", "ESN", "ST","SC")))
    stop(gettextf("family '%s' is not supported", family), domain = NA)
  if (family %in% c("SN","ESN"))
  {  
    if(cp.type == "pseudo") 
      warning("'cp.type=pseudo' makes no sense for SN and ESN families")
#    cp <- msn.dp2cp(dp, aux=aux)
    cp <- msn.dp2cp(dp, limlnk2=limlnk2, aux=aux)
    if(!is.null(upto)) cp <- cp[1:upto]
  }
#  if (family %in% c("SC","ST"))  # Code block turned off by me (PDS) since MAIN.Data does not use the Skew-t or the Skew-Cauchi
#  {
#    if(cp.type=="auto") cp.type <- 
#      if(family == "SC" || dp$nu <= 4) "pseudo" else "proper"
#    if(family == "SC") fixed.nu <- 1
#    cp <- mst.dp2cp(dp, cp.type=cp.type, fixed.nu=fixed.nu, aux=aux, upto=upto)
#    if(is.null(cp))
#    {
#       warning("no CP could be found")
#       return(invisible())
#    }
#  }
  return(cp)
}
  
#msn.dp2cp <- function(dp, aux=FALSE)
msn.dp2cp <- function(dp, limlnk2, aux=FALSE)
{# dp2cp for multivariate SN and ESN 
  alpha <- dp$alpha
  d <- length(alpha)
  Omega <- matrix(dp$Omega, d, d)  
  omega <- sqrt(diag(Omega))
  lot <- delta.etc(alpha, Omega)
  delta <- lot$delta
  delta.star <- lot$delta.star
  alpha.star <- lot$alpha.star
  names(delta) <- names(dp$alpha)
  tau <- if(is.null(dp$tau)) 0 else dp$tau
  mu.z  <- zeta(1, tau) * delta
  sd.z  <- sqrt(1 + zeta(2, tau) * delta^2)
  Sigma <- Omega + zeta(2,tau) * outer(omega*delta, omega*delta)
  gamma1 <- zeta(3, tau) * (delta/sd.z)^3
  if(is.vector(dp[[1]])) { 
    cp <- list(mean=dp[[1]] + mu.z*omega, var.cov=Sigma, gamma1=gamma1)
    }
  else {
    beta <- dp[[1]]  
    beta[1,] <- beta[1,] + mu.z*omega
    cp <- list(beta=beta, var.cov=Sigma, gamma1=gamma1)
  }
  if(!is.null(dp$tau)) cp$tau <- tau
  if(aux){
    lambda <- delta/sqrt(1-delta^2)
    D <- diag(sqrt(1+lambda^2), d, d)
    Ocor <- lot$Omega.cor
    Psi <- D %*% (Ocor-outer(delta,delta)) %*% D
    Psi <- (Psi + t(Psi))/2
#    O.inv <- pdwt.solve(Omega)
    O.inv <- Safepdsolve(Omega,maxlnk2=limlnk2,scale=TRUE)
    O.pcor <- -cov2cor(O.inv) 
    O.pcor[cbind(1:d, 1:d)] <- 1
    R <- force.symmetry(Ocor + zeta(2,tau)*outer(delta,delta))
    ratio2 <- delta.star^2/(1+zeta(2,tau)*delta.star^2)
    mardia <- c(gamma1M=zeta(3,tau)^2*ratio2^3, gamma2M=zeta(4,tau)*ratio2^2)
    # book: (5.74), (5.75) on p.153
    cp$aux <- list(omega=omega, cor=R, Omega.inv=O.inv, Omega.cor=Ocor, 
      Omega.pcor=O.pcor, lambda=lambda, Psi=Psi, delta=delta, lambda=lambda,
      delta.star=delta.star, alpha.star=alpha.star, mardia=mardia)
    }
  return(cp)  
}

delta.etc <- function(alpha, Omega=NULL) 
{ 
  inf <- which(abs(alpha) == Inf)
  if(is.null(Omega)){ # case d=1
    delta <- alpha/sqrt(1+alpha^2)
    delta[inf] <- sign(alpha[inf])
    return(delta)
    }
  else { # d>1
    if(any(dim(Omega) != rep(length(alpha),2))) stop("dimension mismatch")
    Ocor <- cov2cor(Omega)
    if(length(inf) == 0) { # d>1, standard case
      Ocor.alpha <- as.vector(Ocor %*% alpha)
      alpha.sq <- sum(alpha * Ocor.alpha)
      delta <- Ocor.alpha/sqrt(1+alpha.sq)
      alpha. <- sqrt(alpha.sq)
      delta. <- sqrt(alpha.sq/(1+alpha.sq))
      }
     else { # d>1, case with some abs(alpha)=Inf
       if(length(inf) > 1) 
         warning("Several abs(alpha)==Inf, I handle them as 'equal-rate Inf'") 
       k <- rep(0,length(alpha))
       k[inf] <- sign(alpha[inf])
       Ocor.k <- as.vector(Ocor %*% k) 
       delta <- Ocor.k/sqrt(sum(k * Ocor.k))
       delta. <- 1
       alpha. <- Inf
       }
  return(
    list(delta=delta, alpha.star=alpha., delta.star=delta., Omega.cor=Ocor))
  }
}

force.symmetry <- function(x, tol=10*sqrt(.Machine$double.eps)) 
{
  if(!is.matrix(x)) stop("x must be a matrix")
  # err <- abs(x-t(x))
  err <- abs(x-t(x))/(1+abs(x))
  max.err <- max(err/(1+err))
  if(max.err > tol) warning("matrix seems not symmetric")
  if(max.err > 100*tol) stop("this matrix really seems not symmetric")
  return((x + t(x))/2)
}

#mysn.infoMv <- function(dp, x=NULL, y, w, penalty=NULL, norm2.tol=1e-5) 
mysn.infoMv <- function(dp, x=NULL, y, w, limlnk2, penalty=NULL, norm2.tol=1e-5) 
{# computes observed/expected Fisher information matrix for multiv.SN variates
 # using results in Arellano-Valle & Azzalini (JMVA, 2008+erratum)
  type <- if(missing(y)) "expected" else "observed"
  if (type == "observed") { if(!is.matrix(y)) stop("y is not a matrix") }
#  cp <- dp2cpMv(dp, "SN")        #  Original Azzalin code     
#  cp <- mydp2cpMv(dp, "SN")       #  My (PDS) version 
  cp <- mydp2cpMv(dp, limlnk2=limlnk2, "SN")       #  My (PDS) version 
  d <- length(dp$alpha)
  d2 <- d*(d+1)/2
#  if(missing(w)) w <- rep(1, max(NROW(cbind(x,y)),1))       #  Original Azzalin code
  if (is.null(w)) { w <- rep(1, max(NROW(cbind(x,y)),1)) }   #  My (PDS) version
  if (any(w != round(w)) | any(w<0))
    { stop("weights must be non-negative integers") }
  n <- length(w)
  nw <- sum(w)
  if (is.null(x)) {
    p <- 1
    xx <- sum.x <- nw
    x <- matrix(1, nrow=n, ncol=1)
  }  else { 
    p <- NCOL(x)
    # x <- matrix(x, n, p)
    xx <- drop(t(x) %*% (w*x))
    sum.x <- drop(matrix(colSums(w*x)))
  }
  beta <- matrix(dp[[1]],p,d)  ## Note (PDS): Original Azzalini code has as.matrix instead of matrix which seems to be a bug!!
  Omega <- dp$Omega
  omega <- sqrt(diag(Omega))
  alpha <- dp$alpha
  eta   <- alpha/omega
  # vOmega <- Omega[lower.tri(Omega,TRUE)]
  Obar <- cov2cor(Omega)
  Obar.alpha <-  as.vector(Obar %*% alpha)
  alpha.star <- sqrt(sum(alpha * Obar.alpha)) 
  if(alpha.star < 1e-4) {
#    warning("information matrix of multivariate SN not computed near alpha=0")        # Azzalini original code
#   return(NULL)                                                                       # Azzalini original code
    return(list(status="SingInf"))                                                     # My (PDS) version
    }
  # delta.star <- alpha.star/sqrt(1+alpha.star^2)
  c1 <- sqrt(2/pi)/sqrt(1+alpha.star^2)
  c2 <- 1/(pi*sqrt(1+2*alpha.star^2))
  # theta <- c(beta,vOmega,eta)
  D <- duplicationMatrix(d)
  i1 <- 1:prod(dim(beta))
  i2 <- max(i1) + 1:(d*(d+1)/2)
  i3 <- max(i2) + 1:d
  # ind <- list(i1=i1, i2=i2, i3=i3)
#  O.inv <- pdwt.solve(Omega, silent=TRUE)
  O.inv <- Safepdsolve(Omega,maxlnk2=limlnk2,scale=TRUE)
  if (is.null(O.inv)) return(list(status="SingOmega"))                                                      # My (PDS) version     
  if(type == "observed"){ 
    y0 <- y - x %*% beta
    S0 <- t(y0) %*% (w*y0) / nw
    y0.eta <- as.vector(y0 %*% eta)
    z1 <- zeta(1, y0.eta) * w
    z2 <- (-zeta(2, y0.eta) * w)
    # Z2 <- diag(z2, n)
    S1 <- (O.inv %x% t(x)) %*% as.vector(w*y0)- (eta %x% t(x)) %*% z1
    S2 <- (nw/2) * t(D) %*% ((O.inv %x% O.inv) %*% as.vector(S0-Omega))
    S3 <- t(y0) %*% z1
    score <- c(S1,S2,S3)
    u  <- t(x) %*% z1
    U  <- t(x) %*% (z2 * y0)
    V  <- O.inv %*% (2*S0-Omega) %*% O.inv
    # terms as given in the last but one matrix of p.16 
    j11 <- O.inv %x% xx + outer(eta,eta) %x% (t(x) %*% (z2 *x) )
    j12 <- (O.inv %x% (t(x) %*% (w*y0) %*% O.inv))  %*% D
    j13 <- diag(d) %x% u - eta %x% U
    j22 <- (nw/2) * t(D) %*% (O.inv %x% V) %*% D
    j23 <- matrix(0, d*(d+1)/2, d)
    j33 <- t(y0) %*% (z2 * y0)            
    uaA.coef <- NULL
    }
  else { # expected information
    Omega.eta <- omega * Obar.alpha
    mu.c <- Omega.eta/alpha.star^2 
    Omega.c <- Omega - outer(Omega.eta, Omega.eta)/alpha.star^2 
    alpha.bar <- alpha.star/sqrt(1+2*alpha.star^2)
    ginvMills <- function(x, m=0, s=1)  
        # generalized inverse Mills ratio: \phi(x; m, s^2)/\Phi(x)
        exp(-0.5*((x-m)^2/s^2-x^2)+log(zeta(1,x))-log(s))
    fn.u <- function(x, sd, k) x^k * ginvMills(x,0,sd) 
    if(alpha.bar > 0) {
      err<- .Machine$double.eps^0.5
      u0 <- integrate(fn.u, -Inf, Inf, sd=alpha.bar, k=0, rel.tol=err)$value
      u1 <- integrate(fn.u, -Inf, Inf, sd=alpha.bar, k=1, rel.tol=err)$value
      u2 <- integrate(fn.u, -Inf, Inf, sd=alpha.bar, k=2, rel.tol=err)$value }
    else {u0 <- 2; u1<- u2 <- 0}
    a0 <- u0
    a1 <- u1 * mu.c
    A2 <- u2 * outer(mu.c, mu.c) + u0 * Omega.c                    # cfr (19)
    A1 <- (c1*(diag(d)-outer(eta,eta) %*% Omega/(1+alpha.star^2))
           - c2*outer(eta, a1))   # cfr line after (12)
    # terms as given in the last matrix of p.16
    j11 <- (O.inv + c2*a0*outer(eta,eta)) %x% xx
    j12 <- c1*(O.inv %x% outer(sum.x, eta)) %*% D
    j13 <- A1 %x% sum.x
    j22 <- 0.5*nw *t(D) %*% (O.inv %x% O.inv) %*% D
    j23 <- matrix(0, d*(d+1)/2, d)
    j33 <- nw *c2 * A2
    uaA.coef <- list(u0=u0, u1=u1, u2=u2, a1=a1, A1=A1, A2=A2)
    score <- NULL
    }
  I.theta <-rbind(cbind( j11,    j12,   j13),
                  cbind(t(j12),  j22,   j23),
                  cbind(t(j13), t(j23), j33))
#  if(!is.null(penalty))  # Code block turned off by me (PDS) since MAIN.Data does not use the penalized likelihhod version of sn
#  { 
#    # penalization depends on blocks (2,3) of the parameter set only
#    penalty.fn <- if(is.null(penalty)) NULL else get(penalty, inherits=TRUE) 
#    penalty.theta <- function(theta23, penalty, d)
#    {
#      vOmega <- theta23[1:(d*(d+1)/2)]
#      eta <- theta23[(d*(d+1)/2) + (1:d)]
#      Omega <- vech2mat(vOmega)
#      alpha <- eta *sqrt(diag(Omega))
#      penalty(list(alpha=alpha, Omega=Omega))
#    } 
#    i23 <- c(i2,i3)
#    theta23 <- c(Omega[lower.tri(Omega,TRUE)], eta) # beta does not enter here
#    score[i23] <- (score[i23] - 
#      numDeriv::grad(penalty.theta, theta23, penalty=penalty.fn, d=d))
#    jQ <- numDeriv::hessian(penalty.theta, theta23, penalty=penalty.fn, d=d)
#    I.theta[i23, i23] <- I.theta[i23, i23] + jQ
#  }                 
  I.theta <- force.symmetry(I.theta, tol=1e3)
#  inv_I.theta <- pdwt.solve(I.theta, silent=TRUE)
  inv_I.theta <- Safepdsolve(I.theta,maxlnk2=limlnk2,scale=TRUE)
  if(is.null(inv_I.theta)) {
#    warning("numerically unstable information matrix")         # Azzalini original code
#     return(NULL)                                              # Azzalini original code 
      return(list(status="InstInf"))                            # My (PDS) version
     }
  if(type == "observed" ) {
    score.norm2 <- sum(score * as.vector(inv_I.theta %*% score))
#    if(score.norm2/d > norm2.tol) stop("'dp' does not seem to be at MLE")    # Azzalini original code
    if(score.norm2/d > norm2.tol) return(list(status="PositiveScore"))        # My (PDS) version
    }
  D32 <- matrix(0,d, d2)
  tmp32 <- matrix(0,d^2,d^2)
  for(i in 1:d){
    Eii <- matrix(0,d,d)
    Eii[i,i] <- 1
    tmp32 <- tmp32 + Eii %x% Eii
    }
  D32 <- (-0.5)* (t(eta) %x% diag(1/omega^2, d,d)) %*% tmp32 %*% D
  # here we use the expression given in the notes, not in the paper
  Dlow <- cbind(matrix(0,d,d*p), D32, diag(1/omega,d,d))
  Dtheta.dp <- rbind(cbind(diag(d*p+d2), matrix(0,d*p+d2,d)), Dlow)
  I.dp <- t(Dtheta.dp) %*% I.theta %*% Dtheta.dp                     # cfr (14)
  I.dp <- force.symmetry(I.dp, tol=1e3)
  #
  # psi<- c(mu, vSigma, mu0)
  Sigma <- cp$var.cov
  sigma <- sqrt(diag(Sigma))
#  Sigma.inv <- pdwt.solve(Sigma)
  Sigma.inv <- Safepdsolve(Sigma,maxlnk2=limlnk2,scale=TRUE)
  mu0 <- c1* omega * Obar.alpha
  beta0.sq <- as.vector(t(mu0) %*% Sigma.inv %*% mu0)
  beta0 <- sqrt(beta0.sq)
  q1 <- 1/(c1*(1+beta0.sq))
  q2 <- 0.5*q1*(2*c1-q1)
#  Dplus <- pdwt.solve(t(D) %*% D) %*% t(D)
  Dplus <- Safepdsolve(t(D) %*% D,maxlnk2=limlnk2,scale=TRUE) %*% t(D)
  D23 <- Dplus %*% (diag(d) %x% mu0 + mu0 %x% diag(d))
  a <- as.vector(Sigma.inv %*% mu0)
  D32 <- t(-a) %x% (q1 * Sigma.inv - q1*q2*outer(a,a)) %*% D
  D33 <- q1 * Sigma.inv - 2*q1*q2*outer(a,a)
  one00 <- c(1,rep(0,p-1))
  Dtheta.psi <- rbind(
        cbind(diag(p*d),  matrix(0,p*d,d2), -diag(d) %x% one00),
        cbind(matrix(0,d2,p*d),  diag(d2),   D23),
        cbind(matrix(0,d,p*d),    D32,       D33))                # cfr (22a)
  mu0. <- mu0/(sigma*beta0)  # \bar{\mu}_0
  D32. <- matrix(0, d, d2)   # \tilde{D}_{32}
  for(i in 1:d)  {
    Eii <- matrix(0,d,d)
    Eii[i,i] <- 1
    D32. <- D32. + (1/sigma[i])*((t(mu0.) %*% Eii) %x% Eii) %*% D
    }
  D32. <- 0.5* beta0 * D32.
  D33. <- (2/(4-pi)) * diag(sigma/mu0.^2, d, d)/(3*beta0.sq)
  Dpsi.cp <- rbind(cbind(diag(p*d+d2), matrix(0,p*d+d2,d)), 
                   cbind(matrix(0,d,p*d), D32., D33.))            # cfr (22b)
  jacob <- Dtheta.psi %*% Dpsi.cp
  I.cp <- t(jacob) %*% I.theta %*% jacob                          # cfr (17)
  I.cp <- if(any(is.na(I.cp))) NULL else force.symmetry(I.cp)  
#  asyvar.dp <- pdwt.solve(I.dp, silent=TRUE)
  asyvar.dp <- Safepdsolve(I.dp,maxlnk2=limlnk2,scale=TRUE)
  if(is.null(asyvar.dp))  se.dp <- list(NULL) else {
    diags.dp <- sqrt(diag(asyvar.dp))
    se.beta <- matrix(diags.dp[1:(p*d)], p, d)
    se.diagOmega <- diags.dp[p*d + d2 +1 -rev(cumsum(1:d))]
    # se.omega <- se.Omega/(2*omega)
    se.alpha <- diags.dp[p*d +d2 +(1:d)]
    se.dp <- list(beta=se.beta, diagOmega=se.diagOmega, alpha=se.alpha)
    }
#  asyvar.cp <- pdwt.solve(I.cp, silent=TRUE)
  asyvar.cp <- Safepdsolve(I.cp,maxlnk2=limlnk2,scale=TRUE)
  if(is.null(asyvar.cp))  se.cp <- list(NULL) else {
    diags.cp <- sqrt(diag(asyvar.cp))
    se.beta <- matrix(diags.cp[1:(p*d)], p, d)
    se.diagSigma <- diags.cp[p*d + d2 +1 -rev(cumsum(1:d))]
    # se.sigma <- se.Sigma/(2*sigma)
    se.gamma1 <- diags.cp[p*d + d2 +(1:d)]
    se.cp <- list(beta=se.beta, var=se.diagSigma, gamma1=se.gamma1)
    }
  aux <- list(info.theta=I.theta, score.theta=score,
              Dtheta.dp=Dtheta.dp, Dpsi.cp=Dpsi.cp, Dtheta.psi=Dtheta.psi, 
              uaA.coef=uaA.coef)
  list(dp=dp, cp=cp, type=type, info.dp=I.dp, info.cp=I.cp, 
       asyvar.dp=asyvar.dp, asyvar.cp=asyvar.cp, 
#      se.dp=se.dp, se.cp=se.cp, aux=aux)                       # Azzalini original code
       se.dp=se.dp, se.cp=se.cp, aux=aux, status="Regular")     # My (PDS) version
}


