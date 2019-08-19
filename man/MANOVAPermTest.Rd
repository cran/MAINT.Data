\name{MANOVAPermTest}
\alias{MANOVAPermTest}

\title{MANOVA permutation test}
\description{Function MANOVAPermTest performs a MANOVA permutation test allowing for both Gaussian and Skew-Normal distribtuions and homocedastic or heterocedastic setups.}

\usage{

MANOVAPermTest(MANOVAres, Idt, grouping, nrep=200,
    Model=c("Normal","SKNormal","NrmandSKN"), CovCase=1:4,
    SelCrit=c("BIC","AIC"), Mxt=c("Hom","Het","Loc","Gen"), CVtol=1.0e-5, 
    OptCntrl=list(), onerror=c("stop","warning","silentNull"), \dots)

}

\details{Function MANOVAPermTest performs a MANOVA permutation test allowing for both Gaussian and Skew-Normal distribtuions and homocedastic or heterocedastic setups. This test is implemented by simulating the null distribution of the MANOVA likelihood ratio statistc, using many random permutations of the observation group labels. It is intended as an alternative of the classical Chi-squares based MANOVA likelihood ratio tests, when small sample sizes cast doubt on the applicability of the Chi-squared distribution. We note that this test may be computationally intensive, in particular when used for the Skw-Normal model.}


\arguments{
  \item{MANOVAres}{An object representing a MANOVA analysis on interval-valued entities.}

  \item{Idt}{An IData object representing interval-valued entities.}

  \item{grouping}{Factor indicating the group to which each observation belongs to.}

  \item{nrep}{Number of random generated permutations used to approximate the null distribution of the likelihood ratio statistic.}

  \item{Model}{The joint distribution assumed for the MidPoint and LogRanges. Current alternatives are \dQuote{Normal} for Gaussian,
distributions, \dQuote{SKNormal} for Skew-Normal and \dQuote{NrmandSKN} for both Gaussian and Skew-Normal distributions.}

  \item{CovCase}{Configuration of the variance-covariance matrix: a set of integers between 1 and 4.}

   \item{SelCrit}{The model selection criterion.}

   \item{Mxt}{Indicates the type of mixing distributions to be considered. Current alternatives are \dQuote{Hom} (homocedastic) and \dQuote{Het} (hetereocedasic) for Gaussian models, \dQuote{Loc} (location model -- groups differ only on their location parameters) and \dQuote{Gen} \dQuote{Loc} (general model -- groups differ on all parameters) for Skew-Normal models.}

  \item{CVtol}{Tolerance level for absolute value of the coefficient of variation of non-constant variables. When a MidPoint or LogRange has an absolute value within-groups coefficient of variation below CVtol, it is considered to be a constant.}

   \item{OptCntrl}{List of optional control parameters to passed to the optimization routine. See the documentation of RepLOptim for a description of the available options.}

   \item{onerror}{Indicates whether an error in the optimization algorithm should stop the current call, generate a warning, or return silently a NULL object.}

  \item{\dots}{Other named arguments.}
}

\value{the p-value of the MANOVA permutation test.}


\keyword{MANOVA}
\keyword{interval data}
\keyword{likelihood ratio tets}
\keyword{permutation tets}

\seealso{ \code{\link{MANOVA}}, \code{\linkS4class{IdtMANOVA}} }

