\name{phi.test}
\alias{phi.test}
\title{Nonparametric goodness-of-fit test based on phi-divergences}
\description{
Performs a one-sample goodness-of-fit test based on the phi-divergence statistics of Jager and Wellner (2007).
}
\usage{
  phi.test(x, s, fx="punif", ...)
}
\arguments{
  \item{x}{A numeric vector of data values.}
  \item{s}{The parameter value, s, in Sn(s).}
  \item{fx}{A character string naming a cumulative distribution function or a built-in distribution function from \code{R}, such as \code{"pnorm"}.  Default is \code{"punif"}.}
  \item{...}{Parameters of the distribution specified (as a character string) by \code{fx}.  Default for \code{"punif"} is 0 and 1.}
}
\details{
\code{phi.test} performs the generalized goodness-of-fit test described in Jager and Wellner (2007).  The null hypothesis is that \code{x} is drawn from the continuous distribution specified by \code{fx}.  The alternative is that \code{x} is not drawn from this distribution.  One-sided tests are not currently supported.

The test is parameterized by the real-valued parameter \code{s}, which can take values in [-1,2].  Choosing \code{s=1} results in the test statistic of Berk and Jones (1979).  Choosing \code{s=2} results in the supremum form of the Anderson-Darling statistic.  Choosing \code{s=-1} results in the statistic of Jaescke (1979) and Eicker (1979).

Exact p-values for this test are calculated through finite sample methods described in Jager and Wellner (2010).  These methods are an extension of those described in Owen (1995) and rely on a recursion from Khmaladze and Shinjinkashivili (2001). Exact calculations are supported for sample sizes up to n=10,000.  
}

\value{
  A list with class "htest" containing the following components:
  \item{statistic}{The value of the test statistic, Sn(s).}
  \item{p.value}{The p-value of the test.}
  \item{alternative}{A character string describing the alternative hypothesis.  Currently \code{"two-sided"} is all that is supported.} 
  \item{method}{A character string indicating what type of test was performed.}
  \item{data.name}{A character string giving the name of the data.}
}
\references{

  Anderson, T. W. and Darling, D. A. (1952), "Asymptotic theory of a certain 'goodness of fit' criteria based on stochastic processes," Ann. Math. Statist. 23: 193 - 212.

  Berk, R. H. and Jones, D. H. (1979), "Goodness-of-fit test statistics that dominate the Kolmogorov statistics," Z. Wahrsch. Verw. Gebiete 35: 47 - 59.

  Eicker, F. (1979), "The asymptotic distribution of the suprema of the standardized empirical processes," Ann. Statist. 7: 116 - 138.

  Jaeschke, D. (1979), "The asymptotic distribution of the supremum of the standardized empirical distribution function on subintervals," Ann. Statist. 7: 108 - 115.

  Jager, L. and Wellner, J. A. (2007), "Goodness-of-fit tests via phi-divergences," Ann. Statist. 35: 2018 - 2053. 

  Jager, L. and Wellner, J. A. (2010), "Extending the range of application for a family of generalized goodness-of-fit statistics," In preparation. 

  Khmaladze, E. and Shinjinkashivili, E. (2001), "Calculation of non-crossing probabilities for Poisson processes and its corollaries," Adv. in Appl. Probab. 33: 702 - 716. 

  Owen, A. B. (1995), "Nonparametric likelihood confidence bands for a distribution function,"  J. Amer. Statist. Assoc. 90: 516 - 521. 

}

\author{Leah R. Jager \email{jager@usna.edu}}
\seealso{\code{\link{phi.bands}}, \code{\link{plot.phibands}}}
\examples{
  \dontrun{
  ## Create uniform(0,1) data
  ## Run a phi-test with s=.5 to test whether the data is uniform(0,1)
  x = runif(100,0,1) 
  phi.test(x, s=.5)

  ## Create normal(10,5) data
  ## Run a phi-test with s=2 to test whether data is normal(10,5)
  ## Run a phi-test with s=2 to test whether data is normal(10,3)
  y = rnorm(50, 10, 5)
  phi.test(y, s=2, "pnorm", 10, 5) 
  phi.test(y, s=2, "pnorm", 10, 3) 
 }
}

