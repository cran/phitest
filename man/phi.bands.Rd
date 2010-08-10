\name{phi.bands}
\alias{phi.bands}
\title{Confidence bands based on phi-divergences}
\description{
Constructs confidence bands for the cumulative distribution function by inverting the phi-divergence statistics of Jager and Wellner (2007).
}
\usage{
  phi.bands(x, s, alpha=.05)
}
\arguments{
  \item{x}{A numeric vector of data values.}
  \item{s}{The parameter value, s, in Sn(s).}
  \item{alpha}{The confidence level for the confidence bands.  Default is 95\% confidence bands (or alpha = .05)}
}
\details{
\code{phi.bands} inverts the generalized goodness-of-fit test described in Jager and Wellner (2007) to form confidence bands for the underlying distribution function in a single sample of data. Bands are stepwise functions of the form \eqn{Pr(l_i \leq F(X_{(i)}) \leq h_i) = 1 - \alpha}{Pr(l_i <= F(X_{(i)}) <= h_i) = 1 - alpha}. 

The confidence bands are parameterized by the real-valued parameter \code{s}, which can take values in [-1,2].  Choosing \code{s=1} results in inverstion of the test statistic of Berk and Jones (1979).  Choosing \code{s=2} results in inversion of the supremum form of the Anderson-Darling statistic.  Choosing \code{s=-1} results in inversion of the statistic of Jaescke (1979) and Eicker (1979).

Confidence bands are computed exactly through finite sample methods described in Jager and Wellner (2010).  These methods are an extension of those described in Owen (1995) and rely on a recursion from Khmaladze and Shinjinkashivili (2001). Exact calculations are supported for sample sizes up to n=10,000.  However, for extreme values of \code{s} (near -1 and 2), calculations may take several minutes.
}

\value{
  A object of the class "phibands" containing the following components:
  \item{l}{A vector of values for the lower half of the bands.}
  \item{h}{A vector of values for the upper half of the bands.}
  \item{data}{A vector containing the original data.} 
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
\seealso{\code{\link{phi.test}}, \code{\link{plot.phibands}}}
\examples{
  \dontrun{
  ## Create uniform(0,1) data
  ## Create confidence bands
  x = runif(100,0,1) 
  x.phibands = phi.bands(x, .5)
  x.phibands$l
  x.phibands$h

  ## Plot confidence bands
  plot(x.phibands)

  ## Another example
  y = rnorm(100,0,1) 
  y.phibands = phi.bands(y, .5, alpha=.1)
  plot(y.phibands)

 }
}

