\name{plot.phibands}
\alias{plot.phibands}
\alias{plot}
\title{Graphical display of phibands objects}
\description{
Graphical display of phibands objects (confidence bands produced by phi.bands)
}
\usage{
  plot.phibands(x, ...)
}
\arguments{
  \item{x}{A \code{phibands} object, typically the result of \code{phi.bands}.}
  \item{...}{Other parameters to be passed through to plotting functions.}
}
\details{
\code{plot.phibands} creates a plot of the empirical distribution function of a single sample of data, along with corresponding confidence bands produced by the inverting the phi-divergence statistics of Jager and Wellner (2007).  Data set and confidence level are specified when creating the \code{phibands} object using the \code{phi.bands} function.
}

\value{
Nothing of interest.
}
\references{

  Jager, L. and Wellner, J. A. (2007), "Goodness-of-fit tests via phi-divergences," Ann. Statist. 35: 2018 - 2053. 

}

\author{Leah R. Jager \email{jager@usna.edu}}
\seealso{\code{\link{phi.bands}}, \code{\link{phi.test}}}
\examples{
  \dontrun{
  ## Create uniform(0,1) data
  ## Create and plot confidence bands
  x = runif(100,0,1) 
  x.phibands = phi.bands(x, .5)
  plot(x.phibands)

  ## Create normal(10,5) data
  ## Create and plot confidence bands
  y = rnorm(200, 10, 5)
  y.phibands = phi.bands(y, 1, alpha=.1)
  plot(y.phibands)
 }
}

