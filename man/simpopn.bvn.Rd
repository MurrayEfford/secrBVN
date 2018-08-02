\name{simpopn.bvn}
\alias{simpopn.bvn}

\title{
Simulate Population With Bivariate Normal Attributes}

\description{
This wrapper for the \pkg{secr} function sim.popn() adds attributes
's2xy' and 'theta')
}

\usage{
simpopn.bvn(s2xy = c(300, 300), theta = NULL, ...)
}

\arguments{
  \item{s2xy}{numeric vector of s^2(x) and s^2(y), or function to generate these (see Details)}
  \item{theta}{numeric orientation of home range (radians)}
  \item{\dots}{other arguments passed to \code{sim.popn}}
}

\details{
  
  By default all animals receive the same home range (s2xy) attribute values. An alternative is to generate individual values on the fly; this is achieved by providing a function as the argument s2xy. The function should accept a single argument (the number of individuals) and return a 2-column matrix (the values of s2xy, one row per animal).

  NULL \code{theta} causes each animal to be assigned a uniform random
  orientation.
  
  Negative \code{theta} causes a single common random orientation to be
  chosen (all ranges aligned).
  
}

\value{

  Object of class "popn" with additional attributes
  's2xy' and 'theta
}

\seealso{
\code{\link{plotpopn.bvn}}, \code{\link{sim.popn}}
}

\note{

  See \code{\link{plotpopn.bvn}} Examples for variable size (up to
9-fold difference in area in that case).  }

\examples{
tempgrid <- make.grid(detector='proximity')
pop <- simpopn.bvn(s2xy=c(225/2,225*2), core = tempgrid,
    buffer = 100, D = 10)
plot(pop)
plotpopn.bvn(pop)

# random variation
rs2xy <- function(N, scale = 225) {
    matrix(scale * (0.5 + runif(N*2)*2)^2, ncol = 2)
}
pop <- simpopn.bvn(s2xy = rs2xy, core = tempgrid, buffer = 100, D = 10)
plotpopn.bvn(pop)

}

\keyword{ datagen }