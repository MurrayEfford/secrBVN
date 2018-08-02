\name{simcapt.bvn}
\alias{simcapt.bvn}

\title{Simulate Detections from Elliptical Home Ranges }

\description{
This extends the functionality of \code{\link{sim.capthist}} to elliptical home ranges, with some limitations.
}

\usage{

simcapt.bvn(traps, popn, type = c('uniform','BVN'), g0 = 0.2, lambda0 = 0.2, p = 0.95, 
    noccasions = 5, renumber = TRUE)

}

\arguments{
  \item{traps}{object of class 'traps'}
  \item{popn}{object of class 'popn' with attributes 's2xy' and 'theta'}
  \item{type}{uniform ellipse vs bivariate normal}
  \item{g0}{numeric intercept of detection function (type = 'uniform')}
  \item{lambda0}{numeric intercept of hazard function (type = 'BVN')}
  \item{p}{numeric probability corresponding to edge of ellipse (type = 'uniform')}
  \item{noccasions}{integer number of sampling occasions}
  \item{renumber}{logical; if FALSE row names are carried over from
  \code{popn}; otherwise the n detected animals are renumbered 1:n }
}

\details{

  \code{traps} may have detector type `multi' or `proximity'; other types 
  will be coerced to detector type 'proximity'.
  
  \code{popn} will have the necessary attributes if it was
  generated with \code{\link{simpopn.bvn}}.
  
  Two detection functions are offered. In one (type = 'BVN') the
  hazard of detection \eqn{h} by a detector at point x is proportional to
  the probability density of the bivariate normal distribution at x specified
  by \code{s2xy} and \code{theta} (these may vary among individuals
  in \code{popn}). The proportionality is determined by lambda0, which is the
  maximum detection hazard (i.e. when detector is at HR centre). 
  The probability of detection is 1 - exp(-\eqn{h}).

  With the other detection function (type = 'uniform') a uniform
  probability of detection (\code{g0}) applies throughout a home range
  ellipse, with zero probability of detection elsewhere. The boundary is
  determined by both the shape, orientation and location (using
  \code{s2xy} and \code{theta} as before) and the argument \code{p} that
  determines the notional bivariate-normal probability contour to
  use. (The simulated distribution is not bivariate-normal - this is
  just a device to specify an ellipse in a familiar way).
  
  For detector type `multi', the hazard of detection is modelled as bivariate 
  normal or elliptical uniform.
  
}
\value{

capthist object suitable for analysis in \pkg{secr}

}

\note{

  Arguments \code{g0} and \code{p} have no effect when \code{type = 'BVN'}.
  Argument \code{lambda0} has no effect when \code{type = 'uniform'}.

}

\seealso{

\code{\link{simpopn.bvn}}, \code{\link{secr.fit}}

}

\examples{

tempgrid <- make.grid(detector='proximity')
pop <- simpopn.bvn(s2xy=c(225/2,225*2), core = tempgrid,
    buffer = 100, D = 10)
temp <- simcapt.bvn(tempgrid, pop, 'uniform')
plot (tempgrid)
plotpopn.bvn(pop)
plot(temp, add = TRUE)

}

\keyword{ datagen }