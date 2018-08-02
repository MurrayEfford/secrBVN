<<<<<<< HEAD
\name{plotpopn.bvn}
\alias{plotpopn.bvn}

\title{ Plot Population with Elliptical Home Ranges }

\description{

  An \pkg{secr} popn object with attributes 's2xy' and
  'theta' is plotted (over an existing plot if add = TRUE).

}

\usage{
plotpopn.bvn(popn, radius = 2.45, add = FALSE, region = NULL, ...)
}

\arguments{
  \item{popn}{object of class 'popn'}
  \item{radius}{radius corresponding to chosen p values (see examples)}
  \item{add}{logical; if FALSE a new plot is started by calling plot(popn)}
  \item{region}{polygon(s) to which ellipses should be truncated (optional)}
  \item{\dots}{other arguments passed to \code{\link{polygon}}}
  
}

\details{
  
  Coordinates of each ellipse are generated and plotted.
  
  \code{region} may be an sp object that inherits from SpatialPolygons or a 2-column matrix of boundary coordinates. 
  
}

\value{
  
  None. The function is used for its side effect of adding to an
  existing plot.
  
}

\seealso{

\code{\link{simpopn.bvn}}

}

\examples{
tempgrid <- make.grid()
circular.pop <- simpopn.bvn(s2xy = c(625,625), core = tempgrid,
    buffer = 100, D = 1)
elliptical.pop <- simpopn.bvn(s2xy = c(625/3,625*3), core = tempgrid,
    buffer = 100, D = 1)
par(mfrow = c(1,2), par(xpd = TRUE))
plot(circular.pop)   ## this step determines the axes, frame etc.
plotpopn.bvn(circular.pop, add = TRUE)
plotpopn.bvn(elliptical.pop)

## demo aligned ellipses
par(mfrow=c(2,2))
for (i in 1:4) {
    aligned.pop <- simpopn.bvn(s2xy = c(625/3,625*3), theta=-1,
        core = tempgrid, buffer = 100, D = 1)
    plotpopn.bvn(aligned.pop, col='lightblue')
}

## demo variable size (up to 9-fold)
par(mfrow=c(1,1))
pop <- simpopn.bvn(s2xy = c(625,625), 
    core = tempgrid, buffer = 100, D = 1)
attr(pop, 's2xy') <- attr(pop, 's2xy') * (0.5 + runif(nrow(pop)*2))^2
plot(pop, cex=0)
plotpopn.bvn(pop, add = TRUE, col='grey')

}

\keyword{ hplot }
=======
\name{plotpopn.bvn}
\alias{plotpopn.bvn}

\title{ Plot Population with Elliptical Home Ranges }

\description{

  An \pkg{secr} popn object with attributes 's2xy' and
  'theta' is plotted (over an existing plot if add = TRUE).

}

\usage{
plotpopn.bvn(popn, radius = 2.45, add = FALSE, region = NULL, ...)
}

\arguments{
  \item{popn}{object of class 'popn'}
  \item{radius}{radius corresponding to chosen p values (see examples)}
  \item{add}{logical; if FALSE a new plot is started by calling plot(popn)}
  \item{region}{polygon(s) to which ellipses should be truncated (optional)}
  \item{\dots}{other arguments passed to \code{\link{polygon}}}
  
}

\details{
  
  Coordinates of each ellipse are generated and plotted.
  
  \code{region} may be an sp object that inherits from SpatialPolygons or a 2-column matrix of boundary coordinates. 
  
}

\value{
  
  None. The function is used for its side effect of adding to an
  existing plot.
  
}

\seealso{

\code{\link{simpopn.bvn}}

}

\examples{
tempgrid <- make.grid()
circular.pop <- simpopn.bvn(s2xy = c(625,625), core = tempgrid,
    buffer = 100, D = 1)
elliptical.pop <- simpopn.bvn(s2xy = c(625/3,625*3), core = tempgrid,
    buffer = 100, D = 1)
par(mfrow = c(1,2), par(xpd = TRUE))
plot(circular.pop)   ## this step determines the axes, frame etc.
plotpopn.bvn(circular.pop, add = TRUE)
plotpopn.bvn(elliptical.pop)

## demo aligned ellipses
par(mfrow=c(2,2))
for (i in 1:4) {
    aligned.pop <- simpopn.bvn(s2xy = c(625/3,625*3), theta=-1,
        core = tempgrid, buffer = 100, D = 1)
    plotpopn.bvn(aligned.pop, col='lightblue')
}

## demo variable size (up to 9-fold)
par(mfrow=c(1,1))
pop <- simpopn.bvn(s2xy = c(625,625), 
    core = tempgrid, buffer = 100, D = 1)
attr(pop, 's2xy') <- attr(pop, 's2xy') * (0.5 + runif(nrow(pop)*2))^2
plot(pop, cex=0)
plotpopn.bvn(pop, add = TRUE, col='grey')

}

\keyword{ hplot }
>>>>>>> 67bfb5b162dab3c4406b76b4faf48eabee35533f
