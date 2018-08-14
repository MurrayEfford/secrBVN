\name{anisotropic.fit}
\alias{anisotropic.fit}
\alias{predictAniso}

\title{Spatially Explicit Capture--Recapture in Transformed Space}

\description{

Fits SECR model using distances calculated in a space transformed by compression along a single axis. Both the orientation of the axis (\eqn{\psi_A}) and the compression ratio (\eqn{\psi_R}) are treated as parameters to be estimated. The general method was first applied by Murphy et al. (2016).

}

\usage{

anisotropic.fit(..., psiA = pi/4, psiR = 2)

predictAniso(fit, angle = c("degrees", "radians"))

}

\arguments{
  \item{\dots}{arguments passed to ]\code{\link{secr.fit}}}
  \item{psiA}{numeric initial value of phiA}
  \item{psiR}{numeric initial value of phiR}
  
  \item{fit}{fitted model from \code{anisotropic.fit}}
  \item{angle}{units for angle}
}

\details{


The compression ratio \eqn{\psi_R} takes values greater than or equal to 1.0. The corresponding coefficient on the link scale is log(\eqn{\psi_R - 1}).

The function \code{predictAniso} extracts estimates of the transformation parameters from a fitted model and back-transforms psiR.

The estimate of the spatial scale parameter sigma applies in the isotropic compressed space; it may be understood as an estimate along the minor axis of each ellipse.

}

\value{

For \code{anisotropic.fit}, a fitted secr model.

For \code{predictAniso}, a dataframe of two rows (psiA, psiR).

}

\references{

Murphy, S. M., Cox, J. J., Augustine, B. C., Hast, J. T., Guthrie, J. M., Wright, J., McDermott, J., Maehr, S. C. and Plaxico, J. H. (2016) Characterizing recolonization by a reintroduced bear population using genetic spatial capture--recapture. \emph{Journal of Wildlife Management} \bold{80}, 1390--1407.

}

\seealso{

\code{\link{secr.fit}}

}

\examples{

## view the function used internally for distances
secrBVN:::anisodistfn

## simulate data
tr <- make.grid(10, 10, spacing = 25, hollow = TRUE, detector = 'proximity')
pop <- simpopn.bvn(s2xy = (25*c(0.5,2))^2, theta = pi/4, core = tr, buffer = 200, 
                   D = 4, Ndist = 'fixed')
CH <- simcapt.bvn(tr, pop, type = 'BVN', lambda = 0.4, noccasions = 5)

## fit model
fit <- anisotropic.fit(CH, buffer = 200, detectfn = 'HHN', trace = FALSE)

## examine results
coef(fit)
predictAniso(fit)
predict(fit)

}

\keyword{ model }
