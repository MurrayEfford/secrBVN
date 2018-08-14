\name{anisotropic.fit}
\alias{anisotropic.fit}

\title{Fit Spatially Explicit Capture--Recapture Model in Transformed Space}

\description{

Using method of Murphy et al. (2016)
}

\usage{

anisotropic.fit(..., psiA = pi/4, psiR = 2)

}

\arguments{
  \item{\dots}{arguments passed to ]\code{\link{secr.fit}}}
  \item{psiA}{numeric initial value of phiA}
  \item{psiR}{numeric initial value of phiR}
}

\details{

%%  ~~ If necessary, more details than the description above ~~

}

\value{

Fitted secr model

}
\references{

Murphy, S. M., Cox, J. J., Augustine, B. C., Hast, J. T., Guthrie, J. M., Wright, J., McDermott, J., Maehr, S. C. and Plaxico, J. H. (2016) Characterizing recolonization by a reintroduced bear population using genetic spatial capture--recapture. \emph{Journal of Wildlife Management} \bold{80}, 1390--1407.

}

\seealso{

\code{\link{secr.fit}}

}

\examples{

## simulate data
tr <- make.grid(10, 10, spacing = 25, hollow = TRUE, detector = 'proximity')
pop <- simpopn.bvn(s2xy = (25*c(0.5,2))^2, theta = pi/4, core = tr, buffer = 200, 
                   D = 4, Ndist = 'fixed')
CH <- simcapt.bvn(tr, pop, type = 'BVN', lambda = 0.4, noccasions = 5)

## fit model
fit <- anisotropic.fit(CH, buffer = 200, detectfn = 'HHN', trace = FALSE)

## examine coefficients
coef(fit)

}

\keyword{ model }
