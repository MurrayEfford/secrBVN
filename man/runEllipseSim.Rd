\name{runEllipseSim}
\alias{runEllipseSim}
\alias{runCircleSim}

\title{Run SECR simulations with Elliptical or Circular Home Ranges}

\description{
A wrapper for other bvn functions.
}

\usage{

runEllipseSim(nrepl = 100, sigmaX = 25, sigmaY = 25, theta = NULL,
    type = c("uniform","BVN"), g0 = 0.2, lambda0 = 0.4, p = 0.95, traps = NULL,
    noccasions = 5, buffer = 100, D = 10, extractfn = derived, seed = NULL, 
    ncores = NULL, Ndist = "fixed",
    secrfn = c("secr.fit", "anisotropic.fit"), ...) 

}

\arguments{
  \item{nrepl}{integer number of replicates}
  \item{sigmaX}{numeric or function}
  \item{sigmaY}{numeric}
  \item{theta}{real orientation of ellipses (radians)}
  \item{type}{uniform ellipse vs bivariate normal}
  \item{g0}{numeric intercept of detection function (type = "uniform")}
  \item{lambda0}{numeric intercept of hazard function (type = "BVN")}
  \item{p}{numeric truncation probability for uniform elliptical home ranges}
  \item{traps}{secr traps object}
  \item{noccasions}{integer number of sampling occasions}
  \item{buffer}{real width of simulation and estimation buffer metres}
  \item{D}{real population density animals/ha}
  \item{extractfn}{ function used to summarise each simulated model fit}
  \item{seed}{integer random number seed}
  \item{ncores}{integer number of cores to use (used in 1.3.0 to setNumThreads)}
  \item{secrfn}{character; name of function to call or NULL}
  \item{Ndist}{character string for distribution of number of individuals (see \code{\link{sim.popn}})}
  \item{...}{other arguments passed to secrfn}
}
\details{

If \code{traps} is not specified, a default 6 x 6 grid of binary proximity detectors is constructed with \code{\link{make.grid}}.

If \code{sigmaX} is a function it is used as in \code{\link{simpopn.bvn}} to generate individual home ranges (population attribute `s2xy').

If \code{sigmaY = NULL} simulates detections of animals with circular halfnormal and uniform home ranges using \code{\link{sim.capthist}}, rather than the ad hoc code in \code{\link{simcapt.bvn}}. This is technically redundant as the same result may be
achieved with sigmaX = sigmaY, but is a useful check.

}

\value{

A list with one component for each replicate, comprising a list with components

\item{npop}{number of individuals in simulated population}
\item{nCH}{number of individuals detected}
\item{ncapt}{number of detections}

And, if \code{secrfn} not NULL,

\item{fit}{result of extractfn applied to fitted secr model (default, a 2 x 4 matrix from derived}
\item{pred}{result of predict applied to fitted secr model (default, a 2 x 5 data.frame)}

}

\seealso{
  \code{\link{simpopn.bvn}}, 
  \code{\link{simcapt.bvn}}, 
  \code{\link{secr.fit}},
  \code{\link{anisotropic.fit}}
}

\examples{
\dontrun{
nrepl <- 20
tr <- make.grid(6,6, spacing = 50, detector = "proximity")
sigmaX <- 25/2; sigmaY <- 25*2
details <- list(distribution = "binomial")
sims <- runEllipseSim (nrepl, sigmaX, sigmaY, buffer = 200, ncores = 2, 
                       traps = tr, g0 = 0.2, D = 4, detectfngen = 4, CL = TRUE, 
                       details = details) 
}
}

\keyword{ datagen }