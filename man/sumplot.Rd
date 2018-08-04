\name{sumplot}
\alias{sumplot}

\title{Summarise and Plot Ellipse Simulations}

\description{

A list of outputs from \code{\link{runEllipseSim}} is summarised and optionally plotted.

}

\usage{

sumplot(simlist, trueD = 4, plt = TRUE, xval = 1:4, xlim = c(0.7, 4.3), ylim = c(-0.2, 0.2),
    legend = TRUE, pchi = c(21, 16, 22, 23), compact = c("av.nCH", "RB", "RSE", "COV"), 
    dec = 3)

}

\arguments{

  \item{simlist}{list, each component an output from runEllipseSim}
  \item{trueD}{numeric true density (animals / ha)}
  \item{plt}{logical; if TRUE a plot is generated}
  \item{xval}{numeric vector of values}
  \item{xlim}{numeric limits of x-axis}
  \item{ylim}{numeric limits of y-axis}
  \item{legend}{logical; if TRUE a legend is displayed}
  \item{pchi}{numeric vector for symbol to use for each component of \code{simlist}}
  \item{compact}{character vector for summary statistics to include in output (NULL for all)}
  \item{dec}{integer number of decimal places in output}

}

\details{

This function is built to meet the needs of Efford in prep. as shown in the package vignette, and may not be useful for novel applications.

Each component of \code{simlist} should have the same number of components (ellipse scenarios).

}

\value{

List of output matrices (tables), returned invisibly if \code{plt = TRUE}.

}

\seealso{

\code{\link{runEllipseSim}}

}

\examples{

\dontrun{

nrepl <- 20
tr <- make.grid(6,6, spacing = 50, detector = 'proximity')
simrandom <- vector('list')
simrandomBVN <- vector('list')
for (i in 1:4) {
    sigmaX <- 25/i^0.5; sigmaY <- 25*i^0.5
    details <- list(distribution = 'binomial')
    ## uniform
    simrandom[[i]] <- runEllipseSim (nrepl, sigmaX, sigmaY, buffer = 200, ncores = 2,
                                          traps = tr, g0 = 0.2, D = 4, type = 'uniform', 
                                          CL = TRUE, detectfn = 'HHN', details = details)
    ## bvn
    simrandomBVN[[i]] <- runEllipseSim (nrepl, sigmaX, sigmaY, buffer = 200, ncores = 2,
                                             traps = tr, lambda0 = 0.4, D = 4, type = 'BVN', 
                                             CL = TRUE, detectfn = 'HHN', details = details)
    message ('Completed aspect ratio ', i)
}

sumplot(list(Uniform = simrandom, BVN = simrandomBVN), legend = TRUE)

print(sumplot(list(Uniform = simrandom36.2, BVN = simrandomBVN36.2), plt = FALSE))

}

}

\keyword{ hplot }
