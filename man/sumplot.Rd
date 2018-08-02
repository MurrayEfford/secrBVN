\name{sumplot}
\alias{sumplot}

\title{Summarise and Plot Ellipse Simulations}

\description{

A list of outputs from \code{\link{runEllipseSim}} is summarised and optionally plotted.

}

\usage{

sumplot(simlist, trueD = 4, plt = TRUE, xval = 1:4, xlim = c(0.7, 4.3), ylim = c(-0.2, 0.2),
    legend = TRUE, pchi = c(21, 16, 22, 23), compact = c("av.nCH", "RB", "RSE", "COV"), 
    ndec = 3)

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
  \item{ndec}{integer number of decimal places in output}

}

\details{

This function is built to meet the needs of Efford in prep. as shown in the package vignette, and may not be useful for novel applications.

Each component of \code{simlist} should have the same number of components (ellipse scenarios).

}

\value{

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (simlist, trueD = 4, plt = TRUE, xval = 1:4, xlim = c(0.7, 
    4.3), ylim = c(-0.2, 0.2), legend = TRUE, colnames = NULL, 
    pchi = c(21, 16, 22, 23), compact = c("av.nCH", "RB", "RSE", 
        "COV")) 
{
    sumD <- function(x) {
        Dval <- sapply(lapply(x, "[[", "fit"), "[[", "D", "estimate")
        DSE <- sapply(lapply(x, "[[", "fit"), "[[", "D", "SE.estimate")
        Dlcl <- sapply(lapply(x, "[[", "fit"), "[[", "D", "lcl")
        Ducl <- sapply(lapply(x, "[[", "fit"), "[[", "D", "ucl")
        sigval <- sapply(lapply(x, "[[", "pred"), "[[", "sigma", 
            "estimate")
        n <- sum(!is.na(Dval))
        RB <- (Dval - trueD)/trueD
        RSE <- DSE/Dval
        COV <- (trueD >= Dlcl) & (trueD <= Ducl)
        npop <- sapply(x, "[[", "npop")
        nCH <- sapply(x, "[[", "nCH")
        c(av.npop = mean(npop), av.nCH = mean(nCH), nvalid = n, 
            av.Dhat = mean(Dval, na.rm = T), md.Dhat = median(Dval, 
                na.rm = T), sd.Dhat = sd(Dval, na.rm = T), RB = mean(RB, 
                na.rm = T), seRB = sd(Dval, na.rm = T)/trueD/n^0.5, 
            RSE = mean(RSE, na.rm = T), seRSE = sd(RSE, na.rm = T)/n^0.5, 
            COV = mean(COV, na.rm = T), av.sighat = mean(sigval, 
                na.rm = T), se.sighat = sd(sigval, na.rm = T)/sqrt(n))
    }
    tidy <- function(x) {
        tmp <- sapply(x, sumD)
        colnames(tmp) <- colnames
        tmp
    }
    outlist <- lapply(simlist, tidy)
    if (plt) {
        plotone <- function(out, i) {
            segments(xval + offset[i], out["RB", ] - 2 * out["seRB", 
                ], xval + offset[i], out["RB", ] + 2 * out["seRB", 
                ])
            points(xval + offset[i], out["RB", ], pch = pchi[i], 
                bg = "white", cex = 1.2)
        }
        nscen <- length(simlist)
        offset <- seq(-(nscen - 1) * 0.025, (nscen - 1) * 0.025, 
            0.05)
        plot(0, 0, xlim = xlim, ylim = ylim, xlab = "Aspect ratio", 
            ylab = "Relative bias", axes = FALSE, type = "n")
        mapply(plotone, outlist, 1:length(outlist))
        axis(1, at = 1:4)
        axis(2, las = 1)
        abline(h = 0, lty = 2)
        if (legend) 
            legend(par()$usr[2] * 0.7, par()$usr[4] * 0.95, legend = names(simlist), 
                pch = pchi, cex = 0.9, pt.cex = 1.2)
    }
    if (is.null(compact)) {
        lapply(outlist, round, 4)
    }
    else {
        comp <- function(out) {
            round(t(out)[, compact], 3)
        }
        lapply(outlist, comp)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
