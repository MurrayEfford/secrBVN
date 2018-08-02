## random ellipse utility functions

inside.ellipse <- function (xy, centres, theta, s2xy, p = 0.95) {
    xy <- unlist(xy)
    xyc <- centres
    xyc[,] <- NA
    for (i in 1:nrow(centres)) {
        xyc[i,] <- xy - centres[i,]
        xyc[i,] <- rotate (xyc[i,,drop=F], -360/2/pi*theta[i])
    }
    inside <- function (xy) {
        sum(xy^2/s2xy) < circular.r(p)^2
    }
    any(apply(xyc,1,inside))
}

pmvn <- function (xy, centres, theta, s2xy) {
    xy <- unlist(xy)
    xyc <- centres
    xyc[,] <- NA
    maxd <- 1/(2 * pi * prod(s2xy)^0.5)
    d <- logical(nrow(centres))
    for (i in 1:nrow(centres)) {
        xyc[i,] <- xy - centres[i,]
        xyc[i,] <- rotate (xyc[i,,drop=F], -360/2/pi*theta[i])
        d[i] <- runif(1) < (mvtnorm::dmvnorm(xyc[i,], mean=c(0,0),
                                    sigma = matrix(c(s2xy[1],0,0,s2xy[2]), nrow = 2))/maxd)
    }
    any(d)
}

ellipse <- function (center, shape, radius, segments = 51) {
    ## based on ellipse function in 'car' package
    angles <- (0:segments) * 2 * pi/segments
    unit.circle <- cbind(cos(angles), sin(angles))
    Q <- chol(shape, pivot = TRUE)
    order <- order(attr(Q, "pivot"))
    ellipse <- t(center + radius * t(unit.circle %*% Q[, order]))
    colnames(ellipse) <- c("x", "y")
    ellipse
}

random.ellipses <- function(popn, plt = T, s2xy = c(0.5,5), radius = 2.45) {
    popn <- as.matrix(popn)
    N <- nrow(popn)
    cov <- matrix(c(s2xy[1],0,0,s2xy[2]), ncol = 2)
    ## random orientation
    theta <- runif (N, max=2 * pi)
    for (i in 1:N) {
        ellipsi <- ellipse(popn[i,], shape=cov, radius = radius)
        ellipsi <- rotate (ellipsi, theta[i]*360/2/pi, centrexy = popn[i,])
        polygon (ellipsi)
    }

    xy <- cbind (runif(1000, min=par()$usr[1], max=par()$usr[2]),
                 runif(1000, min=par()$usr[3], max=par()$usr[4]))
    points(xy)
    OK <- apply(xy, 1, inside.ellipse, centres=popn, theta=theta, s2xy=s2xy)
    OK <- apply(xy, 1, pmvn, centres=popn, theta=theta, s2xy=s2xy)
    points(xy[OK,], pch=16)
}


