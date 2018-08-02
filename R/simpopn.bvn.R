simpopn.bvn <- function (s2xy = c(300,300), theta = NULL, ...) {
    popn <- sim.popn (...)
    N <- nrow(popn)
    if (is.null(theta))
        theta <- runif(N, max = 2 * pi)
    else {
        ## negative theta to draw random shared theta
        if (theta<0) theta <- runif(1, max=2*pi)
        theta <- rep(theta,N)
    }
    attr(popn, 'theta') <- theta
    if (length(s2xy) == 2)
        attr(popn, 's2xy') <- matrix(rep(s2xy, each = N), ncol = 2)
    else if (is.function(s2xy)) {
        ## for heterogeneity
        attr(popn, 's2xy') <- s2xy(N) # function(N) c(625,625) * (0.5 + runif(N*2))^2
    }
    else
        stop ("require length 2 vector for s2xy)")
    popn
}
