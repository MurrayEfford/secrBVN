simcapt.bvn <- function (traps, popn, type = c('uniform', 'BVN'), g0 = 0.2, lambda0 = 0.2, p = 0.95,
                         noccasions = 5, renumber = TRUE) {
    
    traps.inside.ellipse <- function (traps, centre, theta, s2xy, p) {
        xy <- sweep (traps, MARGIN = 2, STATS = centre, FUN = '-')
        xy <- rotate (xy, -360/2/pi*theta)
        inside <- function (xy) {
            sum(xy^2/s2xy) < circular.r(p, hazard = FALSE)^2   ## 2.45 default
        }
        apply(xy,1,inside)
    }
    traps.pmvn <- function (traps, centre, theta, s2xy) {
        xy <- sweep (traps, MARGIN = 2, STATS = centre, FUN = '-')
        xy <- rotate (xy, -360/2/pi*theta)
        maxd <- 1/(2 * pi * prod(s2xy)^0.5)
        apply (xy, 1, mvtnorm::dmvnorm, mean = c(0,0),
               sigma = matrix(c(s2xy[1],0,0,s2xy[2]), nrow = 2)) / maxd
    }
    
    ## mainline
    if (!(detector(traps) %in% c('multi','proximity')) ){
        warning ("detectors other than multi or proximity will be treated as proximity detectors")
        detector(traps) <- 'proximity'
    }
    type <- match.arg(type)
    N <- nrow(popn)
    K <- nrow(traps)
    w <- array(0, dim = c(N,noccasions,K))
    s2xy <- attr(popn, 's2xy')
    theta <- attr(popn, 'theta')
    popn <- as.matrix(popn)
    for (i in 1:N) {
        if (type == 'uniform')  ## uniform ellipse
            OK <- traps.inside.ellipse( traps = traps, centre = popn[i,],
                                        theta = theta[i], s2xy = s2xy[i,], p = p )
        else                    ## BVN
            OK <- traps.pmvn( traps = traps, centre = popn[i,],
                              theta = theta[i], s2xy = s2xy[i,])
        for (s in 1:noccasions) {
            w[i,s,] <- OK   ## constant over occasions
        }
    }
    ## w remains numeric array dim = c(N,noccasions,K)
    if (detector(traps)=='proximity') {
        if (type == 'uniform') 
            g <- w * g0 ## 2018-07-26, for clarity
        else
            g <- 1 - exp(- w * lambda0)
        w <- 1 * (runif(prod(dim(w))) < g)  
        ## drop null capture histories
        caught <- apply(w,1,sum) > 0
        w <- w[caught,,,drop = FALSE]
    }
    else {
        # multi-catch traps
        # use hazard of detection
        if (type == 'BVN') 
            w <- w * lambda0
        else
            w <- w * -log(1-g0)
        sumhaz <- apply(w, 1:2, sum)
        w <- w / as.numeric(sumhaz)  # relative hazard
        w[is.na(w)] <- 0 # arises with uniform when zero overall Pr(capt)
        caughts <- (1 - exp(-sumhaz)) > runif(prod(dim(sumhaz))) 
        ## drop null capture histories
        caught <- apply(caughts,1,sum) > 0
        w <- w[caught,,,drop = FALSE]
        caughts <- caughts[caught,,drop = FALSE]
        # select one of traps with probability in w 
        # (w provides prob argument of sample.int)
        w <- caughts * apply(w, 1:2, sample.int, n=K, size=1, replace = FALSE)
    }
    
    n <- nrow(w)
    class(w) <- 'capthist'
    traps(w) <- traps
    session(w)           <- '1'           ## dummy session values for now
    if (renumber && (n>0)) rownames(w) <- 1:n
    else {
        rownames(w) <- rownames(popn)[caught]
    }
    w
}
