## modified 2018-07-26 to also report predicted values 
## modified 2022-10-11 for random traps
runEllipseSim <- function(nrepl = 100, sigmaX = 25, sigmaY = 25, theta = NULL, 
    type = c('uniform', 'BVN'), g0 = 0.2, lambda0 = 0.4, p = 0.95, 
    traps = NULL, trapargs = NULL, noccasions = 5, buffer = 100, 
    D = 10, extractfn = derived, seed = NULL, ncores = NULL, 
    Ndist = 'fixed', secrfn = c('secr.fit', 'anisotropic.fit'),...) {
    onereplicate <- function (r) {

        if (is.function(traps)) {
            # function to deliver random trap locations
            # local 'traps' does not alter that argument of outer fn
            traps <- do.call(traps, trapargs)
            mask <- make.mask(traps, buffer = buffer, type = 'trapbuffer')
        }
    
        if (is.null(sigmaY) & !is.function(sigmaX)) {
            ## simple and direct circular using conventional secr::sim.capthist
            pop <- sim.popn(core = traps, buffer = buffer, D = D, Ndist = Ndist)
            CH <- sim.capthist(traps, pop, detectfngen, detectpar, noccasions)            
        }
        else {
            pop <- simpopn.bvn(s2xy = s2xy, theta = theta, core = traps, buffer = buffer, 
                               D = D, Ndist = Ndist)
            CH <- simcapt.bvn(traps, pop, type, g0, lambda0, p, noccasions)
        }
        out <- list(npop = nrow(pop), nCH = nrow(CH), ncapt=sum(CH), fit = nullfit)
        
        if (!is.null(secrfn)) {
            args <- list(...)
            args$capthist <- CH
            args$mask <- mask
            args$trace <- FALSE
            fit <- try(do.call(secrfn, args))
            if (!inherits(fit, "secr")) {
                out$fit <- nullfit
                out$pred <- nullpred
            }
            else {
                out$fit <- extractfn(fit)
                out$pred <- predict(fit)
            }
        }
        cat ('Completed ', r, date(), '\n')
        flush.console()
        out
    }
    secrfn <- match.arg(secrfn)
    if (is.null(traps)) {
        traps <- make.grid(detector = 'proximity')
    }
    else if (!is.function(traps)) {
        mask <- make.mask(traps, buffer = buffer, type = 'trapbuffer')
    }
    columns <- c('estimate','SE.estimate','lcl','ucl')
    
    type <- match.arg(type)
    nullfit <- as.data.frame(matrix(NA, nrow = 2, ncol = 4, 
                                    dimnames = list(c('esa','D'), columns)))
    nullpred <- as.data.frame(matrix(NA, nrow = 2, ncol = 5, 
                                    dimnames = list(c('lambda0','sigma'), c('link',columns))))
    
    ## for simcapt.bvn
    if (is.function(sigmaX))
        s2xy <- sigmaX
    else
        s2xy <- c(sigmaX, sigmaY)^2
    
    ## define detectpar etc. in case needed for sim.capthist
    detectfngen <- if (type=='uniform') 4 else 14
    if (is.null(sigmaY) & !is.function(sigmaX)) {
        if (detectfngen == 4)
            detectpar <- list(g0 = g0, sigma = sigmaX * circular.r(p, hazard = FALSE))
        else
            detectpar <- list(lambda0 = lambda0, sigma = sigmaX)
    }
    
    if (!is.null(ncores)) {
        setNumThreads(ncores)
    }
    
    set.seed(seed)        
    lapply(1:nrepl, onereplicate)
}