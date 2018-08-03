## modified 2018-07-26 to also report predicted values 
runEllipseSim <- function(nrepl = 100, sigmaX = 25, sigmaY = 25, theta = NULL, type = c('uniform', 'BVN'), 
                          g0 = 0.2, lambda0 = 0.4, p = 0.95, traps = NULL, noccasions = 5, buffer = 100, D = 10, 
                          extractfn = derived, seed = NULL, ncores = 1, SECR = TRUE, Ndist = 'fixed', ...) {
    onereplicate <- function (r) {
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
        if (SECR) {
            fit <- try(secr.fit(CH, mask = mask, trace = FALSE, ...))
            if (inherits(fit, 'try-error')) {
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
    
    if (is.null(traps))
        traps <- make.grid(detector = 'proximity')
    mask <- make.mask(traps, buffer = buffer, type = 'trapbuffer')
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
    
    if (ncores == 1) {
        set.seed(seed)        
        output <- lapply(1:nrepl, onereplicate)
    }
    else {
        ## use 'parallel'
        list(...)    ## ensures promises evaluated see parallel vignette 2015-02-02        
        clust <- makeCluster(ncores)
        clusterSetRNGStream(clust, seed)                
        output <- parLapply (clust, 1:nrepl, onereplicate)
        stopCluster(clust)        
    }
    output
}