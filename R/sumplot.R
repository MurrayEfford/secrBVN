simsum <- function (simlist, component = c("fit", "pred"), parm = "D", trueval = 4, 
                    compact = c('av.nCH', 'RB', 'RSE', 'rRMSE', 'COV'), dec = 3) {
    sumD <- function(x) {
        Dval <- sapply(lapply(x, '[[', component), '[[', parm, 'estimate')
        DSE  <- sapply(lapply(x, '[[', component), '[[', parm, 'SE.estimate')
        Dlcl <- sapply(lapply(x, '[[', component), '[[', parm, 'lcl')
        Ducl <- sapply(lapply(x, '[[', component), '[[', parm, 'ucl')
        n <- sum(!is.na(Dval))
        RB <- (Dval-trueval)/trueval
        RSE <- DSE/Dval
        COV <- (trueval>=Dlcl) & (trueval<=Ducl)
        npop <- sapply(x, '[[', 'npop')
        nCH <- sapply(x, '[[', 'nCH')
        c(av.npop = mean(npop), 
          av.nCH = mean(nCH), 
          nvalid = n,
          av.parmhat = mean(Dval, na.rm = T),
          md.parmhat = median(Dval, na.rm = T), 
          sd.parmhat = sd(Dval, na.rm = T),
          RB = mean(RB, na.rm = T), 
          seRB = sd(Dval, na.rm = T)/trueval/n^0.5,
          RSE = mean(RSE, na.rm = T), 
          seRSE = sd(RSE, na.rm = T)/n^0.5,
          rRMSE = sqrt(mean((Dval-trueval)^2)) / trueval,
          COV = mean(COV, na.rm = T))
    }
    tidy <- function (x) {
        tmp <- sapply(x, sumD)
        colnames(tmp) <- names(x)
        tmp
    }
    component <- match.arg(component)
    outlist <- lapply(simlist, tidy)
    if (is.null(compact)) {
        out <- lapply(outlist,round, dec)
    }
    else {
        comp <- function (out) {
            round(t(out)[,compact], dec)
        }
        out <- lapply(outlist, comp)
    }
    out
}

simplot <- function (simlist, component = c("fit", "pred"), parm = "D", trueval = 4, 
                     xval = 1:4, xlim = c(0.7,4.3), ylim = c(-0.2,0.2),
                     legend = TRUE, pchi = c(21, 16, 22, 24), cexi = rep(1.2,4)) {
    plotone <- function (out, i) {
        segments(xval+offset[i], out['RB',]-2*out['seRB',],
                 xval+offset[i], out['RB',]+2*out['seRB',])
        points(xval+offset[i], out['RB',], pch = pchi[i], bg = 'white', cex = cexi[i])
    }
    component <- match.arg(component)
    outlist <- simsum(simlist, component = component, parm = parm, 
                      trueval = trueval, compact = NULL)
    nscen <- length(simlist)
    offset <- seq(-(nscen-1)*0.025, (nscen-1)*0.025, 0.05)
    plot(0,0, xlim=xlim, ylim = ylim, xlab = 'Aspect ratio',
         ylab = 'Relative bias',axes = FALSE, type = 'n')
    mapply(plotone, outlist, 1:length(outlist))
    axis(1, at = 1:4)
    axis(2, las = 1)
    abline(h=0, lty=2)
    if (legend)
        legend ((par()$usr[2] - par()$usr[1]) * 0.1 + par()$usr[1],  
                par()$usr[4]*0.95, legend = names(simlist),
                pch = pchi, cex = par()$cex * 0.8, pt.cex = cexi, adj = 0)
    invisible(outlist)
}