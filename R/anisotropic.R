anisodistfn <- function (xy1, xy2, mask) {
    if (missing(xy1)) return(character(0))
    xy1 <- as.matrix(xy1)
    xy2 <- as.matrix(xy2)
    miscparm <- attr(mask, 'miscparm')
    psiA <- miscparm[1]           # anisotropy angle; identity link
    psiR <- 1 + exp(miscparm[2])  # anisotropy ratio; log link
    aniso.xy1 <- geoR::coords.aniso(xy1, aniso.pars = c(psiA, psiR))
    aniso.xy2 <- geoR::coords.aniso(xy2, aniso.pars = c(psiA, psiR))
    secr::edist(aniso.xy1, aniso.xy2) # nrow(xy1) x nrow(xy2) matrix
}

predictAniso <- function (fit, angle = c("degrees", "radians")) {
    angle <- match.arg(angle)
    co <- coef(fit)
    if (!all( c("psiA", "psiR")  %in% rownames(co)))
        stop ("input is not anisotropic.fit")
    pred <- co[c("psiA", "psiR"), ]
    colnames(pred)[1:2] <- c("estimate", "SE.estimate")
    if (angle == "degrees")
        pred["psiA", ] <- pred["psiA", ] * 360 / 2 / pi
    pred["psiR", ] <- 1 + exp(pred["psiR", ])  
    beta <- co["psiR","beta"]
    sebeta <- co["psiR","SE.beta"]
    pred["psiR", "SE.estimate"] <- exp(beta) * sqrt(exp(sebeta^2)-1)
    pred
}

anisotropic.fit <- function (..., psiA = pi/4, psiR = 2) {
    args <- list(...)
    if (is.null(args$details))
        args$details <- vector('list')
    args$details$userdist <- anisodistfn
    args$details$miscparm <- c(psiA = psiA, psiR = log(psiR-1))
    tmp <- do.call("secr.fit", args)
    tmp$call <- NULL ## drop bulky evaluated call
    tmp
}