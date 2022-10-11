plotpopn.bvn <- function (popn, radius = 2.45, add = FALSE, region = NULL,  ...) {
    if (!is.null(region))
        requireNamespace('sp')
    if (!is.null(region) & !inherits(region, 'SpatialPolygons'))
        region <- SpatialPolygons(list(Polygons(list(Polygon(as.matrix(region))), ID=1)))
    N <- nrow(popn)
    s2xy <- attr(popn, 's2xy')
    theta <- attr(popn, 'theta')
    if (!add) plot(popn)
    popn <- as.matrix(popn)
    for (i in 1:N) {
        cov <- matrix(c(s2xy[i,1],0,0,s2xy[i,2]), ncol = 2)
        ellipsi <- ellipse(popn[i,], shape = cov, radius = radius)
        ellipsi <- rotate (ellipsi, theta[i]*360/2/pi, centrexy = popn[i,])
        ellipsesp <- SpatialPolygons(list(Polygons(list(Polygon(as.matrix(ellipsi))), ID=1)))
        if (!is.null(region)) {
            if (!requireNamespace("rgeos")) stop ("gIntersection uses rgeos that is unavailable")
            ellipsi2 <- gIntersection(ellipsesp, region)
            if (!is.null(ellipsi2))
                plot(ellipsi2, add = TRUE, ...)
            polygon(ellipsi)
        }
        else {
            polygon (ellipsi, ...)
        }
    }
}

