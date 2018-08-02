## temporarily redefine subset.popn to respect s2xy, theta attributes (not needed)
subset.popn <- function (x, subset = NULL, sessions = NULL, poly = NULL,
                         poly.habitat = TRUE, keep.poly = TRUE, renumber = FALSE, ...)
    ## x - popn object
    ## subset - vector (character, logical or integer) to subscript rows (dim(1)) of x
    ## sessions - vector (integer or logical) to subscript sessions
{
    if (ms(x)) {
        if (!is.null(sessions))
            x <- x[sessions]
        out <- vector('list')
        for (i in 1:length(x)) {
            if (is.list(subset))
                sset <- subset[[i]]
            else
                sset <- subset
            if (is.null(subset))  ## 2015-03-17
                out[[i]] <- x[[i]]
            else
                out[[i]] <- subset(x[[i]], subset[[i]], NULL, renumber, ...)
        }
        class(out) <- c('list','popn')
        out
    }
    else {
        #-------------------------
        # default subset is all
        if (is.null(subset))
            subset <- 1:nrow(x)
        #-------------------------
        # restrict to a polygon
        # added 2011-10-20
        
        if (!is.null(poly)) {
            OK <- pointsInPolygon(x, poly)
            if (!poly.habitat)
                OK <- !OK
            subset <- subset[OK]
        }
        #-------------------------
        # apply subsetting
        pop <- x[subset,]
        #-------------------------
        if (renumber)
            rownames(pop) <- 1 : nrow(pop)
        
        class(pop) <- c('popn', 'data.frame')
        attr(pop, 'Ndist') <- NULL     ## no longer known
        attr(pop, 'model2D') <- NULL   ## no longer known
        attr(pop, 'boundingbox') <- attr(x, 'boundingbox')
        if (!is.null(poly) & keep.poly) {
            attr(pop, 'polygon') <- poly
            attr(pop, 'poly.habitat') <- poly.habitat
        }
        if (!is.null(covariates(x))) {
            covariates(pop) <- covariates(x)[subset,,drop=FALSE]
        }
        attr(pop, 's2xy') <- attr(x, 's2xy')[subset, ,drop = FALSE]
        attr(pop, 'theta') <- attr(x, 'theta')[subset]
        pop
    }
}
