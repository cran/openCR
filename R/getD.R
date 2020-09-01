## uses untransform()
getD <- function (designD, beta, mask, parindx, link, fixed,
                  grouplevels, sessionlevels, parameter = 'D') {
    if ((is.null(designD) | nrow(designD)==0) & (is.null(fixed[[parameter]]))) return(NULL)
    if (ms(mask))
        nmask <- max(sapply(mask, nrow))
    else
        nmask <- nrow(mask)
    ngroup <- length(grouplevels)
    nsession <- length(sessionlevels)
    D <- array(dim = c(nmask, ngroup, nsession))
    dimnames(D) <- list(1:nrow(D), grouplevels, sessionlevels)
    if (!is.null(fixed[[parameter]])) {
        D[,,] <- fixed[[parameter]]
    }
    else {
        D[,,] <- designD %*% beta[parindx[[parameter]]]   # linear predictor
        D[,,] <- untransform (D, link[[parameter]])
        # silently truncate D at zero
        if (parameter %in% c('D')) D[D<0] <- 0
    }
    D
}
###############################################################################

## from secrloglik

    # D.modelled <- !CL & is.null(fixedpar$D)
    # if (!CL ) {
    #     grplevels <- if (length(grp)>0) levels(grp[[1]]) else 1
    #   D <- getD (designD, beta, mask, parindx, link, fixedpar,
    #              grplevels, sessionlevels, parameter = 'D')
    # 
    #   if (!is.na(sumD <- sum(D)))
    #       if (sumD <= 0)
    #           warning ("invalid density <= 0")
    # }

