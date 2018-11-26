################################################################################
## package 'openCR'
## openCR.esa.R
## 2011-12-30, 2013-01-20
################################################################################


# CJSsecr = 6
# JSSAsecrf = 7
# JSSAsecrD = 8
# JSSAsecrfCL = 9
# JSSAsecrlCL = 10
# JSSAsecrbCL = 11
# JSSAsecrl = 12
# JSSAsecrb = 13
# JSSAsecrB = 14
# JSSAsecrg = 24
# JSSAsecrgCL = 25

openCR.esa <- function (object, bysession = FALSE) {

    if (!(grepl('secr',object$type) & inherits(object,'openCR')) )
        stop ("requires fitted openCR secr model")
    # Return the esa
    beta     <- object$fit$par
    parindx  <- object$parindx
    link     <- object$link
    fixed    <- object$fixed
    capthist <- object$capthist
    mask     <- object$mask
    details  <- object$details
    type     <- object$type
    intervals <- object$intervals
    cumss    <- getcumss(object$capthist)
    movementmodel <- object$movementmodel
    binomN <- object$binomN
    if (is.null(binomN)) binomN <- 1
    ## openCR >= 1.2.0
    individual <- object$design0$individual
    ## openCR < 1.2.0
    if (is.null(individual)) individual <- individualcovariates(object$design0$PIA)

    ###### movement kernel and related
    cellsize <- mqarray <- 0
    kernel <- data.frame(rownames=NULL)
    if (movementmodel %in% c('normal','exponential','user')) {
        k2 <- details$kernelradius
        cellsize <- attr(mask,'area')^0.5 * 100   ## metres, equal mask cellsize
        kernel <- expand.grid(x = -k2:k2, y = -k2:k2)
        kernel <- kernel[(kernel$x^2 + kernel$y^2) <= (k2+0.5)^2, ]
        mqarray <- mqsetup (mask, kernel, cellsize)
    }
    movemodel <- switch (movementmodel, static = 0, uncorrelated = 1, normal = 2, 
                         exponential = 3, user = 4) 

    #--------------------------------------------------------------------
    # Fixed beta
    fb <- details$fixedbeta
    if (!is.null(fb)) {
        fb[is.na(fb)] <- beta
        beta <- fb    ## complete
    }
    #--------------------------------------------------------------------
    # Real parameters

    realparval0 <- makerealparameters (object$design0, beta, parindx, link, fixed)
    nc  <- nrow(capthist)
    J <- length(intervals) + 1

    type <- switch(type, CJSsecr = 6, JSSAsecrfCL = 9, JSSAsecrlCL = 10, JSSAsecrbCL = 11, JSSAsecrgCL = 25,
                   JSSAsecrf = 7, JSSAsecrl = 12, JSSAsecrb = 13, JSSAsecrB = 14, JSSAsecrg = 24,
                   JSSAsecrD = 8)

    trps <- traps(capthist)
    k <- nrow(trps)
    m <- nrow(mask)
    if (!is.null(mask)) area <- attr(mask,'area')
    else area <- 0
    distrib <- switch (object$distribution, poisson = 0, binomial = 1)

    usge <- usage(traps(capthist))
    if (is.null(usge)) usge <- matrix(1, nrow = k, ncol = ncol(capthist))

    if (object$ncores>1)
        temp <- makegkParallelcpp (as.integer(object$detectfn), 
                               as.integer(.openCRstuff$sigmai[type]),
                               as.integer(details$grain),
                               as.matrix(realparval0),
                               as.matrix(trps),
                               as.matrix(object$mask))
    else
        temp <- makegkcpp(
            as.integer(nrow(realparval0)),     # number of rows in lookup table
            as.integer(k),                     # detectors
            as.integer(m),                     # mask points
            as.integer(object$detectfn),
            as.integer(.openCRstuff$sigmai[type]),                # column index of sigma in realparval
            as.double(realparval0),
            as.double(unlist(trps)),
            as.double(unlist(mask)))
        
    
    gk0 <- temp[[1]]

    ## mixture proportions
    if (details$nmix > 1) {
        # temp <- fillpmix3cpp(   CHECK 2018-02-02
        temp <- fillpmix2(nc, details$nmix, object$design0$PIA, realparval0)
        pmix <- matrix(temp, ncol = nc)
    }
    else {
        pmix <- matrix(1, nrow = details$nmix, ncol = nc)
    }
    onea <- function (x) {
        if (object$ncores == 1) {
            pch1 <-  PCH1secrcpp(
                as.integer(type),
                as.integer(x-1),
                as.integer(nc),
                as.integer(J),
                as.integer(cumss),
                as.integer(k),
                as.integer(m),
                as.matrix(realparval0),
                as.integer(object$design0$PIA),
                as.integer(object$design$PIAJ),
                as.double (gk0),
                as.integer(binomN),
                as.matrix (usge),
                as.double (intervals),
                as.integer(object$moveargsi),
                as.integer(movemodel),
                as.character(object$usermodel),
                as.matrix(kernel),
                as.matrix(mqarray),
                as.double(cellsize))
        }
        else {
            pch1 <-  PCH1secrparallelcpp(
                as.integer(x-1),
                as.integer(type),
                as.integer(object$details$grain),
                as.logical(individual),
                as.integer(J),
                as.integer(m),
                as.integer(nc),
                as.integer(cumss),
                as.matrix (realparval0),
                as.integer(object$design0$PIA),
                as.integer(object$design0$PIAJ),
                as.double (gk0),
                as.integer(binomN),
                as.matrix (usge),
                as.double (intervals),
                as.integer(object$moveargsi),
                as.integer(movemodel),
                as.matrix(kernel),
                as.matrix(mqarray),
                as.double (cellsize))
        }
        pmix[x,] * pch1
    }
    
    oneabysession <- function (x) {
        a0 <- PCH0secrjcpp (
            as.integer(type),
            as.integer(x-1),
            as.integer(nc),
            as.integer(J),
            as.integer(cumss),
            as.integer(k),
            as.integer(m),
            as.integer(nrow(realparval0)),
            as.integer(object$design0$PIA),
            as.double (gk0),
            as.integer(binomN),
            as.matrix (usge)) 
        a0 <- matrix(a0, nrow = nc, ncol = J)
        1-sweep(a0, MARGIN=1, STATS=pmix[x,], FUN="*")
    }
    if (bysession) {
        a <- sapply(1:details$nmix, oneabysession, simplify = FALSE)
        a <- t(apply(abind(a, along=3), 1:2, sum))   ## session * ch, summed over mixture classes
    }
    else {
        a <- sapply(1:details$nmix, onea)
        a <- apply(a, 1, sum)    ## sum over latent classes
    }
    a * area * m
}