################################################################################
## package 'openCR'
## openCR.esa.R
## 2019-05-17
################################################################################


# CJSsecr = 6
# JSSAsecrf = 7
# JSSAsecrD = 8
# JSSAsecrfCL = PLBsecrf = 9
# JSSAsecrlCL = PLBsecrl = 10
# JSSAsecrbCL = PLBsecrb = 11
# JSSAsecrl = 12
# JSSAsecrb = 13
# JSSAsecrB = 14
# JSSAsecrg = 24
# JSSAsecrgCL = PLBsecrg = 25

# secrCL = 30
# secrD = 31

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
    binomN <- object$binomN
    if (is.null(binomN)) binomN <- 1
    ## openCR >= 1.2.0
    individual <- object$design0$individual
    ## openCR < 1.2.0
    if (is.null(individual)) individual <- individualcovariates(object$design0$PIA)

    ###### movement kernel and related
    cellsize <- mqarray <- 0
    kernel <- data.frame(rownames = NULL)
    movementcode <- movecode(object$movementmodel)
    edgecode <- edgemethodcode(object$edgemethod)
    if (object$movementmodel %in% c('normal','exponential','user','t2D')) {
        k2 <- details$kernelradius
        cellsize <- attr(mask,'area')^0.5 * 100   ## metres, equal mask cellsize
        kernel <- expand.grid(x = -k2:k2, y = -k2:k2)
        kernel <- kernel[(kernel$x^2 + kernel$y^2) <= (k2+0.5)^2, ]
        mqarray <- mqsetup (mask, kernel, cellsize, edgecode)   
    }

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

    type <- typecode(type)  # convert to integer
    trps <- traps(capthist)
    k <- nrow(trps)
    m <- nrow(mask)
    if (!is.null(mask)) area <- attr(mask,'area')
    else area <- 0
    distrib <- switch (object$distribution, poisson = 0, binomial = 1)

    detectr <- detector(trps)[1]
    if (detectr == 'count') {
        detectr <- if (object$binomN == 0) "poissoncount" else "binomialcount"
    }

    usge <- usage(trps)
    if (is.null(usge)) usge <- matrix(1, nrow = k, ncol = ncol(capthist))

    temp <- makegkParallelcpp (as.integer(object$detectfn), 
                               as.integer(.openCRstuff$sigmai[type]),
                               as.integer(details$grain),
                               as.matrix(realparval0),
                               as.matrix(trps),
                               as.matrix(object$mask))
    gk0 <- array(temp[[1]], dim=c(nrow(realparval0), k, m)) # 2020-10-28 as array
    hk0 <- array(temp[[2]], dim=c(nrow(realparval0), k, m)) # 2020-10-28 as array

    ## mixture proportions
    if (details$nmix > 1) {
        temp <- fillpmix2(nc, details$nmix, object$design0$PIA, realparval0)
        pmix <- matrix(temp, ncol = nc)
    }
    else {
        pmix <- matrix(1, nrow = details$nmix, ncol = nc)
    }
    onea <- function (x) {
        if (object$details$R) {
            pch1 <-  PCH1secr(
                type,
                as.logical(object$design0$individual),
                x,
                nc,
                J,
                cumss,
                k,
                m,
                realparval0,
                object$design0$PIA,
                object$design0$PIAJ,
                if (detectr %in% c("multi", "poissoncount")) hk0 else gk0,
                binomN,
                usge,
                intervals,
                object$moveargsi,
                movementcode,                     # what about edgecode? 2020-12-12
                object$usermodel,
                kernel,
                mqarray,
                cellsize)
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
                as.double (if (detectr == "poissoncount") hk0 else gk0),  ## 2019-05-08, 17
                as.integer(binomN),
                as.matrix (usge),
                as.double (intervals),
                as.integer(object$moveargsi),
                as.integer(movementcode),
                as.integer(edgecode),
                as.character(object$usermodel),
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
        a <- t(apply(abind(a, along=3), 1:2, sum))   ## session * ch, summed over latent classes
        a <- a * area * m
        sess <- primarysessions(intervals(object$capthist)) ## differs from 'intervals'
        OK <- apply(object$capthist, 1, by, sess, sum)>0
        freq <- sweep(OK, MARGIN=2, STATS=covariates(object$capthist)$freq, FUN = "*")
        a <- lapply(1:J, function(j) rep(a[j,], freq[j,]))
    }
    else {
        a <- sapply(1:details$nmix, onea)
        a <- apply(a, 1, sum)    ## sum over latent classes
        a <- rep(a, covariates(object$capthist)$freq)
        a <- a * area * m
    }
    ## adjusts for covariates(object$capthist)$freq - no need to do this downstream
    a 
}
