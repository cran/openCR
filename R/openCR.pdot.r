################################################################################
## package 'openCR'
## openCR.pdot.R
## 2011-12-30, 2013-01-20
################################################################################

openCR.pdot <- function (object, bysession = FALSE) {

    # Return pdot for open population model
    beta     <- object$fit$par
    parindx  <- object$parindx
    link     <- object$link
    fixed    <- object$fixed
    design0  <- object$design    # object$design0  stopgap
    capthist <- object$capthist
    details  <- object$details
    type     <- object$type
    intervals <- object$intervals
    cumss    <- getcumss(object$capthist)
    
    #--------------------------------------------------------------------
    # Fixed beta
    fb <- details$fixedbeta
    if (!is.null(fb)) {
        fb[is.na(fb)] <- beta
        beta <- fb    ## complete
    }
    #--------------------------------------------------------------------
    # Real parameters
    realparval0 <- makerealparameters (design0, beta, parindx, link, fixed)
    nc  <- nrow(capthist)
    J <- length(intervals) + 1
    type <- switch(type, CJS = 1, JSSAb = 2, JSSAl = 3, JSSAf = 4, JSSAg = 22, JSSAgCL = 23,
        JSSAfCL = 15, JSSAlCL = 16, JSSAbCL = 17, JSSAB = 18, JSSAN = 19, JSSARET = 21,
        Pradel = 20, Pradelg = 26, JSSAk = 28, JSSAkCL = 29)
    distrib <- switch (object$distribution, poisson = 0, binomial = 1)
    binomN <- details$binomN
    
    ## mixture proportions
    if (details$nmix > 1) {
        pmix <- fillpmix2(nc, details$nmix, design0$PIA, realparval0)
    }
    else {
        pmix <- matrix(1, nrow = details$nmix, ncol = nc)
    }
    onep <- function (x) {
        if (details$R) {
            pch1 <- PCH1(
                type,
                x,
                nc,
                cumss,
                details$nmix,
                realparval0,
                design0$PIA,
                design0$PIAJ,
                intervals)
        }
        else {
            pch1 <-  PCH1cpp(
            as.integer(type),
            as.integer(x-1),
            as.integer(nc),
            as.integer(J),
            as.integer(cumss),
            as.integer(details$nmix),
            as.matrix(realparval0),
            as.integer(design0$PIA),
            as.integer(design0$PIAJ),
            as.double(intervals))
        }
        pmix[x,] * pch1
    }
    onepbysession <- function (x) {
        one <- function(n) {    
            p <- getp (n, x, realparval0, design0$PIA)
            sessp <- function (j) {
                s <- (cumss[j]+1) : cumss[j+1]
                1-prod(1 - p[s])       ## Pr detected
            }
            sapply(1:J, sessp)
        }
        p1 <- sapply(1:nc, one)   # or rep(one(1), nc) if all the same
        sweep(p1, MARGIN=2, STATS=pmix[x,], FUN="*")
    }
    if (bysession) {
        p <- sapply(1:details$nmix, onepbysession, simplify = FALSE)
        apply(abind(p, along=3), 1:2, sum)   ## session * ch, summed over mixture classes
    }
    else {
        p <- sapply(1:details$nmix, onep)
        apply(p,1,sum)
    }
}
############################################################################################


