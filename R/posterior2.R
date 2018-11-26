###############################################################################
# posterior2.R
## 2018-11-19 draft of classMembership as standalone function
##            cf posterior.allocation previously called in openCR.fit
###############################################################################

matchch <- function (x, sq) {
    df <- as.data.frame(matrix(x, nrow = nrow(x)))
    if (!is.null(covariates(x))) {
        if (nrow(covariates(x))>0)
            df <- cbind(df, covariates(x))
    }
    
    dfsq <- as.data.frame(matrix(sq, nrow = nrow(sq)))
    if (!is.null(covariates(sq))) {
        if (nrow(covariates(sq))>0)
            dfsq <- cbind(dfsq, covariates(sq))
    }
    df$freq <- NULL
    dfsq$freq <- NULL
    
    if (!all(colnames(df) == colnames(dfsq)))
        stop ("sq, x do not match")
    cx <- do.call("paste", c(df[, , drop = FALSE], sep = "\r"))
    csq <- do.call("paste", c(dfsq[, , drop = FALSE], sep = "\r"))
    match(cx, csq, nomatch = NA)
}

classMembership <- function (object, ...) UseMethod("classMembership")

classMembership.openCR <- function (object, fullCH = NULL, ...) {
    
    # Return the probability of membership in latent classes for h2 and h3 models
    
    nmix <- object$details$nmix
    if(object$movementmodel != 'static') {
        stop ("posterior allocation not ready for movement models")
    }
    capthist <- object$capthist   # squeezed!
    if (nmix == 1) return (rep(1,nrow(capthist)))
    secr <- grepl('secr', object$type)
    
    #--------------------------------------------------------------------
    # Real parameters
    realparval  <- makerealparameters (object$design, complete.beta(object), 
                                       object$parindx, object$link, object$fixed)
    # Parameter Index Array
    PIA <- object$design$PIA
    PIAJ <- object$design$PIAJ
    #--------------------------------------------------------------------
    
    cumss <- getcumss(capthist)
    J <- length(cumss)-1 
    nc <- nrow(capthist)   ## beware freq!
    intervals <- intervals(capthist)
    primarysession <- primarysessions(intervals)
    
    
    ##--------------------------------------------------------------
    ## get fi,li, re-form capthist as CH
    lost <- which(apply(capthist,1,min, drop = FALSE)<0)
    twoD <- apply(abs(capthist), 1:2, sum, drop = FALSE)
    CH <- twoD
    if (J==1)
        twoD <- as.matrix(apply(twoD, 1, function(x) tapply(x,primarysession,max)))
    else
        twoD <- t(apply(twoD, 1, function(x) tapply(x,primarysession,max)))  # in terms of primary sessions
    fi <- apply(twoD, 1, function(x) min(which(x>0)))
    li <- apply(twoD, 1, function(x) max(which(x>0)))
    twoD[cbind(lost, li[lost])] <- -1
    li[lost] <- -li[lost]
    covariates(CH) <- covariates(capthist)
    covariates(twoD) <- covariates(capthist)
    JScounts <- unlist(JS.counts(twoD))
    ##--------------------------------------------------------------
    
    if (secr) {
        type <- switch(object$type, CJSsecr = 6, JSSAsecrf = 7, JSSAsecrD = 8,
                       JSSAsecrfCL = 9, JSSAsecrlCL = 10, JSSAsecrbCL = 11,
                       JSSAsecrl = 12, JSSAsecrb = 13, JSSAsecrB = 14,
                       JSSAsecrg = 24, JSSAsecrgCL = 25, secrCL = 30, secrD = 31, -1)
        if (type < 0) stop ("Invalid likelihood type for posterior allocation")
        
        trps <- traps(object$capthist)
        detectr <- detector(trps)[1]
        k <- nrow(trps)
        m <- nrow(object$mask)
        
        if (!is.null(object$mask)) area <- attr(object$mask,'area')
        else area <- 0
        binomN <- switch (detectr, multi = 1, proximity = 1, count = object$binomN, -1)
        
        usge <- usage(traps(capthist))
        if (is.null(usge) | object$details$ignoreusage) 
            usge <- matrix(1, nrow=k, ncol= cumss[J+1])  # in terms of secondary sessions
        ## 2017-11-26 collapse data from exclusive detectors; modified 2018-01-17
        CH <- capthist
        if (detectr == 'multi') {
            CH <- abs(capthist)
            CH <- apply(CH,1:2, which.max) *  (apply(CH,1:2, max)>0)
            lost <- apply(capthist,1:2, min)<0
            CH[lost] <- -CH[lost]
            class (CH) <- 'capthist'
            traps(CH) <- traps(capthist)
        }
        ##--------------------------------------------------------------
        
        if (detectr == 'multi') {
            onehistory <- function (n, pmix) {
                px <- numeric(nmix)
                if ((type == 6) & fi[n] == J) rep(NA, nmix)
                else {
                    for (x in 1:nrow(pmix)) {
                        hx <- matrix(haztemp$h[x,,], nrow = m)
                        px[x] <- prwisecrcpp(
                            as.integer(type),
                            as.integer(0),                       ## always just one n
                            as.integer(x-1),
                            as.integer(1),                       ## nc one at a time
                            as.integer(J),
                            as.integer(k),
                            as.integer(m),
                            as.integer(nmix),
                            as.integer(cumss),
                            as.integer(CH[n,,drop = FALSE]),       ## 2-D CH
                            as.integer(fi[n]),
                            as.integer(li[n]),
                            as.double (hk),                      ## hazard instead of probability
                            as.matrix (realparval),
                            as.integer(PIA[n,,,]),
                            as.integer(PIAJ[n,,]),
                            as.integer(binomN),
                            as.matrix (usge),
                            as.double (intervals),
                            as.integer(c(-2,-2)),
                            as.matrix (hx),                      ## lookup sum_k (hazard)
                            as.matrix (haztemp$hindex[n,,drop=FALSE]),      ## index to h
                            as.integer(object$details$CJSp1),
                            as.integer(0),                        ## object$movemodel
                            as.character(""),                     ## object$usermodel
                            as.matrix(-1),                        ## object$kernel
                            as.matrix(-1),                        ## mqarray
                            as.double (-1))                       ## data$cellsize
                    }
                    px/sum(px)
                }
            }
        }
        else{
            onehistory <- function (n, pmix) {
                nmix <- nrow(pmix)
                if ((type == 6) & fi[n] == J) rep(NA, nmix)
                else {
                    px <- numeric(nmix)
                    for (x in 1:nmix) {
                        px[x]<- prwisecrcpp(
                            as.integer(type),
                            as.integer(0),    # n-1
                            as.integer(x-1),
                            as.integer(1),
                            as.integer(J),
                            as.integer(k),
                            as.integer(m),
                            as.integer(nmix),
                            as.integer(cumss),
                            as.integer(CH[n,,, drop = FALSE]),   ## 3-D CH
                            as.integer(fi[n]),
                            as.integer(li[n]),
                            as.double (gk),                ## precomputed probability
                            as.matrix (realparval),
                            as.integer(PIA[n,,,, drop = FALSE]),
                            as.integer(PIAJ[n,,]),
                            as.integer(binomN),
                            as.matrix (usge),
                            as.double (intervals),
                            as.integer(c(-2,-2)),           ## moveargsi  
                            as.matrix(-1),                  ## not multicatch
                            as.matrix(-1),                  ## not multicatch
                            as.integer(object$details$CJSp1),
                            as.integer(0),
                            as.character(""),
                            as.matrix(-1),
                            as.matrix(-1),
                            as.double (-1)
                        )
                    }
                    px/sum(px)
                }
            }
        }
        
        temp <- makegkParallelcpp (as.integer(object$detectfn), 
                                   as.integer(.openCRstuff$sigmai[type]),
                                   as.integer(object$details$grain),
                                   as.matrix(realparval),
                                   as.matrix(trps),
                                   as.matrix(object$mask))
        gk <- temp[[1]]
        hk <- temp[[2]]
        pmix <- fillpmix2(nc, nmix, PIA, realparval)
        
        if (detectr=='multi') {
            haztemp <- gethcpp(
                as.integer(nc),
                as.integer(nrow(realparval)),
                as.integer(nmix),
                as.integer(k),
                as.integer(cumss[J+1]),
                as.integer(m),
                as.integer(PIA),
                as.matrix(usge),
                as.double(hk))
            haztemp$h <- array(haztemp$h, dim = c(nmix, m, max(haztemp$hindex)+1))
        }
        else {
            haztemp <- list(h = array(-1, dim=c(nmix,1,1)), hindex = matrix(-1))
        }
    }
    else {
        ## NON-SPATIAL
        type <- switch(object$type, CJS = 1, JSSAb = 2, JSSAl = 3, JSSAf = 4,
                       JSSAfCL = 15, JSSAlCL = 16, JSSAbCL = 17, JSSAB = 18, JSSAN = 19,
                       Pradel = 20, JSSARET = 21, JSSAg = 22, JSSAgCL = 23, Pradelg = 26, JSSAfgCL = 27,
                       -1)
        if (type<0) stop ("Invalid likelihood type")
        if (type %in% c(20,26)) {
            stop("mixtures not expected in Pradel models")
        }
        
        onehistory <- function (n, pmix) {
            nmix <- nrow(pmix)
            if ((type == 1) & fi[n] == J) rep(NA, nmix)
            else {
                px <- numeric(nmix)
                for (x in 1:nmix) {
                    px[x] <- prwicpp (
                        as.integer(type),
                        as.integer(0),
                        as.integer(x-1),
                        as.integer(1),
                        as.integer(J),
                        as.integer(cumss),
                        as.integer(object$details$nmix),
                        as.integer(CH[n,]),
                        as.integer(fi[n]),
                        as.integer(li[n]),              # may be negative if censored 2018-01-17
                        as.matrix (realparval),
                        as.integer(PIA[n,,,]),
                        as.integer(PIAJ[n,,]),
                        as.double (intervals),
                        as.integer(object$details$CJSp1)
                    )
                }
                px/sum(px)
            }
        }
        if (!is.null(fullCH)) {
            fullCH <- reduce(fullCH, outputdetector = 'nonspatial', verify = FALSE, dropunused=FALSE)
            # FAILS TO COMPRESS LAST DIM (bug in reduce.capthist 2018-11-20), SO REPEAT 
            fullCH <- reduce(fullCH, verify = FALSE, dropunused=FALSE)
        }
    }
    pmix <- fillpmix2(nc, nmix, PIA, realparval)
    out <- data.frame(t(sapply(1:nc, onehistory, pmix = pmix)))
    names(out) <- paste0('class', 1:nmix)
    out$maxclass <- max.col(out)
    
    if (is.null(fullCH)) {   ## untested
        fullCH <- unsqueeze(object$capthist)
    }
    else {
        if (ms(fullCH)) fullCH <- join(fullCH, drop.sites = !secr)
        nf <- sum(covariates(capthist)$freq)
        if (nrow(fullCH) != nf)
            stop ("fullCH has differing number of histories (", nc, " vs ", nf, ")")
    }
    i <- matchch(fullCH, capthist)
    out <- out[i,]
    rownames(out) <- rownames(fullCH)  
    out
}

######################################################################################################
