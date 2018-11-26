###############################################################################
# posterior.R
## 2018-04-20 openCR 1.2.0
## 2018-11-11 debug 'not a matrix'
## 2018-11-20 replaced by classMembership method (see posterior2.R)
###############################################################################

posterior.allocation <- function (beta, data)
    
    # Return the probability of membership in latent classes for h2 and h3 models
    
{
    nmix <- data$details$nmix
    if(data$movemodel>1) {
        warning ("posterior allocation not ready for movement models")
        nmix <- 1
    }
    if (nmix == 1) return (rep(1,nrow(data$capthist)))
    secr <- grepl('secr', data$type)
    
    #--------------------------------------------------------------------
    # Fixed beta
    fb <- data$details$fixedbeta
    if (!is.null(fb)) {
        fb[is.na(fb)] <- beta
        beta <- fb    ## complete
    }
    
    #--------------------------------------------------------------------
    # Real parameters
    realparval  <- makerealparameters (data$design, beta, data$parindx, data$link, data$fixed)
    # Parameter Index Array
    PIA <- data$design$PIA
    PIAJ <- data$design$PIAJ
    #--------------------------------------------------------------------
    
    if (secr) {
        type <- switch(data$type, CJSsecr = 6, JSSAsecrf = 7, JSSAsecrD = 8,
                       JSSAsecrfCL = 9, JSSAsecrlCL = 10, JSSAsecrbCL = 11,
                       JSSAsecrl = 12, JSSAsecrb = 13, JSSAsecrB = 14,
                       JSSAsecrg = 24, JSSAsecrgCL = 25, secrCL = 30, secrD = 31, -1)
        if (type < 0) stop ("Invalid likelihood type for posterior allocation")
        trps <- traps(data$capthist)
        if (!is.null(data$mask)) area <- attr(data$mask,'area')
        else area <- 0
        binomN <- switch (detector(trps)[1], multi = 1, proximity = 1, count = data$binomN, -1)
        
        if (data$multi) {
            onehistory <- function (n, pmix) {
                px <- numeric(nmix)
                if ((type == 6) & data$fi[n] == data$J) rep(NA, nmix)
                else {
                    for (x in 1:nrow(pmix)) {
                        hx <- matrix(haztemp$h[x,,], nrow = data$m)
                        px[x] <- prwisecrcpp(
                            as.integer(type),
                            as.integer(0),                       ## always just one n
                            as.integer(x-1),
                            as.integer(1),                       ## nc one at a time
                            as.integer(data$J),
                            as.integer(data$k),
                            as.integer(data$m),
                            as.integer(data$details$nmix),
                            as.integer(data$cumss),
                            as.integer(data$capthist[n,,drop = FALSE]),       ## 2-D CH
                            as.integer(data$fi[n]),
                            as.integer(data$li[n]),
                            as.double (hk),                      ## hazard instead of probability
                            as.matrix (realparval),
                            as.integer(PIA[n,,,]),
                            as.integer(PIAJ[n,,]),
                            as.integer(binomN),
                            as.matrix (data$usge),
                            as.double (data$intervals),
                            as.integer(data$moveargsi),
                            as.matrix (hx),                      ## lookup sum_k (hazard)
                            as.matrix (haztemp$hindex[n,,drop=FALSE]),      ## index to h
                            as.integer(data$details$CJSp1),
                            as.integer(data$movemodel),
                            as.character(data$usermodel),
                            as.matrix(data$kernel),
                            as.matrix(data$mqarray),
                            as.double (data$cellsize))
                    }
                    px/sum(px)
                }
            }
        }
        else{
            onehistory <- function (n, pmix) {
                nmix <- nrow(pmix)
                if ((type == 6) & data$fi[n] == data$J) rep(NA, nmix)
                else {
                    px <- numeric(nmix)
                    for (x in 1:nmix) {
                        px[x]<- prwisecrcpp(
                            as.integer(type),
                            as.integer(0),    # n-1
                            as.integer(x-1),
                            as.integer(1),
                            as.integer(data$J),
                            as.integer(data$k),
                            as.integer(data$m),
                            as.integer(data$details$nmix),
                            as.integer(data$cumss),
                            as.integer(data$capthist[n,,, drop = FALSE]),   ## 3-D CH
                            as.integer(data$fi[n]),
                            as.integer(data$li[n]),
                            as.double (gk),                ## precomputed probability
                            as.matrix (realparval),
                            as.integer(PIA[n,,,, drop = FALSE]),
                            as.integer(data$design$PIAJ[n,,]),
                            as.integer(binomN),
                            as.matrix (data$usge),
                            as.double (data$intervals),
                            as.integer(data$moveargsi),
                            as.matrix(-1),                  ## not multicatch
                            as.matrix(-1),                  ## not multicatch
                            as.integer(data$details$CJSp1),
                            as.integer(data$movemodel),
                            as.character(data$usermodel),
                            as.matrix(data$kernel),
                            as.matrix(data$mqarray),
                            as.double (data$cellsize)
                        )
                    }
                    px/sum(px)
                }
            }
        }
        
        
        temp <- makegkParallelcpp (as.integer(data$detectfn), 
                                   as.integer(.openCRstuff$sigmai[type]),
                                   as.integer(data$details$grain),
                                   as.matrix(realparval),
                                   as.matrix(trps),
                                   as.matrix(data$mask))
        gk <- temp[[1]]
        hk <- temp[[2]]
        pmix <- fillpmix2(data$nc, nmix, PIA, realparval)
        
        if (data$multi) {
            haztemp <- gethcpp(
                as.integer(data$nc),
                as.integer(nrow(realparval)),
                as.integer(data$details$nmix),
                as.integer(data$k),
                as.integer(data$cumss[data$J+1]),
                as.integer(data$m),
                as.integer(PIA),
                as.matrix(data$usge),
                as.double(hk))
            haztemp$h <- array(haztemp$h, dim = c(data$details$nmix, data$m, max(haztemp$hindex)+1))
        }
        else {
            haztemp <- list(h = array(-1, dim=c(data$details$nmix,1,1)), hindex = matrix(-1))
        }
        
    }
    else {
        ## NON-SPATIAL
        type <- switch(data$type, CJS = 1, JSSAb = 2, JSSAl = 3, JSSAf = 4,
                       JSSAfCL = 15, JSSAlCL = 16, JSSAbCL = 17, JSSAB = 18, JSSAN = 19,
                       Pradel = 20, JSSARET = 21, JSSAg = 22, JSSAgCL = 23, Pradelg = 26, JSSAfgCL = 27,
                       -1)
        if (type<0) stop ("Invalid likelihood type")
        if (type %in% c(20,26)) {
            stop("mixtures not expected in Pradel models")
        }
        
        onehistory <- function (n, pmix) {
            nmix <- nrow(pmix)
            if ((type == 1) & data$fi[n] == data$J) rep(NA, nmix)
            else {
                px <- numeric(nmix)
                for (x in 1:nmix) {
                    px[x] <- prwicpp (
                        as.integer(type),
                        as.integer(0),
                        as.integer(x-1),
                        as.integer(1),
                        as.integer(data$J),
                        as.integer(data$cumss),
                        as.integer(data$details$nmix),
                        as.integer(data$capthist[n,]),
                        as.integer(data$fi[n]),
                        as.integer(data$li[n]),              # may be negative if censored 2018-01-17
                        as.matrix (realparval),
                        as.integer(PIA[n,,,]),
                        as.integer(data$design$PIAJ[n,,]),
                        as.double (data$intervals),
                        as.integer(data$details$CJSp1)
                    )
                }
                px/sum(px)
            }
        }
    }
    pmix <- fillpmix2(data$nc, nmix, PIA, realparval)
    out <- data.frame(t(sapply(1:data$nc, onehistory, pmix = pmix)))
    names(out) <- paste0('class', 1:nmix)
    out$maxclass <- max.col(out)
    out
}

######################################################################################################
