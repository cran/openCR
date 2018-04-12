# 2017-05-17 not ready for mixtures
# 2017-05-17 pass just one animal to onehistory
# 2017-05-17 PCH0 called with single representative history
# 2017-05-17 capthist 2D for non-spatial
# 2017-05-22 use data argument (an environment)
# 2018-02-06 drop 'par.' prefix from function names
# 2018-03-26 switch pch0 to pch1

# types

# CJS = 1
# JSSAb = 2
# JSSAl = 3
# JSSAf = 4
# JSSAfCL = 15
# JSSAlCL = 16
# JSSAbCL = 17
# JSSAB = 18
# JSSAN = 19
# Pradel = 20
# JSSARET = 21
# JSSAg = 22
# JSSAgCL = 23
# Pradelg = 26
# JSSAfgCL = 27
# JSSAk = 28
# JSSAkCL = 29

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

#---------------------------------------------------------

open.loglikfn <- function (beta, dig = 3, betaw = 8, cluster = NULL, oneeval = FALSE, data)

    # Return the negative log likelihood
    # Transformed parameter values are passed in the vector 'beta'
    # details$trace=T sends a one-line report to the screen

{
    #--------------------------------------------------------------------
    # Fixed beta
    fb <- data$details$fixedbeta
    if (!is.null(fb)) {
        fb[is.na(fb)] <- beta
        beta <- fb    ## complete
    }
    if (data$details$debug>0) {
        print(beta)
    }
    freq <- covariates(data$capthist)$freq

    if (is.null(freq)) freq <- rep(1, data$nc)
    if (length(freq) == 1) freq <- rep(freq, data$nc)
    ncf <- sum(freq)

    #--------------------------------------------------------------------
    # Real parameters
    realparval  <- makerealparameters (data$design, beta, data$parindx, data$link, data$fixed)
    if (data$learnedresponse)
    realparval0 <- makerealparameters (data$design0, beta, data$parindx, data$link, data$fixed)
    else realparval0 <- realparval

    #-----------------------------------------
    # check valid parameter values
    if (!all(is.finite(realparval))) {
        cat ('beta vector :', beta, '\n')
        cat ('real vector :', realparval, '\n')
        warning ("extreme 'beta' in 'openloglikfn' ",
                 "(try smaller stepmax in nlm Newton-Raphson?)")
        return (1e10)
    }
    type <- typecode(data$type)
    if (type<0) stop ("Invalid likelihood type")
    if (type %in% 28:29 & any(data$intervals !=1)) stop ("kappa parameterisation available only if all intervals = 1")
    
    PIA <- data$design$PIA
    PIAJ <- data$design$PIAJ

    if (data$details$debug>1) {
        message('Type ', type)
        message('J    ', data$J)
        message('nmix ', data$details$nmix)
        print(realparval)
        print (table(PIA))
        print(data$intervals)
        flush.console()
    }

    if (data$details$debug>2) browser()

    if (type %in% c(20,26)) {
        # Pradel model
        if (data$details$R) {
            comp <- pradelloglik(type, data$JScounts, realparval,  PIAJ, data$intervals)
        }
        else {
            comp <- pradelloglikcpp(
                as.integer(type),
                as.integer(data$JScounts),
                as.integer(data$nc),         ## needed for nrows of PIAJ
                as.integer(data$J),
                as.integer(data$details$nmix),
                as.double(realparval),
                as.integer(nrow(realparval)),       # number of rows in lookup table
                as.integer(PIAJ),             # index of nc,S,mix to rows
                as.double(data$intervals))                # number of interval == J-1
        }
    }
    else if (data$details$R & (type %in% c(28,29))) {
        comp <- kappaloglik (type, realparval,  PIA, PIAJ, data) 
    }
    else {
        onehistory <- function (n, pmix) {
                sump <- 0
                if (data$details$R) {
                    for (x in 1:nrow(pmix)) {
                        temp <- prwi(
                            type,
                            1,   # n
                            x,
                            data$J,
                            data$cumss,
                            data$details$nmix,
                            data$capthist[n,, drop = FALSE],
                            data$fi[n],
                            data$li[n],
                            realparval,
                            PIA[n,,,,drop = FALSE],
                            PIAJ[n,,,drop = FALSE],
                            data$intervals,
                            data$details$CJSp1)
                        sump <- sump + pmix[x,n] * temp
                    }
                }
                else {
                    for (x in 1:nrow(pmix)) {
                        temp <- prwicpp (
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
                            as.double (realparval),
                            as.integer(nrow(realparval)),
                            as.integer(PIA[n,,,]),
                            as.integer(PIAJ[n,,]),
                            as.double (data$intervals),
                            as.integer(data$details$CJSp1))
                        sump <- sump + pmix[x,n] * temp
                    }
                }
                if (any(sump<=0)) {
                    -1e10
                }
                else
                    freq[n] * log(sump)
        }

        comp <- numeric(5)
        pmix <- fillpmix2(data$nc, data$details$nmix, PIA, realparval)

        #####################################################################
        # Component 1: Probability of observed histories - all models
        if (is.null(cluster))
            temp <- sapply(1:data$nc, onehistory, pmix = pmix)
        else
            temp <- parSapply(cluster, 1:data$nc, onehistory, pmix = pmix, USE.NAMES = FALSE )

        comp[1] <- sum(temp)
        #####################################################################
        # Component 2: Probability of missed animals (all-zero histories)

        if (type %in% c(2:4,15:19, 21, 22, 23, 27, 28, 29)) {
            pdot <- rep(0, data$nc)
            if (data$learnedresponse)
                PIA0 <- data$design0$PIA
            else
                PIA0 <- PIA
            for (x in 1:data$details$nmix) {   # loop over latent classes
                if (data$details$R) {
                    pch1 <- PCH1(
                        type,
                        x,
                        data$nc,
                        data$cumss,
                        data$details$nmix,
                        realparval0,
                        PIA0,
                        PIAJ,
                        data$intervals)
                }
                else {
                pch1 <-  PCH1cpp(
                            as.integer(type),
                            as.integer(x-1),
                            as.integer(data$nc),
                            as.integer(data$J),
                            as.integer(data$cumss),
                            as.integer(data$details$nmix),
                            as.double(unlist(realparval0)),
                                 as.integer(nrow(realparval0)),
                                 as.integer(PIA0),
                                 as.integer(PIAJ),
                            as.double(data$intervals))
                }
                pdot <- pdot + pmix[x] * pch1
            }
            comp[2] <- - sum(freq * log(pdot))
        }

        #####################################################################
        # Component 3: Probability of observing nc animals
        if (type %in% c(2:4,18,19,21, 22, 28)) {
            if (type %in% c(2,3,4,22,28)) {
                superN <- realparval[nrow(realparval)*3+1] # Nsuper direct
            }
            else {
                # if (data$details$R) {
                    superN <- getN(type, ncf, data$J, data$details$nmix, pmix, realparval, PIAJ, data$intervals)
                # }
                # else {
                #     temp <- getNcpp(
                #         as.integer(type),
                #         as.integer(data$nc),
                #         as.integer(ncf),
                #         as.integer(data$J),
                #         as.integer(data$details$nmix),
                #         as.double(pmix),
                #         as.double(data$intervals),
                #         as.double(realparval),
                #         as.integer(nrow(realparval)),
                #         as.integer(PIAJ))
                #     
                #     superN <- temp[1]
                # }
            }
            meanpdot <- ncf / sum(1/rep(pdot,freq))  ## cf CLmeanesa in 'secr'
            comp[3] <- switch (data$distrib+1,
                               dpois(ncf, superN * meanpdot, log = TRUE),
                               ## lnbinomial (ncf, superN, meanpdot),
                               lnbinomial (ncf, superN + ncf, meanpdot),
                               NA)
        }
        
    }

    ## optional multinomial term
    if (data$details$multinom & (type %in% c(2,3,4,15,16,17,18,19,21,22,23,27))) {
        nh <- table(rep(apply(data$capthist, 1, paste, collapse=''), freq))
        comp[5] <- lgamma(ncf+1) - sum(lgamma(nh+1))
    }

    ## log-likelihood as sum of components 1-3 and 5
    loglik <- sum(comp)

    ## debug
    if (data$details$debug>=1) {
        cat("Likelihood components (comp) ", format(comp, digits=10), "\n")
        cat("Total ", format(loglik, digits = 10), "\n")
        browser()
    }

    .openCRstuff$iter <- .openCRstuff$iter + 1   ## moved outside loop 2011-09-28
    if (data$details$trace) {
        if (!is.null(data$details$fixedbeta))
            beta <- beta[is.na(data$details$fixedbeta)]
        cat(format(.openCRstuff$iter, width=4),
            formatC(round(loglik,dig), format='f', digits=dig, width=10),
            formatC(beta, format='f', digits=dig+1, width=betaw),
            '\n', sep = " ")

        flush.console()
    }

    if (oneeval) {
        c(loglik, beta)
    }
    else {
        if (is.finite(loglik)) -loglik   # return the negative loglikelihood
        else 1e10
    }

}

open.secr.loglikfn <- function (beta, dig = 3, betaw = 8, cluster = NULL, oneeval = FALSE, data)

    # Return the negative log likelihood
    # Transformed parameter values are passed in the vector 'beta'
    # details$trace=T sends a one-line report to the screen

    # an existing cluster is used

{
    #--------------------------------------------------------------------
    # Fixed beta
    fb <- data$details$fixedbeta
    if (!is.null(fb)) {
        fb[is.na(fb)] <- beta
        beta <- fb    ## complete
    }
    freq <- covariates(data$capthist)$freq
    if (is.null(freq)) freq <- rep(1, data$nc)
    if (length(freq) == 1) freq <- rep(freq, data$nc)  # shoudn't be needed
    ncf <- sum(freq)

    #--------------------------------------------------------------------
    # Real parameters
    realparval  <- makerealparameters (data$design, beta, data$parindx, data$link, data$fixed)
    if (data$learnedresponse)
        realparval0 <- makerealparameters (data$design0, beta, data$parindx, data$link, data$fixed)
    else realparval0 <- realparval

    #-----------------------------------------
    # check valid parameter values
    if (!all(is.finite(realparval))) {
        return (1e10)
    }

    #-----------------------------------------
    type <- typecode(data$type)
    if (type < 0) stop ("Invalid likelihood type")
    
    trps <- traps(data$capthist)
    if (!is.null(data$mask)) area <- attr(data$mask,'area')
    else area <- 0
    
    binomN <- switch (detector(trps)[1], multi = 1, proximity = 1, count = data$binomN, -1)
    if (binomN < 0)
        stop("open-population secr requires multi, proximity or count detector type")
    
    if (data$details$debug>1) {
        print(type)
        print(summary(updateCH(data$capthist)))
        print(data$nc)
        print(data$J)
        print(data$k)
        print(data$m)
        print(data$cumss)
        print(data$details)
        print(summary(trps))
        print(summary(data$mask))
        print(area)
        print(data$design)
        print(realparval)
        print(nrow(realparval))
        print(data$intervals)
        print(data$movemodel)
        print(binomN)
    }
    
    if (data$details$debug>2) browser()
    
    PIA <- data$design$PIA
    PIAJ <- data$design$PIAJ
    if (data$learnedresponse)
        PIA0 <- data$design0$PIA
    else
        PIA0 <- PIA
    
    onehistory <- function (n, pmix, gk) {
        sump <- 0
        if (data$details$R) {
            for (x in 1:nrow(pmix)) {
                temp <- prwisecr(
                    type,
                    1,   # n
                    x,
                    1,
                    data$J,
                    data$k,
                    data$m,
                    data$details$nmix,
                    data$cumss,
                    data$capthist[n,,, drop = FALSE],   ## 3-D CH
                    data$fi[n],
                    data$li[n],
                    gk,                ## precomputed probability
                    realparval,
                    PIA[n,,,,drop = FALSE],
                    data$design$PIAJ[n,,,drop = FALSE],
                    binomN,
                    data$usge,
                    data$intervals,
                    data$details$CJSp1,
                    data$moveargsi,
                    data$movemodel,
                    get(data$usermodel),
                    data$kernel,
                    data$mqarray,
                    data$cellsize)
                sump <- sump + pmix[x,n] * temp
            }
        }
        else {
            for (x in 1:nrow(pmix)) {
                temp <- prwisecrcpp(
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
                    as.double (realparval),
                    as.integer(nrow(realparval)),
                    as.integer(PIA[n,,,, drop = FALSE]),
                    as.integer(data$design$PIAJ[n,,]),
                    as.integer(binomN),
                    as.double (data$usge),
                    as.double (data$intervals),
                    as.integer(data$moveargsi),
                    as.integer(data$details$CJSp1),
                    as.integer(data$movemodel),
                    as.character(data$usermodel),
                    as.integer(nrow(data$kernel)),
                    as.integer(as.matrix(data$kernel)),
                    as.integer(data$mqarray),
                    as.double (data$cellsize))
                sump <- sump + pmix[x,n] * temp
            }
        }
        freq[n] * log(sump)
    }
    onehistorymulti <- function (n, pmix, hk) {
        sump <- 0
        for (x in 1:nrow(pmix)) {
            if (data$details$R) {
                temp <- prwisecrmulti (
                    type,
                    1, x, 
                    data$J,
                    data$m, 
                    data$cumss, 
                    data$capthist[n,,drop = FALSE], 
                    data$fi[n], 
                    data$li[n], 
                    hk, 
                    realparval,
                    PIA[n,,,,drop = FALSE],
                    PIAJ[n,,,drop = FALSE],
                    binomN, 
                    data$usge, 
                    data$intervals,  
                    data$moveargsi, 
                    haztemp$h, 
                    haztemp$hindex[n,,drop = FALSE]+1, 
                    data$details$CJSp1, 
                    data$movemodel,
                    data$usermodel, 
                    data$kernel, 
                    data$mqarray, 
                    data$cellsize)
            }
            else {
                temp <- prwisecrmulticpp(
                    as.integer(type),
                    as.integer(0),                       ## always just one n
                    as.integer(x-1),
                    as.integer(1),                       ## nc one at a time
                    as.integer(data$J),
                    as.integer(data$cumss),
                    as.integer(data$k),
                    as.integer(data$m),
                    as.integer(data$details$nmix),
                    as.integer(data$capthist[n,]),       ## 2-D CH
                    as.integer(data$fi[n]),
                    as.integer(data$li[n]),
                    as.double (hk),                      ## hazard instead of probability
                    as.double (realparval),
                    as.integer(nrow(realparval)),
                    as.integer(PIA[n,,,]),
                    as.integer(PIAJ[n,,]),
                    as.integer(binomN),
                    as.double (data$usge),
                    as.double (data$intervals),
                    as.integer(data$moveargsi),
                    as.double (haztemp$h),               ## lookup sum_k (hazard)
                    as.integer(haztemp$hindex[n,]),     ## index to h
                    as.integer(data$details$CJSp1),
                    as.integer(data$movemodel),
                    as.character(data$usermodel),
                    as.integer(nrow(data$kernel)),
                    as.integer(as.matrix(data$kernel)),
                    as.integer(data$mqarray),
                    as.double (data$cellsize))
            }
            # cat('n ', n, ' LL ', temp, '\n')
            sump <- sump + pmix[x,n] * temp
        }
        freq[n] * log(sump)
        }
    temp <- makegkcpp(
        as.integer(nrow(realparval)),      # number of rows in lookup table
        as.integer(data$k),                # detectors
        as.integer(data$m),                # mask points
        as.integer(data$detectfn),
        as.integer(.openCRstuff$sigmai[type]),                # column index of sigma in realparval
        as.double(realparval),
        as.double(unlist(trps)),
        as.double(unlist(data$mask)))
    gk <- array(temp[[1]], dim=c(nrow(realparval), data$k, data$m))  # array form for R use
    hk <- array(temp[[2]], dim=c(nrow(realparval), data$k, data$m))  # array form for R use

    # ch <- split.default(capthist, slice.index(capthist,1))
    # pia <- split.default(PIA, slice.index(PIA, 1))
    
    comp <- numeric(5)
    
    pmix <- fillpmix2(data$nc, data$details$nmix, PIA, realparval)
    S <- ncol(data$capthist)
    
    haztemp <- gethcpp(
        as.integer(data$nc),
        as.integer(nrow(realparval)),
        as.integer(data$details$nmix),
        as.integer(data$k),
        as.integer(data$J),
        as.integer(data$m),
        as.integer(PIA),
        as.integer(data$cumss),
        as.double(data$usge),
        as.double(hk))
    names(haztemp) <- c('hc0','h','hindex')
    haztemp$hindex <- matrix(haztemp$hindex, nrow = data$nc)
    ## discard surplus h for speed in C   DISCARD INTERNALLY
    haztemp$h <- haztemp$h[1:(data$m * data$details$nmix * (max(haztemp$hindex)+1))]
    if (data$details$R) {
        haztemp$h <- array(haztemp$h, dim = c(data$details$nmix, data$m, max(haztemp$hindex)+1))
    }
    #####################################################################
    # Component 1: Probability of observed histories
    
    if (is.null(cluster)) {
        if (data$multi) {
            ## temp <- allhistoriesmulti(pmix, hk)  ## little speed gain, or slower
            temp <- sapply(1:data$nc, onehistorymulti, pmix = pmix, hk = hk, USE.NAMES = FALSE )
        }
        else {
            ## temp <- allhistories(pmix, gk)  ## little speed gain, or slower
            temp <- sapply(1:data$nc, onehistory, pmix = pmix, gk = gk, USE.NAMES = FALSE)
        }
    }
    else {
        if (data$multi) {
            clusterExport(cluster, c("pmix", "hk", "haztemp"), environment())
            temp <- parSapply(cluster, 1:data$nc, onehistorymulti, pmix = pmix, hk = hk,
                              USE.NAMES = FALSE )
        }
        else {
            clusterExport(cluster, c("pmix", "gk"), environment())
            temp <- parSapply(cluster, 1:data$nc, onehistory, pmix = pmix, gk = gk,
                              USE.NAMES = FALSE )
        }
    }
    comp[1] <- sum(temp)
    #####################################################################
    # Component 2: Probability of unobserved histories
    
    if (type %in% c(7:14,24,25, 37:44)) {
        if (data$learnedresponse) {   ## use model for naive animal
            temp <- makegkcpp(
                as.integer(nrow(realparval0)), # number of rows in lookup table
                as.integer(data$k),            # detectors
                as.integer(data$m),            # mask points
                as.integer(data$detectfn),
                as.integer(.openCRstuff$sigmai[type]),   # column index of sigma in realparval
                as.double(realparval0),
                as.double(unlist(trps)),
                as.double(unlist(data$mask)))
            gk <- temp[[1]]
            hk <- temp[[2]]
        }
        pdot <- rep(0, data$nc)
        for (x in 1:data$details$nmix) {   # loop over latent classes
            
            if (data$details$R) {
                pch1 <-  PCH1secr(
                    type, 
                    x, 
                    1,            # data$nc, 
                    data$J,
                    data$cumss,
                    data$k, 
                    data$m, 
                    realparval0,
                    PIA0[1,,,,drop = FALSE], 
                    PIAJ[1,,,drop = FALSE], 
                    gk, 
                    binomN, 
                    data$usge,
                    data$intervals, 
                    data$moveargsi, 
                    data$movemodel,
                    data$kernel, 
                    data$mqarray, 
                    data$cellsize) 
            }
            else {
                pch1 <-  PCH1secrcpp(
                    as.integer(type),
                    as.integer(x-1),
                    as.integer(1),
                    as.integer(data$J),
                    as.integer(data$cumss),
                    as.integer(data$k),
                    as.integer(data$m),
                    as.double (unlist(realparval0)),
                    as.integer(nrow(realparval0)),
                    as.integer(PIA0[1,,,]),
                    as.integer(PIAJ[1,,]),
                    as.double (gk),             ## gk0
                    as.integer(binomN),
                    as.double (data$usge),
                    as.double (data$intervals),
                    as.integer(data$moveargsi),
                    as.integer(data$movemodel),
                    as.character(data$usermodel),
                    as.integer(nrow(data$kernel)),
                    as.integer(as.matrix(data$kernel)),
                    as.integer(data$mqarray),
                    as.double (data$cellsize))
            }
            pch1 <- rep(pch1, data$nc)      ## let's ignore individual variation!!!
            pdot <- pdot + pmix[x] * pch1
        }
        comp[2] <- - sum(freq * log(pdot))
    }
    #####################################################################
    # Component 3: Probability of observing nc animals
    
    if (type %in% c(7,8,12,13,14,24)) {
        if (type %in% c(7, 12, 13, 24))
            Dsuper <- realparval[nrow(realparval)*3+1] # Dsuper direct
        else  {        # type %in% c(8, 14))
            Dsuper <- getD(type, data$J, data$details$nmix, pmix,
                           realparval, data$design$PIAJ,
                           data$intervals)
        }
        A <- maskarea(data$mask)
        N <- Dsuper * A
        # impose constraint: return with invalid result code if not possible
        if (N < ncf) return;
        
        meanpdot <- data$nc / sum(1/pdot)
        ## cf CLmeanesa in 'secr'
        
        comp[3] <- switch (data$distrib+1,
                           dpois(ncf, N * meanpdot, log = TRUE),
                           lnbinomial (ncf, N, meanpdot),
                           NA)
    }
    

    ## optional multinomial term (not if CJS)
    if (data$details$multinom & !(type %in% c(6,36))) {
        # nh <- table(rep(apply(data$capthist, 1, paste, collapse=''), freq))
        comp[5] <- lgamma(ncf + 1) - sum(lgamma(freq + 1))
    }
    ## log-likelihood as sum of components
    loglik <- sum(comp)
    ## debug
    if (data$details$debug>=1) {
        cat("Likelihood components (comp) ", format(comp, digits=10), "\n")
        cat("Total ", format(loglik, digits = 10), "\n")
        browser()
    }

    .openCRstuff$iter <- .openCRstuff$iter + 1
    if (data$details$trace) {
        if (!is.null(data$details$fixedbeta))
            beta <- beta[is.na(data$details$fixedbeta)]
        cat(format(.openCRstuff$iter, width=4),
            formatC(round(loglik,dig), format='f', digits=dig, width=betaw),
            formatC(beta, format='f', digits=dig+1, width=betaw),
            '\n', sep = " ")
        flush.console()
    }

    if (oneeval) {
        c(loglik, beta)
    }
    else {
        if (is.finite(loglik)) -loglik   # return the negative loglikelihood
        else 1e10
    }

}
############################################################################################

