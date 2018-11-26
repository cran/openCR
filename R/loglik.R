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

# secr = 5
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

# secrCL = 30
# secrD = 31

#---------------------------------------------------------

## open.loglikfn <- function (beta, dig = 3, betaw = 8, cluster = NULL, oneeval = FALSE, data)
open.loglikfn <- function (beta, dig = 3, betaw = 8, oneeval = FALSE, data)

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
                as.matrix(realparval),
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
                            as.matrix (realparval),
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
        allhistparallel <- function () {
            sump <- numeric(data$nc)
            for (x in 1:nrow(pmix)) {
                temp <-  allhistparallelcpp(
                    as.integer(x-1),
                    as.integer(type),
                    as.integer(data$nc),
                    as.integer(data$details$CJSp1),
                    as.integer(data$details$grain),
                    ## as.matrix (pmix),
                    as.double (data$intervals),
                    as.integer(data$cumss),
                    as.integer(data$capthist),
                    as.integer(data$fi),
                    as.integer(data$li),
                    as.matrix (realparval),
                    as.integer(PIA),
                    as.integer(data$design$PIAJ))

                ## sump <- sump + freq * log(temp)
                sump <- sump + pmix[x,] * temp
            }
            sump  ## return vector of individual LL contributions
            freq * log(sump)  ## return vector of individual LL contributions
        }


        comp <- numeric(5)
        pmix <- fillpmix2(data$nc, data$details$nmix, PIA, realparval)
        #####################################################################
        # Component 1: Probability of observed histories - all models
        # if (is.null(cluster)) {
            if (data$ncores>1)
                temp <- allhistparallel()
            else
                temp <- sapply(1:data$nc, onehistory, pmix = pmix)
        # }
        # else {
        #         temp <- parSapply(cluster, 1:data$nc, onehistory, pmix = pmix, USE.NAMES = FALSE )
        # }
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
                            as.matrix(realparval0),
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
                superN <- getN(type, ncf, data$J, data$details$nmix, pmix, realparval, PIAJ, data$intervals)
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

open.secr.loglikfn <- function (beta, dig = 3, betaw = 8, oneeval = FALSE, data)

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

    binomN <- switch (detector(trps)[1], multi = -1, proximity = -1, count = data$binomN, -9)
    if (binomN < -2)
        stop("open-population secr requires multi, proximity or count detector type")

    if (data$details$debug>2) browser()

    PIA <- data$design$PIA
    PIAJ <- data$design$PIAJ
    if (data$learnedresponse)
        PIA0 <- data$design0$PIA
    else
        PIA0 <- PIA
    #-----------------------------------------

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
                hx <- matrix(haztemp$h[x,,], nrow = data$m)
                temp <- prwisecrcpp(
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
                    as.matrix (hx),               ## lookup sum_k (hazard)
                    as.matrix(haztemp$hindex[n,,drop = FALSE]),      ## index to h
                    as.integer(data$details$CJSp1),
                    as.integer(data$movemodel),
                    as.character(data$usermodel),
                    as.matrix(data$kernel),
                    as.matrix(data$mqarray),
                    as.double (data$cellsize))
            }
            sump <- sump + pmix[x,n] * temp
        }
        freq[n] * log(sump)
    }
    allhistsecrparallel <- function () {
        sump <- numeric(data$nc)
        for (x in 1:nrow(pmix)) {
            hx <- if (data$multi) matrix(haztemp$h[x,,], nrow = data$m) else -1 ## lookup sum_k (hazard)
            hi <- if (data$multi) haztemp$hindex else -1                        ## index to hx
            temp <-  allhistsecrparallelcpp(
                as.integer(x-1),
                as.integer(type),
                as.integer(data$m),
                as.integer(data$nc),
                as.integer(binomN),
                as.integer(data$details$CJSp1),
                as.integer(data$details$grain),
                as.double (data$intervals),
                as.integer(data$cumss),
                as.matrix (data$capthist),     
                as.integer(data$fi),
                as.integer(data$li),
                as.double (if (data$multi) hk else gk),   ## precomputed probability or hazard
                as.matrix (realparval),
                as.integer(PIA),
                as.integer(data$design$PIAJ),
                as.matrix (data$usge),
                as.matrix (hx),                
                as.matrix (hi),      
                as.integer(data$movemodel),
                as.integer(data$moveargsi),
                as.matrix (data$kernel),
                as.matrix (data$mqarray),
                as.double (data$cellsize))
            sump <- sump + pmix[x,] * temp
        }
        freq * log(sump)
    }

    #-----------------------------------------

    ## number of threads was set in openCR.fit
    temp <- makegkParallelcpp (as.integer(data$detectfn),
                               as.integer(.openCRstuff$sigmai[type]),
                               as.integer(data$details$grain),
                               as.matrix(realparval),
                               as.matrix(trps),
                               as.matrix(data$mask))
    gk <- array(temp[[1]], dim=c(nrow(realparval), data$k, data$m))  # array form for R use
    hk <- array(temp[[2]], dim=c(nrow(realparval), data$k, data$m))  # array form for R use
    pmix <- fillpmix2(data$nc, data$details$nmix, PIA, realparval)
    S <- ncol(data$capthist)

    #-----------------------------------------

    if (data$multi) {
        ## R alternative
        if (data$details$R)
            haztemp <- gethR(data$m, PIA, data$usge, hk)
        else
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
    #####################################################################
    ## Vector to store components of log likelihood
    comp <- numeric(5)
    #####################################################################
    # Component 1: Probability of observed histories
    if (data$ncores>1) {
        temp <- allhistsecrparallel()
    }
    else {
        if (data$multi) {
            temp <- sapply(1:data$nc, onehistorymulti, pmix = pmix, hk = hk, USE.NAMES = FALSE )
        }
        else {
            temp <- sapply(1:data$nc, onehistory, pmix = pmix, gk = gk, USE.NAMES = FALSE)
        }
    }
    comp[1] <- sum(temp)
    #####################################################################
    # Component 2: Probability of unobserved histories

    if (type %in% c(7:14,24,25, 30, 31)) {
        if (data$learnedresponse) {   ## use model for naive animal
            gk <- makegkParallelcpp (as.integer(data$detectfn),
                                       as.integer(.openCRstuff$sigmai[type]),
                                       as.integer(data$details$grain),
                                       as.matrix(realparval0),
                                       as.matrix(trps),
                                       as.matrix(data$mask))[[1]]
        }
        ## else use gk as gk0

        pdot <- rep(0, data$nc)
        for (x in 1:data$details$nmix) {   # loop over latent classes
            if (data$details$R) {
                
                pch1 <-  PCH1secr(
                    type,
                    as.logical(data$design0$individual),
                    x,
                    data$nc,
                    data$J,
                    data$cumss,
                    data$k,
                    data$m,
                    realparval0,
                    PIA0,
                    PIAJ,
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
                if (data$ncores==1) {
                    pch1 <-  PCH1secrcpp(
                        as.integer(type),
                        as.logical(data$design0$individual),
                        as.integer(x-1),
                        as.integer(data$nc),
                        as.integer(data$J),
                        as.integer(data$cumss),
                        as.integer(data$k),
                        as.integer(data$m),
                        as.matrix (realparval0),
                        as.integer(PIA0),
                        as.integer(PIAJ),
                        as.double (gk),
                        as.integer(binomN),
                        as.matrix (data$usge),
                        as.double (data$intervals),
                        as.integer(data$moveargsi),
                        as.integer(data$movemodel),
                        as.character(data$usermodel),
                        as.matrix(data$kernel),
                        as.matrix(data$mqarray),
                        as.double (data$cellsize))
                }
                else {
                    pch1 <-  PCH1secrparallelcpp(
                        as.integer(x-1),
                        as.integer(type),
                        as.integer(data$details$grain),
                        as.logical(data$design0$individual),
                        as.integer(data$J),
                        as.integer(data$m),
                        as.integer(data$nc),
                        as.integer(data$cumss),
                        as.matrix (realparval0),
                        as.integer(PIA0),
                        as.integer(PIAJ),
                        as.double (gk),
                        as.integer(binomN),
                        as.matrix (data$usge),
                        as.double (data$intervals),
                        as.integer(data$moveargsi),
                        as.integer(data$movemodel),
                        as.matrix(data$kernel),
                        as.matrix(data$mqarray),
                        as.double (data$cellsize))
                }
            }
            pdot <- pdot + pmix[x] * pch1
        }
        comp[2] <- - sum(freq * log(pdot))
    }
    #####################################################################
    # Component 3: Probability of observing nc animals

    if (type %in% c(7,8,12,13,14,24, 31)) {
        if (type %in% c(7, 12, 13, 24, 31))
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
    #####################################################################

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
    ## optionally display message on console for this iteration
    .openCRstuff$iter <- .openCRstuff$iter + 1
    if (data$details$trace) {
        if ((.openCRstuff$iter %% data$details$trace) == 0) {
            if (!is.null(data$details$fixedbeta))
                beta <- beta[is.na(data$details$fixedbeta)]
            message(format(.openCRstuff$iter, width=4), "   ",
                    formatC(round(loglik,dig), format='f', digits=dig, width=betaw+2),
                    formatC(beta, format='f', digits=dig+1, width=betaw+1),
                    sep = " ")
            flush.console()
        }
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

