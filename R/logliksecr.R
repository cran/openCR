# spatial likelihood

# 2019-04-09 split from loglik.R
# 2019-04-09 explicit treatment of count detector; dropuse of data$multi
# 2019-04-14 R option failed with movement in PCH1secr because argument usermodel omitted
# 2019-04-23 removed type 5 'secr' (unused)
# 2019-05-06 1.4.0
# 2019-06-19 onehistory modified for single call to prwisecr (merged prwimulti)
# 2020-09-01 changed return; to return(1e10) in open.secr.loglikfn component 3
# 2020-10-25 fixed bug in open.secr.loglikfn component 3

# types

# CJSsecr = 6

# JSSAsecrf 7
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

#---------------------------------------------------------

open.secr.loglikfn <- function (beta, dig = 3, betaw = 8, oneeval = FALSE, data)

    # Return the negative log likelihood
    # Transformed parameter values are passed in the vector 'beta'
    # details$trace=T sends a one-line report to the screen

{
    ####################################################################
    # functions for one history at a time are used mostly for debugging
    ####################################################################

    onehistory <- function (n) {
        sump <- 0
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
                if (detectr == "multi") data$capthist[n,, drop = FALSE] else data$capthist[n,,, drop = FALSE],   ## 3-D CH
                data$fi[n],
                data$li[n],
                if (detectr %in% c("poissoncount", "multi")) hk else gk, ## precomputed probability or hazard
                realparval,
                PIA[n,,,,drop = FALSE],
                data$design$PIAJ[n,,,drop = FALSE],
                binomN,
                data$usge,
                data$intervals,
                haztemp$h,
                haztemp$hindex[n,,drop = FALSE]+1,
                data$details$CJSp1,
                data$moveargsi,
                data$movementcode,
                data$edgecode,
                get(data$usermodel),
                data$kernel,
                data$mqarray,
                data$cellsize)
            
            sump <- sump + pmix[x,n] * temp
        }
        ## message('n ', n, ' sump ', sump)
        if (sump<=0) NA else freq[n] * log(sump)
    }
    
    ##############################################
    # this is the function that is generally used
    ##############################################
    
    allhistsecrparallel <- function () {
        sump <- numeric(data$nc)
        for (x in 1:nrow(pmix)) {
            hx <- if (detectr == "multi") matrix(haztemp$h[x,,], nrow = data$m) else -1 ## lookup sum_k (hazard)
            hi <- if (detectr == "multi") haztemp$hindex else -1                        ## index to hx
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
                as.double (if (detectr %in% c("multi", "poissoncount")) hk else gk), ## precomputed probability or hazard
                as.matrix (realparval),
                as.integer(PIA),
                as.integer(data$design$PIAJ),
                as.matrix (data$usge),
                as.matrix (hx),                
                as.matrix (hi),      
                as.integer(data$movementcode),
                as.integer(data$edgecode),
                as.character(data$usermodel),  ## 2019-05-07
                as.integer(data$moveargsi),
                as.matrix (data$kernel),
                as.matrix (data$mqarray),
                as.double (data$cellsize))
            sump <- sump + pmix[x,] * temp
        }
        if (any(is.na(sump)) || any(sump<=0)) NA else freq * log(sump)
    }
    
    #####################################################################
    # main line
    
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
    if (data$details$debug>0) print(realparval)    
    
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

    detectr <- detector(trps)[1]
    if (detectr == 'count') {
        detectr <- if (data$binomN == 0) "poissoncount" else "binomialcount"
    }
    binomN <- switch (detectr, multi = -2, proximity = -1, 
                      poissoncount = 0, binomialcount = data$binomN, -9)
    if (binomN < -2)
        stop("open-population secr requires multi, proximity or count detector type")

    if (data$details$debug>1) browser()
    PIA <- data$design$PIA
    PIAJ <- data$design$PIAJ
    if (data$learnedresponse)
        PIA0 <- data$design0$PIA
    else
        PIA0 <- PIA
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
    if (sum(hk)==0) {
        return(1e10)
    }

    if (data$details$debug>0) message ("sum(gk) = ", sum(gk))

    pmix <- fillpmix2(data$nc, data$details$nmix, PIA, realparval)
    S <- ncol(data$capthist)
    #-----------------------------------------

    if (detectr == "multi") {
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
    if (!data$details$R) {
        ## this is the streamlined option; always uses C code
        prwi <- allhistsecrparallel()
    }
    else {
        ## clunky option using R for debugging
        prwi <- sapply(1:data$nc, onehistory, USE.NAMES = FALSE )
    }
    comp[1] <- sum(prwi)

    #####################################################################
    # Component 2: Probability of unobserved histories a^{-n} in likelihood for uniform-D
    
    if ((type %in% c(9,10,11,25,30)) ## CL and require global pdot for component 3
        | (type %in% c(7,8,12,13,14,24,31))) {                   ## all other
        
        if (data$learnedresponse) {   ## overwrite gk with model for naive animal
            temp <- makegkParallelcpp (as.integer(data$detectfn),
                                       as.integer(.openCRstuff$sigmai[type]),
                                       as.integer(data$details$grain),
                                       as.matrix(realparval0),
                                       as.matrix(trps),
                                       as.matrix(data$mask))
            gk <- array(temp[[1]], dim=c(nrow(realparval), data$k, data$m))  # array form for R use
            hk <- array(temp[[2]], dim=c(nrow(realparval), data$k, data$m))  # array form for R use
        }
        ## else use gk as gk0, hk as hk0
        pdot <- rep(0, data$nc)  # vector: one value for each unique observed history
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
                    if (detectr %in% c("multi", "poissoncount")) hk else gk,  ## 2019-05-19
                    binomN,
                    data$usge,
                    data$intervals,
                    data$moveargsi,
                    data$movementcode,
                    data$edgecode,
                    get(data$usermodel),   # bug fixed 2019-04-14, 2019-05-07
                    data$kernel,
                    data$mqarray,
                    data$cellsize)
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
                    as.double (if (detectr %in% c("multi", "poissoncount")) hk else gk),  ## 2019-05-19
                    as.integer(binomN),
                    as.matrix (data$usge),
                    as.double (data$intervals),
                    as.integer(data$moveargsi),
                    as.integer(data$movementcode),
                    as.integer(data$edgecode),
                    as.character(data$usermodel),
                    as.matrix(data$kernel),
                    as.matrix(data$mqarray),
                    as.double (data$cellsize))
            }
            pdot <- pdot + pmix[x] * pch1
        }
        pdot <- rep(pdot, freq)
        comp[2] <- - sum(log(pdot))    ## log(1 / a_i)
    }
    #####################################################################
    # Component 3: Probability of observing nc animals (non-CL types)
    if (type %in% c(7,8,12,13,14,24, 31)) {
        if (type %in% c(7, 12, 13, 24, 31))
            Dsuper <- realparval[nrow(realparval)*3+1] # Dsuper direct
        else  {        # type %in% c(8, 14)) D or B parameterisation
            Dsuper <- getD(type, data$J, data$details$nmix, pmix,
                           realparval, data$design$PIAJ,
                           data$intervals)
        }
        A <- maskarea(data$mask)
        N <- Dsuper * A
        # impose constraint: return with invalid result code if not possible
        if (N < ncf) return (1e10);
        ## meanpdot <- data$nc / sum(1/pdot)
        ## possible bug <1.4.0 did not use freq
        ## meanpdot <- data$nc / sum(freq * 1/pdot)
        ## 2020-10-25 use expanded n (ncf)
        # meanpdot <- ncf / sum(rep(1/pdot,freq))
        meanpdot <- ncf / sum(1/pdot)
        ## cf CLmeanesa in 'secr'

        comp[3] <- switch (data$distrib+1,
                           dpois(ncf, N * meanpdot, log = TRUE),
                           lnbinomial (ncf, N, meanpdot),
                           NA)
    }
    
    #####################################################################
    ## optional multinomial term (not if CJS)
    if (data$details$multinom & !(type %in% c(6))) {
        # nh <- table(rep(apply(data$capthist, 1, paste, collapse=''), freq))
        comp[5] <- lgamma(ncf + 1) - sum(lgamma(freq + 1))
    }
    #####################################################################
    ## log-likelihood as sum of components
    loglik <- sum(comp)
    ## debug
    if (data$details$debug>=1) {
        message("Likelihood components (comp) ", paste(format(comp, digits=10), collapse = ' '))
        message("Total ", format(loglik, digits = 10))
        if (data$details$debug>1) browser()
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
        out <- c(loglik, beta)
        attr(out, 'components') <- comp
        out
    }
    else {
        if (is.finite(loglik)) -loglik   # return the negative loglikelihood
        else 1e10
    }

}
############################################################################################

