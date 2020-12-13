# non-spatial likelihood

# 2017-05-17 not ready for mixtures
# 2017-05-17 pass just one animal to onehistory
# 2017-05-17 PCH0 called with single representative history
# 2017-05-17 capthist 2D for non-spatial
# 2017-05-22 use data argument (an environment)
# 2018-02-06 drop 'par.' prefix from function names
# 2018-03-26 switch pch0 to pch1
# 2019-04-09 ejected open.secr.loglikfn to logliksecr.R
# 2019-04-09 explicit treatment of count detector; dropuse of data$multi
# 2020-12-07 CJSmte (Markovian Temporary Emigration) trial - not completed

# types

# CJS = 1
# JSSAb = 2
# JSSAl = 3
# JSSAf = 4
# CJSmte = 5
# JSSAfCL = PLBf = 15
# JSSAlCL = PLBl = 16
# JSSAbCL = PLBb = 17
# JSSAB = 18
# JSSAN = 19
# Pradel = 20
# JSSARET = 21
# JSSAg = 22
# JSSAgCL = PLBg = 23
# Pradelg = 26
# JSSAfgCL = 27
# JSSAk = 28
# JSSAkCL = PLBk = 29

#---------------------------------------------------------

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
        message ('beta vector : ', paste(beta, collapse=', '))
        message ('real vector : ', paste(realparval, collapse=','))
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
        message('realparval')
        print(realparval)
        message('table(PIA)')
        print (table(PIA))
        message('data$intervals')
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
    # optional code using m.array
    # else if (data$details$R & (type %in% c(1))) {
    #     if (data$details$nmix>1) stop ("R CJS does not use mixture")
    #     comp <- CJSloglik (
    #         type,
    #         x = 1,
    #         data$capthist,
    #         data$marray,
    #         data$J,
    #         data$cumss,
    #         nmix = 1,
    #         realparval,
    #         PIA[1,,,,drop=FALSE],
    #         PIAJ[1,,,drop=FALSE],
    #         data$intervals)
    # }
    else {
        onehistory <- function (n, pmix) {
            sump <- 0
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
                    data$details$CJSp1
                    # , data$moveargsi
                )
                sump <- sump + pmix[x,n] * temp
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
                    as.double (data$intervals),
                    as.integer(data$cumss),
                    as.integer(data$capthist),
                    as.integer(data$fi),
                    as.integer(data$li),
                    as.matrix (realparval),
                    as.integer(PIA),
                    as.integer(data$design$PIAJ))
                sump <- sump + pmix[x,] * temp
            }
            freq * log(sump)  ## return vector of individual LL contributions
        }


        comp <- numeric(5)
        pmix <- fillpmix2(data$nc, data$details$nmix, PIA, realparval)
        
        #####################################################################
        # Component 1: Probability of observed histories - all models
        # if (is.null(cluster)) {
        if (data$details$R)
            temp <- sapply(1:data$nc, onehistory, pmix = pmix)
        else
            temp <- allhistparallel()
        comp[1] <- sum(temp)
        
        #####################################################################
        # Component 2: Probability of missed animals (all-zero histories)
        # not CJS, CJSmte
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
        # not CJS, CJSmte
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
        message("Likelihood components (comp) ", format(comp, digits=10))
        message("Total ", format(loglik, digits = 10))
        if (data$details$debug>1) browser()
    }

    .openCRstuff$iter <- .openCRstuff$iter + 1   ## moved outside loop 2011-09-28
    if (data$details$trace) {
        if (!is.null(data$details$fixedbeta))
            beta <- beta[is.na(data$details$fixedbeta)]
        # cat(format(.openCRstuff$iter, width=4),
        #     formatC(round(loglik,dig), format='f', digits=dig, width=10),
        #     formatC(beta, format='f', digits=dig+1, width=betaw),
        #     '\n', sep = " ")
        message(format(.openCRstuff$iter, width=4), ' ',
            formatC(round(loglik,dig), format='f', digits=dig, width=10), ' ',
            paste(formatC(beta, format='f', digits=dig+1, width=betaw), 
                collapse = ' '))
        
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

