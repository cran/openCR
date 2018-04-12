## package 'openCR'
## derived.R
## 2011-12-30, 2013-01-20, 2017-11-21, 2017-12-21, 2018-02-15
## does not deal with mixtures
## does not allow specified newdata or all levels
################################################################################

derived.openCRlist <- function (object, Dscale = 1,  HTbysession = FALSE, ...) {
    lapply(object, derived, Dscale = Dscale, HTbysession = HTbysession, ...)
}

derived.openCR <- function (object, Dscale = 1, HTbysession = FALSE, ...) {

    allvars <- unlist(sapply(object$model, all.vars))
    if ('h2' %in% allvars | 'h3' %in% allvars)
        warning ("derived.openCR does not handle finite mixtures")

    beta <- object$fit$par
    parindx <- object$parindx
    link <- object$link
    fixed <- object$fixed
    design0 <- object$design0
    capthist <- object$capthist
    if (!grepl("secr", object$type)) capthist <- apply(capthist,1:2,sum)
    mask <- object$mask
    details <- object$details
    type <- object$type
    totaln <- sum(covariates(object$capthist)$freq)
    intervals <- object$intervals
    getreal <- function (par) {
        j <- if (grepl("super",par)) 1.0 else J
        fxd <- object$fixed
        if (!is.null(fxd[[par]]))
            rep(fxd[[par]], j)
        else {
            predict(object)[[par]][1:J,'estimate']
        }
    }
    getfbeta <- function (beta, phij) {
        # return fj for inputs phij and beta
        # phij is per session phi (adjusted for intervals)
        d <- beta  # for beta[1]
        for (j in J1) d[j+1] <- d[j] * phij[j] + beta[j+1]
        c(beta[-1]/d[-J], NA)
    }
    getkappa <- function(p, phi, f) {
        J <- length(p)
        tau <- kappa <- numeric(J)
        tau[1] <- 1/p[1]
        kappa[1] <- 1
        for (j in 1:(J-1)) {
            kappa[j+1] <- kappa[j]/p[j] * (1-p[j]) * phi[j] * p[j+1] + tau[j] * f[j] * p[j+1]
            # if(j==(J-1)) browser()
            tau[j+1] <- tau[1] * prod((phi+f)[1:j])
        }
        kappa[1] <- NA
        kappa
    }

    ## alternatively: makerealparameters (design0, beta, parindx, link, fixed)
    J <- length(intervals) + 1
    J1 <- 1:(J-1)
    nc <- nrow(capthist)
    persession <- function (x) c(x[-J]^intervals, NA)
    stdrate <- function (x) c(x[-J]^(1/intervals), NA)
    # note b, B always on per interval basis - no scaling

    phi <- getreal('phi')
    phij <- persession(phi)

    fj <- NULL
    out <- NULL

    if (type %in% c('CJS', 'CJSsecr')) {
        warning ("derived.openCR is not implemented for type ", type)
        return(NULL)
    }

    if (type %in% "JSSAN") {
        Nj <- getreal('N')
        superN <- Nj[1] + sum (Nj[-1] - (Nj * phij)[-J])
        B <- c(Nj[1], Nj[-1] - (Nj*phij)[-J])
        b <- B/sum(B)
    }
    else {
        if (type %in% "JSSAsecrD") {
            Dj <- getreal('D')
            B <- c(Dj[1], Dj[-1] - (Dj*phij)[-J])
            superD <- sum(B)
            b <- B/sum(B)
            fj <- c(B[-1]/Dj[-J], NA)
        }
        ## otherwise extract session-specific b
        else if (type %in% c('JSSAb','JSSAbCL','JSSAsecrb','JSSAsecrbCL','JSSARET')) {
            b <- getreal('b') # predicted$b[,'estimate']
        }
        else if (type %in% c('JSSAf','JSSAfCL','JSSAsecrf','JSSAsecrfCL',
                             'JSSAg','JSSAgCL','JSSAsecrg','JSSAsecrgCL', 'Pradelg',
                             'JSSAk', 'JSSAkCL',
                             'JSSAl','JSSAlCL','JSSAsecrl','JSSAsecrlCL', 'Pradel')) {

            if (type %in% c('JSSAl','JSSAlCL','JSSAsecrl','JSSAsecrlCL', 'Pradel')) {
                lambda <- getreal('lambda')
                f <- lambda-phi
            }
            else if (type %in% c('JSSAg','JSSAgCL','JSSAsecrg','JSSAsecrgCL', 'Pradelg')) {
                gamma <- getreal('gamma')
                gam1 <- c(gamma[-1],NA)
                f <- phi * (1/gam1 - 1)
            }
            else if (type %in% c('JSSAk','JSSAkCL')) {
                kappa <- getreal('kappa')
                kappa <- c(1, kappa[-1])
                f <- tau <- numeric(J)
                p <- getreal('p')
                tau[1] <- 1/p[1]
                for (j in 1:(J-1)) {
                    f[j] <- (kappa[j+1] - kappa[j]/p[j] * (1 - p[j]) * phi[j] * p[j+1]) / (tau[j] * p[j+1])
                    tau[j+1] <- tau[1] * prod(phi[1:j] + f[1:j])
                }
            }
            else {
                f <- getreal('f')
            }
            lambdaj <- persession(phi + f)
            fj <- persession(f)

            d <- numeric(J)
            d[1] <- 1;
            for (j in 2:J) d[j] <- d[j-1] * lambdaj[j-1]

            b <- numeric(J)
            b[1] = 1;
            for (j in 2:J) b[j] <- fj[j-1] * d[j-1]
            b <- b / sum(b)

        }
        else if (type %in% "JSSAB") {
            BN <- getreal('BN')
            b <-  BN / sum(BN)
        }
        else if (type %in% "JSSAsecrB") {
            BD <- getreal('BD')
            b <-  BD / sum(BD)
        }
        else stop("unrecognized model type")
    }

    if (is.null(fj)) fj <- getfbeta(b, phij)

    df <- data.frame(session = 1:J, JS.counts(object$capthist),
                     time = cumsum(c(0,intervals)))
    parnames <- names(predict(object))
    for (i in parnames) {
        df[,i] <- getreal(i)
    }
    #
    # if (grepl('secr', type)) {
    #     df$lambda0 <- getreal('lambda0')
    #     df$sigma <- getreal('sigma')
    # }
    # else {
    #     df$p <- getreal('p')
    # }
    # df$phi <- phi


    df$b <- b
    df$lambda <- stdrate(phij+ fj)
    df$f <- df$lambda - df$phi
    df$gamma <- c(NA, (df$phi / (df$f + df$phi))[-J])
    df$kappa <- getkappa (df$p, df$phi, df$f)
    ## Nonspatial
    ## superN and time-specific Nj
    if (type %in% c('JSSAN',
                    'JSSAf','JSSAfCL','JSSAl','JSSAlCL', 'JSSAg','JSSAgCL',
                    'JSSAk', 'JSSAkCL',
                    'Pradel','Pradelg',
                    'JSSAb','JSSAbCL','JSSARET', 'JSSAB')) {

        if (type %in% "JSSAN") {
            ## as before...
            Nj <- getreal('N')
            superN <- Nj[1] + sum (Nj[-1] - (Nj * phij)[-J])
        }
        else {
            if (type %in% c('JSSAf','JSSAb','JSSAl','JSSAk','JSSAg','JSSARET'))
                superN <- getreal('superN')[1]
            else if (type %in% "JSSAB" )
                superN <- sum(BN)
            else {
                p <- openCR.pdot(object, bysession = FALSE)
                superN <- sum(1 / rep(p, covariates(object$capthist)$freq))
            }
            # assume no mixture for now but see Pledger et al 2010 p885
            Nj <- numeric(J)
            if (HTbysession) {
                p <- openCR.pdot(object, bysession = TRUE)  ## session x ch
                sess <- primarysessions(intervals(object$capthist)) ## differs from 'intervals'
                OK <- apply(capthist,1, by, sess, sum)>0
                freq <- sweep(OK, MARGIN=2, STATS=covariates(object$capthist)$freq, FUN = "*")
                p1 <- lapply(1:J, function(j) 1/ rep(p[j,], freq[j,]))
                Nj <- sapply(p1, sum)
            }
            else {
                Nj[1] <- superN * b[1]
                for (j in 2:J) Nj[j] <- Nj[j-1] * phij[j-1] + superN*b[j]
            }
        }
        df$N <- Nj
        out <- list(totalobserved = totaln, parameters = parnames, superN = superN, estimates = df)
    }

    ## Spatial
    ## superD and time-specific Dj
    if (type %in% c('JSSAsecrf','JSSAsecrfCL','JSSAsecrl','JSSAsecrlCL',
                    'JSSAsecrg','JSSAsecrgCL',
                    'JSSAsecrb','JSSAsecrbCL', 'JSSAsecrB',
                    'JSSAsecrD')) {
        if (type %in% c('JSSAsecrf','JSSAsecrb','JSSAsecrl','JSSAsecrg'))
            superD <- getreal('superD')[1]
        else if (type %in% "JSSAsecrB" )
            superD <- sum(BD)
        else {
            a <- openCR.esa (object, bysession = FALSE)
            superD <- sum(1/ rep(a, covariates(object$capthist)$freq))
        }
        # assume no mixture for now but see Pledger et al 2010 p885
        Dj <- numeric(J)
        if (HTbysession) {
            a <- openCR.esa(object, bysession = TRUE)  ## session x ch
            sess <- primarysessions(intervals(object$capthist)) ## differs from 'intervals'
            OK <- apply(capthist, 1, by, sess, sum)>0
            freq <- sweep(OK, MARGIN=2, STATS=covariates(object$capthist)$freq, FUN = "*")
            a1 <- lapply(1:J, function(j) 1/ rep(a[j,], freq[j,]))
            Dj <- sapply(a1, sum)
        }
        else {
            Dj[1] <- superD * b[1]
            for (j in 2:J) Dj[j] <- Dj[j-1] * phij[j-1] + superD*b[j]
        }

        df$D <- Dj

        out <- list(totalobserved = totaln, parameters = parnames, superD = superD,
                    estimates = df, Dscale = Dscale)
    }
    class (out) <- c("derivedopenCR","list")
    out
}

print.derivedopenCR <- function (x, Dscale = NULL, legend = FALSE, ...) {

    args <- list(...)
    args$x <- x$estimates
    if (!"row.names" %in% names(args))
        args$row.names <- FALSE
    if (!"digits" %in% names(args))
        args$digits <- 4

    if (is.null(Dscale)) Dscale <- x$Dscale

    D.units <- switch(paste0('D',Dscale),
                      D1 = 'per ha',
                      D100 = 'per km^2',
                      D10000 = 'per 100 km^2',
                      D100000 = 'per 1000 km^2',
                      paste0('per ', Dscale, ' ha'))

    fields <- c(
        '-------@-----------------------------------------',
        'session@primary session',
        'n@number observed',
        'R@number released',
        'm@number already marked',
        'r@number recaptured in later session',
        'z@number known alive but not caught',
        'time@accumulated time since start',
        'p@detection probability per secondary session',
        'lambda0@intercept of detection function per secondary session',
        'sigma@spatial scale of detection function in metres',
        'p@detection probability per secondary session',
        'phi@apparent survival per unit time',
        'b@entry probabilities',
        'lambda@population growth rate per unit time',
        'f@per capita recruitment per unit time',
        'gamma@seniority (cf reverse-time phi)',
        'kappa@recruitment parameter of Link and Barker (2005)',
        'BN@number of recruits',
        'BD@density of recruits',
        'N@population size',
        paste0('D@density  ', D.units))
    leg <- data.frame(do.call(rbind, strsplit(fields, '@')))
    names(leg) <- c('Field', 'Definition')

    cat ("Total number observed", x$totalobserved, "\n")
    cat ("Parameters in model", paste(x$parameters, collapse=', '), "\n")
    if (!is.null(x$superN))
        cat ("Superpopulation size", x$superN, "\n")
    if (!is.null(x$superD)) {
        cat ("Superpopulation density", x$superD  * Dscale, D.units, "\n")
        x$estimates$D <- x$estimates$D * Dscale
    }
    cat ("Session-specific counts and estimates:\n\n")
    do.call(print.data.frame, args)
    cat ("\n")

    if (legend) {
        fields <- names(x$estimates)
        ord <- match(fields, leg[,'Field'])
        leg <- leg[c(1,ord),]
        OK <- leg$Field %in% names(x$estimates) | leg$Field == '-------'
        print(leg[OK ,,drop=FALSE], right=FALSE, row.names=FALSE)
        cat("\n")
    }


    }
