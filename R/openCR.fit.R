################################################################################
## package 'openCR'
## openCR.fit.R
## 2011-12-30, 2013-01-21
## 2013-01-21 modified to balance with terminal beta, tau
## 2015-01-30 removed make.lookup (see secr)
## 2015-01-30 moved miscellaneous functions to utility.r
## 2015-01-31 deleted utility functions not needed or can be called from secr
## 2015-02-06 reconciled this current version with forked 1.2.0
## 2015-02-06 removed pdot, esa, derived
## 2017-05-15 reconciled versions; revision in progress
## 2017-05-18 2.2.0 ditched old openCR.fit; renamed openCR.MCfit
## 2017-11-20 general revision 2.2.1
## 2017-11-20 refined start options
## 2017-11-20 method default Newton-Raphson
## 2018-01-20 remember compileAttributes('d:/open populations/openCR')
## 2018-01-25 single coerced to multi
## 2018-02-02 detectfn HHR etc.
## 2018-05-01 intermediate variable allbetanames to fix problem with fixedbeta
## 2018-10-29 CJSp1 argument for openCR.design
## 2018-11-20 dropped posterior (see classMembership method)
## 2019-04-07 check whether require autoini
## 2019-04-09 removed data$multi
## 2020-10-19 added agecov
## 2020-11-02 movemodel renamed movementcode
## 2020-12-07 CJSmte experimental - trial abandoned 2020-12-12

################################################################################

openCR.fit <- function (capthist, type = "CJS", model = list(p~1, phi~1, sigma~1),
  distribution = c("poisson", "binomial"), mask = NULL, 
  detectfn = c('HHN','HHR','HEX','HAN','HCG','HVP'), binomN = 0, 
  movementmodel = c('static','uncorrelated','normal','exponential', 't2D', 'uniform'),
  edgemethod = c('truncate', 'wrap', 'none'), start = NULL, link = list(), 
  fixed = list(), timecov = NULL, sessioncov = NULL, agecov = NULL, 
  dframe = NULL, dframe0 = NULL, details = list(), method = 'Newton-Raphson', 
  trace = NULL, ncores = NULL, ...)
  
{
    # Fit open population capture recapture model
    #
    # Some arguments:
    #
    #  capthist   -  capture history object (includes traps object as an attribute)
    #  model      -  formulae for real parameters in terms of effects and covariates
    #  start      -  start values for maximization (numeric vector link scale);
    #  link       -  list of parameter-specific link function names 'log', 'logit', 'loglog',
    #                'identity', 'sin', 'neglog', 'mlogit'
    #  fixed      -  list of fixed values for named parameters
    #  sessioncov -  dataframe of session-level covariates
    #  mask
    #  detectfn
    #  dframe     -  optional data frame of design data for detection model (tricky & untested)
    #  details    -  list with several additional settings, mostly of special interest
    #  method     -  optimization method (indirectly chooses
    #  trace      -  logical; if TRUE output each likelihood as it is calculated
    #  ...        -  other arguments passed to join()

    #########################################################################
    ## Use input 'details' to override various defaults
  defaultdetails <- list(
    autoini = NULL, 
    CJSp1 = FALSE, 
    contrasts = NULL, 
    control = list(),
    debug = 0, 
    grain = 1,
    hessian = 'auto', 
    ignoreusage = FALSE, 
    initialage = 0, 
    LLonly = FALSE,
    kernelradius = 10,
    minimumage = 0, 
    maximumage = 1,
    multinom = FALSE,
    R = FALSE, 
    squeeze = TRUE, 
    trace = FALSE
    )
  
  if (is.logical(details$hessian))
    details$hessian <- ifelse(details$hessian, 'auto', 'none')
  details <- replace (defaultdetails, names(details), details)
    if (!is.null(trace)) details$trace <- trace
    if (details$LLonly)  details$trace <- FALSE
    if (details$R) ncores <- 1    ## force 2018-11-12
    if (!is.null(ncores) && (ncores == 1)) details$grain <- -1
    #########################################################################

    distribution <- match.arg(distribution)
    distrib <- switch (distribution, poisson = 0, binomial = 1)
    
    ##############################################
    # Multithread option 2018-04-11, 2020-11-02
    ##############################################

    secr::setNumThreads(ncores, stackSize = "auto")  # change to match secr 
    
    
    if (is.character(detectfn)) {
        detectfn <- match.arg(detectfn)
        detectfn <- secr:::detectionfunctionnumber(detectfn)
    }
    if (is.character(dframe)) {
        dframename <- dframe; rm(dframe)
        dframe <- get(dframename, pos=-1)
    }
    if (is.character(dframe0)) {
        dframename <- dframe0; rm(dframe0)
        dframe0 <- get(dframename, pos=-1)
    }
    if (is.character(capthist)) {
        capthistname <- capthist; rm(capthist)
        capthist <- get(capthistname, pos=-1)
    }
    ##############################################
    ## Standard form for capthist
    ##############################################

    marray <- m.array(capthist)
    capthist <- stdcapthist(capthist, type, details$nclone, details$squeeze,  ...)
    inputcapthist <- capthist  ## PROCESSED
    ##############################################
    ## check type argument
    ##############################################

    if (type %in% .openCRstuff$suspendedtypes)
        stop (type, " not currently available")
    secr <- grepl('secr', type)
    if (secr) {
        if (is.null(mask))
            stop("requires valid mask")
        if (ms(mask)) {
            mask <- mask[[1]]
            warning("multi-session mask provided; using first")
        }
        if (is.character(mask)) {
            maskname <- mask; rm(mask)
            mask <- get(maskname, pos=-1)
        }
        
        if (is.function (movementmodel)) {
            moveargs <- formalArgs(movementmodel)
            usermodel <- as.character(substitute(movementmodel))
            movementmodel <- "user"
        }
        else {
            usermodel <- ""
            movementmodel <- match.arg(movementmodel)
        }
        ## integer code for movement model
        movementcode <- movecode(movementmodel)
        edgemethod <- match.arg(edgemethod)
        edgecode <- edgemethodcode(edgemethod)  # 0 none, 1 wrap, 2 truncate
        if (!is.null(mask) && attr(mask, 'type') != 'traprect' && 
                movementcode > 1 && edgemethod == 'wrap') {
            stop("edgemethod = 'wrap' requires mask of type 'traprect'")        
        }
        if (movementcode > 1 && edgemethod == 'none') {
            warning("specify edgemethod 'wrap' or 'truncate' to avoid ", 
              "bias in movement models")
        }
    }
    else {
        if (!is.null(mask)) warning("mask not used in non-spatial analysis")
        movementcode <- -1
        edgecode <- -1
        usermodel <- ""
    }

    ##############################################
    ## Remember start time and call
    ##############################################

    ptm  <- proc.time()
    starttime <- format(Sys.time(), "%H:%M:%S %d %b %Y")
    cl   <- match.call(expand.dots = TRUE)

    if (type %in% c("secrCL","secrD"))
        intervals(capthist) <- rep(0, ncol(capthist)-1)
    intervals <- intervals(capthist)
    intervals <- intervals[intervals>0]          ## primary intervals only
    sessnames <- sessionlabels(capthist)
    cumss <- getcumss(capthist)                  ## cumulative secondary sessions per primary session
    if (is.null(sessnames))
        sessnames <- 1:(length(cumss)-1)
    nc <- nrow(capthist)
    if (nc == 0) warning ("no detection histories")
    J <- length(cumss)-1                         ## number of primary sessions
    primarysession <- primarysessions(intervals(capthist)) # rep(1:J, diff(cumss))      ## map secondary to primary
    k <- nrow(traps(capthist))                   ## number of detectors (secr only)
    m <- if (is.null(mask)) 0 else nrow(mask)
    marea <- if (is.null(mask)) NA else maskarea(mask)

    ##############################################
    ## Use input formula to override defaults
    ##############################################

    if ('formula' %in% class(model)) model <- list(model)
    model <- secr:::stdform (model)  ## named, no LHS
    defaultmodel <- list(p = ~1, lambda0 = ~1, phi = ~1, b = ~1, f = ~1, lambda = ~1, g = ~1,
                         gamma = ~1, kappa = ~1, BN = ~1, BD = ~1, N=~1, D = ~1, superN = ~1,
                         superD = ~1, sigma = ~1, z = ~1, move.a = ~1, move.b = ~1, tau = ~1)
    model <- replace (defaultmodel, names(model), model)
    
    pnames <- switch (type,
        CJS = c('p', 'phi'),                                      # 1
        CJSmte = c('p', 'phi', 'move.a', 'move.b'),               # 5
        JSSAb = c('p', 'phi','b','superN'),                       # 2
        JSSAl = c('p', 'phi','lambda','superN'),                  # 3
        JSSAf = c('p', 'phi','f','superN'),                       # 4
        JSSAg = c('p', 'phi','gamma','superN'),                   # 22
        JSSAk = c('p', 'phi','kappa','superN'),                   # 28

        JSSAfCL = c('p', 'phi','f'),                              # 15
        JSSAlCL = c('p', 'phi','lambda'),                         # 16
        JSSAbCL = c('p', 'phi','b'),                              # 17
        JSSAgCL = c('p', 'phi','gamma'),                          # 23
        JSSAkCL = c('p', 'phi','kappa'),                          # 29
        
        PLBf = c('p', 'phi','f'),                                 # 15
        PLBl = c('p', 'phi','lambda'),                            # 16
        PLBb = c('p', 'phi','b'),                                 # 17
        PLBg = c('p', 'phi','gamma'),                             # 23
        PLBk = c('p', 'phi','kappa'),                             # 29
        
        JSSAB = c('p', 'phi','BN'),                               # 18
        JSSAN = c('p', 'phi','N'),                                # 19
        Pradel = c('p', 'phi','lambda'),                          # 20
        Pradelg = c('p', 'phi','gamma'),                          # 26
        JSSARET = c('p', 'phi','b','superN','tau'),               # 21
        
        JSSAfgCL = c('p', 'phi','f','g'),                         # 27    # experimental temporary emigration
        
        CJSsecr = c('lambda0', 'phi','sigma'),                    # 6

        JSSAsecrfCL = c('lambda0', 'phi','f','sigma'),            # 9
        JSSAsecrlCL = c('lambda0', 'phi','lambda','sigma'),       # 10
        JSSAsecrbCL = c('lambda0', 'phi','b','sigma'),            # 11
        JSSAsecrgCL = c('lambda0', 'phi','gamma','sigma'),        # 25
        
        PLBsecrf = c('lambda0', 'phi','f','sigma'),               # 9
        PLBsecrl = c('lambda0', 'phi','lambda','sigma'),          # 10
        PLBsecrb = c('lambda0', 'phi','b','sigma'),               # 11
        PLBsecrg = c('lambda0', 'phi','gamma','sigma'),           # 25
        
        JSSAsecrf = c('lambda0', 'phi','f','superD','sigma'),       # 7
        JSSAsecrl = c('lambda0', 'phi','lambda','superD','sigma'),  # 12
        JSSAsecrb = c('lambda0', 'phi','b','superD','sigma'),       # 13
        JSSAsecrg = c('lambda0', 'phi','gamma','superD','sigma'),   # 24
        JSSAsecrB = c('lambda0', 'phi','BD','sigma'),               # 14
        JSSAsecrD = c('lambda0', 'phi','D','sigma'),                # 8
        secrCL = c('lambda0', 'phi', 'b','sigma'),                # 30
        secrD = c('lambda0', 'phi', 'b', 'superD', 'sigma'),      # 31
        
        "unrecognised type")
    
    moveargsi <- c(-2,-2)
    if (secr) {
        if (movementmodel %in% c('normal','exponential')) {
            pnames <- c(pnames, 'move.a')
            moveargsi[1] <- .openCRstuff$sigmai[typecode(type)] + 1 + (detectfn %in% c(15,17,18,19))
        }
        else if (movementmodel == 't2D') {
                pnames <- c(pnames, 'move.a', 'move.b')
                moveargsi[1] <- .openCRstuff$sigmai[typecode(type)] + 1 + (detectfn %in% c(15,17,18,19))
                moveargsi[2] <- moveargsi[1]+1
        }
        else if (movementmodel == 'user') {
            if (! ("r" == moveargs[1]))
                stop ("user-supplied movement model must have r as first argument")
            if ("a" %in% moveargs) {
                pnames <- c(pnames, 'move.a')
                moveargsi[1] <- .openCRstuff$sigmai[typecode(type)] + 1 + (detectfn %in% c(15,17,18,19))
                if ("b" %in% moveargs) {
                    pnames <- c(pnames, 'move.b')
                    moveargsi[2] <- moveargsi[1] + 1
                }
            }
        }
        else if (movementmodel == 'uniform') {
            ## no parameters, no action needed
        }
        if (type %in% c("secrCL","secrD")) {
            ## closed population
            ## fix survival and recruitment
            fixed <- replace(list(phi = 1.0, b = 1.0), names(fixed), fixed)
        }
    }

    if (any(pnames == 'unrecognised type'))
        stop ("'type' not recognised")
    if (detectfn %in% c(15,17:19)) pnames <- c(pnames, 'z')
    
    ########################################
    # Finite mixtures
    ########################################
    nmix <- secr:::get.nmix(model, capthist, NULL)
    if ((nmix>1) & (nmix<4)) {
        if (type %in% c('Pradel', 'Pradelg')) stop ("Mixture models not implemented for Pradel models")
        model$pmix <- as.formula(paste('~h', nmix, sep=''))
        if (!all(all.vars(model$pmix) %in% c('session','g','h2','h3')))
            stop ("formula for pmix may include only 'session', 'g' or '1'")
        pnames <- c(pnames, 'pmix')
    }
    details$nmix <- nmix
    
    if (type == 'CJSmte') {
      moveargsi[1] <- 1 + nmix
      moveargsi[2] <- moveargsi[1] + 1
    }
    
    #################################
    # Link functions (model-specific)
    #################################
    defaultlink <- list(p = 'logit', lambda0 = 'log', phi = 'logit', b = 'mlogit', f = 'log',
                        gamma = 'logit', kappa = 'log', g = 'logit',
                        lambda = 'log', BN = 'log', BD = 'log', D = 'log', N = 'log',
                        superN = 'log', superD = 'log', sigma = 'log', z = 'log', pmix='mlogit',
                        move.a = 'log', move.b = 'log', tau = 'mlogit')
    link <- replace (defaultlink, names(link), link)
    link[!(names(link) %in% pnames)] <- NULL
    if (details$nmix==1) link$pmix <- NULL

    pnamesR <- pnames[!(pnames %in% names(fixed))]
    model[!(names(model) %in% pnamesR)] <- NULL
    if ((length(model) == 0) & (length(fixed)>0))
        stop ("all parameters fixed")   ## assume want only LL
    vars <-  unlist(lapply(model, all.vars))

    ##############################################
    # Prepare detection design matrices and lookup
    ##############################################
    memo ('Preparing design matrices', details$trace)
    design <- openCR.design (capthist, model, type,
      timecov = timecov,
      sessioncov = sessioncov,
      agecov = agecov,
      dframe = dframe,
      naive = FALSE,
      contrasts = details$contrasts,
      initialage = details$initialage,
      minimumage = details$minimumage,
      maximumage = details$maximumage,
      CJSp1 = details$CJSp1)
    allvars <- unlist(lapply(model, all.vars))
    learnedresponse <- any(.openCRstuff$learnedresponses %in% allvars) || !is.null(dframe)
    mixturemodel <- "h2" %in% allvars | "h3" %in% allvars
    design0 <- if (learnedresponse) {
        if (is.null(dframe0)) dframe0 <- dframe
        openCR.design (capthist, model, type,
          timecov = timecov,
          sessioncov = sessioncov,
          agecov = agecov,
          dframe = dframe0,
          naive = TRUE,
          contrasts = details$contrasts,
          initialage = details$initialage,
          minimumage = details$minimumage,
          maximumage = details$maximumage,
          CJSp1 = details$CJSp1)
    }
    else {
        design
    }
    
    ############################
    # Parameter mapping (general)
    #############################
    np <- sapply(design$designMatrices, ncol)
    NP <- sum(np)
    parindx <- split(1:NP, rep(1:length(np), np))
    names(parindx) <- names(np)

    ##########################
    # Movement kernel
    ##########################

    cellsize <- mqarray <- 0
    kernel <- mqarray <- matrix(0,1,2)  ## default
    if (secr && (movementmodel %in%  
            c('normal', 'exponential', 'user', 't2D', 'uniform'))) {
        ## movement kernel
        k2 <- details$kernelradius
        cellsize <- attr(mask,'area')^0.5 * 100   ## metres, equal mask cellsize
        kernel <- expand.grid(x = -k2:k2, y = -k2:k2)
        kernel <- kernel[(kernel$x^2 + kernel$y^2) <= (k2+0.5)^2, ]
        mqarray <- mqsetup (mask, kernel, cellsize, edgecode)  
    }

    ###########################################
    # Choose likelihood function 
    ###########################################

    if (secr)
        loglikefn <- open.secr.loglikfn  # see logliksecr.R
    else
        loglikefn <- open.loglikfn       # see loglik.R

    ##########################
    # Variable names (general)
    ##########################

    allbetanames <- unlist(sapply(design$designMatrices, colnames))
    names(allbetanames) <- NULL
    realnames <- names(model)
    allbetanames <- sub('..(Intercept))','',allbetanames)
    ## allow for fixed beta parameters 
    if (!is.null(details$fixedbeta))
        betanames <- allbetanames[is.na(details$fixedbeta)]
    else
        betanames <- allbetanames
    betaw <- max(c(nchar(betanames),8))  # for 'trace' formatting

    ###################################################
    # Option to generate start values from previous fit
    ###################################################
    if (inherits(start, 'secr') | inherits(start, 'openCR')) {
        start <- mapbeta(start$parindx, parindx, coef(start)$beta, NULL)
    }
    else if (is.list(start) & (inherits(start[[1]], 'secr') | inherits(start[[1]], 'openCR')) ) {
        start2 <- if (length(start)>1) mapbeta(start[[2]]$parindx, parindx, coef(start[[2]])$beta, NULL) else NULL
        start <- mapbeta(start[[1]]$parindx, parindx, coef(start[[1]])$beta, NULL)
        if (!is.null(start2)) {
            start[is.na(start)] <- start2[is.na(start)]  ## use second as needed
        }
    }
    else if (is.numeric(start) & !is.null(names(start))) {
        ## optionally reorder and subset beta values by name
        OK <- allbetanames %in% names(start)
        if (!all(OK))
            stop ("beta names not in start : ", paste(allbetanames[!OK], collapse=', '))
        start <- start[allbetanames]
    }
    ###############################
    # Start values (model-specific)
    ###############################

    if (is.null(start)) start <- rep(NA, NP)
    freq <- covariates(capthist)$freq
    ncf <- if (!is.null(freq)) sum(freq) else nc

    if (any(is.na(start)) | is.list(start)) {
        rpsv <- if(secr) RPSV(capthist, CC = TRUE) else NA
        ## assemble start vector
  
        default <- list(
            p = 0.6,
            lambda0 = 0.6,
            phi = 0.7,
            gamma = 0.7,
            kappa = 2,
            b = 0.1,
            f = 0.3,
            lambda = 1.0,
            g = 0.2,   # random temporary emigration parameter
            # tau = 1/(details$M+1),
            BN = 20,
            BD = (ncf + 1) / marea,
            D = (ncf + 1) / marea,
            N = ncf + 1,
            # superN = ncf + 20,
            superN = ncf*(1-distrib) + 20,   ## use N-n for binomial 2018-03-12
            superD = (ncf + 20) / marea,
            # superD = (ncf*(1-distrib) + 20) / marea,  ## not a good idea 2018-05-28
            sigma =  rpsv,
            z = 2,
            move.a = if (secr) rpsv/2 else 0.6,
            move.b = if (secr) 1 else 0.2,
            pmix = 0.25
        )

        getdefault <- function (par) transform (default[[par]], link[[par]])
        defaultstart <- rep(0, NP)
        for ( i in 1:length(parindx) ) {
            defaultstart[parindx[[i]][1]] <- getdefault (names(model)[i])
        }

        if(details$nmix>1) {
            ## scaled by mlogit.untransform
            defaultstart[parindx[['pmix']]] <- (2:details$nmix)/(details$nmix+1)
        }
        
        if('b' %in% names(parindx)) {
            ## scaled by mlogit.untransform
            defaultstart[parindx[['b']]] <- 1/J
        }

        # if('tau' %in% names(parindx))
        #     ## scaled by mlogit.untransform
        #     defaultstart[parindx[['tau']]] <- 1/(details$M+1)

        if (secr & !(type %in% c('CJSsecr'))) {
            
            requireautoini <- (is.null(start) | !all(names(parindx) %in% names(start)))   
            if (requireautoini) {   ## condition added 2019-04-07
                if (!is.null(details$autoini))
                    start3 <- autoini (subset(capthist, occasions = primarysession==details$autoini), mask)
                else
                    start3 <- autoini (capthist, mask)  ## 2019-05-07
                    
                if (any(is.na(unlist(start3))))
                    warning ("initial values not found")
                defaultstart[parindx[['lambda0']][1]] <- transform (-log(1-start3[['g0']]), link[['lambda0']])
                defaultstart[parindx[['sigma']][1]] <- transform (start3[['sigma']], link[['sigma']])
                if (type == 'JSSAsecrD')
                    defaultstart[parindx[['D']][1]] <- transform (start3[['D']], link[['D']])
                else if (type == 'JSSAsecrB')
                    defaultstart[parindx[['BD']][1]] <- transform (start3[['D']]/J, link[['BD']])
                else if (type %in% c('JSSAsecrf','JSSAsecrl','JSSAsecrb', 'JSSAsecrg'))
                    defaultstart[parindx[['superD']][1]] <- transform (start3[['D']], link[['superD']])
            }
            # CL types do not need density
        }
    }
    tmp <- start
    if (is.null(start) | is.list(start)) start <- rep(NA, NP)
    if (any(is.na(start))) {
        start[is.na(start)] <- defaultstart[is.na(start)]
    }
    if (is.list(tmp)) {
        # 2020-10-31 protect against bad start list
        ok <- names(tmp) %in% names(link)
        if (any(!ok)) {
            warning("ignoring parameter(s) in start not in model: ", 
                paste(names(tmp)[!ok], collapse = ', '))
            tmp <- tmp[ok]
        }
        for (i in names(tmp)) {
            if (i == 'b') {
                start[parindx[[i]][1]] <- tmp[[i]]
            }
            else {
                start[parindx[[i]][1]] <- transform (tmp[[i]], link[[i]])
            }
        }
    }
    
    ##########################
    # Fixed beta parameters
    ##########################
    fb <- details$fixedbeta
    if (!is.null(fb)) {
        if (!(length(fb)== NP))
            stop ("invalid fixed beta - require NP-vector")
        if (sum(is.na(fb))==0)
            stop ("cannot fix all beta parameters")
        start <- start[is.na(fb)]  ## drop unwanted betas; remember later to adjust parameter count
    }
    #########################
    # capthist statistics
    #########################
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
    if (secr) {
        usge <- usage(traps(capthist))
        if (is.null(usge) | details$ignoreusage) 
            usge <- matrix(1, nrow=k, ncol= cumss[J+1])  # in terms of secondary sessions
        ## 2017-11-26 collapse data from exclusive detectors; modified 2018-01-17
        CH <- capthist
        if (detector(traps(capthist))[1] == "multi") {
            CH <- abs(capthist)
            CH <- apply(CH,1:2, which.max) *  (apply(CH,1:2, max)>0)
            lost <- apply(capthist,1:2, min)<0
            CH[lost] <- -CH[lost]
            class (CH) <- 'capthist'
            traps(CH) <- traps(capthist)
        }
    }
    else {
        usge <- NULL
    }
    data <- new.env(parent = emptyenv())
    assign("capthist", CH,       pos = data)
    assign("type",     type,     pos = data)
    assign("mask",     mask,     pos = data)
    assign("detectfn", detectfn, pos = data)
    assign("distrib",  distrib,  pos = data)
    assign("binomN",   binomN,   pos = data)
    assign("link",     link,     pos = data)
    assign("fixed",    fixed,    pos = data)
    assign("details",  details,  pos = data)
    assign("ncores",   ncores,   pos = data)
    assign("design",   design,   pos = data)
    assign("design0",  design0,  pos = data)
    assign("parindx",  parindx,  pos = data)
    assign("intervals", intervals, pos = data)
    assign("nc",       nc,       pos = data)
    assign("J",        J,        pos = data)
    assign("cumss",    cumss,    pos = data)
    assign("k",        k,        pos = data)
    assign("m",        m,        pos = data)
    assign("betaw",    betaw,    pos = data)
    assign("fi",       fi,       pos = data)
    assign("li",       li,       pos = data)
    assign("JScounts", JScounts, pos = data)
    assign("marray",   marray,   pos = data)    # 2020-12-08
    assign("usge",     usge,     pos = data)
    assign("moveargsi",       moveargsi,       pos = data)
    assign("movementcode",    movementcode,    pos = data)
    assign("edgecode",        edgecode,        pos = data)
    assign("usermodel",       usermodel,       pos = data)
    assign("kernel",          kernel,          pos = data)
    assign("cellsize",        cellsize,        pos = data)
    assign("mqarray",         mqarray,         pos = data)
    assign("learnedresponse", learnedresponse, pos = data)
    assign("mixturemodel",    mixturemodel,    pos = data)
    # assign("PIA0njx",  PIA0njx,  pos = data)

    #############################
    # Single evaluation option
    #############################
    .openCRstuff$iter <- 0
    if (details$LLonly) {
        if (is.null(start))
            stop ("provide transformed parameter values in 'start'")
        args <- list(beta = start,
                     oneeval = TRUE,
                     data = data)
        LL <- do.call(loglikefn, args)
        names(LL) <- c('logLik', betanames)
        attr(LL, 'parindx') <- parindx
        return(LL)
    }

    #####################
    # Maximize likelihood
    #####################

    ## modified 2017-05-16 to assume most data are in the environment, not needing to be passed
    memo('Maximizing likelihood...', details$trace)
    if (details$trace) {
      message('Eval       Loglik', paste(str_pad(betanames, width = betaw), collapse = ' '))
    }

    if (tolower(method) %in% c('newton-raphson', 'nr')) {
        args <- list (p        = start,
                      f        = loglikefn,
                      data     = data,   # environment(),
                      betaw    = betaw,
                      hessian  = tolower(details$hessian)=='auto',
                      stepmax  = 10)
                      ## cluster  = cluster)
        this.fit <- do.call (nlm, args)
        this.fit$par <- this.fit$estimate     # copy for uniformity
        this.fit$value <- this.fit$minimum    # copy for uniformity
        if (this.fit$code > 2)
            warning ("possible maximization error: nlm returned code ",
                     this.fit$code, ". See ?nlm")
    }
    else if (tolower(method) %in% c('none')) {
        # Hessian-only
        memo ('Computing Hessian with fdHess in nlme', details$trace)
        loglikfn <- function (beta) {
            ## args <- list(beta = beta, data = data, cluster = cluster)
            args <- list(beta = beta, data = data)
            do.call(loglikefn, args)
        }
        grad.Hess <- nlme::fdHess(start, fun = loglikfn, .relStep = 0.001, minAbsPar=0.1)
        this.fit <- list (value = loglikfn(start), par = start,
                          gradient = grad.Hess$gradient,
                          hessian = grad.Hess$Hessian)
    }
    else {

        args <- list(par     = start,
                     fn      = loglikefn,
                     data    = data,
                     hessian = tolower(details$hessian)=='auto',
                     control = details$control,
                     method  = method)
                     # cluster = cluster)

        this.fit <- do.call (optim, args)
        # default method = 'BFGS', control=list(parscale=c(1,0.1,5))
        if (this.fit$convergence != 0)
            warning ("probable maximization error: optim returned convergence ",
                     this.fit$convergence, ". See ?optim")
    }

    this.fit$method <- method         ## remember what method we used...
    covar <- NULL
    if (this.fit$value > 1e9) {     ## failed
        this.fit$beta[] <- NA
        eigH <- NA
    }
    else {

        ############################
        # Variance-covariance matrix
        ############################

        if (tolower(details$hessian)=='fdhess') {
            memo ('Computing Hessian with fdHess in nlme', details$trace)

            loglikfn <- function (beta) {
                args <- list (beta    = beta,
                              parindx    = parindx,
                              env     = data) # environment(),
                              ## cluster = cluster)
                -do.call(loglikefn, args)
            }
            grad.Hess <- nlme::fdHess(this.fit$par, fun = loglikfn, .relStep = 0.001, minAbsPar=0.1)
            this.fit$hessian <- -grad.Hess$Hessian
        }

        hess <- this.fit$hessian
        eigH <- NA
        NP <- length(betanames)
        covar <- matrix(nrow = NP, ncol = NP)
        if (!is.null(hess)) {
            eigH <- eigen(this.fit$hessian)$values
            ## eigH <- eigH/max(eigH)
            eigH <- abs(eigH)/max(abs(eigH))   ## 2020-05-28
            covar <- try(MASS::ginv(hess))
            if (inherits(covar, "try-error")) {
                warning ("could not invert Hessian to compute ",
                         "variance-covariance matrix")
                covar <- matrix(nrow = NP, ncol = NP)
            }
            else if (any(diag(covar)<0)) {
                warning ("variance calculation failed for ",
                         "some beta parameters; confounding likely")
            }
        }
        dimnames(covar) <- list(betanames, betanames)
    }

    desc <- packageDescription("openCR")  ## for version number
    temp <- list (call = cl,
                  capthist = inputcapthist,
                  type = type,
                  model = model,
                  distribution = distribution,
                  mask = mask,
                  detectfn = detectfn,
                  binomN = binomN,
                  movementmodel = movementmodel,
                  edgemethod = edgemethod,
                  usermodel = usermodel,
                  moveargsi = moveargsi,
                  start = start,
                  link = link,
                  fixed = fixed,
                  timecov = timecov,
                  sessioncov = sessioncov,
                  agecov = agecov,
                  dframe = dframe,
                  dframe0 = dframe0,
                  details = details,
                  method = method,
                  ncores = ncores,
                  design = design,
                  design0 = design0,
                  parindx = parindx,
                  intervals = intervals,
                  vars = vars,
                  betanames = betanames,
                  realnames = realnames,
                  sessionlabels = sessnames,
                  fit = this.fit,
                  beta.vcv = covar,
                  eigH = eigH,
                  version = desc$Version,
                  starttime = starttime,
                  proctime = proc.time()[3] - ptm[3]
    )
    if (secr) temp <- c(temp, list(mask=mask))
    attr (temp, 'class') <- 'openCR'

    ###############################################
    ## if (!is.null(cluster)) stopCluster(cluster)
    ###############################################

    memo(paste('Completed in ', round(temp$proctime,2), ' seconds at ',
                      format(Sys.time(), "%H:%M:%S %d %b %Y"),
                      sep=''), details$trace)
    temp

}

################################################################################

