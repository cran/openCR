###############################################################################
## package 'openCR'
## openCR.design.R
## 2011 12 29
## 2012-06-08 robust design (J from intervals, not ncol(capthist))
## 2012-12-20 JSSARET
## 2017-11-19 robust design intervals
## 2017-12-10 expanded for secondary-session effects
## 2018-01-23 timecov added; does sessioncov work?
## 2018-02-23 Bsession added, bsession corrected
## 2018-10-29 CJSp1 argument
## 2018-11-22 bk, Bk added; bsession etc. redefined
## 2018-11-23 tidy up
## 2018-12-21 trap covariates
## 2020-10-19 agecov for recoding age effects
################################################################################

openCR.design <- function (capthist, models, type, naive = FALSE, timecov = NULL, 
    sessioncov = NULL, agecov = NULL, dframe = NULL, 
    contrasts = NULL, initialage = 0, minimumage = 0, maximumage = 1, 
    CJSp1 = FALSE, ...) {
    
## Generate design matrix, reduced parameter array, and parameter index array (PIA)
## for each parameter
## 'capthist' must be of class 'capthist' or 'list'

    findvars <- function (cov, vars, dimcov, scov = FALSE) {
        ## function to add covariates to a design data frame 'dframe'
        ## uses pad1 and insertdim from utility.R
        if (is.null(cov) | (length(cov)==0) | (length(vars)==0)) return()
        else {
            found <- names(cov) %in% vars
            if (is.data.frame(cov) & any(found)) {
                found <- names(cov)[found]
                values <- as.data.frame(cov[,found])
                values <- secr:::stringsAsFactors(values)  ## updated 2020-12-04
                names(values) <- found
                if (length(values)>0) {
                    for (variable in found) {
                        vals <- values[,variable]
                        if (scov) vals <- rep(vals, secondarysessions)
                        dframe[,variable] <<- secr::insertdim (vals, dimcov, dims)
                    }
                }
                vars <<- vars[!(vars %in% found)]
            }
        }
    }

    findvars.covtime <- function (covindices, vars) {
        ## function to add time-specific covariates to a design data frame 'dframe'
        ## covindices should be a list
        dimcov <- c(1,2)   ## animal, secondarysession
        ## covindices is list of numeric or character index vectors, one component per session
        if (length(covindices[[1]]) != J)
            stop ("require one index per primary session")
        covnames <- names(covindices)
        found <- covnames[covnames %in% vars]
        vars <<- vars[!(vars %in% found)]
        for (variable in found) {
            firstcol <- zcov[,covindices[[1]][1]]
            factorlevels <- NULL
            if (is.factor(firstcol)) {
                ## all must have same levels!!
                factorlevels <- levels(firstcol)
            }
            getvals <- function (indices, zcov) {
                notOK <- is.na(zcov[,indices])
                if (any(notOK)) {
                    warning ("covariate missing values set to -1")
                    zcov[,indices][notOK] <- -1
                }
                mat <- as.matrix(zcov[,indices]) ## detectors x occasions
            }
            vals <- getvals(covindices[[variable]], zcov)
            vals <- vals[,primarysessions(intervals)]
            vals <- unlist(vals)  
            if (!is.null(factorlevels))
                vals <- factor(vals, factorlevels)
            dframe[,variable] <<- secr::insertdim (vals, dimcov, dims)
        }
    }
    #--------------------------------------------------------------------------------
    npar     <- length(models)                # real parameters
    parnames <- names(models)                 # c('p','phi') for CJS
    vars     <- unique (unlist(sapply (models, all.vars)))
    trps     <- traps(capthist)
    trapcov  <- covariates(trps)
    used     <- usage(trps)>0
    nmix     <- secr:::get.nmix(models, capthist, NULL)
    zcov     <- covariates(capthist)          # individual covariates
    n        <- nrow(capthist)
    S        <- ncol(capthist)
    K        <- if (grepl('secr', type)) nrow(trps) else 1
    intervals         <- attr(capthist, 'intervals')
    if (is.null(intervals)) intervals <- rep(1, ncol(capthist)-1)
    J                 <- sum(intervals>0) + 1 
    primarysession    <- primarysessions(intervals)
    secondarysessions <- tabulate(primarysession)
    firstofsession    <- match(1:J, primarysession)
    validlevels       <- getvalidlevels (type, parnames, J, CJSp1)
    
    # primary.secondary <- 
    
    #--------------------------------------------------------------------------

    if (!grepl('secr', type) & any(.openCRstuff$traplearnedresponses %in% vars))
        stop ("cannot use detector-specific predictor with non-spatial model")
        
    if (sum(.openCRstuff$learnedresponses %in% vars) > 1)
        stop ("model should not use more than one type of behavioural response")
    
    if (any(.openCRstuff$learnedresponses %in% vars) &
        packageDescription("openCR")$Version<"1.4.0")  # sunset < 1.4.0 for this msg
        warning ("learned response models were re-defined in version 1.3 - check vignette")
    
    #--------------------------------------------------------------------------
    # session covariates   (primary sessions)
    if (!is.null(sessioncov)) {
        scov <- sessioncov  # trick to name a vector 'scov'
        sessioncov <- as.data.frame(scov)
        sessioncov <- secr:::stringsAsFactors(sessioncov)  ## updated 2020-12-04
        if (nrow(sessioncov) != J)
            stop("number of rows in 'sessioncov' should equal ",
                 "number of primary sessions")
    }
    #--------------------------------------------------------------------------
    # time covariates   (secondary sessions)
    if (!is.null(timecov)) {
        tcov <- timecov  # trick to name a vector 'tcov'
        timecov <- as.data.frame(tcov)
        timecov <- secr:::stringsAsFactors(timecov)  ## updated 2020-12-04
        if (nrow(timecov) != S)
            stop("number of rows in 'timecov' should equal ",
                 "number of occasions (secondary sessions")
    }
    #--------------------------------------------------------------------------
    # age covariates   (primary sessions)
    if (!is.null(agecov)) {
        acov <- agecov  # trick to name a vector 'acov'
        agecov <- as.data.frame(acov)
        agecov <- secr:::stringsAsFactors(agecov)
        if (nrow(agecov) != (maximumage - minimumage + 1))
            stop("number of rows in 'agecov' should equal ",
                "number of possible ages")
    }
    
    #--------------------------------------------------------------------------
    dims <- c(n,S,K,nmix)       # virtual dimensions
    dframenrow <- prod(dims)    # number of rows
    autovars <- c(.openCRstuff$learnedresponses, 'session', 't', 'tt', 'Session', 
                  'h2', 'h3', 'age', 'Age', 'Age2')
    #--------------------------------------------------------------------------
    # user-specified dframe
    if (is.null(dframe)) {
        dframe <- data.frame(dummy=rep(1, dframenrow))
        dframevars <- ""
    }
    else {
        if (nrow(dframe) !=  dframenrow )
            stop ("dframe should have ", n*S*K*nmix, " rows ( n*S*K*nmix )")
        dframevars <- names(dframe)
    }
    #--------------------------------------------------------------------------

    dframe$session <- factor(secr::insertdim (rep(1:J, secondarysessions), 2, dims))
    dframe$Session <- secr::insertdim (rep(0:(J-1), secondarysessions), 2, dims)
    dframe$tt <- factor(secr::insertdim (1:S, 2, dims))
    
    ## t as synonym of session
    if ('t' %in% vars) {
        dframe$t <- factor(secr::insertdim (rep(1:J, secondarysessions), 2, dims))
    }
    ## firstage <- as.numeric(grepl('CJS', type))
    # firstage <- 0 # replaced by minimumage
    
    # rearranged 2020-10-19
    if (any(c('age','Age','Age2', names(agecov)) %in% vars)) {
        age <- age.matrix(capthist, initialage, minimumage, maximumage)
    }
    if ('age' %in% vars) {
        dframe$age <- factor(secr::insertdim (factor(age), 1:2, dims))
    } 
    if ('Age' %in% vars) {
        dframe$Age <- secr::insertdim (age, 1:2, dims)
    } 
    if ('Age2' %in% vars) {
        dframe$Age2 <- secr::insertdim (age^2, 1:2, dims)
    } 
    for (i in names(agecov)) {
        if (i %in% vars) {
            agecovi <- agecov[,i][age-minimumage+1]
            dframe[,i] <- secr::insertdim (agecovi, 1:2, dims)
        }
    }
    
    #--------------------------------------------------------------------------

    ## behavioural response fields
    
    makeb <- function (caphist) {      ## global response
        temp0 <- apply(abs(caphist), 1:2, sum)
        ## convert to n x (S+1) CH
        t(apply(temp0, 1, prevcapt))
    }
    ## by primary session
    makebJ <- function (caphist) {      ## global response
        temp0 <- apply(abs(caphist), 1:2, sum)
        temp0 <- t(apply(temp0, 1, tapply, primarysession, sum))
        temp0 <- t(apply(temp0, 1, prevcapt))
        ## convert to n x (S+1) CH
        t(apply(temp0, 1, rep, secondarysessions))
    }
    #--------------------------------------------------------------------------
    if ('bsession' %in% vars) {
        if (naive) dframe$bsession <- rep(FALSE, dframenrow)
        else {
            prevcapt <- function(x) {
                #prevcapt1 <- function(x) c(FALSE, cumsum(x[-S])>0)  BUG 2018-11-22
                prevcapt1 <- function(x) c(FALSE, cumsum(x[-length(x)])>0)
                ## apply within each primary session
                unlist(lapply( split(x, primarysession), prevcapt1))
            }
            temp <- makeb(capthist)
            dframe$bsession <- secr::insertdim (as.vector(temp), c(1,2), dims)
        }
    }
    
    #------------------------------------------------
    if ('Bsession' %in% vars) {
        if (naive) dframe$Bsession <- rep(FALSE, dframenrow)
        else {
            prevcapt <- function(x) {
                prevcapt1 <- function(x) c(FALSE, x[-length(x)]>0)
                ## apply within each primary session
                unlist(lapply( split(x, primarysession), prevcapt1))
            }
            temp <- makeb(capthist)
            dframe$Bsession <- secr::insertdim (as.vector(unlist(temp)), c(1,2), dims)
        }
    }
    
    #------------------------------------------------
    if ('b' %in% vars) {
        if (naive) dframe$b <- rep(FALSE, dframenrow)
        else {
            prevcapt <- function(x) c(FALSE, cumsum(x[-length(x)])>0)
            temp <- makeb(capthist)
            dframe$b <- secr::insertdim (as.vector(temp), c(1,2), dims)
        }
    }
    
    #------------------------------------------------
    if ('B' %in% vars) {
        if (naive) dframe$B <- rep(FALSE, dframenrow)
        else {
            prevcapt <- function(x) {
                prevcapt1 <- function(x) c(FALSE, cumsum(x[-length(x)])>0)
                ## apply within each primary session
                capt1 <- unlist(lapply( split(x, primarysession), prevcapt1))
                pcapt <- unique(primarysession[x>0])
                ## EITHER same primary as ct and later OR next primary
                (capt1 & (primarysession %in% pcapt)) | (primarysession %in% (pcapt+1))
            }
            temp <- makeb(capthist)
            dframe$B <- secr::insertdim (as.vector(temp), c(1,2), dims)
        }
    }
    # old Bsession
    # if ('B' %in% vars) {
    #     if (naive) dframe$B <- rep(FALSE, dframenrow)
    #     else {
    #         prevcapt <- function(x) c(FALSE, x[-J]>0)
    #         temp <- makebJ(capthist)
    #         dframe$B <- secr::insertdim (as.vector(temp), c(1,2), dims)
    #     }
    # }
    #------------------------------------------------
    
        
    #------------------------------------------------
    ## individual trap-specific responses
    
    makebk <- function (caphist) {     
        # condition added 2016-10-01
        if (nrow(caphist)==0) 
            array(dim = c(0,S,K))
        else {
            temp <- apply(abs(caphist), c(1,3), prevcapt)
            aperm(temp, c(2,1,3))
        }
    }
    
    #------------------------------------------------
    if ('bksession' %in% vars) {
        if (naive) dframe$bksession <- rep(FALSE, dframenrow)
        else {
            prevcapt <- function(x) {
                prevcapt1 <- function(x) c(FALSE, cumsum(x[-length(x)])>0)
                ## apply within each primary session
                unlist(lapply( split(x, primarysession), prevcapt1))
            }
            temp <- makebk(capthist) 
            dframe$bksession <- secr::insertdim(temp, 1:3, dims)  
        }
    }
    
    #------------------------------------------------
    if ('Bksession' %in% vars) {
        if (naive) dframe$Bksession <- rep(FALSE, dframenrow)
        else {
            prevcapt <- function(x) c(FALSE, x[-S]>0)
            temp <- makebk(capthist)  # one session
            dframe$Bksession <- secr::insertdim(temp, 1:3, dims)
        }
    }
    
    #------------------------------------------------
    if ('bk' %in% vars) {
        if (naive) dframe$bk <- rep(FALSE, dframenrow)
        else {
            prevcapt <- function(x) c(FALSE, cumsum(x[-length(x)])>0)
            temp <- makebk(capthist) 
            dframe$bk <- secr::insertdim(temp, 1:3, dims)  
        }
    }
    
    #------------------------------------------------

    ## 2018-11-22
    ## Detector response (k) is an approximation because
    ## "naive" state refers to animals not detectors --
    ## undocumented for now (k matches secr)
    
    makek <- function (caphist) {      ## trap responds to capture of any animal
        temp <- apply(abs(caphist), c(2,3), sum) # occasion x trap
        apply(temp, 2, prevcapt)
    }

    #------------------------------------------------
    if ('k' %in% vars) {
        if (naive) dframe$k <- rep(FALSE, dframenrow)
        else {
            prevcapt <- function(x) c(FALSE, cumsum(x[-length(x)])>0)
            temp <- makek(capthist)
            dframe$k <- secr::insertdim(temp, 2:3, dims)
        }
    }

    #------------------------------------------------
    if ('ksession' %in% vars) {
        if (naive) dframe$ksession <- rep(FALSE, dframenrow)
        else {
            prevcapt <- function(x) {
                prevcapt1 <- function(x) c(FALSE, cumsum(x[-length(x)])>0)
                ## apply within each primary session
                unlist(lapply( split(x, primarysession), prevcapt1))
            }
            temp <- makek(capthist)
            dframe$ksession <- secr::insertdim(temp, 2:3, dims)
        }
    }

    #------------------------------------------------
    if ('Ksession' %in% vars) {
        if (naive) dframe$Ksession <- rep(FALSE, dframenrow)
        else {
            prevcapt <- function(x) {
                prevcapt1 <- function(x) c(FALSE, x[-length(x)]>0)
                ## apply within each primary session
                unlist(lapply( split(x, primarysession), prevcapt1))
            }
            temp <- makek(capthist)
            dframe$Ksession <- secr::insertdim(temp, 2:3, dims)
        }
    }

    #------------------------------------------------
    ## h2 or h3
    if (nmix > 1) {
        mixture <- paste('h',nmix,sep='')
        dframe[,mixture] <- secr::insertdim(factor(1:nmix), 4, dims)
    }

    #--------------------------------------------------------------------------
    ## all autovars should have now been dealt with
    vars <- vars[!(vars %in% c(autovars, names(agecov), dframevars))]

    #--------------------------------------------------------------------------
    # add zcov, sessioncov etc.

    findvars (zcov, vars, 1)
    findvars (timecov, vars, 2)
    findvars (trapcov, vars, 3)  ## 2018-12-21
    findvars (sessioncov, vars, 2, scov = TRUE)  ## expands correctly 2018-05-07
    findvars (agecov, vars, 1, scov = TRUE)      ## new 2020-10-19
    tvc <- timevaryingcov(capthist)
    if (!is.null(tvc) & (length(vars)>0)) {
        findvars.covtime (tvc, vars)
    }

    if (length(vars)>0) {
        if (!is.null(zcov)) {
            if (is.data.frame(zcov))
                znames <- names(zcov)
            else
                znames <- unlist(lapply(zcov, names))
        }
        stop ("covariate(s) ", paste(vars,collapse=","), " not found")
    }

    make.designmatrix <- function (formula, prefix, ...) {
     # combine formula and dframe to generate design matrix
        if (is.null(formula)) {
            list (model = NULL, index = rep(1,dframenrow))
        }
        else {
            # adjust for unidentifiable parameters
            dframe <- adjustlevels(prefix, dframe, validlevels)
            tempmat <- model.matrix(formula, data = dframe, contrasts.arg = contrasts, ...)
            ## drop pmix beta0 column from design matrix
            if (prefix=='pmix') tempmat <- tempmat[,-1,drop=FALSE]
            ## temp <- secr::make.lookup (tempmat)   # retain unique rows
            temp <- makelookupcpp (tempmat)   # retain unique rows   ## 2018-11-06
            list (model=temp$lookup, index=temp$index)
        }
    }
    dframe[is.na(dframe)] <- 0
    # list with one component per real parameter
    # each of these is a list with components 'model' and 'index'
    designMatrices <- sapply (1:length(models), simplify=FALSE,
        function (x) make.designmatrix(models[[x]], names(models[x])))
    names(designMatrices) <- names(models)
    
    ## dim(indices) = c(n*S*K*nmix, npar)
    indices <- sapply (designMatrices, function(x) x$index)
    indices <- matrix(unlist(indices), ncol = npar)

    # retain just the 'model' components of 'designMatrices'
    designMatrices <- lapply (designMatrices, function(x)x$model )

    # prefix column names in 'designMatrices' with parameter name
    for (i in 1:npar) {
        colnames(designMatrices[[i]]) <- paste (parnames[i], '.',
            colnames(designMatrices[[i]]), sep='')
    }

    # repackage indices to define unique combinations of parameters
    ##indices2 <- secr::make.lookup(indices)
    indices2 <- makelookupcpp(indices)
    
    #--------------------------------------------------------------------
    # PIA = Parameter Index Array
    #       index to row of parameterTable for a given n,s,nmix
    # dim(parameterTable) = c(uniqueparcomb, npar)
    #       index to row of designMatrix for each real parameter
    #--------------------------------------------------------------------
    PIA <- array(indices2$index, dim = dims)
    PIAJ <- njx(PIA, primarysession)
    
    parameterTable <- indices2$lookup
    colnames(parameterTable) <- parnames

    #--------------------------------------------------------------------
    # Zero the index of trap+time pairs that were 'not set'
    # the external C code checks for this and sets p(detection) to zero
    #--------------------------------------------------------------------

    if (grepl('secr', type)) {
        if ((!is.null(unlist(used))) & (length(used)>0)) {
            allused <- unlist(used)
            if (!is.null(allused)) {
                if (any(!allused)) {
                    PIA[ , , ,] <- PIA[ , , ,] * rep(rep(t(used),rep(n,S*K)),nmix)
                }
            }
        }
    }
    individual <- individualcovariates(PIA)
    #--------------------------------------------------------------------
    list(designMatrices = designMatrices, parameterTable = parameterTable, PIA = PIA,
         PIAJ = PIAJ, validlevels = validlevels, individual = individual)
}
############################################################################################

## PIA has dim = c(nc, ss, kk, xx)
## primarysession gives primary session (j) for each secondary session (s)
# returns list with components 
#    lookup - matrix of unique rows, each row ss * kk long, containing indices as in PIA to 
#    the rows of realparval0
#    index  vector of length nc * jj * xx; values are indices to rows of lookup

# this could also be used in secr to speed up esa

# BUT NOT WORKING 2018-02-13
njxlookup <- function (PIA, primarysession) {
    dims <- lapply(dim(PIA), seq, from=1)
    names(dims) <- c('n','s','k','x')
    df <- data.frame(pia = as.numeric(PIA), do.call(expand.grid, dims))
    ss <- dim(PIA)[2]
    kk <- dim(PIA)[3]
    names(df) <- c('pia', 'n','s','k','x')
    df$j <- formatC(primarysession[df$s], width=3, flag="0")
    df$n <- formatC(df$n, width=4, flag="0")
    # splitter <- apply(df[,c('n','j','x')], 1, paste, collapse='.')
    splitter <- apply(df[,c('x','j','n')], 1, paste, collapse='.')
    splitdf <- split(df[,c('pia','s')], splitter)
    fixpiask <- function (x) {
        pia <- matrix(0,ss,kk)
        pia[x$s,] <- x$pia
        as.numeric(pia)
    }
    piask <- lapply(splitdf, fixpiask)
    njxIA <- do.call(rbind, piask)
    ## lookup <- secr::make.lookup(njxIA)
    lookup <- makelookupcpp(njxIA)
    lookup
}

# direct index to session PIA in njx array 2018-02-10
# slow - must be a better way 2018-04-11
njx <- function (PIA, primarysession) {
    n <- dim(PIA)[1]
    J <- max(primarysession)
    xx <- dim(PIA)[4]
    # why all this complexity? 2018-04-11
    # dims <- lapply(dim(PIA), seq, from=1)
    # names(dims) <- c('n','s','k','x')
    # df <- data.frame(pia = as.numeric(PIA), do.call(expand.grid, dims))
    # names(df) <- c('pia', 'n','s','k','x')
    # df$j <- formatC(primarysession[df$s], width=3, flag="0")
    # df$n <- formatC(df$n, width=4, flag="0")
    # splitter <- apply(df[,c('x', 'j', 'n')], 1, paste, collapse='.')  # deliberately reverse order
    # splitdf <- split(df[,c('pia','s')], splitter)
    # piaval <- sapply(splitdf, function(x) x$pia[which.min(!(x$pia != 0))])
    # array(piaval, dim = c(n, J, xx))

    # instead 2018-04-11   
    s1 <- match(1:J, primarysession)
    piaj <- PIA[,s1,1,]
    array(piaj, dim = c(n,J,xx))
}

