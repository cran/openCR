############################################################################################
## package 'openCR'
## openCR.make.newdata.R
## 2011 12 09
## Create (neutral) design data suitable for 'predict'
## 2015-02-06 reconciled this current version with forked 1.2.0:
## 2017-12 revamped
## 2018-04-12 allow single session
## 2018-11-22 new learned responses
############################################################################################

openCR.make.newdata <- function (object, all.levels = FALSE) {
    # 'Session', 't' are handled separately at end
    autovars <- c(.openCRstuff$learnedresponses, 'session','h2','h3')
    capthist <- object$capthist
    mask <- object$mask
    vars <- object$vars
    timecov <- object$timecov
    sessioncov <- object$sessioncov
    dframe <- object$dframe
    nmix <- object$details$nmix
    if(is.null(nmix)) nmix <- 1
    J <- sum(object$intervals>0) + 1
    sessions <- 1:J
    ngrp <- 1  # no groups for now

    #############################################################
    findvars <- function (basevars, cov) {
        ## function to add covariates to a list
        if (is.null(cov) | (length(cov)==0) | (length(vars)==0))
            return(basevars)
        else {
            found <- ''
            for (v in vars) {
                if (v %in% names(cov)) {
                    vals <- cov[,v]
                    if (is.character(vals)) vals <- factor(vals)
                    basevars[[v]] <- if (is.factor(vals))
                        factor(levels(vals), levels = levels(vals))
                    else
                        unique(vals)

                    found <- c(found, v)
                }
            }
            vars <<- vars[!(vars %in% found)]
            return(basevars)
        }
    }
    basevars <- list(session = factor(sessions))

    ## if (ngrp>1) basevars$g <- factor(grouplevels)
    mixvar <- 'h2'   ## stop gap 2018-01-22
    if (nmix>1) basevars[mixvar] <- list(as.character(1:nmix))
    for (v in vars) {
        if (v=='T')  basevars$T <- 0
        
        # superceded 2018-11-15            
        # if (v=='b')  basevars$b <- factor(0:1)
        # if (v=='k')  basevars$k <- factor(0:1)
        # if (v=='bk') basevars$bk <- factor(0:1)
        # if (v=='bsession')  basevars$bsession <- factor(0:1)
        # if (v=='ksession')  basevars$ksession <- factor(0:1)
        # if (v=='bksession') basevars$bksession <- factor(0:1)
        # if (v=='Bsession')  basevars$Bsession <- factor(0:1)
        # if (v=='Ksession')  basevars$Ksession <- factor(0:1)
        # if (v=='Bksession')  basevars$Bksession <- factor(0:1)
        for (i in .openCRstuff$learnedresponses) {
            if (v == i) basevars[,i] <- factor(0,1)
        }
            
        # firstage <- as.numeric(grepl('CJS', object$type))
        firstage <- 0
        if (v=='age')  basevars$age <- factor(firstage:object$details$maximumage)
        if (v=='Age')  basevars$Age <- firstage:object$details$maximumage
    }
    ## all autovars should now have been dealt with
    vars <- vars[!vars %in% autovars]
    if (ngrp==1)
        basevars <- findvars (basevars, covariates(capthist))  ## individual covariates
    
    if (!is.null(timecov)) {
        tcov <- timecov
        timecov <- as.data.frame(tcov)
    }
    if (!is.null(sessioncov)) {
        scov <- sessioncov
        sessioncov <- as.data.frame(scov)
    }
    basevars <- findvars (basevars, timecov)
    basevars <- findvars (basevars, covariates(traps(capthist)))
    if (!is.null(mask))
        basevars <- findvars (basevars, covariates(mask))
    if (!is.null(dframe))
        basevars <- findvars (basevars, dframe)
    
    ## revert to first level (original default)
    if (length(v)>0) {
        if (!all.levels & !(v %in% c('session', 'g', 'h2','h3'))) {
            basevars[[v]] <- basevars[[v]][1]
        }
    }
    newdata <- expand.grid(basevars)
    if ('Session' %in% vars)
        ## based on sequence in capthist
        newdata$Session <- as.numeric(newdata$session) - 1
    if (!is.null(sessioncov) & any(names(sessioncov) %in% vars))
        newdata <- cbind(newdata, sessioncov[newdata$session,,drop=FALSE])
    if ('t' %in% vars) { ## synonym 
        newdata$t <- newdata$session
    }
    newdata
}
############################################################################################

