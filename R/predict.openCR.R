###############################################################################
## package 'openCR'
## predict.openCR.R
## moved from methods 2017-12-25
## 2019-04-11 ignores contrasts arg of lpredictor?
## 2019-04-13 failed with single parameter when others fixed
###############################################################################

predict.openCRlist <- function (object, newdata = NULL, se.fit = TRUE,
                                alpha = 0.05, savenew = FALSE, ...) {
    lapply(object, predict, newdata, se.fit, alpha, savenew, ...)
}

predict.openCR <- function (object, newdata = NULL, se.fit = TRUE, alpha = 0.05,
                            savenew = FALSE, ...) {
    if (is.null(object$fit)) {
        warning ("empty (NULL) object")
        return(NULL)
    }
    if (is.null(newdata)) newdata <- openCR.make.newdata (object, ...)
    beta <- complete.beta(object)
    beta.vcv <- complete.beta.vcv(object)
    if (is.null(beta.vcv)) se.fit <- FALSE
    getfield <- function (x) {
        lpredictor (newdata = newdata, model = object$model[[x]],
                    indx = object$parindx[[x]], beta = beta, field = x,
                    beta.vcv = beta.vcv, validlevels = object$design$validlevels)
    }
    predict <- sapply (object$realnames, getfield, simplify = FALSE)

    if(se.fit) {
        realvcv <- vcov(object, realnames = names(predict), newdata = newdata, byrow = TRUE)
        nreal <- length(predict)
        realvcv <- array(unlist(realvcv), dim=c(nreal, nreal, nrow(newdata)))
        realSE <- apply(realvcv,3, function(x) suppressWarnings(sqrt(diag(x))))
        # bug fix 2019-04-13
        # if single parameter, apply() returns vector instead of matrix
        if (is.null(dim(realSE))) realSE <- matrix(realSE, nrow = 1)
        ####################
        rownames(realSE) <- names(predict)
    }
    if (!is.null(predict$pmix)) {
        nmix <- object$details$nmix
        # assuming mixture is always last dimension...
        temp <- matrix(predict$pmix$estimate, ncol = nmix)
        #        temp2 <- apply(temp, 1, function(est) logit(mlogit.untransform(est, 1:nmix)))
        temp2 <- apply(temp, 1, function(est) logit(invmlogit(est)))
        predict$pmix$estimate <- as.numeric(t(temp2))
        ######################
        predict$pmix$se <- NA    ## uncertain
    }
    nsess <- sum(object$intervals>0) + 1
    z <- abs(qnorm(1-alpha/2))
    for (i in names(predict)) {
        nc <- ncol(predict[[i]])
        if (i %in% c('superN','superD')) {
            predict[[i]] <- predict[[i]][1,]
            predict[[i]]$session <- NULL
            predict[[i]]$h2 <- NULL
            predict[[i]]$h3 <- NULL
        }
        else if (i %in% c('tau')) {
            M <- object$details$M
            predict[[i]] <- predict[[i]][c(1:M, M),]
            predict[[i]][M+1,] <- rep(NA,3)
            predict[[i]]$session <- 1:(object$details$M+1)
            rownames(predict[[i]]) <- 1:(object$details$M+1)
        }
        else if (i %in% c('pmix')) {
            predict[[i]] <- predict[[i]][predict[[i]]$session==1,,drop=FALSE]
            predict[[i]]$session <- NULL
        }
        else {
            # replace factor session with integer codes
            if (!is.null(predict[[i]]$session)) predict[[i]]$session <- 1:nsess   ## repeats as necessary
            predict[[i]]$t <- NULL
            # drop spurious mixtures
            if (('h2' %in% names(predict[[i]])) & (!('h2' %in% all.vars(object$model[[i]])))) {
                predict[[i]] <- predict[[i]][predict[[i]]$h2 == 1,]
                predict[[i]]$h2 <- NULL
            }
            if (('h3' %in% names(predict[[i]])) & (!('h3' %in% all.vars(object$model[[i]])))) {
                predict[[i]] <- predict[[i]][predict[[i]]$h3 == 1,]
                predict[[i]]$h3 <- NULL
            }
        }
        lpred  <- predict[[i]][,'estimate']
        if (i %in% c('b')) {
            tempmat <- matrix(lpred, nrow=nsess)
            predict[[i]]$estimate <- as.numeric(apply(tempmat,2,invmlogit))
        }
        else
            if (i %in% c('tau')) {
                tempmat <- matrix(lpred, nrow=object$details$M+1)
                predict[[i]]$estimate <- as.numeric(apply(tempmat,2,invmlogit))
            }
        else
            if (i %in% c('pmix')) {
                predict[[i]]$estimate <- invmlogit(lpred)
            }
        else
            predict[[i]]$estimate <- untransform(lpred, object$link[[i]])
        if (se.fit) {
            selpred <- predict[[i]][,'se']
            if (i %in% c('b', 'tau')) {
                predict[[i]]$SE.estimate <- rep(NA,length(selpred))
                ## doubtful
                ## warning("doubtful confidence limits for 'b' or 'tau' model in predict() - to be tested")
                tempmat <- matrix(lpred-z*selpred, nrow=nsess)
                predict[[i]]$lcl <- as.numeric(apply(tempmat,2,invmlogit))
                tempmat <- matrix(lpred+z*selpred, nrow=nsess)
                predict[[i]]$ucl <- as.numeric(apply(tempmat,2,invmlogit))
            }
            else {
                ## 2018-02-25 adapted to use vcov for real parameters from vcov.openCR 
                ##            rather than se.untransform (lpred, selpred, object$link[[i]])
                se <- realSE[i,]
                se[is.na(predict[[i]]$estimate)] <- NA
                predict[[i]]$SE.estimate <- se[1:nrow(predict[[i]])]   ## only 1 for superN, superD
                predict[[i]]$lcl <- untransform(lpred-z*selpred, object$link[[i]])
                predict[[i]]$ucl <- untransform(lpred+z*selpred, object$link[[i]])
            }
            predict[[i]][is.na(predict[[i]])] <- NA
        }
        else {
            predict[[i]]$SE.estimate <- rep(NA, nrow(predict[[i]])) 
            predict[[i]]$lcl <- rep(NA, nrow(predict[[i]])) 
            predict[[i]]$ucl <- rep(NA, nrow(predict[[i]])) 
        }
        predict[[i]]$se <- NULL
  
        freq <- covariates(object$capthist)$freq
        if (i == 'superN') {
            if (object$distribution == 'binomial') {
                ncf <- sum(freq)
                predict[[i]][,c('estimate','lcl','ucl')] <- predict[[i]][,c('estimate','lcl','ucl')] + ncf
            }
        }
        if (i == 'N1') {
            # not ready
            # n1 <- sum(freq[abs(object$capthist[,1] > 0)])
            # predict[[i]][,c('estimate','lcl','ucl')] <- predict[[i]][,c('estimate','lcl','ucl')]+ n1
        }
        sessnames <- object$sessionlabels
        if (!is.null(sessnames) ) {
            if (length(sessnames) == nrow(predict[[i]])) {
                row.names(predict[[i]]) <- sessnames
            }
            else if (nrow(predict[[i]])==1) {
                row.names(predict[[i]]) <- sessnames[length(sessnames) %/% 2 + 1]
            }
        }
    }
    if (savenew) attr(predict, 'newdata') <- newdata
    predict
}
############################################################################################

