###############################################################################
# openCR
## 2018-04-18 openCR 1.2.0
## 2018-04-20 m.array removed to m.array.R
## 2018-04-20 JS.counts removed to JS.counts.R
## 2018-04-20 posterior.allocation removed to posterior.R
## 2018-11-22 n.unique.rows drop = FALSE 
## 2018-11-22 learnedresponses global vector
## 2020-11-02 PLB aliases
## 2020-12-07 CJSmte type experimental - to be added?
## 2020-12-08 primaryonly tweaked to give 3-D output (single NULL traps)
## 2020-12-12 mqsetup moved from prwisecr.R
###############################################################################

typecode <- function (type) {
    type <- switch(type, 
        CJS = 1, CJSmte = 5,
        JSSAf = 4, JSSAl = 3, JSSAb = 2, JSSAg = 22, JSSAk = 28, 
        JSSAfCL = 15, JSSAlCL = 16, JSSAbCL = 17, JSSAgCL = 23, JSSAkCL = 29,
        JSSAB = 18, JSSAN = 19,
        
        CJSsecr = 6, 
        JSSAsecrfCL = 9, JSSAsecrlCL = 10, JSSAsecrbCL = 11, JSSAsecrgCL = 25, 
        JSSAsecrf = 7, JSSAsecrl = 12, JSSAsecrb = 13, JSSAsecrg = 24, 
        JSSAsecrB = 14, JSSAsecrD = 8,
        
        Pradel = 20, JSSARET = 21, Pradelg = 26, JSSAfgCL = 27,
        secrCL = 30, secrD = 31, 
        
        # added aliases 2020-11-01
        PLBf = 15, PLBl = 16, PLBb = 17, PLBg = 23, PLBk = 29,
        PLBsecrf = 9, PLBsecrl = 10, PLBsecrb = 11, PLBsecrg = 25, 
        
        -1)
    type
}

movecode <- function (movementmodel) {
    switch (movementmodel, static = 0, uncorrelated = 1, normal = 2, 
        exponential = 3, user = 4, t2D = 5, uniform = 6)
}

edgemethodcode <- function (edgemethod) {
    if (is.null(edgemethod)) edgemethod <- 'none' 
    if (!edgemethod %in% c('none','wrap','truncate')) {
        stop ("unrecognised edge method - should be 'none','wrap' or 'truncate'")
    }
    switch (edgemethod, none = 0, wrap = 1, truncate = 2, 0)
}

.openCRstuff <- new.env()
#.openCRstuff$packageType <- ' pre-release'
.openCRstuff$packageType <- ''
.openCRstuff$iter <- 0
.openCRstuff$suspendedtypes <- c('JSSAfgCL', 'JSSAfg')

# location of 'sigma' in parameter vector
# zero-based, so take care to add 1 when used as index in R 2019-04-21
# moveargs (move.a etc) follow sigma
.openCRstuff$sigmai <- rep(3, 50)   ## default 3
.openCRstuff$sigmai[c(6,15)] <- 2
.openCRstuff$sigmai[c(7,12,13,24,31)] <- 4

.openCRstuff$learnedresponses <- c('b','bk', 'B', 'bsession',
                                   'k', 'ksession', 'Ksession',
                                   'bksession', 'Bsession','Bksession')
.openCRstuff$traplearnedresponses <- c('bk', 'Bk', 'k', 'ksession', 'Ksession',
                                   'bksession', 'Bksession')

replacedefaults <- function (default, user) replace(default, names(user), user)

discreteN <- function (n, N) {
    tN <- trunc(N)
    if (N != tN) tN + sample (x = c(1,0), prob = c(N-tN, 1-(N-tN)),
        replace = T, size = n)
    else rep(tN,n)
}

memo <- function (text, trace) {
    ## could use message(text), but does not immediately flush console
    if (trace) { cat (text, '\n')
    flush.console() }
}

pad1 <- function (x, n) {
## pad x to length n with dummy (first value)
    if (is.factor(x)) {
        xc <- as.character(x)
        xNA <- c(xc, rep(xc[1], n-length(xc)))
        out <- factor(xNA, levels=levels(x))
    }
    else out <- c(x, rep(x[1], n-length(x)))
    out
}

padarray <- function (x, dims) {
    temp <- array(dim=dims)
    dimx <- dim(x)
    if (length(dimx)<2 | length(dimx)>3)
        stop ("invalid array")
    if (length(dimx)>2) temp[1:dimx[1], 1:dimx[2], 1:dimx[3]] <- x
    else temp[1:dimx[1], 1:dimx[2]] <- x
    temp
}

## regularize a list of formulae
## added 2009 08 05
LHS <- function (form) {
    trms <- as.character (form)
    if (length(trms)==2) '' else trms[2]
}
RHS <- function (form) {
    trms <- as.character (form)
    if (length(trms)==3) as.formula(paste(trms[c(1,3)])) else form
}
stdform <- function (flist) {
    lhs <- sapply(flist, LHS)
    temp <- lapply(flist, RHS)
    if (is.null(names(flist))) names(temp) <- lhs
    else names(temp) <- ifelse(names(flist) == '', lhs, names(flist))
    temp
}

## Start of miscellaneous functions

sine     <- function (x) asin (x*2-1)
invsine  <- function (y) (sin(y)+1) / 2
odds     <- function (x) x / (1-x)
invodds  <- function (y) y / (1+y)
############################################################################################

lnbinomial <- function (x,size,prob) {
  lgamma (size+1) - lgamma (size-x+1) - lgamma (x+1) +
      x * log(prob) + (size-x) * log (1-prob)
}

############################################################################################

var.in.model <- function(v,m) v %in% unlist(lapply(m, all.vars))
############################################################################################

get.nmix <- function (model) {
    model$D <- NULL  ## ignore density model
    nmix <- 1
    if (any(var.in.model('h2', model))) {
        nmix <- 2
        if (any(var.in.model('h3', model)))
            stop ("do not combine h2 and h3")
    }
    if (any(var.in.model('h3', model)))
        nmix <- 3
##    if ((nmix == 3))
##        stop ('3-part mixtures not yet implemented')
    nmix
}

add.cl <- function (df, alpha, loginterval, lowerbound = 0) {

## add lognormal or standard Wald interval to dataframe with columns
## 'estimate' and 'SE.estimate'
## lowerbound added 2011-07-15
    z <- abs(qnorm(1-alpha/2))
    if (loginterval) {
        delta <- df$estimate - lowerbound
        df$lcl <- delta / exp(z * sqrt(log(1 + (df$SE.estimate /
                        delta)^2))) + lowerbound
        df$ucl <- delta * exp(z * sqrt(log(1 + (df$SE.estimate /
                        delta)^2))) + lowerbound
    }
    else {
        df$lcl <- pmax(lowerbound, df$estimate - z * df$SE.estimate)
        df$ucl <- df$estimate + z * df$SE.estimate
    }
    df
}

###############################################################################

model.string <- function (model, userDfn) {
    temp <- paste (names(model), as.character(model), collapse=' ', sep='')
    temp
}
fixed.string <- function (fixed) {
    if (is.null(fixed) | length(fixed)==0) 'none'
    else paste (names(fixed), as.character(fixed), collapse=', ', sep=' = ')
}
############################################################################################


## lifted from secr score.test

mapbeta <- function (parindx0, parindx1, beta0, betaindex)

## Extend beta vector from simple model (beta0) to a more complex (i.e. general)
## model, inserting neutral values (zero) as required.
## For each real parameter, a 1:1 match is assumed between
## beta values until all beta values from the simpler model are
## used up. THIS ASSUMPTION MAY NOT BE JUSTIFIED.
## betaindex is a user-controlled alternative.

{
    ## list of zeroed vectors, one per real parameter
    beta1 <- lapply(parindx1, function (x) {x[]<-0; x})
    if (!is.null(betaindex)) {
        beta1 <- unlist(beta1)
        if (sum(betaindex>0) != length(beta0))
            stop ("invalid 'betaindex'")
        beta1[betaindex] <- beta0
        beta1
    }
    else {
        indx <- lapply(parindx0, function(x) x-x[1]+1)
        # if (!all(names(beta1) %in% names(beta0)))
        #     stop ("incompatible parameters in old model used for start")
        # for (j in 1:length(beta1))
        #     beta1[[j]][indx[[j]]] <- beta0[parindx0[[j]]]
        for (j in names(beta1)) {
            if (j %in% names(parindx0)) {
                newi <- indx[[j]]  ## 2017-11-21
                if (length(parindx1[[j]]) < length(newi))
                    newi <- newi[1:length(parindx1[[j]])]   # clip tail
                beta1[[j]][newi] <- beta0[parindx0[[j]]][newi]
            }
            else
                beta1[[j]][] <- NA
        }
        unlist(beta1)
    }
}
############################################################################################

first <- function(y) {
    indf <- min(which(y>0))
    y[] <- 0
    y[indf] <- 1
    y
}
############################################################################################

last <- function(y) {
    indl <- max(which(y>0))
    y[] <- 0
    y[indl] <- 1
    y
}
############################################################################################

primarysessions <- function(intervals) {
    primarysession <- cumsum(c(0,intervals))
    match(primarysession, unique(primarysession))
}

secondarysessions <- function(intervals) {
    primary <- primarysessions(intervals)
    unname(unlist(sapply(table(primary), seq_len)))  
}

primaryonly <- function(object) {
    intervals <- intervals(object)
    object <- unsqueeze(object)    ## 2018-02-06
    if (!is.null(intervals)) {
        if (any(intervals == 0)) {
            primarysession <- primarysessions(intervals)
            # reduce not ready for not for non-spatial
            # so use more elaborate alternative in primaryonly: 
            if (!is.null(traps(object))) {
                newocc <- split(1:ncol(object), primarysession)
                object <- suppressWarnings(reduce(object, newoccasions = newocc, verify = FALSE))
            }
            else {
                lost <- which(apply(object,1,min, drop = FALSE)<0)
                twoD <- apply(abs(object), 1:2, sum)
                twoD <- t(apply(twoD, 1, function(x) tapply(x,primarysession,max)))  # in terms of primary sessions
                li <- apply(twoD, 1, function(x) max(which(x>0)))
                twoD[cbind(lost, li[lost])] <- -1
                dim(twoD) <- c(dim(twoD),1)   ## 2020-12-08
                class(twoD) <- 'capthist'
                object <- twoD
            }
        }
    }
    object
}

# expandfreq <- function (CH, freq) {
#     f <- rep(1:nrow(CH), freq)    
#     out <- CH[f,,]
#     out
# }

##### constrain to 0 < sum(p) <= 1.0
invmlogit <- function (x, fillin = TRUE) {
    if (!any(is.na(x))) x[1] <- NA
    x2 <- sapply(x,exp) / (1 + sum(exp(x), na.rm=T))
    if (fillin) {
        x2[is.na(x)] <- 1-sum(x2, na.rm=T)
    }
    x2
}
############################################################################################

transform <- function (x, link) {
    switch (link,
            identity = x,
            log = log(x),
            neglog = log(-x),
            logit = logit(x),
            loglog = -log(-log(x)),
            # mlogit = x,
            mlogit = logit(x),   ## 2019-04-20
            odds = odds(x),
            sin = sine(x)
    )
}
############################################################################################

untransform <- function (beta, link) {
    switch (link,
            identity = beta,
            log = exp(beta),
            neglog = -exp(beta),
            logit = invlogit(beta),
            loglog = exp(-exp(-beta)),
            mlogit = beta,   ## delay
            odds = invodds(beta),
            sin = invsine(beta))
}
############################################################################################

se.untransform <- function (beta, sebeta, link) {
    switch (link,
            identity = sebeta,
            log = exp(beta) * sqrt(exp(sebeta^2)-1),
            neglog = exp(beta) * sqrt(exp(sebeta^2)-1),
            logit = invlogit(beta) * (1-invlogit(beta)) * sebeta,
            loglog = -exp(-exp(-beta)) * sebeta * -exp(-beta),
            mlogit = NA,      ## don't know how
            sin = NA)         ####!!!!
}
############################################################################################

# vector version of transform()
Xtransform <- function (real, linkfn, varnames) {
    out <- real
    for (i in 1:length(real)) {
        vn <- varnames[i]
        out[i] <- switch (linkfn[[vn]],
                          identity = real[i],
                          log = log(real[i]),
                          neglog = log(-real[i]),
                          logit = logit(real[i]),
                          loglog = -log(-log(real[i])),
                          mlogit = real[i],      ## delay
                          odds = odds(real[i]),
                          sin = sine(real[i]))
    }
    out
}
############################################################################################

se.Xtransform <- function (real, sereal, linkfn, varnames) {
    out <- real
    for (i in 1:length(real)) {
        vn <- varnames[i]
        out[i] <- switch (linkfn[[vn]],
                          identity = sereal[i],
                          log = log((sereal[i]/real[i])^2 + 1)^0.5,
                          neglog = log((sereal[i]/-real[i])^2 + 1)^0.5,
                          logit = sereal[i] / real[i] / (1 - real[i]),
                          loglog = -real[i] * sereal[i] * log(real[i]),
                          mlogit = NA,      ## don't know how
                          sin = NA)
    }
    out
}
############################################################################################

# vector version of untransform()
Xuntransform <- function (beta, linkfn, varnames) {
    out <- beta
    for (i in 1:length(beta)) {
        vn <- varnames[i]
        out[i] <- switch (linkfn[[vn]],
                          identity = beta[i],
                          log = exp(beta[i]),
                          neglog = -exp(beta[i]),
                          logit = invlogit(beta[i]),
                          loglog = exp(-exp(-beta[i])),
                          mlogit = beta[i],   # delay
                          odds = invodds(beta[i]),
                          sin = invsine(beta[i]))
    }
    out
}
############################################################################################

se.Xuntransform <- function (beta, sebeta, linkfn, varnames)
    # Approximate translation of SE to untransformed scale
    # Delta method cf Lebreton et al 1992 p 77
{
    out <- beta
    if (length(beta)!=length(sebeta))
        stop ("'beta' and 'sebeta' do not match")
    if (!all(varnames %in% names(linkfn)))
        stop ("'linkfn' component missing for at least one real variable")
    for (i in 1:length(beta)) {
        vn <- varnames[i]
        out[i] <- switch (linkfn[[vn]],
                          identity = sebeta[i],
                          log = exp(beta[i]) * sqrt(exp(sebeta[i]^2)-1),
                          neglog = exp(beta[i]) * sqrt(exp(sebeta[i]^2)-1),
                          logit = invlogit(beta[i]) * (1-invlogit(beta[i])) * sebeta[i],
                          loglog = NA, 
                          mlogit = NA,    ## don't know how
                          sin = NA)         ####!!!!
    }
    out
}
############################################################################################

adjustlevels <- function (field, dframe, validlevels ) {
    OK <- validlevels[field,]
    if ('session' %in% names(dframe)) {
        levels(dframe$session)[!OK] <- levels(dframe$session)[match(TRUE,OK)]
    }
    if ('t' %in% names(dframe)) {   ## synonym for 'session'
        levels(dframe$t)[!OK] <- levels(dframe$t)[match(TRUE,OK)]
    }
    dframe
}

############################################################################################
# getvalidlevels <- function (type, parnames, J) {
# 2018-10-29
getvalidlevels <- function (type, parnames, J, CJSp1) {
        
    # flag impossible parameters by primary session

    validlevels <- matrix(TRUE, nrow = length(parnames), ncol = J)
    dimnames(validlevels) <- list(parnames, 1:J)
    
    # all types
    if ('phi' %in% parnames) {
        validlevels['phi',J] <- FALSE
    }
    if ('move.a' %in% parnames) {
        validlevels['move.a',J] <- FALSE
    }
    if ('move.b' %in% parnames) {
        validlevels['move.b',J] <- FALSE
    }
    if (type %in% c('CJS', 'CJSmte')) {
        # if ('p' %in% parnames)
        # 2018-10-29
        if ('p' %in% parnames & !CJSp1)
            validlevels['p',1] <- FALSE
    }
    if (type %in% c('JSSAf','JSSAfCL', 'JSSAsecrf','JSSAsecrfCL', 'PLBf', 'PLBsecrf')) {
        if ('f' %in% parnames)
            validlevels['f',J] <- FALSE
    }
    if (type %in% c('JSSAg','JSSAgCL', 'JSSAsecrg','JSSAsecrgCL','Pradelg', 'PLBg', 'PLBsecrg')) {
        if ('gamma' %in% parnames)
            validlevels['gamma',1] <- FALSE
    }
    if (type %in% c('JSSAk','JSSAkCL', 'PLBk')) {
        if ('kappa' %in% parnames)
            validlevels['kappa',1] <- FALSE
    }
    if (type %in% c('JSSAl','JSSAlCL', 'JSSAsecrl','JSSAsecrlCL', 'PLBl', 'PLBsecrl',
                    'Pradel')) {
        if ('lambda' %in% parnames)
            validlevels['lambda',J] <- FALSE
    }
    if (type %in% c('JSSAb', 'JSSAbCL', 'JSSAsecrb','JSSAsecrbCL', 'PLBb', 'PLBsecrb')) {
        if ('b' %in% parnames)
            validlevels['b',1] <- FALSE  
    }
    if (type %in% c('CJSsecr')) {
        if (!CJSp1) {              ## 2018-10-29
            if ('lambda0' %in% parnames)
                validlevels['lambda0',1] <- FALSE
            if ('sigma' %in% parnames)
                validlevels['sigma',1] <- FALSE
        }
        # if ('move.a' %in% parnames)
        #     validlevels['move.a',1] <- FALSE
        # if ('move.b' %in% parnames)
        #     validlevels['move.b',1] <- FALSE
    }
    validlevels
}
############################################################################################
lpredictor <- function (model, newdata, indx, beta, field, beta.vcv=NULL, validlevels, contrasts = NULL) {
    vars <- all.vars(model)
    if (any(!(vars %in% names(newdata))))
        stop ("one or more model covariates not found in 'newdata'")
    newdata <- as.data.frame(newdata)

    lpred <- matrix(ncol = 2, nrow = nrow(newdata),dimnames=list(NULL,c('estimate','se')))

    # adjust for CJS nonidentifiable parameters
    # see also make.designmatrix in openCR.design.R
    if (!is.null(newdata$session))
        notOK <- !validlevels[field,newdata$session]   ## BEFORE adjust levels
    else
        notOK <- !validlevels[field,newdata$t]   ## BEFORE adjust levels
    newdata <- adjustlevels(field, newdata, validlevels)
    mat <- model.matrix(model, data=newdata, contrasts.arg = contrasts)

    ## drop pmix beta0 column from design matrix (always zero)
    if (field=='pmix') {
        mat <- mat[,-1,drop=FALSE]
    }
    lpred[,1] <- mat %*% beta[indx]
    if ('session' %in% names(newdata) | 't' %in% names(newdata))
        lpred[notOK,1] <- NA
    if (is.null(beta.vcv)) return ( cbind(newdata,lpred) )
    else {
        vcv <- beta.vcv[indx,indx]    ## OR maybe all betas?
        nrw <- nrow(mat)
        vcv <- apply(expand.grid(1:nrw, 1:nrw), 1, function(ij)
            mat[ij[1],, drop=F] %*% vcv %*% t(mat[ij[2],, drop=F]))  # link scale
        vcv <- matrix (vcv, nrow = nrw)
        lpred[,2] <- diag(vcv)^0.5
        temp <- cbind(newdata,lpred)
        attr(temp, 'vcv') <- vcv
        return(temp)
    }
}

############################################################################################

makerealparameters <- function (design, beta, parindx, link, fixed) {
    modelfn <- function(i) {
        ## linear predictor for real parameter i
        Yp <- design$designMatrices[[i]] %*% beta[parindx[[i]]]
        # delay untransform until inside C code
        if (i == "b") {
            link[["b"]] <- "identity"
        }
        if (i == "pmix") {
            Yp <- invmlogit(Yp)
        }
        else {
            Yp <- untransform(Yp, link[[i]])
        }
        Yp[design$parameterTable[,i]]   ## replicate as required
    }

    ## construct matrix of parameter values
    nrealpar  <- length(design$designMatrices)
    pnames <- names(link)  ## should be complete!

    if (length(fixed)>0)
        for (a in names(fixed))     ## bug fixed by adding this line 2011-09-28
            link[[a]] <- NULL
    if (length(link) != nrealpar)
        stop ("number of links does not match design matrices")

    temp <- sapply (names(parindx), modelfn)
    if (nrow(design$parameterTable)==1) temp <- t(temp)
    nrw <- nrow(temp)
    ## make new matrix and insert columns in right place
    if (is.null(nrw)) stop ("bug in makerealparameters")
    temp2 <- as.data.frame(matrix(nrow = nrw, ncol = length(parindx)))
    names(temp2) <- names(parindx)
    temp2[ , names(design$designMatrices)] <- temp       ## modelled
    if (!is.null(fixed) & length(fixed)>0)
        temp2[ , names(fixed)] <- sapply(fixed, rep, nrw)    ## fixed
    as.matrix(temp2[,pnames])   ## pnames to get order right 2011-12-13
}
## End of miscellaneous functions
############################################################################################

# from Shirley P
# see d:/density users/will rayment

# Logit function:
# ---------------

# Convert a vector to logits:

logit.fn <- function(invec)
{
    y <- log(invec/(1-invec))
    y[invec<=0] <- -10
    y[invec>=1] <- 10
    y
}

# Expit function:
# ---------------

# Inverse of logit:

expit.fn <- function(invec)
    1/(1+exp(-invec))

# PLG function:
# -------------

# Convert probs (summing to 1) to logit function of conditional probs.
# (one fewer element). e.g. (a1,a2,a3,a4) goes to
# (logit(a1),logit(a2/(1-a1)),logit(a3/(1-a1-a2)))

plg.fn <- function(invect)
   {
   len <- length(invect)
   if (len==2) newvect <- invect[1]
   if (len>2)
      {
      # Cumulative sum:
      csin <- cumsum(invect)
      # Ratio vector:
      newvect <- invect[1:(len-1)]/(1-c(0,csin[1:(len-2)]))
      }
   # return logit:
   logit.fn(newvect)
   }


# LGP function:
# -------------

# Inverse of plg.fn.
# Logit ratios to probs (sum to 1) function, one more element.

lgp.fn <- function(invect)
   {
   len <- length(invect)
   if (len==1) p.v <- expit.fn(invect)
   if (len>1)
      {
      # Get ratios vector:
      ratvec <- expit.fn(invect)
      # Convert to beta vector: need cum product of 1-ratios
      cp <- cumprod(1-ratvec)
      # Return probability vector:
      p.v <- ratvec[1:len]*c(1,cp[1:(len-1)])
      }
   c(p.v,1-sum(p.v))
   }
############################################################################################

## vector of cumulative number in each primary session
getcumss <- function(capthist) {
    intervals <- intervals(capthist)
    if (is.null(intervals)) intervals <- rep(1, ncol(capthist)-1)
    primarysession <- primarysessions(intervals)
    cumsum(c(0,table(primarysession)))  ## secondary sessions per primary session
}
############################################################################################

############################################################################################
fillpmix2 <- function (nc, nmix, PIA, realparval) {
    pmix <- matrix(1, nrow = nc, ncol = nmix) 
    if (nmix>1) {        
        c <- as.numeric(PIA[,1,1,])
        pmix[] <- realparval[c, 'pmix']   # assumes dim 2 of realparval has name
    }
    t(pmix)
}
############################################################################################

age.matrix <- function (capthist, initialage = 0, minimumage = 0, maximumage = 1, collapse = FALSE) #, bysession = FALSE) 
{    
    makeage <- function (n) {
        ch <- capthist[n,,,drop=FALSE]
        temp0 <- apply(abs(ch), 1:2,sum)  ## sum over traps
        first <- match(TRUE, cumsum(temp0)>0)
        firstsession <- primarysession[first]
        age <- (primarysession - firstsession) + initialage[n]
        # age <- rep(initialage[n], S)
        # recruited <- primarysession>=firstsession
        # age[recruited] <- initialage[n] + cumsum(c(0, intervals[recruited[-S]]))
        age <- pmax(age, minimumage)
        age <- pmin(age, maximumage)
        age
    }
    if (ms(capthist)) stop("age.matrix requires single-session capthist")
    n <- nrow(capthist)
    S <- ncol(capthist)
    primarysession <- primarysessions(intervals(capthist))
    if (is.character(initialage)) {
        if (!initialage %in% names(covariates(capthist)))
            stop ("initialage covariate ", initialage, " not found")
        initialage <- covariates(capthist)[[initialage]]
    }
    else
        initialage <- rep(initialage, length.out = n)
    out <- t(sapply(1:n, makeage))
    dimnames(out) <- dimnames(capthist)[1:2]
    if (collapse)
        out <- as.matrix(apply(out,1,paste,collapse=''))
    out
} 
############################################################################################

squeeze <- function(x) {
    if (ms(x)) {
        out <- lapply(x, squeeze)
        names(out) <- names(x)
        intervals(out) <- intervals(x)
        sessionlabels(out) <- sessionlabels(x)
        class(out) <- c('capthist','list')
        out
    }
    else {
        freq <- covariates(x)$freq
        if (is.null(freq) | length(unique(freq))==1) {
            dim <- dim(x)
            SK <- prod(dim[2:3])
            df <- as.data.frame(matrix(x, nrow = nrow(x)))
            if (!is.null(covariates(x)))
                if (nrow(covariates(x))>0)
                    df <- cbind(df, covariates(x))
            df <- plyr::count(df)
            ch <- array(unlist(df[,1:SK]), dim=c(nrow(df), dim[2:3]))
            class(ch) <- 'capthist'
            if (!is.null(traps(x)))
                traps(ch) <- traps(x)
            rownames(ch) <- 1:nrow(ch)
            intervals(ch) <- intervals(x)
            sessionlabels(ch) <- sessionlabels(x)
            timevaryingcov(ch) <- timevaryingcov(x)
            covariates(ch) <- df[, (SK+1):ncol(df), drop = FALSE]
            if (!is.null(freq))
                covariates(ch)$freq <- covariates(ch)$freq * freq[1]
            ch
        } 
        else {
            # warning("capthist already has variable freq covariate; not squeezed")
            x
        }
    }
}
############################################################################################

unsqueeze <- function(x) {
    if (ms(x)) {
        out <- lapply(x, unsqueeze)
        names(out) <- names(x)
        intervals(out) <- intervals(x)
        sessionlabels(out) <- sessionlabels(out)
        class(out) <- class(x)
        out
    }
    else {
        freq <- covariates(x)$freq
        n <- nrow(x)
        if (!is.null(freq) | (length(freq)==n)) {
            i <- rep(1:n, freq)
            if (length(dim(x))==2)
                ch <- x[i,, drop = FALSE]
            else if (length(dim(x))==3)
                ch <- x[i,,, drop = FALSE]
            class(ch) <- class(x)
            if (!is.null(traps(x)))
                traps(ch) <- traps(x)
            rownames(ch) <- 1:nrow(ch)
            intervals(ch) <- intervals(x)
            sessionlabels(ch) <- sessionlabels(x)
            timevaryingcov(ch) <- timevaryingcov(x)
            covariates(ch) <- covariates(x)[i,, drop = FALSE]
            covariates(ch)$freq[] <- 1
            ch
        } 
        else {
            # warning("capthist has no freq covariate; not unsqueezed")
            x
        }
    }
}
############################################################################################


## code used in openCR.fit to put capthist object in a standard form
## 2018-02-11

stdcapthist <- function (capthist, type, nclone, squ, ...) {
    if (!inherits(capthist, 'capthist'))
        stop ("requires 'capthist' object")
    if (!ms(capthist) & is.null(intervals(capthist))) {
        # warning ("intervals for single-session capthist set to 1")
        intervals(capthist) <- rep(1, ncol(capthist)-1)
    }
    
    if (ms(capthist)) {    ## collapse ms to single
        capthist <- join(capthist, drop.sites = !grepl('secr',type), ...)
    }
    interv <- intervals(capthist)
    sessnames <- sessionlabels(capthist)
    timevarcov <- timevaryingcov(capthist)
    if (grepl('secr', type)) {
        if (detector(traps(capthist))[1] == 'single') {
            capthist <- reduce(capthist, output = 'multi', dropunused = FALSE)
            warning ("capthist coerced to 'multi' detector type")
        }
        else if (!(detector(traps(capthist))[1] %in% c('proximity','count','multi'))) {
            capthist <- reduce(capthist, output = 'proximity', dropunused = FALSE)
            warning ("capthist coerced to 'proximity' detector type")
        }
    }
    intervals(capthist) <- interv   ## restore if reduce.capthist has lost them
    sessionlabels(capthist) <- sessnames
    timevaryingcov(capthist) <- timevarcov
    
    # no cloning multiplier
    if (all(nclone == 1)) {
        if (squ) capthist <- squeeze(capthist)
    }
    # cloning multiplier
    else {   
        if (is.null(covariates(capthist)$freq)) {
            if (is.null(covariates(capthist))) {
                covariates(capthist) <- data.frame(freq = rep(1,nrow(capthist)))
            }
            else {
                covariates(capthist)$freq <- rep(1,nrow(capthist))
            }
        }
        covariates(capthist)$freq <- covariates(capthist)$freq * nclone
    }
    capthist
}

########################################################################################

complete.beta <- function (object) {
    fb <- object$details$fixedbeta
    if (!is.null(fb)) {
        names(fb) <- unlist(sapply(object$design$designMatrices, colnames))
        fb[is.na(fb)] <- object$fit$par
        beta <- fb
    }
    else {
        beta <- object$fit$par
    }
    beta
}
###############################################################################

complete.beta.vcv <- function (object) {
    fb <- object$details$fixedbeta
    if (!is.null(fb)) {
        names(fb) <- unlist(sapply(object$design$designMatrices, colnames))
        nbeta <- length(fb)
        beta.vcv <- matrix(0, nrow = nbeta, ncol = nbeta, dimnames=list(names(fb), names(fb)))
        if (!is.null(object$beta.vcv))
            beta.vcv[is.na(fb[row(beta.vcv)]) & is.na(fb[col(beta.vcv)])] <- object$beta.vcv
    }
    else {
        beta.vcv <- object$beta.vcv
    }
    beta.vcv
}
###############################################################################

## Based on Charles C. Berry on R-help 2008-01-13
## drop = FALSE 2018-11-22
n.unique.rows <- function(x) {
    order.x <- do.call(order, as.data.frame(x))
    equal.to.previous <- rowSums(x[tail(order.x,-1),,drop = FALSE] != x[head(order.x,-1),,drop = FALSE])==0 
    1 + sum(!equal.to.previous)
}

individualcovariates <- function (PIA) {
    pia <- matrix(PIA, nrow = nrow(PIA))
    n.unique.rows(pia) > 1
}
###############################################################################

# moved from prwisecr.R 2020-12-12

# The mqarray is a lookup array giving the pixel in the output mask
# that corresponds to a particular m in the input mask and q in the
# kernel [q * mm + m].
#
# Destinations that lie outside the mask receive a value of -1.
#
# mqsetup() initialises mqarray for a particular mask and kernel.

rectwrap <- function(oldx, oldy, newx, newy, concat = TRUE) {
  # improved wrapping 2020-10-29
  # assumes integer coordinates
  
  # rbind(-5:5, -5:5 %% 3)
  #      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11]
  # [1,]   -5   -4   -3   -2   -1    0    1    2    3     4     5
  # [2,]    1    2    0    1    2    0    1    2    0     1     2
  
  xmin <- min(oldx)
  ymin <- min(oldy)
  xmax <- max(oldx)
  ymax <- max(oldy)
  newx <- newx %% (xmax-xmin+1)
  newy <- newy %% (ymax-ymin+1)
  if (concat) 
    paste(newx, newy)
  else 
    list(newx=newx, newy=newy)
}

mqsetup <- function (
  mask,      ## x,y points on mask (first x, then y)
  kernel,    ## list of integer dx,dy points on kernel with p(move|dx,dy) mask (first x, then y)
  cellsize,   ## side of grid cell (m) for each mask point
  edgecode
)
{
  ## assuming cells of mask and kernel are same size
  ## and kernel takes integer values centred on current mask point
  
  oldx <- round((mask$x-min(mask$x))/cellsize)
  oldy <- round((mask$y-min(mask$y))/cellsize)
  
  newx <- as.integer(outer(oldx, kernel$x, "+"))
  newy <- as.integer(outer(oldy, kernel$y, "+"))
  
  # mqarray shared with C++ so indices are zero-based
  if (edgecode == 1)    # "wrap"
    newxy <- rectwrap(oldx,oldy,newx,newy)
  else                  # "truncate", "none"
    newxy <- paste(newx, newy)
  
  i <- match(newxy, paste(oldx,oldy)) - 1
  i[is.na(i)] <- -1
  matrix(i, nrow = nrow(mask), ncol = nrow(kernel))
  
}

###############################################################################
