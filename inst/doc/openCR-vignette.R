## ----settings, echo = FALSE---------------------------------------------------
mycache <- FALSE
Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")  ## to dodge CRAN ASAN issue

## ----setup, message = FALSE, eval = TRUE, warning = FALSE-------------------------------
library(openCR)                   # also loads secr
options(digits = 4, width = 90)   # for more readable output

## ----dipper, eval = TRUE, message = FALSE, warning = FALSE------------------------------
m.array(dipperCH, never.recap = T)   # compare Lebreton et al. 1992 Table 10

## ----dipperCJS, fig.width = 4.5, fig.height = 3.3---------------------------------------
dipper.phi.t <- openCR.fit(dipperCH, type = 'CJS', model = phi~t)
predict(dipper.phi.t)
plot(dipper.phi.t, par = 'phi', ylim = c(0,1), pch = 16, col = 'red')

## ----figfunc, ref.label = "figurefunctions", echo = FALSE-------------------------------
onemulti <- function(st = c(0,6,11,15), le = c(5,4,3,5), yb = 7, col=col1, outer = TRUE) {
    col <- rep(col, le)
    xl <- unlist(mapply(":",st,le+st-1))
    yb <- rep(yb,length(xl))
    xr <- xl + width
    yt <- yb + height
    rect(xl,yb,xr,yt,col=col)
    text(xl+width/2, yb+height/2, unlist(mapply(":", 1, le)))
    
    xl <- st - margin
    yb <- rep(yb[1], length(xl)) - margin
    xr <- st+le-1+width+margin
    yt <- yb+height+2*margin
    rect(xl,yb,xr,yt)
    text(st+le/2, rep(yb[1]+2*margin,length(st))+height+0.5, paste('session',1:length(st)))
    if (outer) {
        rect(st[1]-3*margin, yb[1]-2*margin, tail(st+le-1,1)+width+3*margin, 
          yb[1]+height+8*margin)
    }
}

onejoined <- function(offset = 1.5, le = c(5,4,3,5), yb = 2.2, col=col1, intervals = TRUE,
                      intlabel = 'intervals', leftlabel = '', outer = TRUE) {
    col <- rep(col, le)
    xl <- 0:(sum(le)-1)+offset
    yb <- rep(yb,length(xl))
    xr <- xl + width
    yt <- yb + height
    rect(xl,yb,xr,yt,col=col)
    text(xl+width/2, yb+height/2, c(1:length(xl)))
    if (intervals) {
        xi <- offset + (1:(length(xl)-1)) - (1-width)/2
        xip <- cumsum(le)[-length(le)]   # intermediate between primary sessions
        intervals <- rep(0,length(xi))
        intervals[xip] <- 1
        text(xi, yb [-1]-0.8, intervals)
        text(-0.2, yb[1]-0.8, intlabel)
        segments(xi[xip], rep(yb[1]-0.4,length(xip)), xi[xip], rep(yb[1]+0.4,
          length(xip))+height)
    }
    text (0.4, yb[1]+height/2, leftlabel, adj = c(1,0.5))
    if (outer) {
        rect(offset-2*margin, yb[1]-2*margin, sum(le)-1+offset+width+2*margin, 
          yb[1]+height+2*margin)
    }
}

## ----fig1plot, ref.label = "fig1", fig.width=8, fig.height=4, echo = FALSE--------------
# Fig. 1 Single-stratum data
par(cex=1, xpd = TRUE, mfrow = c(1,1), mar=c(1,4,1,4))
width <- 0.85
height <- 1.1
margin <- 0.15
col1 <- c('salmon','pink','brown', 'red')
col2 <- c('green','lightgreen','darkgreen', 'lightblue')
MASS::eqscplot(0,0,xlim=c(0,20), ylim=c(0,8), type='n', axes=F,xlab='',ylab='')
onemulti(col = col1)
text(9, 5.2, 'join()', cex=1.1)
arrows (10.7,6.2,10.7,4.2)
onejoined(leftlabel='')

## ----fig2plot, ref.label = "fig2", fig.width=8, fig.height=4.5, echo = FALSE------------
# Fig. 2 Multi-stratum data
par(cex = 0.9, xpd = TRUE, mfrow = c(1,1), mar = c(1,4,1,4))
MASS::eqscplot(0,0,xlim=c(-3,20), ylim=c(-2,8), type='n', axes=FALSE, xlab = '',ylab='')
onejoined(leftlabel='stratum 1', yb = 6.5, intlabel='')
onejoined(leftlabel='stratum 2', yb = 3, intlabel='')
onejoined(leftlabel='stratum 3', yb = -0.5, le = c(4,3,4,4), intlabel='', col = col2)
rect(-3, -2, 19.3, 8.7)

## ----compare, cache = mycache, warning = FALSE------------------------------------------
msk <- make.mask(traps(captdata), buffer = 100, type = 'trapbuffer')

fit_secr <- secr.fit(captdata, detectfn = 'HHN', mask = msk, trace = FALSE)
fit_openCR <- openCR.fit(captdata, detectfn = 'HHN', mask = msk, type = 'secrD')

# massage the predict.openCR results to the same format as predict.secr
pred_openCR <- plyr::rbind.fill(predict(fit_openCR))

pred_openCR <- pred_openCR[c(2,1,3), !(names(pred_openCR) %in% c('stratum','session'))]
rownames(pred_openCR) <- fit_secr$realnames

# compare estimates
predict(fit_secr)[,-1]
pred_openCR

## ----timing-----------------------------------------------------------------------------
# compare timings in seconds
c(secr = fit_secr$proctime, openCR = fit_openCR$proctime)

## ----multinom---------------------------------------------------------------------------
# compare maximised log likelihoods
c(secr.logLik = logLik(fit_secr), openCR.logLik = logLik(fit_openCR) + logmultinom(captdata))

## ----makedf, cache = mycache------------------------------------------------------------
makedf.b <- function (ch, spatial = FALSE, nmix = 1, naive = FALSE) {
  R <- 1 # assume single stratum
  ch <- squeeze(ch)
  # Construct matrix of logical values TRUE iff caught before 
  detected <- apply(abs(ch),1:2,sum)>0
  detected <- t(apply(detected, 1, cumsum)>0)
  if (naive)
    b <- rep(FALSE, prod(dim(ch)[1:2]))
  else
    b <- t(apply(detected, 1, function(x) {x[which.max(x)] <- FALSE; x}))
  # For a simple non-spatial case: data.frame(customb = as.vector(b))  
  # More generally:
  n <- nrow(ch)
  S <- ncol(ch)
  K <- if (spatial) dim(ch)[3] else 1
  data.frame(customb = insertdim(b, c(2,3,1), c(R,n,S,K,nmix)))  
}

## ----customb, cache = mycache, warning = FALSE------------------------------------------
ovenj <- join(ovenCH)
fitb <- openCR.fit(ovenj, model = p ~ b)
fitbc <- openCR.fit(ovenj, model = p ~ customb, dframe = makedf.b(ovenj))
AIC(fitb, fitbc)

## ----customb2, cache = mycache, warning = FALSE-----------------------------------------
fitb2 <- openCR.fit(ovenj, model = p ~ b, type = 'JSSAfCL', start = fitb)
fitbc2 <- openCR.fit(ovenj, model = p ~ customb,  type = 'JSSAfCL', 
                    dframe = makedf.b(ovenj), dframe0 = makedf.b(ovenj, naive = TRUE))
AIC(fitb2, fitbc2)

## ----transient, cache = mycache---------------------------------------------------------
makedf.resident <- function (ch, spatial = FALSE, nmix = 1) {
  nstrata <- 1 # assume single stratum
    ch <- squeeze(ch)
    n <- nrow(ch)
    S <- ncol(ch)
    K <- if (spatial) dim(ch)[3] else 1
    primary <- primarysessions(intervals(ch))
    detected <- apply(abs(ch),1:2,sum)>0
    nprimary <- apply(detected, 1, function(x) length(unique(primary[x])))
    data.frame(resident = insertdim(nprimary>1, 1, c(nstrata, n, S, K, nmix)))  
}

## ----transient2, cache = mycache--------------------------------------------------------
addresidentcov <- function (ch) {
    primary <- primarysessions(intervals(ch))
    detected <- apply(abs(ch), 1:2, sum)>0
    nprimary <- apply(detected, 1, function(x) length(unique(primary[x])))
    covariates(ch) <- data.frame(residentcov =  nprimary>1)
    ch
}

## ----transient3, cache = mycache--------------------------------------------------------
ovenj <- join(ovenCH)
ovenj <- addresidentcov(ovenj)
fitnull <- openCR.fit(ovenj, model = phi ~ 1)
fitcov  <- openCR.fit(ovenj, model = phi ~ residentcov)
fitdf   <- openCR.fit(ovenj, model = phi ~ resident, dframe = makedf.resident(ovenj))
fits <- openCRlist(fitnull, fitcov, fitdf)
AIC(fits)
pred <- predict(fits, newdata = data.frame(resident = TRUE, residentcov = TRUE))
do.call(rbind, lapply(pred, '[[', 'phi'))

## ----transient4, eval = FALSE-----------------------------------------------------------
#  addresidentcov2 <- function (ch, d = 1) {
#      primary <- primarysessions(intervals(ch))
#      secondary <- secondarysessions(intervals(ch))
#      detected <- apply(abs(ch), 1:2, sum)>0
#      nprimary <- apply(detected, 1, function(x) length(unique(primary[x])))
#      dsecondary <- apply(detected, 1, function(x)
#          max(by(secondary[x], primary[x], function(y) diff(range(y)))))
#      covariates(ch) <- data.frame(residentcov1 = nprimary>1,
#                                   residentcov2 = nprimary>1 | dsecondary>=d)
#      ch
#  }

## ----dummy, eval = TRUE, cache = mycache------------------------------------------------
fit0 <- openCR.fit(ovenCH, model = p~t)
fitd <- openCR.fit(ovenCH, model = p ~ -1+t)
coef(fit0)
coef(fitd)

## ----dummy2, eval = TRUE, cache = mycache-----------------------------------------------
fitd2 <- openCR.fit(ovenCH, model = p~t, details = list(dummyvariablecoding = 't'))
coef(fitd2)

## ----contrsum, eval = TRUE, cache = mycache---------------------------------------------
fit <- openCR.fit(dipperCH, model = phi~t, details = list(contrasts = list(t = contr.sum)))
invlogit(coef(fit)['phi',c('beta','lcl','ucl')])

## ----sparsekernel, fig.width=8, fig.height=3.5------------------------------------------
par(mar = c(3,1,4,5))
k <- make.kernel(movementmodel = 'BVN', kernelradius = 10, spacing = 10, move.a = 40, 
  sparse = TRUE, clip = TRUE)
plot(k)
symbols(0,0, add = TRUE, circles = 100, inches = FALSE)

## ----plotkernel, fig.width = 8, fig.height = 3.5----------------------------------------
par (mar = c(3,3,4,6), cex = 0.9)
k <- make.kernel (movementmodel = 'BVN', spacing = 10, move.a = 40, clip = TRUE)
plot(k, contour = TRUE)
summary(k)

## ----settlementexample, eval = FALSE----------------------------------------------------
#  ovenCHb <- reduce(ovenCHp, by = 'all', outputdetector = 'count')
#  msk <- make.mask(traps(ovenCHp[[1]]), buffer = 500, spacing = 40, type = 'trapbuffer')
#  # uniform settlement
#  fit0 <- openCR.fit(ovenCHb, type = 'PLBsecrf', mask = msk, binomN = 1,
#    movementmodel = 'BVN', details = list(settlemodel = FALSE))
#  # logarithmic N-S gradient in settlement
#  fit1 <- openCR.fit(ovenCHb, type = 'PLBsecrf', mask = msk, binomN = 1,
#    movementmodel = 'BVN', details = list(settlemodel = TRUE), model = settle~y)

## ----loadsettlementresults, echo = FALSE------------------------------------------------
load('settlement.RData')

## ----settlementexampleresults-----------------------------------------------------------
AIC(fit0, fit1)[,-6]
coef(fit1)

## ----derived, cache = mycache-----------------------------------------------------------
dipperCL <- openCR.fit(dipperCH, type = 'JSSAlCL', 
            model = list(lambda~t, phi~t))
# only these parameters are in the model and estimated directly,
names(predict(dipperCL))
# but we can derive b, f, gamma and N, as well as the super-population N
d <- derived(dipperCL)
print(d, digits = 3, legend = TRUE)

## ----fitnm, cache = mycache-------------------------------------------------------------
fitnr <- openCR.fit(ovenCH, type = 'JSSAlCL', model = list(phi ~ t, lambda~t))
fitnm <- openCR.fit(ovenCH, type = 'JSSAlCL', model = list(phi ~ t, lambda~t),
                    method = "Nelder-Mead", details = list(control = list(maxit = 5000)))

## ----aic, cache = mycache---------------------------------------------------------------
AIC(fitnm,fitnr)

## ----fitnm2, cache = mycache------------------------------------------------------------
fitnm <- openCR.fit(ovenCH, type = 'JSSAlCL', model = list(phi ~ t, lambda~t),
                    method = "Nelder-Mead", details = list(control = list(maxit = 2000)),
                    start = fitnr)
AIC(fitnm,fitnr)

## ----ncores, eval = FALSE---------------------------------------------------------------
#  # RCPP_PARALLEL_NUM_THREADS
#  # recommended for quad-core Windows PC
#  setNumThreads(7)

## ----ucare------------------------------------------------------------------------------
if (requireNamespace("R2ucare"))
    ucare.cjs(dipperCH, verbose = FALSE, by = 'sex')

## ----figurefunctions, eval = FALSE------------------------------------------------------
#  onemulti <- function(st = c(0,6,11,15), le = c(5,4,3,5), yb = 7, col=col1, outer = TRUE) {
#      col <- rep(col, le)
#      xl <- unlist(mapply(":",st,le+st-1))
#      yb <- rep(yb,length(xl))
#      xr <- xl + width
#      yt <- yb + height
#      rect(xl,yb,xr,yt,col=col)
#      text(xl+width/2, yb+height/2, unlist(mapply(":", 1, le)))
#  
#      xl <- st - margin
#      yb <- rep(yb[1], length(xl)) - margin
#      xr <- st+le-1+width+margin
#      yt <- yb+height+2*margin
#      rect(xl,yb,xr,yt)
#      text(st+le/2, rep(yb[1]+2*margin,length(st))+height+0.5, paste('session',1:length(st)))
#      if (outer) {
#          rect(st[1]-3*margin, yb[1]-2*margin, tail(st+le-1,1)+width+3*margin,
#            yb[1]+height+8*margin)
#      }
#  }
#  
#  onejoined <- function(offset = 1.5, le = c(5,4,3,5), yb = 2.2, col=col1, intervals = TRUE,
#                        intlabel = 'intervals', leftlabel = '', outer = TRUE) {
#      col <- rep(col, le)
#      xl <- 0:(sum(le)-1)+offset
#      yb <- rep(yb,length(xl))
#      xr <- xl + width
#      yt <- yb + height
#      rect(xl,yb,xr,yt,col=col)
#      text(xl+width/2, yb+height/2, c(1:length(xl)))
#      if (intervals) {
#          xi <- offset + (1:(length(xl)-1)) - (1-width)/2
#          xip <- cumsum(le)[-length(le)]   # intermediate between primary sessions
#          intervals <- rep(0,length(xi))
#          intervals[xip] <- 1
#          text(xi, yb [-1]-0.8, intervals)
#          text(-0.2, yb[1]-0.8, intlabel)
#          segments(xi[xip], rep(yb[1]-0.4,length(xip)), xi[xip], rep(yb[1]+0.4,
#            length(xip))+height)
#      }
#      text (0.4, yb[1]+height/2, leftlabel, adj = c(1,0.5))
#      if (outer) {
#          rect(offset-2*margin, yb[1]-2*margin, sum(le)-1+offset+width+2*margin,
#            yb[1]+height+2*margin)
#      }
#  }

## ----fig1, eval = FALSE-----------------------------------------------------------------
#  # Fig. 1 Single-stratum data
#  par(cex=1, xpd = TRUE, mfrow = c(1,1), mar=c(1,4,1,4))
#  width <- 0.85
#  height <- 1.1
#  margin <- 0.15
#  col1 <- c('salmon','pink','brown', 'red')
#  col2 <- c('green','lightgreen','darkgreen', 'lightblue')
#  MASS::eqscplot(0,0,xlim=c(0,20), ylim=c(0,8), type='n', axes=F,xlab='',ylab='')
#  onemulti(col = col1)
#  text(9, 5.2, 'join()', cex=1.1)
#  arrows (10.7,6.2,10.7,4.2)
#  onejoined(leftlabel='')

## ----fig2, fig.width=8, fig.height=4.5, eval = FALSE------------------------------------
#  # Fig. 2 Multi-stratum data
#  par(cex = 0.9, xpd = TRUE, mfrow = c(1,1), mar = c(1,4,1,4))
#  MASS::eqscplot(0,0,xlim=c(-3,20), ylim=c(-2,8), type='n', axes=FALSE, xlab = '',ylab='')
#  onejoined(leftlabel='stratum 1', yb = 6.5, intlabel='')
#  onejoined(leftlabel='stratum 2', yb = 3, intlabel='')
#  onejoined(leftlabel='stratum 3', yb = -0.5, le = c(4,3,4,4), intlabel='', col = col2)
#  rect(-3, -2, 19.3, 8.7)

