## R functions to avoid C for prwi (spatial and nonspatial)

## prwi
## prwisecr

## 2018-02-12, 2019-05-19

## openval has
## column 1 p or lambda0
## column 2 phi
## column 3 f, g, l  (recruitment parameter)
## column 4 sigma (pmix if non-spatial)
## column 5 pmix

#########################################################################################

 # probability of count for session s, detector k, animal i
 # The argument 'g' is understood to be a cumulative hazard if binomN=0,
 # a probability otherwise

pski <- function(binomN, count, Tski, g) {

    result <- 1.0

    if (binomN == -1) {                              ## binary proximity detectors : Bernoulli
        if (any(abs(Tski-1) > 1e-10)) {              ## effort not unity; adjust g
            g <- 1 - (1 - g)^Tski
        }
        if (count>0)
            result <- g
        else
            result <- 1 - g
    }
    else if (binomN == 0) {                          ## count detectors : Poisson
        if (count == 0)
            result <- exp(-Tski * g)                 ## routinely apply Tsk adjustment to cum. hazard
        else
            result <- dpois(count, Tski * g, FALSE)
    }
    else if (binomN == 1) {                          ## count detectors : Binomial, size from Tsk
        result <- dbinom (count, round(Tski), g, FALSE)
    }
    else if (binomN > 1) {                           ## count detectors : Binomial, specified size
        if (abs(Tski-1) > 1e-10) {                   ## effort not unity, adjust g
            g <- 1 - (1 - g)^Tski
        }
        result <- dbinom (count, binomN, g, FALSE)
    }
    else stop("binomN < -1 not allowed")
    result
}

#########################################################################################

# The mqarray is a lookup array giving the pixel in the output mask
# that corresponds to a particular m in the input mask and q in the
# kernel [q * mm + m].
#
# Destinations that lie outside the mask receive a value of -1.
#
# mqsetup() initialises mqarray for a particular mask and kernel.

mqsetup <- function (
    mask,      ## x,y points on mask (first x, then y)
    kernel,    ## list of dx,dy points on kernel with p(move|dx,dy) mask (first x, then y)
    cellsize   ## side of grid cell (m) for each mask point
)
{
    ## assuming cells of mask and kernel are same size
    ## and kernel takes integer values centred on current mask point

    oldx <- floor(mask$x/cellsize)
    oldy <- floor(mask$y/cellsize)

    newx <- as.numeric(outer(oldx, kernel$x, "+"))
    newy <- as.numeric(outer(oldy, kernel$y, "+"))

    i <- match(paste(newx,newy), paste(oldx,oldy)) - 1
    i[is.na(i)] <- -1

    matrix(i, nrow = nrow(mask), ncol = nrow(kernel))

}

convolvemq <- function (
    j,         ## session number 1..jj
    kernelp,   ## p(move|dx,dy) for points in kernel
    mqarray,   ## input
    pjm
)
{
    mm <- nrow(mqarray)
    kn <- ncol(mqarray)
    workpjm <- numeric(mm)

    ## convolve movement kernel and pjm...
    for (m in 1:mm) {
        for (q in 1:kn) {              ## over movement kernel
            mq <- mqarray[m,q] + 1     ## which new point corresponds to kernel point q relative to current mask point m
            if (mq >= 0) {             ## post-dispersal site is within mask
                if (mq>mm) stop("mq > mm")
                workpjm[mq] <- workpjm[mq] + pjm[m] * kernelp[q,j]   ## probability of this move CHECK DIM KERNELP
            }
        }
    }
    workpjm
}

###################################################################################
prwi <- function (type, n, x, jj, cumss, nmix, w, fi, li, openval, PIA, PIAJ,
                      intervals, CJSp1) {
    # get session-specific real parameter values
    p <- getp (n, x, openval, PIA)
    phij <- getphij (n, x, openval, PIAJ, intervals)
    if (type == 1) {
        minb <- fi[n]
        cjs <- 1-CJSp1
    }
    else {
        minb <- 1
        cjs <- 0
        beta <- getbeta (type, n, x, openval, PIAJ, intervals, phij)
    }
    maxb <- fi[n]
    mind <- abs(li[n])
    maxd <- jj
    if (li[n] < 0) maxd <- mind     # possible censoring

    pdt <- 0

    # loop over possible birth and death times
    for (b in minb:maxb) {
        for (d in mind:maxd) {
            pbd <- if (type==1) 1 else beta[b]
            if (b<d) pbd <- pbd * prod(phij[b:(d-1)])
            if (li[n]>0)    # not censored
                pbd <- pbd * (1-phij[d])
            # pj now accounts for birth at b and survival to d
            # next multiply by conditional probability of observed CH
            pjt <- 1
            if ((b+cjs) <= d) {
                for (j in (b+cjs):d) {
                    # detection probability for each secondary session
                    # in primary session j
                    s <- (cumss[j]+1) : cumss[j+1]
                    counts <- abs(w[n, s])
                    pjt <- pjt * prod(ifelse(counts>0, p[s], 1-p[s]))
                    dead <- any(w[n, s]<0)
                }
            }
            pdt <- pdt + pbd * pjt
        }
    }
    pdt
}

###################################################################################

prw <- function (n, j, x, nc, jj, kk, binomN, cumss, PIA, gk, Tsk, w, pjm) {
    dead <- FALSE
    for (s in (cumss[j]+1):cumss[j+1]) {
        for (k in 1:kk) {
            c <- PIA[n,s,k,x]
            if (c >= 1) {    # drops unset traps
                count <- w[n,s,k]
                if (count<0) {count <- -count; dead <- TRUE }
                pjm <- pjm * pski(binomN,count,Tsk[k,s], gk[c,k,])
            }
        }
        if (dead) break;   # after processing all traps on this occasion
    }
    pjm
}
###################################################################################

prwisecr <- function (type, n, x, nc, jj, kk, mm, nmix, cumss, w, fi, li, gk,
                      openval, PIA, PIAJ, binomN, Tsk, intervals,
                      CJSp1, calcpdotbd, moveargsi, movemodel,
                      usermodel, kernel = NULL, mqarray = NULL, cellsize = NULL) {
    phij <- getphij (n, x, openval, PIAJ, intervals)

    if (movemodel > 1) {
        moveargsi <- pmax(moveargsi,0)
        moveargs <- getmoveargs (n, x, openval, PIAJ, intervals, moveargsi)
        kernelp <- fillkernelp (jj, movemodel-2, kernel, usermodel, cellsize, moveargsi, moveargs, normalize = TRUE)
    }

    if(type==6) {
        minb <- fi[n]
        cjs <- 1 - CJSp1
    }
    else {
        minb <- 1
        cjs <- 0
        beta <- getbeta (type, n, x, openval, PIAJ, intervals, phij)
    }
    if (calcpdotbd) {
        ## Precompute session-specific Pr for all unique parameter combinations 2019-05-05
        ## Does not allow for learned response
        pjmat <- pr0njmx(n, x, cumss, jj, mm, binomN, PIA, gk, Tsk)  ## PCH1.R
    }
    maxb <- fi[n]
    mind <- abs(li[n])
    maxd <- jj
    if (li[n] < 0) maxd <- mind     # possible censoring
    pdt <- 0
    pdotbd <- 1.0
    for (b in minb:maxb) {
        for (d in mind:maxd) {
            if (calcpdotbd) {
                # compute pdot for this b,d
                if (movemodel == 1) {   ## uncorrelated; product over primary sessions
                    if (d>(b+cjs)) {
                        pdotbd <- 1 - prod(apply(pjmat[(b+cjs):d,]/mm,2,sum))
                    }
                    else {
                        pdotbd <- 1
                    }
                }
                else {
                    pr0 <- pjmat[b+cjs,] / mm
                    if (d>(b+cjs)) {
                        for (j in (b+cjs+1):d) {
                            if (movemodel>1) pr0 <- convolvemq(j-1, kernelp, mqarray, pr0)
                            pr0 <- pr0 * pjmat[j,]
                        }
                    }
                    pdotbd <- 1-sum(pr0)
                }
            }
            if (type==6) pbd <- 1 else pbd <- beta[b]
            if (b<d) pbd <- pbd * prod(phij[b:(d-1)])
            if (li[n]>0)    # not censored
                pbd <- pbd * (1-phij[d])
            if ((b+cjs) <= d) {
                if (movemodel == 0) {
                    alpha <- rep(1/mm,mm)
                    for (j in (b+cjs):d) {
                        alpha <- prw (n, j, x, nc, jj, kk, binomN, cumss, PIA,
                                    gk, Tsk, w, alpha)
                    }
                    prwi <- sum(alpha) / pdotbd
                }
                else if ( movemodel == 1) { # uncorrelated; product over primary sessions
                    prwi <- 1.0
                    for (j in (b+cjs):d) {
                        alpha <- prw(n, j, x, nc, jj, kk, binomN, cumss, PIA,
                                   gk, Tsk, w, rep(1/mm,mm))
                        prwi <- prwi * sum(alpha)
                    }
                    prwi <- prwi / pdotbd
                }
                else { # movemodel>1
                    alpha <- prw (n, b+cjs, x, nc, jj, kk, binomN, cumss, PIA,
                                      gk, Tsk, w, rep(1/mm, mm))
                    if (d>(b+cjs)) {
                        for (j in (b+cjs+1):d) {
                            alpha <- convolvemq(j-1, kernelp, mqarray, alpha)
                            alpha <- prw (n, j, x, nc, jj, kk, binomN, cumss, PIA,
                                          gk, Tsk, w, alpha)
                        }
                    }
                    prwi <- sum(alpha) / pdotbd
                    # cat("n ", n, " b ", b, " d ", d, " prwi ", prwi, "\n")
                }
                pdt <- pdt + pbd * prwi
            }
        }
    }
    pdt
}

#==============================================================================

prwmulti <- function (n, j, x, cumss, w, PIA, hk, Tsk, h, p0, hindex, pjm) {

    dead <- FALSE
    ## over secondary sessions (occasions) in this primary session
    for (s in (cumss[j]+1):cumss[j+1]) {
        wi <- w[n, s]
        if (wi < 0) dead <- TRUE;
        k <- abs(wi)         ## trap number 1..kk k = 0 if not caught

        ## Not captured in any trap on occasion s
        if (k < 1) {
            OK <- h[x,,hindex[n,s]] > 1e-8
            pjm[OK] <- pjm[OK] * p0[x,OK,hindex[n,s]]
        }
        ## Captured in trap k
        else {
            c <- PIA[n, s, k, x]
            if (c >= 1) {    ## drops unset traps
                pjm <- pjm * Tsk[k,s] * (1-p0[x,,hindex[n,s]]) *  hk[c, k, ] / h[x,,hindex[n,s]]
            }
        }
        if (dead) break   ## out of s loop
    }
    pjm
}

prwisecrmulti <- function (type, n, x, jj, mm, cumss, w, fi, li, hk, openval,
                      PIA, PIAJ, binomN, Tsk, intervals,  moveargsi, h, hindex,
                      CJSp1, calcpdotbd, movemodel, usermodel, kernel = NULL,
                      mqarray = NULL, cellsize = NULL) {
    ## precompute p0 to save time
    p0 <- exp(-h)
    phij <- getphij (n, x, openval, PIAJ, intervals)
    if (movemodel > 1) {
        moveargsi <- pmax(moveargsi,0)
        moveargs <- getmoveargs (n, x, openval, PIAJ, intervals, moveargsi)
        kernelp <- fillkernelp (jj, movemodel-2, kernel, usermodel, cellsize, moveargsi, moveargs)
    }
    if(type==6) {
        minb <- fi[n]
        cjs <- 1-CJSp1
    }
    else {
        minb <- 1
        cjs <- 0
        beta <- getbeta (type, n, x, openval, PIAJ, intervals, phij)
    }
    if (calcpdotbd) {
        pjmat <- pr0njmx(n, x, cumss, jj, mm, binomN, PIA, hk, Tsk)  ## PCH1.R
    }
    maxb <- fi[n]
    mind <- abs(li[n])
    maxd <- jj
    if (li[n] < 0) maxd <- mind     # possible censoring

    pdt <- 0
    pdotbd <- 1
    for (b in minb:maxb) {
        for (d in mind:maxd) {
            if (calcpdotbd) {
                # compute pdot for this b,d
                if (movemodel == 1) {   ## uncorrelated; product over primary sessions
                    if (d>(b+cjs)) {
                        pdotbd <- 1 - prod(apply(pjmat[(b+cjs):d,]/mm,2,sum))
                    }
                    else {
                        pdotbd <- 1
                    }
                }
                else {
                    pr0 <- pjmat[b+cjs,] / mm
                    if (d>(b+cjs)) {
                        for (j in (b+cjs+1):d) {
                            if (movemodel>1) pr0 <- convolvemq(j-1, kernelp, mqarray, pr0)
                            pr0 <- pr0 * pjmat[j,]
                        }
                    }
                    pdotbd <- 1 - sum(pr0)
                }
            }
            if (type==6) pbd <- 1 else pbd <- beta[b]
            if (b<d) pbd <- pbd * prod(phij[b:(d-1)])
            if (li[n]>0)    # not censored
                pbd <- pbd * (1-phij[d])
            prwi <- 1
            if ((b+cjs) <= d) {
                if (movemodel == 0) {
                    alpha <- rep(1/mm,mm)
                    for (j in (b + cjs):d) {
                        alpha <- prwmulti (n, j, x, cumss, w, PIA, hk, Tsk, h, p0, hindex, alpha)
                    }
                    prwi <- sum(alpha) / pdotbd
                }
                else if (movemodel == 1) {
                    prwi <- 1.0
                    for (j in (b+cjs):d) {
                        alpha <- prwmulti (n, j, x, cumss, w, PIA, hk, Tsk, h, p0, hindex, rep(1/mm,mm))
                        prwi <- prwi * sum(alpha);   # product over primary sessions
                    }
                    prwi <- prwi / pdotbd
                }
                else { # movemodel>1
                    alpha <- prwmulti (n, b+cjs, x, cumss, w, PIA, hk, Tsk, h, p0, hindex, rep(1/mm,mm))
                    if (d>(b+cjs)) {
                        for (j in (b+cjs+1):d) {
                            alpha <- convolvemq(j-1, kernelp, mqarray, alpha)
                            alpha <- prwmulti (n, j, x, cumss, w, PIA, hk, Tsk, h, p0, hindex, alpha)
                        }
                    }
                    prwi <- sum(alpha) / pdotbd
                }
                pdt <- pdt + pbd * prwi
            }
        }
    }
    pdt
}

#==============================================================================
