## R functions to avoid C for prwi (spatial)

## prwisecr

## 2018-02-12, 2019-05-19,
## 2019-06-19 removed nonspatial to prwi.R
## 2019-06-19 merged multicatch into general prwisecr

## openval has
## column 1 lambda0
## column 2 phi
## column 3 f, g, l  (recruitment parameter)
## column 4 sigma
## column 5 pmix

#########################################################################################

# Probability of count for session s, detector k, animal i
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

prw <- function (n, j, x, kk, binomN, cumss, w, PIA, Tsk, gk, h, p0, hindex, pjm) {
    # gk assumed to be hazard if multi (binomN=-2) or Poisson count (binomN=0)
    dead <- FALSE
    for (s in (cumss[j]+1):cumss[j+1]) {

        if (binomN == -2) {      ## multi-catch traps, 2-D w
            wi <- w[n, s]
            if (wi < 0) dead <- TRUE
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
                    pjm <- pjm * Tsk[k,s] * (1-p0[x,,hindex[n,s]]) *  gk[c, k, ] / h[x,,hindex[n,s]]
                }
            }
        }
        else {
            for (k in 1:kk) {
                c <- PIA[n,s,k,x]
                if (c >= 1) {    # drops unset traps
                    count <- w[n,s,k]
                    if (count<0) {count <- -count; dead <- TRUE }
                    pjm <- pjm * pski(binomN,count,Tsk[k,s], gk[c,k,])
                }
            }
        }

        if (dead) break;   # after processing all traps on this occasion
    }
    ## cat("j ", j, " pjm ", sum(pjm), '\n')
    pjm
}
###################################################################################

prwisecr <- function (type, n, x, nc, jj, kk, mm, nmix, cumss, w, fi, li, gk,
                      openval, PIA, PIAJ, binomN, Tsk, intervals, h, hindex,
                      CJSp1, calcpdotbd, moveargsi, movemodel,
                      usermodel, kernel = NULL, mqarray = NULL, cellsize = NULL) {

    ## precompute p0 to save time (multicatch only)
    p0 <- if (binomN == -2) exp(-h) else 1

    phij <- getphij (n, x, openval, PIAJ, intervals)

    if (movemodel > 1) {
        moveargsi <- pmax(moveargsi,0)
        moveargs <- getmoveargs (n, x, openval, PIAJ, intervals, moveargsi)
        kernelp <- fillkernelp (jj, movemodel-2, kernel, usermodel, cellsize,
                                moveargsi, moveargs, normalize = TRUE)
    }

    if(type==6) {
        minb <- fi[n]
        cjs <- 1 - CJSp1
    }
    else {
        minb <- 1
        cjs <- 0
        beta <- getbeta (type, n, x, openval, PIAJ, intervals, phij)
        if (calcpdotbd) {
            ## Precompute session-specific Pr for all unique parameter combinations 2019-05-05
            ## Does not allow for learned response
            pjmat <- pr0njmx(n, x, cumss, jj, mm, binomN, PIA, gk, Tsk)  ## PCH1.R
        }
    }
    maxb <- fi[n]
    mind <- abs(li[n])
    maxd <- jj
    if (li[n] < 0) maxd <- mind     # possible censoring
    pdt <- 0
    pdotbd <- 1.0
    for (b in minb:maxb) {
        for (d in mind:maxd) {
            if (type==6) {     ## CJS
                pbd <- 1
            }
            else {
                pbd <- beta[b]
                pdotbd <- 1
                if (calcpdotbd) {
                    # compute pdot for this b,d
                    if (movemodel == 1) {   ## uncorrelated; product over primary sessions
                        if (d>(b+cjs)) {
                            pdotbd <- 1 - prod(apply(pjmat[(b+cjs):d,]/mm,2,sum))
                        }
                        else {
                            pdotbd <- 0
                        }
                    }
                    else {
                        if (d>=(b+cjs)) {
                            pr0 <- pjmat[b+cjs,] / mm
                            if (d>(b+cjs)) {
                                for (j in (b+cjs+1):d) {
                                    if (movemodel>1) pr0 <- convolvemq(j-1, kernelp, mqarray, pr0)
                                    pr0 <- pr0 * pjmat[j,]
                                }
                            }
                            pdotbd <- 1-sum(pr0)
                        }
                        else {
                            pdotbd <- 0
                        }
                    }
                }
            }
            if (b<d) pbd <- pbd * prod(phij[b:(d-1)])
            if ((li[n]>0) & (d<jj))    # if not censored, died
                pbd <- pbd * (1-phij[d])

            prwi <- 1.0
            if (d >= (b+cjs)) {
                if (movemodel == 0) {

                    alpha <- rep(1.0/mm,mm)
                    for (j in (b+cjs):d) {
                        alpha <- prw(n, j, x, kk, binomN, cumss, w, PIA, Tsk, gk, h, p0, hindex, alpha)
                    }
                    prwi <- sum(alpha)
                }
                else if ( movemodel == 1) { # uncorrelated; product over primary sessions
                    prwi <- 1.0
                    for (j in (b+cjs):d) {
                        alpha <- rep(1.0/mm,mm)
                        alpha <- prw(n, j, x, kk, binomN, cumss, w, PIA, Tsk, gk, h, p0, hindex, alpha)
                        prwi <- prwi * sum(alpha)
                    }
                }
                else { # movemodel>1
                    alpha <- rep(1.0/mm, mm)
                    alpha <- prw (n, b+cjs, x, kk, binomN, cumss, PIA, gk, Tsk, w, alpha)
                    if (d>(b+cjs)) {
                        for (j in (b+cjs+1):d) {
                            alpha <- convolvemq(j-1, kernelp, mqarray, alpha)
                            alpha <- prw(n, j, x, kk, binomN, cumss, w, PIA, Tsk, gk, h, p0, hindex, alpha)
                        }
                    }
                    prwi <- sum(alpha)
                }
                # cat("n ", n, " b ", b, " d ", d, " pbd ", pbd, " prwi ", prwi, " pdotbd ", pdotbd, "\n")
            }
            pdt <- pdt + pbd * prwi / pdotbd
        }
    }
    # cat("n ", n, " pdt ", pdt, "\n")
    pdt
}

#==============================================================================

