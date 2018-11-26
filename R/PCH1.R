###############################################################################
## openCR
## pch0.R
## 2018-02-26 openCR 1.0.0
###############################################################################

##-----------------------------------------------------------
## return JSSA probability animal n detected at least once   
## Pledger et al. 2010 Eqn (3) (mixtures outside)            
##-----------------------------------------------------------
PCH1 <- function (type, x, nc, cumss, nmix, openval0, PIA0, PIAJ, intervals) {
    J <- length(intervals)+1
    one <- function(n) {    
        p <- getp (n, x, openval0, PIA0)
        phij <- getphij (n, x, openval0, PIAJ, intervals)
        # getg (type, n, x, nc, jj, openval0, cc0, PIAJ, g);
        beta <- getbeta (type, n, x, openval0, PIAJ, intervals, phij)
        pdt <- 0
        for (b in 1:J) {
            for (d in b:J) {
                pbd <- beta[b]                      ## entered at b
                if (d>b)
                    pbd <- pbd * prod (phij[b:(d-1)]) ## survived
                pbd <- pbd * (1-phij[d]);             ## departed at d
                ptmp <- 1
                for (j in b:d) {
                    s <- (cumss[j]+1) : cumss[j+1]
                    ptmp <- ptmp * prod(1 - p[s])       ## not detected
                }
                pdt = pdt + pbd * (1 - ptmp)
            }
        }
        pdt
    }
    sapply(1:nc, one)   # or rep(one(1), nc) if all the same
}

#==============================================================================
# 
# prepare matrix n x j x m of session-specific Pr
pr0njmx <- function (n, x, cumss, jj, kk, mm, binomN, PIA0, gk0, Tsk) {
    pjm <- array(1, dim = c(jj, mm))
    for (j in 1:jj) {
        # browser()        
        # s is vector of indices to secondary sessions in this primary session
        s <- (cumss[j]+1) : cumss[j+1]
        S <- length(s)
        csk <- PIA0[n,s, ,x, drop = FALSE]
        cski <- rep(as.numeric(csk),mm)
        i <- cbind(cski, rep(rep(1:kk, each=S), mm), rep(1:mm, each=S*kk))
        gsk <- array(0, dim=c(S, kk, mm))
        gsk[cski>0] <- gk0[i]
        size <- t(Tsk[,s])      
        ## if (binomN != 1) size[] <- binomN
        pjm[j, ] <- if (binomN == 0) 
            apply(1-gsk, 3, function(x) prod(exp(-size * -log(x))))  ## Poisson
        else if (all(size==1))
            apply(1-gsk,3, prod)
        else {
            apply(1-gsk, 3, function(x) prod (x^size))  ## Binomial or Bernoulli
        }
    }
    pjm
}

PCH1secr <- function (type, individual, x, nc, jj, cumss, kk, mm, openval0, PIA0, PIAJ, 
                      gk0, binomN, Tsk, intervals, moveargsi, movemodel, usermodel,
                      kernel, mqarray, cellsize) {
    
    One <- function (n) {
        ## precompute session-specific Pr for all unique parameter combinations
        pjm <- pr0njmx(n, x, cumss, jj, kk, mm, binomN, PIA0, gk0, Tsk)
        phij <- getphij (n, x, openval0, PIAJ, intervals)
        beta <- getbeta (type, n, x, openval0, PIAJ, intervals, phij)
        if (movemodel>1) {
            moveargsi <- pmax(moveargsi,0)
            moveargs <- getmoveargs (type, n, x, openval0, PIAJ, intervals, moveargsi)
            kernelp <- fillkernelp (jj, movemodel-2, kernel, usermodel, cellsize, 
                                    moveargsi, moveargs)
        }
        pdt <- 0
        for (b in 1:jj) {
            for (d in b:jj) {
                pbd <- beta[b]                        ## entered at b
                if (d>b)
                    pbd <- pbd * prod (phij[b:(d-1)]) ## survived
                pbd <- pbd * (1-phij[d]);             ## departed at d
                
                # static home ranges: take sum over M of product over J
                if (movemodel==0) {
                    prodpj <- apply(pjm[b:d,, drop = FALSE], 2, prod)
                    prw0 <- sum(prodpj) / mm
                }
                else if (movemodel==1) {
                    ## over primary sessions in which may have been alive
                    ## centers allowed to differ between primary sessions
                    ## mobile home ranges: take product over J of sum over M
                    sumpm <- apply(pjm[,,drop=FALSE], 1, sum )
                    prw0 <- prod(sumpm[b:d]) / mm
                }
                else {  ## normal, exponential etc.
                    # use dummy (w=0) capthist until can re-write pr0njmx 2018-03-14
                    pjm <- rep(1,mm)
                    for (j in b:d) {
                        pjm <- prw (n, j, x, nc, jj, kk, binomN, cumss, PIA0,
                                     gk0, Tsk, w, pjm)
                        ## replace pjm with version updated for movement 
                        ## pi_1(X|omega_i1) = Pr(omega_i1|X_i1) pi_0(X) / Pr (omega_i1) 
                        ## the denominator Pr(omega_i) is omitted as it is constant across the mask 
                        if (j<d) pjm <- convolvemq(j, kernelp, mqarray, pjm)
                    }
                    prw0 <- sum(pjm) / mm
                }
                pdt <- pdt + pbd * (1 - prw0)
            }
        }
        pdt
    }
    if (movemodel>1) w <- array(0, dim=c(nc,cumss[jj+1],kk)) else w <- NA
    if (individual) {
        sapply(1:nc, One)
    }
    else {
        rep(One(1), nc)
    }
}

