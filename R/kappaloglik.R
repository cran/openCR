dbinomraw <- function(x, n, p, q)
{
    lgamma(n+1) - lgamma (x+1) - lgamma(n-x+1) + x*log(p) + (n-x) * log(q)
}

kappaloglik <- function (type, openval,  PIA, PIAJ, data) {
    
    onehistory <- function (n) {
        prw <- prwi(
            1,   # type = 'CJS'
            1,   # n
            1,
            data$J,
            data$cumss,
            data$details$nmix,
            data$capthist[n,, drop = FALSE],
            data$fi[n],
            data$li[n],
            openval,
            PIA[n,,,,drop = FALSE],
            PIAJ[n,,,drop = FALSE],
            data$intervals,
            data$details$CJSp1)
        if (prw<=0) {
            -1e10
        }
        else
            freq[n] * log(prw)
    }
    
    J <- length(data$intervals) + 1
    ch <- data$capthist
    w <- matrix(data$JScounts, nrow=J)
    ni <- w[,1]       # number viewed at i
    u <- ni - w[,3]   # number viewed for first time at i
    n <- sum(u)
    p <- getpj   (1, 1, openval, PIAJ)
    phij <- getphij (1, 1, openval, PIAJ, data$intervals)
    kapj <- getkapj (1, 1, openval, PIAJ)
    beta <- getbetak (1, 1, openval, PIAJ, phij, data$intervals)  ## needed for beta0
    xi <- kapj/sum(kapj)
    ## f0
    if (type == 28) {
        ## superN <- openval[PIAJ[1,1,1], 4]
        pi <- beta[1] * p[1] * sum(kapj)
        ## distribution 0 = Poisson, 1 = binomial
        if (data$distrib == 1) {
            superN <- openval[PIAJ[1,1,1], 4] + n
            ## superN <- openval[PIAJ[1,1,1], 4]
            if (n < superN) {
                f0 <- dbinomraw (n, superN, pi, 1-pi)
                if (is.na(f0)) browser()
            }
            else {
                 f0 <- -1e10
            }
        }
        else {
            superN <- openval[PIAJ[1,1,1], 4]
            f0 <- dpois (n, pi * superN, log = TRUE)
        }
    }
    else f0 <- 0   ## CL
    
    ## f1
    f1 <- lgamma(n+1) - sum(lgamma(u+1)) + sum(u * log(xi))

    ## f2
    f2 <- 0   ## ignore losses
    
    ## f3
    freq <- covariates(data$capthist)$freq
    f3 <- sum(sapply(1:data$nc, onehistory))
    return (c(f0, f1, f2, f3))   # return log likelihood components
    
}
