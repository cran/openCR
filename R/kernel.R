################################################################################
## package 'openCR'
## kernel.R
## 2019-06-09, 11
## 2020-10-06 tweak to ed, bR
## 2021-02-23 annular etc.
## 2021-02-25 functions not exported: empirical
## 2021-03-01 functions exported pkernel, matchscale
## 2021-03-16 annular2
## 2021-07-18 many changes - gkernel(), frE, frG, frL
## 2021-07-29 BVNzi, BVEzi, uniformzi, frEzi  (11,12,13,14)

################################################################################

## fillkernelp is used in prwisecr.R
## otherwise, these functions are used for exploration of kernel characteristics
## kernel is coded independently within likelihood functions

circleintersectsline <- function (x1, x2, y1, y2, R) {
    # from https://mathworld.wolfram.com/Circle-LineIntersection.html
    # assume origin at 0,0
    dx <- x2-x1
    dy <- y2-y1
    dr2 <- dx^2 + dy^2
    D <- x1*y2 - x2 * y1
    inc <- R^2*dr2 - D^2
    sgndy <- if(dy<0) -1 else 1
    if (inc > 0 ) {
        data.frame(
            x = c(
                (D*dy+sgndy*dx*sqrt(inc)) / dr2,
                (D*dy-sgndy*dx*sqrt(inc)) / dr2
            ),
            y = c(
                (D*dx + abs(dy)*sqrt(inc)) / dr2,
                (D*dx - abs(dy)*sqrt(inc)) / dr2
            )
            
        )
    }
    else {
        NULL
    }
}
#-----------------------------------------------------------------------
## fill cells of kernel with probability of movement from central point

fillkernelC <- function (
    J,              ## number of sessions (J-1 potential movements)
    kerneltype,     ## 0 BVN, 1 BVE, 2 usermodel, 3 BVT, etc.
    sparsekernel,   ## TRUE iff sparse
    kernel,
    usermodel = "",
    cellsize,
    moveargsi,
    moveargs,       ## J x npar matrix
    normalize)
{
    kernel <- as.matrix(kernel)   # convert from mask dataframe to matrix
    kernel <- sweep(kernel, FUN="-", STATS=apply(kernel,2,mean), MARGIN=2)
    kernel[,] <- as.integer(round(kernel / cellsize))
    kernelp <- fillkernelcpp (
       as.matrix   (kernel), 
       as.integer  (kerneltype), 
       as.logical  (sparsekernel),
       as.double   (cellsize),
       as.integer  (J),
       as.character(usermodel),
       as.integer  (moveargsi), 
       as.double   (moveargs),
       as.logical  (normalize)    
    )   
    matrix(kernelp, ncol = J-1)
}

make.kernel <- function (
    movementmodel = c('BVN', 'BVE', 'BVT','frE', 'frG', 'frL', 'uniform'),
    kernelradius = 10, 
    spacing, 
    move.a, 
    move.b, 
    sparsekernel = FALSE, 
    clip = FALSE,
    normalize = TRUE,
    stat = c('estimate','lcl', 'ucl'),
    ...) 
{
    if (inherits(movementmodel, 'openCR')) {
        fit <- movementmodel
        if (fit$version < '2.0.0') stop ("model fitted with openCR < 2.0.0")
        if (fit$movementmodel %in% .openCRstuff$movementmodels) {
            stat <- match.arg(stat)
            pred <- predict(fit, ...)
            kernel <- make.kernel(
                movementmodel = fit$movementmodel, 
                kernelradius  = fit$kernelradius, 
                spacing       = secr::spacing(fit$mask), 
                move.a        = pred$move.a[1, stat],   # gather lcl, ucl as well 2021-08-08
                move.b        = pred$move.b[1, 'estimate'],   # gather lcl, ucl as well 2021-08-08
                sparsekernel  = fit$sparsekernel, 
                clip          = TRUE,
                normalize     = TRUE)
            if (!is.null(fit$kernel) && nrow(kernel) != nrow (fit$kernel)) 
                stop ('bad match to kernel')
            kernel
        }
        else {
            NULL
        }
    }
    else {
        if (is.function (movementmodel)) {
            moveargs <- formalArgs(movementmodel)
            usermodel <- as.character(substitute(movementmodel))
            movementmodel <- "user"
        }
        else {
            usermodel <- ""   ## changed from NULL 2021-07-27
            movementmodel <- match.arg(movementmodel[1],
                choices = .openCRstuff$movementmodels) 
        }

        if (missing(move.a) && movementmodel != "uniform" && movementmodel %in% .openCRstuff$movementmodels) { 
            stop ("move.a required for movementmodel ", movementmodel)
        }
        if (missing(move.b) && movementmodel %in% c('t2D','annular2','annularR', 
            'BVT','frL','frG', 'BVNzi', 'BVEzi','frEzi')) {
            stop ("move.b required for movementmodel ", movementmodel)
        }
        movementcode <- movecode(movementmodel)
        moveargsi <- c(0,0)
        if (missing(move.a) | (movementmodel == 'uniform')) {
            pars <- move.a <- move.b <- NULL
        }
        else {
            pars <- move.a[1]
            moveargsi[1] <- 1
        }
        if (missing(move.b)) {
            move.b <- NULL
        }
        else {
            pars <- c(pars, move.b[1])
            moveargsi[2] <- 2
        }
        if (is.null(pars)) pars <- c(0,0)
        moveargs <- matrix(pars, nrow = 1)   # J-row matrix for fillkernelp
        k2 <- kernelradius
        kernel <- make.mask(type = 'rectangular',
            spacing = spacing,
            buffer = 0, nx = 2 * k2+1, ny = 2 * k2+1)
        ## centre
        kernel[,] <- sweep(kernel, MARGIN=2, FUN = "-", STATS = rep((k2+0.5)*spacing,2))
        
        if (movementmodel %in% c('annular', 'annular2')) {
            r <- (kernel$x^2 + kernel$y^2)^0.5
            origin <- r==0
            ring2 <- (r >= (k2-0.5) * spacing) & (r<(k2+0.5) * spacing)
            if (movementmodel == 'annular') {
                ok <- (origin | ring2)
            }
            else {
                ring1 <- (r >= (k2/2-0.5) * spacing) & (r<(k2/2+0.5) * spacing)
                ok <- (origin | ring1 | ring2)
            }
            kernel <- subset(kernel, ok)
        }
        
        if (sparsekernel) {
            ok <- kernel$x==0 | kernel$y == 0 | kernel$x == kernel$y | kernel$x == -kernel$y
            kernel <- kernel[ok,]
        }
        
        ## optional clipping 
        if (clip) {
            outside <- (kernel$x^2 + kernel$y^2) > ((k2+0.5)*spacing)^2
            kernel <- subset(kernel, !outside)
        }
        
        # call wrapper for C++ function fillkernelcpp 
        kernelp <- fillkernelC (
            2,              ## number of sessions+1
            movementcode-2, ## kerneltype 0 BVN, 1 BVE, 2 usermodel, 3 BVT, 4 uniform, 
                            ## 5 annular, 6 annular2, 7 annularR,
                            ## 8 frE, 9 frG, 10 frL,
                            ## 11 BVNzi, 12 BVEzi, 13 uniformzi, 14 frEzi 
            sparsekernel,   ## TRUE iff sparse
            kernel,         ## coordinates  
            usermodel,      ## name of user function (grain = 0 only) 
            cellsize = spacing,      
            moveargsi,      ## which column has move.a, move.b
            t(matrix(moveargs,2,2)),
            normalize = normalize   ## delay normalization?
        )[,1]

        covariates(kernel) <- data.frame(kernelp = kernelp)
        
        ## 2021-03-01 cumulative distribution function
        r <- apply(kernel^2,1,sum)^0.5
        p <- covariates(kernel)$kernelp[order(r)]
        attr(kernel, 'distribution') <- data.frame(r = sort(r), cumprob = cumsum(p))
        
        attr(kernel, 'movementmodel') <- movementmodel
        attr(kernel, 'sparsekernel') <- sparsekernel
        attr(kernel, 'k2') <- k2
        attr(kernel, 'move.a') <- move.a
        attr(kernel, 'move.b') <- move.b
        class(kernel) <- c('kernel','mask','data.frame')
        kernel
    }
}

getpstring <- function (movementmodel, move.a, move.b, sep = ' ') {
    npar <- nparmove(movementmodel)
    pstring <- ''
    if (npar>0) pstring <- paste0('move.a = ', move.a)
    if (npar>1) pstring <- paste(c(pstring, paste0(' move.b = ', move.b)), collapse = sep)
    pstring
}

plot.kernel <- function (x, contour = FALSE, levels = NULL, text = FALSE, 
    title = NULL, ...) {
    spacing <- spacing(x)
    k2 <- attr(x, 'k2')
    move.a <- attr(x, 'move.a')
    move.b <- attr(x, 'move.b')
    movementmodel <- attr(x, 'movementmodel')
    npar <- nparmove(movementmodel)
    pstring <- getpstring(movementmodel, round(move.a,3), round(move.b,3))
    kernelp <- covariates(x)$kernelp
    dots <- list(...)
    if ('border' %in% names(dots)) {
        border <- dots$border
        dots$border <- NULL
    }
    else border <- 1
    msk <- x
    class(msk) <- c('mask','data.frame')

    if (sum(!is.na(kernelp))==0) {
        warning ("no non-missing kernel probabilities")
        plot(msk, dots = FALSE, meshcol = 'white', border = border, ...)
    }
    else {
        plot(msk, dots = FALSE, meshcol = 'white', covariate = 'kernelp', border = border, ...)
    }
        
    centrecell <- subset(msk, (x$x^2 + x$y^2)<1e-6)
    # 2021-03-17 modify plot call to allow add in ...
    dots$x <- centrecell
    dots$dots <- FALSE
    dots$meshcol <- 'black'
    dots$col <- NA
    dots$add = TRUE
    do.call(plot, dots)

    if (movementmodel %in% c('annular','annular2')) {
        rad <- k2 * spacing
        symbols(0, 0, circles = rad, inches = FALSE, add = TRUE)
        if (movementmodel == 'annular2') {
            rad <- k2 * spacing
            symbols(0, 0, circles = rad, inches = FALSE, add = TRUE)
        }
    }
    
    if (movementmodel %in% c('annularR')) {
        rad <- move.b * k2 * spacing
        symbols(0, 0, circles = rad, inches = FALSE, add = TRUE)
    }
    
    if (contour) {
        kp <- kernelp
        z <- matrix(nrow = 2*k2+1, ncol = 2*k2+1)
        kxy <- as.matrix(x) / spacing + k2 + 1
        z[kxy] <- kp
        if (is.null(levels))
            levels <- pretty(c(0, max(kp, na.rm = TRUE)), 10)
        contour(add = TRUE, (-k2:k2)*spacing, (-k2:k2) * spacing, z, levels = levels)
    }
    if (text) text(x$x, x$y, round(kernelp,3), cex=0.6)
    if (is.null(title)) {
        title <-  paste0('kernel = ', movementmodel, ', spacing = ', spacing, 
            ', kernelradius = ', k2, ', ', pstring, ', ncells = ', nrow(x))
    }
    mtext (side=3, line = 1, title, cex = par()$cex.main)

    # invisible(data.frame(x, kernelp=kernelp[!is.na(kernelp)]))  suppress 2021-07-28
    invisible(data.frame(x, kernelp=kernelp))
}

ed <- function(movementmodel, move.a, move.b, radius) {
    if (movementmodel=='user') return (NA)  # 2020-10-06
    movementmodel <- match.arg(movementmodel[1], choices = 
            .openCRstuff$movementmodels)
    if (movementmodel %in% c('normal', 'BVN'))
        move.a * (pi/2)^0.5         ## Nathan et al 2012 Table 15.1
    else if (movementmodel %in% c('exponential', 'BVE'))
        move.a * 2              ## Nathan et al 2012 Table 15.1
    else if (movementmodel %in% c('t2D', 'BVT')) {
        if (missing(move.b)) stop ("t2D model has two parameters")
        a <- move.a
        b <- move.b + 1
        if (b<3/2)
            Inf
        else
            a * pi^0.5 / 2 * exp(lgamma(b-3/2) - lgamma(b-1)) 
    }
    else if (movementmodel == 'annular') {
        (1-move.a) * radius
    }
    else if (movementmodel == 'annular2') {
        move.b * radius/2 + (1-move.a-move.b) * radius
    }
    else if (movementmodel == 'annularR') {
        (1-move.a) * move.b
    }
    else if (movementmodel == 'frE') {
        move.a
    }
    else if (movementmodel == 'frG') {
        move.a * move.b
    }
    else if (movementmodel == 'frL') {
        mu <- log(move.a)
        s <- sqrt(log(1 + 1/move.b))
        exp(mu + s^2/2)         ## Cousens et al 2008 Table 5.1; 
                                ## cf Nathan et al 2012 Table 15.1 use a = exp(mu)
    }
    else if (movementmodel == 'BVNzi') {
        move.a * (pi/2)^0.5 * (1-move.b)
    }
    else if (movementmodel == 'BVEzi') {
        move.a * 2 * (1-move.b)
    }
    else if (movementmodel == 'uniformzi') {
        NA * (1-move.b)
    }
    else if (movementmodel == 'frEzi') {
        move.a * (1-move.b)
    }
    else NA
}

bR <- function(R, movementmodel, move.a, move.b) {
    if (movementmodel=='user') return (NA)  # 2020-10-06
    1 - pkernel(R, movementmodel, move.a, move.b)
}

summary.kernel <- function (object, ...) {
    spacing <- spacing(object)
    k2 <- attr(object, 'k2')
    r <- sqrt(apply(object^2,1,sum))
    movementmodel <- attr(object, 'movementmodel')
    move.a <- attr(object, 'move.a')
    move.b <- attr(object, 'move.b')
    kernelp <- covariates(object)$kernelp
    kernelp <- kernelp / sum(kernelp)   # force normalization for E(r)
    
    result <- list(k2 = k2,
        spacing = spacing,
        ncells = nrow(object),
        movementmodel = movementmodel,
        move.a = move.a[1],
        move.b = move.b[1],
        expectedmove = ed(movementmodel, move.a, move.b, k2*spacing),
        expectedmovetr = sum(r * kernelp),
        ptruncated = bR((k2+0.5)*spacing, movementmodel, move.a, move.b),
        expectedq50 = qkernel(0.5, movementmodel, move.a, move.b),
        expectedq90 = qkernel(0.9, movementmodel, move.a, move.b),
        expectedq50tr = qkernel(0.5, movementmodel, move.a, move.b, truncate = k2*spacing),
        expectedq90tr = qkernel(0.9, movementmodel, move.a, move.b, truncate = k2*spacing)
    )
    class(result) <- 'summary.kernel'
    result
}

print.summary.kernel <- function(x,...) {
    cat('Kernel radius (cells)     : ', x$k2, '\n')
    cat('Spacing (side of cell)    : ', x$spacing, ' (m)\n')
    cat('Number of cells           : ', x$ncells, '\n')
    if (is.function(x$movementmodel))
        cat('Movement model            : user function\n')
    else
        cat('Movement model            : ', x$movementmodel, '\n')
    cat('Parameter(s)              : ', getpstring(x$movementmodel, x$move.a, x$move.b, sep=','), '\n')
    cat('Proportion truncated      : ', x$ptruncated, '\n')
    cat('Movement as truncated at edge of kernel\n')
    cat('Expected distance         : ', x$expectedmovetr, ' (m)\n')
    cat('50th percentile (median)  : ', x$expectedq50tr, ' (m)\n')
    cat('90th percentile           : ', x$expectedq90tr, ' (m)\n')
    cat('Movement, untruncated kernel\n')
    cat('Expected distance         : ', x$expectedmove, ' (m)\n')
    cat('50th percentile (median)  : ', x$expectedq50, ' (m)\n')
    cat('90th percentile           : ', x$expectedq90, ' (m)\n')
}

# useful functions not exported 2021-02-25
# see kerneltest.R

empirical <- function (kernel) {
    r <- apply(kernel^2,1,sum)^0.5
    p <- covariates(kernel)$kernelp
    p <- p[order(r)]
    r <- sort(r)
    out <- cbind(r, cumsum(p))
}

gkernel <- function (r, movementmodel = c('BVN', 'BVE', 'BVT','frE', 'frG', 'frL'),
    move.a, move.b, truncate = Inf) {
    movementmodel <- try(match.arg(movementmodel, choices = 
            .openCRstuff$movementmodels), silent = TRUE)
    if (inherits(movementmodel, 'try-error')) {
        warning("gkernel defined only for movement models BVN, BVE, BVT, frE, frG, frL")
        rep(NA, length(r))
    }
    else if (movementmodel %in% c('BVNzi', 'BVEzi','uniformzi','frEzi')) {
        warning("gkernel not defined for zero-inflated models BVNzi, BVEzi, uniformzi, frEzi")
        rep(NA, length(r))
    }
    else {
        if (movementmodel %in% c('normal', 'BVN')) {
            g <- 1/move.a^2/2/pi * exp(-r^2/2/move.a^2)
        }    
        else if (movementmodel %in% c('exponential', 'BVE')) {
            g <- 1/move.a^2/2/pi * exp(-r/move.a)
        }
        else if (movementmodel %in% c('t2D', 'BVT')) {
            g <- move.b / pi/ move.a^2 * (1 + r^2/move.a^2)^-(move.b+1)
        }
        else if (movementmodel %in% c('frE')) {
            g <- exp(-r/move.a) / move.a / pi / 2 / r
        }
        else if (movementmodel %in% c('frG')) {
            g <- exp(-r/move.a) / gamma(move.b) / move.a^move.b / pi / 2 * r^(move.b-2)
        }
        else if (movementmodel %in% c('frL')) {
            mu <- log(move.a)
            s <- sqrt(log(1 + 1/move.b))
            g <- dlnorm(r, mu, s) / 2 / pi / r 
        }
        else if (movementmodel %in% c('uniform')) {
            g <- rep(1,length(r)) / truncate^2 /pi
        }
        else {
            g <- rep(NA, length(r))
            warning('pdf not available for movement kernel "', 
                movementmodel, '"')
        }
        if (truncate != Inf) {
            if (movementmodel == 'uniform')
                ptrunc <- 1
            else
                ptrunc <- pkernel(truncate, movementmodel, move.a, move.b, Inf, TRUE)
            g <- g / ptrunc
            g[r>truncate] <- 0
        }
        g
    }
}

# Distribution function for distance moved
dkernel <- function (r, movementmodel = c('BVN', 'BVE', 'BVT','frE', 'frG', 'frL'),
    move.a, move.b, truncate = Inf) {
    movementmodel <- try(match.arg(movementmodel, choices = 
            .openCRstuff$movementmodels), silent = TRUE)
    if (inherits(movementmodel, 'try-error')) {
        warning("dkernel defined only for movement models BVN, BVE, BVT, frE, frG, frL, frZ, uniform")
        rep(NA, length(r))
    }
    else if (movementmodel %in% c('BVNzi', 'BVEzi','uniformzi','frEzi')) {
        warning("dkernel not defined for zero-inflated models BVNzi, BVEzi, uniformzi, frEzi")
        rep(NA, length(r))
    }
    else {
        if (movementmodel %in% c('normal', 'BVN')) {
            d <- r/move.a^2 * exp(-r^2/2/move.a^2)
        }    
        else if (movementmodel %in% c('exponential', 'BVE')) {
            d <- r/move.a^2 * exp(-r/move.a)
        }
        else if (movementmodel %in% c('t2D', 'BVT')) {
            d <- 2 * move.b * r / move.a^2 * (1 + r^2/move.a^2)^-(move.b+1)
        }
        else if (movementmodel %in% c('frE')) {
            d <- exp(-r/move.a) / move.a 
        }
        else if (movementmodel %in% c('frG')) {
            d <- exp(-r/move.a) / gamma(move.b) / move.a^move.b * r^(move.b-1)
        }
        else if (movementmodel %in% c('frL')) {
            mu <- log(move.a)
            s <- sqrt(log(1 + 1/move.b))
            d <- dlnorm(r, mu, s) 
        }
        else if (movementmodel %in% c('uniform')) {
            A <- pi * truncate^2
            d <- 2 * pi * r / A
        }
        else {
            d <- rep(NA, length(r))
            warning('pdf not available for movement kernel "', 
                movementmodel, '"')
        }
        if (truncate != Inf) {
            if (movementmodel == 'uniform')
                ptrunc <- 1
            else
                ptrunc <- pkernel(truncate, movementmodel, move.a, move.b, Inf, TRUE)
            d <- d / ptrunc
            d[r>truncate] <- 0
        }
        d
    }
}

# Distribution function for distance moved
pkernel <- function (q, movementmodel = c('BVN', 'BVE', 'BVT','frE', 'frG', 'frL'),
    move.a, move.b, truncate = Inf, lower.tail = TRUE) {
    movementmodel <- try(match.arg(movementmodel, choices = 
            .openCRstuff$movementmodels), silent = TRUE)
    zeroinflated <- grepl('zi', movementmodel)
    if (zeroinflated) movementmodel <- gsub("zi","",movementmodel)
    if (inherits(movementmodel, 'try-error')) {
        warning("pkernel defined only for movement models BVN, BVE, BVT, frE, frG, frL, uniform")
        rep(NA, length(q))
    }
    else {
        if (movementmodel %in% c('normal', 'BVN')) {
            p <- exp(-q^2/2/move.a^2) 
        }    
        else if (movementmodel %in% c('exponential', 'BVE')) {
            p <- (q/move.a+1) * exp(-q/move.a)
        }
        else if (movementmodel %in% c('t2D', 'BVT')) {
            p <- (move.a^2/(move.a^2+q^2))^move.b
        }
        else if (movementmodel %in% c('frE')) {
            p <- exp(-q/move.a)
        }
        else if (movementmodel %in% c('frG')) {
            p <- 1 - pgamma(q, move.b, scale = move.a)
        }
        else if (movementmodel %in% c('frL')) {
            mu <- log(move.a)
            s <- sqrt(log(1 + 1/move.b))
            p <- 1 - pnorm((log(q) - mu) / s)
        }
        else if (movementmodel %in% c('uniform')) {
            p <- ifelse (q>truncate, 0, 1 - (q/truncate)^2)
        }
        else {
            p <- rep(NA, length(q))
            warning('probability not available for movement kernel "', 
                movementmodel, '"')
        }
        if (lower.tail) p <- 1-p
        if (truncate != Inf) {
            if (!lower.tail) stop ("truncation incompatible with !lower.tail")
            p[q>truncate] <- 1.0
            if (movementmodel != 'uniform') {
                p <- p / pkernel(truncate, movementmodel, move.a, move.b, Inf, TRUE)
            }
        }
        if (zeroinflated) {
            bz <- if (movementmodel %in% c('uniform')) move.a else move.b
            if (bz>1) stop ("requested zero-inflation outside range 0-1")
            if (lower.tail)
                p <- bz + (1-bz) * p
            else
                p <- 1 - (bz + (1-bz)*(1-p))
        }
        p
    }
}

qkernel <- function(p, movementmodel = c('BVN', 'BVE', 'BVT','frE', 'frG', 'frL'),
    move.a, move.b, truncate = Inf, lower.tail = TRUE) {
    movementmodel <- try(match.arg(movementmodel, choices = 
            .openCRstuff$movementmodels), silent = TRUE)
    zeroinflated <- grepl('zi', movementmodel)
    if (zeroinflated) movementmodel <- gsub("zi","",movementmodel)
    if (inherits(movementmodel, 'try-error')) {
        warning("qkernel defined only for movement models BVN, BVE, BVT, frE, frG, frL, uniform")
        rep(NA, length(p))
    }
    else {
        if (zeroinflated) {
            bz <- if (movementmodel %in% c('uniform')) move.a else move.b
            p0 <- p
            if (lower.tail)
                p <- (p-bz)/(1-bz)
            else
                p <- 1 - ((1-p)-bz)/(1-bz)
        }
        else {
            bz <- 0
            p0 <- 1
        }
        
        if (truncate != Inf) {
            if (!lower.tail) stop ("cannot combine truncation and upper tail")
            ptrunc <- pkernel(truncate, movementmodel, move.a, move.b, Inf, TRUE)
            p <- p * ptrunc
        }
        if (lower.tail) {
            p <- 1-p
        }
        if (movementmodel %in% c('normal','BVN')) {
            q <- sqrt(-log(p)*2*move.a^2)
        }    
        else if (movementmodel %in% c('exponential', 'BVE')) {
            # p <- (q/move.a+1) * exp(-q/move.a)
            # Doesn't work W <- VGAM::lambertW(-p/exp(1))
            # q <- (-W - 1) * move.a
            onep <- function(p) {
                onem <- function(m) {
                    # fn <- function (q, p) (q/m+1) * exp(-q/m) - p
                    fn <- function (q, p) 1 - (q/m+1) * exp(-q/m) - p
                    out <- try(uniroot(f = fn, interval = c(m*1e-3, m*1e3), p = p)$root, silent = TRUE)
                    if (inherits(out, 'try-error')) NA else out
                }
                sapply(move.a, onem)
            }
            q <- sapply(1-p, onep)   # adjusted for tail 2021-07-13
        }
        else if (movementmodel %in% c('t2D', 'BVT')) {
            q <- move.a * sqrt(p^(-1/move.b) - 1)
        }
        else if (movementmodel %in% c('frE')) {
            q <- qexp(p, 1/move.a, lower.tail = FALSE)
        }
        else if (movementmodel %in% c('frG')) {
            q <- qgamma (p, shape = move.b, scale = move.a, lower.tail = FALSE)
        }
        else if (movementmodel %in% c('frL')) {
            mu <- log(move.a)
            s <- sqrt(log(1 + 1/move.b))
            q <- qlnorm(p, meanlog = mu, sdlog = s, lower.tail = FALSE)
        }
        else if (movementmodel %in% c('uniform')) {
            q <- truncate * sqrt(p)
        }
        else {
            q <- rep(NA, length(p))
            warning('probability not available for movement kernel "', 
                movementmodel, '"')
        }
        q[q>truncate] <- NA
        q[bz>p0] <- 1-lower.tail
        q    
    }
}

matchscale <- function(movementmodel, q = 40, p = 0.5, lower = 1e-5, upper = 1e5, 
    move.b = 1, truncate = Inf) {
    mfn <- function(move.a) pkernel (q, movementmodel, move.a, move.b, truncate)-p
    if (movementmodel %in% c(
        'normal', 'BVN', 'BVNzi', 
        'exponential','BVE','BVEzi',
        't2D','BVT',
        'frE','frEzi',
        'frG',
        'frL')) {
        out <- try(uniroot(mfn, c(lower,upper))$root, silent = TRUE)
    }    
    else if (movementmodel %in% c('uniformzi')) {
        upper <- min(upper, 1.0)  # do not exceed feasible limit
        if (is.finite(truncate))
        out <- try(uniroot(mfn, c(lower,upper))$root, silent = TRUE)
        else
            stop("uniformzi kernel must be truncated at finite radius")
    }
    else if (movementmodel %in% c('uniform')) {
        mfnt <- function(truncate) pkernel (q, 'uniform', truncate = truncate)-p
        if (is.finite(truncate))
            out <- try(uniroot(mfnt, c(lower,upper))$root, silent = TRUE)
        else
            stop("uniform kernel must be truncated at finite radius")
    }
    else {
        stop(movementmodel, " not available in matchscale()")
    }
    if (inherits(out, 'try-error')) out <- NA
    out
}

# dkernelMasked <- function (..., mask) {
#     dk <- dkernel(...)  ## unconstrained
#     realize <- function (xy) {
#         sweep(mask, MARGIN = 2, STATS = xy, FUN="+")
#     }
#     apply(mask,1,realize)
# }
