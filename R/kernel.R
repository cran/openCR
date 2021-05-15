################################################################################
## package 'openCR'
## kernel.R
## 2019-06-09, 11
## 2020-10-06 tweak to ed, bR
## 2021-02-23 annular etc.
## 2021-02-25 functions not exported: empirical
## 2021-03-01 functions exported pkernel, matchscale
## 2021-03-16 annular2
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

fillkernelp <- function (
    J,              ## number of sessions
    kerneltype,     ## 0 normal, 1 exponential, 2 usermodel, 3 2Dt, 4 uniform, 
                    ## 5 annular, 6 annular2, 7 annularR
    sparsekernel,   ## TRUE iff sparse
    # zeroinflated, 
    kernel,
    usermodel,
    cellsize,
    moveargsi,
    moveargs,
    normalize = TRUE)
{
    kn <- nrow(kernel)
    kernelp <- matrix(0, nrow = kn, ncol = J)
    r2 <- (kernel$x^2 + kernel$y^2) * cellsize^2
    r <- sqrt(r2)
    ## diag is sqrt(2) adjustment for diagonal spokes of sparse kernels
    diag <- ifelse(abs(kernel$x) == abs(kernel$y) & (r>0), sqrt(2), 1)
    zero <- kernel$x == 0 & kernel$y == 0
    for (j in 1:J) {
        if (kerneltype == 0) {        ## normal (Gaussian) kernel
            sigma <- moveargs[j]
            kernelp[,j] <- exp(-r2 / 2 / sigma^2) / 2 / pi / sigma^2
        }
        else if (kerneltype == 1) {  ## negative exponential kernel
            sigma <- moveargs[j]
            kernelp[,j] <- exp(-r / sigma) / 2 / pi / sigma^2
        }
        else if (kerneltype == 2) {
            if (!is.function(usermodel))
                stop ("fillkernelp requires valid user model for movement")
            if (moveargsi[2]>0)
                kernelp[,j] <- usermodel(r, moveargs[j,1], moveargs[j,2])  ##  2 parameters
            else if (moveargsi[1]>0)
                kernelp[,j] <- usermodel(r, moveargs[j,1])
            else
                kernelp[,j] <- usermodel(r)   ## no parameters!
        }
        else if (kerneltype == 3) {   ## 2-D t
            a2 <- moveargs[j,1]^2
            b <- moveargs[j,2] + 1
            kernelp[,j] <-  (b-1) / pi / a2 / (1 + r^2/a2)^b
        }
        else if (kerneltype == 4) {   ## uniform
            kernelp[,j] <-  1 / nrow(kernelp)
        }
        else if (kerneltype == 5) {   ## annular
            p0 <- moveargs[j,1]
            kernelp[,j] <-  (1-p0)/(kn-1)  # assume intermediate r already removed
            kernelp[r==0,j] <-  p0
        }
        else if (kerneltype == 6) {   ## double annular
            p0 <- moveargs[j,1]
            p1 <- moveargs[j,2]
            ri <- round(r/spacing(kernel))
            origin <- ri==0
            K2 <- max(ri)
            ring1 <- ri>0 & (ri<(K2-0.5))
            ring2 <- ri > (K2-0.5)
            kernelp[origin,j] <-  p0
            kernelp[ring1,j] <-  p1/sum(ring1)
            kernelp[ring2,j] <-  (1-(p0+p1))/sum(ring2)
        }
        else if (kerneltype == 7) {   ## annularR
            # annular with variable radius (move.b) as fraction of kernel radius
            # and cells weighted according to relative length of intersecting arc (dtheta)
            move.a <- moveargs[j,1]
            move.b <- moveargs[j,2]
            rge <- range(kernel$x) + spacing(kernel) * c(-0.5,0.5)
            x <- y <- seq(rge[1], rge[2], spacing(kernel))
            minxy <- min(x)*1.1
            maxxy <- max(x)*1.1
            rad <- move.b * max(kernel$x)
            vx <- function(x) circleintersectsline(x,x, minxy, maxxy, rad)
            hy <- function(y) circleintersectsline(minxy, maxxy, y, y, rad)
            xi <- do.call(rbind, lapply(x, vx))
            yi <- do.call(rbind, lapply(y, hy))
            pts <- rbind(xi,yi)
            df <- data.frame(pts, theta = atan2(pts$y, pts$x))
            df <- df[order(df$theta),]
            df <- rbind(df, df[1,])
            dtheta <- diff(df$theta)
            dtheta[dtheta<0] <- dtheta[dtheta<0] + 2*pi
            df$dtheta <- c(dtheta,NA)
            df$cx <- df$x + c(diff(df$x)/2, NA)
            df$cy <- df$y + c(diff(df$y)/2, NA)
            # associate each arc with cell in which it falls
            incell <- nearesttrap(df[,c('cx','cy')], kernel)
            incell <- incell[incell>0]
            kernelp[,j] <- rep(0, nrow(kernel))
            kernelp[trunc(nrow(kernel)/2)+1,j] <- move.a   ## centre
            kernelp[incell,j] <- (1-move.a) * df$dtheta[-nrow(df)] / sum(df$dtheta[-nrow(df)])
        }
        else stop("unrecognised kerneltype")
        if (sparsekernel) kernelp[,j] <- 2 * pi * r * diag * kernelp[,j]
        # if (zeroinflated) {
        #     if (!(kerneltype %in% c (0,1))) stop ("Zero inflation only an option for normal, exponential")
        #     pzero <- moveargs[j,2]
        #     kernelp[,j] <- kernelp[,j] * (1-pzero) * spacing(kernel)^2 + zero * pzero
        # }
    }
    ## normalise
    if (normalize) {
        sumj <- apply(kernelp,2,sum, na.rm = TRUE)
        kernelp <- sweep(kernelp, MARGIN=2, STATS=sumj, FUN="/")
    }
    kernelp
}

make.kernel <- function (
    movementmodel = c('normal','exponential','t2D','uniform'),
    kernelradius = 10, spacing, move.a, move.b, 
    sparsekernel = FALSE, # zeroinflated = FALSE, 
    clip = FALSE,
    normalize = TRUE,
    ...) 
{
    if (inherits(movementmodel, 'openCR')) {
        fit <- movementmodel
        if (fit$movementmodel %in% c('normal','exponential','t2D','uniform','annular')) {
            pred <- predict(fit, ...)
            kernel <- make.kernel(
                movementmodel = fit$movementmodel, 
                kernelradius  = fit$kernelradius, 
                spacing       = secr::spacing(fit$mask), 
                move.a        = round(pred$move.a[1,'estimate'], 3),
                move.b        = pred$move.b[1,'estimate'], 
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
            usermodel <- movementmodel
            movementmodel <- "user"
        }
        else {
            usermodel <- NULL
            movementmodel <- match.arg(movementmodel[1],
                choices = c('normal','exponential','t2D','uniform','annular',
                    'annular2','annularR')) 
        }
        if (missing(move.a) && movementmodel %in% c('normal','exponential','t2D','annular',
            'annular2', 'annularR')) { 
            stop ("move.a required for movementmodel ", movementmodel)
        }
        if (missing(move.b) && movementmodel %in% c('t2D','annular2','annularR')) {
            stop ("move.b required for movementmodel ", movementmodel)
        }
        # if (missing(move.b) && zeroinflated) {
        #     stop ("move.b required for zero inflation")
        # }
        
        movementcode <- movecode(movementmodel)
        moveargsi <- c(0,0)
        if (missing(move.a) | (movementmodel == 'uniform')) {
            pars <- move.a <- move.b <- NULL
        }
        else {
            pars <- move.a
            moveargsi[1] <- 1
        }
        if (missing(move.b)) {
            move.b <- NULL
        }
        else {
            pars <- c(pars, move.b)
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
        # cellsize = 1 because already inflated
        kernelp <- fillkernelp (1, movementcode-2, sparsekernel, # zeroinflated, 
            kernel, usermodel, cellsize = 1, moveargsi, moveargs, FALSE)[,1]
        ## optional clipping 
        if (clip) {
            outside <- (kernel$x^2 + kernel$y^2) > ((k2+0.5)*spacing)^2
            kernel <- subset(kernel, !outside)
            kernelp[outside] <- NA
        }
        
        if (normalize) {
            kernelp <- kernelp/ sum(kernelp[!is.na(kernelp)])
        }
        covariates(kernel) <- data.frame(kernelp = kernelp[!is.na(kernelp)])
        
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
    npar <- switch(movementmodel, static = 0, uncorrelated = 0, 
        normal = 1, exponential = 1, t2D = 2, 
        uniform = 0, annular = 1, annular2 = 2, annularR = 2, 
        0)
    pstring <- ''
    if (npar>0) pstring <- paste0('move.a = ', move.a)
    if (npar>1) pstring <- paste(c(pstring, paste0(' move.b = ', move.b)), collapse = sep)
    pstring
}

plot.kernel <- function (x, contour = FALSE, levels = NULL, text = FALSE, title = NULL, ...) {
    spacing <- spacing(x)
    k2 <- attr(x, 'k2')
    move.a <- attr(x, 'move.a')
    move.b <- attr(x, 'move.b')
    movementmodel <- attr(x, 'movementmodel')
    npar <- switch(movementmodel, static = 0, uncorrelated = 0, normal = 1, 
        exponential = 1, t2D = 2, uniform = 0, annular = 1, annular2 = 2, 
        annularR = 2)
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
    plot(msk, dots = FALSE, meshcol = 'white', covariate = 'kernelp', border = border, ...)
    
    centrecell <- subset(msk, (x$x^2 + x$y^2)<1e-6)
    # 2021-03-17 modify plot call to allow add in ...
    dots$x <- centrecell
    dots$dots <- FALSE
    dots$meshcol <- 'black'
    dots$col <- NA
    dots$add = TRUE
    do.call(plot, dots)
    #plot(centrecell, dots = FALSE, meshcol = 'black', col = NA, add = T, ...)

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
    mtext (side=3, line = 1, title, cex=0.9)

    invisible(data.frame(x, kernelp=kernelp[!is.na(kernelp)]))
}

ed <- function(movementmodel = c('normal', 'exponential', 't2D', 'uniform', 
    'annular', 'annular2', 'annularR'), move.a, move.b, radius) {
    if (movementmodel=='user') return (NA)  # 2020-10-06
    movementmodel <- match.arg(movementmodel[1], choices = c('normal',
        'exponential', 't2D', 'uniform', 'annular', 'annular2', 'annularR'))
    if (movementmodel %in% c('normal'))
        move.a * (pi/2)^0.5         ## Nathan et al 2012 Table 15.1
    else if (movementmodel %in% c('exponential'))
        move.a * 2              ## Nathan et al 2012 Table 15.1
    else if (movementmodel %in% c('t2D')) {
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
    else NA
}

bR <- function(R, movementmodel = c('normal', 'exponential', 't2D', 'uniform',
    'annular', 'annular2', 'annularR'), move.a, move.b) {
    if (movementmodel=='user') return (NA)  # 2020-10-06
    movementmodel <- match.arg(movementmodel[1], choices = c('normal',
        'exponential','t2D','uniform','annular', 'annular2', 'annularR')) 
    if (movementmodel %in% c('normal')) {
        alpha <- move.a
        exp(-R^2/alpha^2/2)
    }
    else if (movementmodel %in% c('exponential')) {
        (R/move.a + 1) * exp(-R/move.a)   # fixed 2021-02-22
    }
    else if (movementmodel %in% c('t2D')) {
        if (missing(move.b)) stop ("t2D model has two parameters")
        a <- move.a
        b <- move.b + 1
        (a^2 / (a^2 + R^2)) ^b
    }
    else if (movementmodel %in% c('uniform', 'annular', 'annular2', 'annularR')) {
        NA
    }
    else NA
}

summary.kernel <- function (object, ...) {
    spacing <- spacing(object)
    k2 <- attr(object, 'k2')
    movementmodel <- attr(object, 'movementmodel')
    move.a <- attr(object, 'move.a')
    move.b <- attr(object, 'move.b')
    kernelp <- covariates(object)$kernelp
    kernelp <- kernelp / sum(kernelp)   # force normalization for E(r)
    r <- sqrt(apply(object^2,1,sum))
    result <- list(k2 = k2,
        spacing = spacing,
        ncells = nrow(object),
        movementmodel = movementmodel,
        move.a = move.a,
        move.b = move.b,
        expectedmove = ed(movementmodel, move.a, move.b, k2 * spacing),
        expectedmovetr = sum(r * kernelp),
        ptruncated = bR((k2+0.5)*spacing, movementmodel, move.a, move.b),
        expectedq50 = qkernel(0.5, movementmodel, move.a, move.b),
        expectedq90 = qkernel(0.9, movementmodel, move.a, move.b)
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
    cat('Expected movement (trunc) : ', x$expectedmovetr, ' (m)\n')
    cat('Expected movement (full)  : ', x$expectedmove, ' (m)\n')
    cat('Median movement (full)    : ', x$expectedq50, ' (m)\n')
    cat('90th percentile (full)    : ', x$expectedq90, ' (m)\n')
}

# useful functions not exported 2021-02-25
# see kerneltest.R

quantile.kernel <- function (x, probs = seq(0, 1, 0.25), na.rm = FALSE,
    names = TRUE, type = 7, truncate = Inf, ...) 
{
    movementmodel <- attr(x, 'movementmodel')
    move.a <- attr(x, 'move.a')
    move.b <- attr(x, 'move.b')
    if (movementmodel == 'normal') {
        pfn <- function(R, q) (1-exp(-R^2/2/move.a^2))-q
    }    
    else if (movementmodel == 'exponential') {
        pfn <- function(R, q) (1-(R/move.a+1) * exp(-R/move.a))-q
    }
    else if (movementmodel == 't2D') {
        pfn <- function(R, q) (1-(move.a^2/(move.a^2+R^2))^move.b)-q
    }
    upper <- 1e5*move.a
    qfn <- function(q) {
        out <- try(uniroot(pfn, c(0,upper), q=q)$root, silent = TRUE)
        if (inherits(out, 'try-error')) NA else out
    }
    
    # optionally adjust for truncation at edge
    if (truncate != Inf) probs <- probs * pfn(truncate,0)
    
    q <- sapply(probs, qfn)
    q[q==upper] <- Inf
    q
}

empirical <- function (kernel) {
    r <- apply(kernel^2,1,sum)^0.5
    p <- covariates(kernel)$kernelp
    p <- p[order(r)]
    r <- sort(r)
    out <- cbind(r, cumsum(p))
}

# Distribution function for distance moved
dkernel <- function (r, movementmodel = c('normal','exponential','t2D'),
    move.a, move.b, truncate = Inf) {
    movementmodel <- match.arg(movementmodel)
    movementmodel <- try(match.arg(movementmodel), silent = TRUE)
    if (inherits(movementmodel, 'try-error')) {
        warning("dkernel defined only for movement models normal, exponential, t2D")
        rep(NA, length(r))
    }
    else {
        if (movementmodel == 'normal') {
            d <- r/move.a^2 * exp(-r^2/2/move.a^2)
        }    
        else if (movementmodel == 'exponential') {
            d <- r/move.a^2 * exp(-r/move.a)
        }
        else if (movementmodel == 't2D') {
            d <- 2 * move.b * r / move.a^2 * (1 + r^2/move.a^2)^-(move.b+1)
        }
        else {
            d <- rep(NA, length(r))
            warning('pdf not available for movement kernel "', 
                movementmodel, '"')
        }
        if (truncate != Inf) {
            ptrunc <- pkernel(truncate, movementmodel, move.a, move.b, Inf, TRUE)
            d <- d / ptrunc
            d[r>truncate] <- 0
        }
        d
    }
}

# Distribution function for distance moved
pkernel <- function (q, movementmodel = c('normal','exponential','t2D'),
    move.a, move.b, truncate = Inf, lower.tail = TRUE) {
    movementmodel <- try(match.arg(movementmodel), silent = TRUE)
    if (inherits(movementmodel, 'try-error')) {
        warning("pkernel defined only for movement models normal, exponential, t2D")
        rep(NA, length(q))
    }
    else {
        if (movementmodel == 'normal') {
            p <- exp(-q^2/2/move.a^2)
        }    
        else if (movementmodel == 'exponential') {
            p <- (q/move.a+1) * exp(-q/move.a)
        }
        else if (movementmodel == 't2D') {
            p <- (move.a^2/(move.a^2+q^2))^move.b
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
            p <- p / pkernel(truncate, movementmodel, move.a, move.b, Inf, TRUE)
        }
        p
    }
}

qkernel <- function(p, movementmodel = c('normal','exponential','t2D'),
    move.a, move.b, truncate = Inf, lower.tail = TRUE) {
    movementmodel <- try(match.arg(movementmodel), silent = TRUE)
    if (inherits(movementmodel, 'try-error')) {
        warning("qkernel defined only for movement models normal, exponential, t2D")
        rep(NA, length(p))
    }
    else {
        if (truncate != Inf) {
            if (!lower.tail) stop ("cannot combine truncation and upper tail")
            ptrunc <- pkernel(truncate, movementmodel, move.a, move.b, Inf, TRUE)
            p <- p * ptrunc
        }
        if (lower.tail) {
            p <- 1-p
        }
        if (movementmodel == 'normal') {
            q <- sqrt(-log(p)*2*move.a^2)
        }    
        else if (movementmodel == 'exponential') {
            # p <- (q/move.a+1) * exp(-q/move.a)
            # Doesn't work W <- VGAM::lambertW(-p/exp(1))
            # q <- (-W - 1) * move.a
            onep <- function(p) {
                onem <- function(m) {
                    fn <- function (q, p) (q/m+1) * exp(-q/m) - p
                    out <- try(uniroot(f = fn, interval = c(m*1e-3, m*1e3), p = p)$root, silent = TRUE)
                    if (inherits(out, 'try-error')) NA else out
                }
                sapply(move.a, onem)
            }
            q <- sapply(p, onep)
        }
        else if (movementmodel == 't2D') {
            q <- move.a * sqrt(p^(-1/move.b) - 1)
        }
        else {
            q <- rep(NA, length(p))
            warning('probability not available for movement kernel "', 
                movementmodel, '"')
        }
        q[q>truncate] <- NA
        q    
    }
}

# t2Dratio <- function (b) sqrt((0.9^(-1/b) - 1) / (0.5^(-1/b) - 1))

## test
# matchscale('normal', q=100, p = pkernel(100, 'normal', 33.97))
# [1] 33.97
# pkernel(100, 'uniform', 33.97, 1)

matchscale <- function(movementmodel, q = 40, p = 0.5, lower = 1e-5, upper = 1e5, move.b = 1) {
    if (movementmodel == 'normal') {
        mfn <- function(move.a) (1-exp(-q^2/2/move.a^2))-p
    }    
    else if (movementmodel == 'exponential') {
        mfn <- function(move.a) (1-(q/move.a+1) * exp(-q/move.a))-p
    }
    else if (movementmodel == 't2D') {
        mfn <- function(move.a) (1-(move.a^2/(move.a^2+q^2))^move.b)-p
    }
    out <- try(uniroot(mfn, c(lower,upper))$root, silent = TRUE)
    if (inherits(out, 'try-error')) out <- NA
    out
}
