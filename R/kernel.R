################################################################################
## package 'openCR'
## kernel.R
## 2019-06-09, 11
## 2020-10-06 tweak to ed, bR
################################################################################

## fillkernelp is used in prwisecr.R
## otherwise, these functions are used for exploration of kernel characteristics
## kernel is coded independently within likelihood functions

#-----------------------------------------------------------------------
## fill cells of kernel with probability of movement from central point

fillkernelp <- function (
    J,              ## number of sessions
    kerneltype,     ## 0 normal, 1 exponential, 2 usermodel, 3 2Dt, 4 uniform
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
        else stop("unrecognised kerneltype")
    }
    ## normalise
    if (normalize) {
        sumj <- apply(kernelp,2,sum)
        kernelp <- sweep(kernelp, MARGIN=2, STATS=sumj, FUN="/")
    }
    kernelp
}

make.kernel <- function (movementmodel = c('normal','exponential','t2D','uniform'),
                        kernelradius = 10, spacing, 
                        # pars, clip = FALSE,
                        move.a, move.b, clip = FALSE,
                        normalize = TRUE) {
    if (is.function (movementmodel)) {
        moveargs <- formalArgs(movementmodel)
        usermodel <- movementmodel
        movementmodel <- "user"
    }
    else {
        usermodel <- NULL
        movementmodel <- match.arg(movementmodel)
    }
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

    # cellsize = 1 because already inflated
    kernelp <- fillkernelp (1, movementcode-2, kernel, usermodel, cellsize=1, moveargsi, moveargs, FALSE)[,1]

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
    attr(kernel, 'movementmodel') <- movementmodel
    attr(kernel, 'k2') <- k2
    attr(kernel, 'move.a') <- move.a
    attr(kernel, 'move.b') <- move.b
    class(kernel) <- c('kernel','mask','data.frame')
    kernel
}

getpstring <- function (movementcode, move.a, move.b, sep = ' ') {
    npar <- switch(movementcode, static = 0, uncorrelated = 0, normal = 1, exponential = 1, t2D = 2, 0)
    pstring <- ''
    if (npar>0) pstring <- paste0('move.a = ', move.a)
    if (npar>1) pstring <- paste(c(pstring, paste0(' move.b = ', move.b)), collapse = sep)
    pstring
}

plot.kernel <- function (x, contour = FALSE, levels = NULL, text = FALSE, ...) {
    spacing <- spacing(x)
    k2 <- attr(x, 'k2')
    move.a <- attr(x, 'move.a')
    move.b <- attr(x, 'move.b')
    movementmodel <- attr(x, 'movementmodel')
    npar <- switch(movementmodel, static = 0, uncorrelated = 0, normal = 1, exponential = 1, t2D = 2)
    pstring <- getpstring(movementmodel, move.a, move.b)
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
    plot(centrecell, dots = FALSE, meshcol = 'black', col = NA, add = T, ...)

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
    title <-  paste('kernel = ', movementmodel, ' spacing =', spacing, ' kernelradius =', k2, ' ', pstring)
    mtext (side=3, line = 1, title, cex=0.9)

    invisible(data.frame(x, kernelp=kernelp[!is.na(kernelp)]))
}

ed <- function(movementmodel = c('normal', 'exponential', 't2D'), pars) {
    if (movementmodel=='user') return (NA)  # 2020-10-06
    movementmodel <- match.arg(movementmodel)
    if (movementmodel == 'normal')
        pars[1] * (pi/2)^0.5         ## Nathan et al 2012 Table 15.1
    else if (movementmodel == 'exponential')
        pars[1] * 2              ## Nathan et al 2012 Table 15.1
    else if (movementmodel == 't2D') {
        if (length(pars)<2) stop ("t2D model has two parameters")
        a <- pars[1]
        b <- pars[2] + 1
        if (b<3/2)
            Inf
        else
            a * pi^0.5 / 2 * exp(lgamma(b-3/2) - lgamma(b-1)) 
    }
    else NA
}

bR <- function(R, movementmodel = c('normal', 'exponential', 't2D', 'uniform'), pars) {
    if (movementmodel=='user') return (NA)  # 2020-10-06
    movementmodel <- match.arg(movementmodel)
    if (movementmodel == 'normal') {
        alpha <- pars[1]
        exp(-R^2/alpha^2/2)
    }
    else if (movementmodel == 'exponential') {
        exp(-R/pars[1])
    }
    else if (movementmodel == 't2D') {
        if (length(pars)<2) stop ("t2D model has two parameters")
        a <- pars[1]
        b <- pars[2] + 1
        (a^2 / (a^2 + R^2)) ^b
    }
    else if (movementmodel == 'uniform') {
        NA  ## un defined
    }
    else NA
}

summary.kernel <- function (object, ...) {
    spacing <- spacing(object)
    k2 <- attr(object, 'k2')
    move.a <- attr(object, 'move.a')
    move.b <- attr(object, 'move.b')
    pars <- c(move.a, move.b)
    movementmodel <- attr(object, 'movementmodel')
    kernelp <- covariates(object)$kernelp
    kernelp <- kernelp / sum(kernelp)   # force normalization for E(r)
    r <- sqrt(apply(object^2,1,sum))
    result <- list(k2 = k2,
         spacing = spacing,
         ncells = nrow(object),
         movementmodel = movementmodel,
         move.a = move.a,
         move.b = move.b,
         expectedmove = ed(movementmodel, pars),
         expectedmovetr = sum(r * kernelp),
         ptruncated = bR((k2+0.5)*spacing, movementmodel, pars)
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
}

