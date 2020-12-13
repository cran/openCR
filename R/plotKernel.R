plotKernel <- function (movementmodel = c('normal','exponential','t2D','uniform'),
                        kernelradius = 10, spacing, pars, clip = FALSE,
                        plt = TRUE, contour = FALSE, levels = NULL,
                        text = FALSE, normalize = TRUE, ...) {
    
    warning('plotkernel() is deprecated in openCR>=1.4.0. Use make.kernel() and plot instead')
    
    if (is.function (movementmodel)) {
        moveargs <- formalArgs(movementmodel)
        usermodel <- movementmodel
        movementmodel <- "user"
    }
    else {
        usermodel <- NULL
        movementmodel <- match.arg(movementmodel)
    }
    movementcode <- movecode (movementmodel)
    moveargsi <- c(0,0)
    if (missing(pars)) pars <- NULL
    title <-  paste('spacing =', spacing,
                    ' kernelradius =', kernelradius)
    if (length(pars)>0) title <- paste(title, ' pars =', pars[1])
    if (length(pars)>1) title <- paste(title, pars[2], sep=', ')
    if (length(pars)>0) moveargsi[1] <- 1
    if (length(pars)>1) moveargsi[2] <- 2
    if (is.null(pars)) pars <- c(0,0)
    moveargs <- matrix(pars, nrow = 1)   # J-row matrix for fillkernelp
    k2 <- kernelradius
    kernel <- make.mask(type = 'rectangular',
                        spacing = spacing,
                        buffer = 0, nx = 2 * k2+1, ny = 2 * k2+1)
    ## centre
    kernel[,] <- sweep(kernel, MARGIN=2, FUN = "-", STATS = rep((k2+0.5)*spacing,2))

    # cellsize = 1 because already inflated
    kernelp <- fillkernelp (1, movementcode-2, kernel, usermodel, cellsize=1,
                            moveargsi, moveargs, normalize)[,1]

    ## optional clipping (incompatible with contour)
    if (clip) {
        outside <- (kernel$x^2 + kernel$y^2) > ((k2+0.5)*spacing)^2
        kernel <- subset(kernel, !outside)
        kernelp[outside] <- NA
        if (normalize)
            kernelp <- kernelp/ sum(kernelp[!is.na(kernelp)])
    }

    covariates(kernel) <- data.frame(kernelp = kernelp[!is.na(kernelp)])

    if (plt) {
        centrecell <- subset(kernel, (kernel$x^2 + kernel$y^2)<1e-6)
        plot(kernel, dots = FALSE, meshcol = 'white', covariate = 'kernelp', border = 0, ...)
        plot(centrecell, dots = FALSE, meshcol = 'black', col = NA, add = T)
        if (contour) {
            if (is.null(levels))
                levels <- pretty(c(0, max(kernelp, na.rm = TRUE)), 10)
            contour(add = TRUE, (-k2:k2)*spacing, (-k2:k2) * spacing,
                    matrix(kernelp, nrow = 2*k2+1), levels = levels)
        }
        if (text)
            text(kernel$x, kernel$y, round(kernelp,3), cex=0.6)
        mtext (side=3, line = 1, title)

    }
    invisible(data.frame(kernel, kernelp=kernelp[!is.na(kernelp)]))
}

