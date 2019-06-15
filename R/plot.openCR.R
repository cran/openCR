plot.openCR <- function(x, par = 'phi', newdata = NULL, add = FALSE, xoffset = 0, ylim = NULL,
                        useintervals = TRUE, CI = TRUE, intermediate.x = TRUE,  alpha = 0.05, ...) {
    if (useintervals)
        xv <- cumsum(c(0, x$intervals))
    else
        xv <- 0:length(x$intervals)
    sessionxv <- xv
    if (intermediate.x & (par %in% c('phi', 'f', 'lambda', 'b', 'BN','BD')))
        xv <- (xv + c(xv[-1],NA))/2
    xv <- xv + xoffset

    sessnames <- sessionlabels(x$capthist)
    if (is.null(sessnames)) sessnames <- 1:length(xv)

    pred <- predict(x, newdata=newdata, alpha = alpha)[[par]][1:length(xv),]

    xl <- range(sessionxv)
    yl <- ylim
    if (is.null(yl)) {
        yl <- c(min(pred$lcl, na.rm=TRUE)*0.8, max(pred$ucl, na.rm=TRUE)*1.05)
        if (yl[1]<0.05) yl[1] <- 0
        if (yl[2]>0.95 & yl[2]<1) yl[2] <- 1
    }
    args <- list(...)
    argsbase <- args[names(args) %in% c('xlab','ylab','axes','col','cex')]
    defaultbase <- list(x = sessionxv, y = pred$estimate, xlab = 'Session', ylab = par,
                        ylim = yl, xlim = xl, type = 'n', yaxs = 'i')
    base <- replace(defaultbase, names(argsbase), argsbase)
    base$axes <- FALSE
    argspt <- args[names(args) %in% c('pch','cex', 'col', 'fg', 'bg', 'type', 'xpd', 'xaxs', 'yaxs')]
    argspt$x <- xv
    argspt$y <- pred$estimate
    argsseg <- args[names(args) %in% c('lty','lwd', 'col', 'xpd')]
    argsseg$x0 <- xv
    argsseg$y0 <- pred$lcl
    argsseg$x1 <- xv
    argsseg$y1 <- pred$ucl
    if (!add) {
        #plot (sessionxv, pred$estimate, type = 'n', xlab = 'Session', ylab = par,
        #      ylim = yl, xlim = xl, axes = FALSE, yaxs = 'i')
        do.call(plot, base)
        if (is.null(argsbase$axes) || argsbase$axes) {
            axis(1, at = sessionxv, labels = sessnames)
            axis (2, at = pretty(yl), las=1)
            box(bty = 'l')
        }
    }
        
    if (CI) {
        do.call(segments, argsseg)
        #segments(xv, pred$lcl, xv, pred$ucl)
    }
    # points (xv, pred$estimate, ...)
    do.call(points, argspt)
    invisible(xv)
}
