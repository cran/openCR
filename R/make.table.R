###############################################################################
## openCR
## make.table.R
## 2018-02-26 openCR 1.0.0
## 2018-03-27 openCR 1.1.0 select first nsess[i] rows from pred[i]
## 2021-04-21 stratification
###############################################################################

make.table <- function (fits, parm = 'phi', fields = 'estimate', strata = 1, ...) {
    if (!inherits(fits, 'openCRlist')) fits <- openCRlist(fits)
    if (is.null(names(fits))) names(fits) <- paste0('fit',1:length(fits))
    if (length(strata)>1) {
        lapply (strata, make.table, fits=fits, parm = parm, fields = fields)
    }
    else {
        # strata of length 1
        nsess <- sapply(fits, function(x) length(primaryintervals(x)[[strata]])+1) 
        pred <- predict(fits, ...)
        rown <- names(fits)
        labellist <- lapply(fits, '[[', 'sessionlabels')
        # for backward compatibility 2021-04-26
        if (!is.list(labellist[[1]])) labellist <- list(labellist)
        coln <- unique(as.character(unlist(labellist)))
        tab <- matrix (nrow = length(rown), ncol = length(coln), 
            dimnames=list(model=rown, session = coln))
        tab <- as.table(tab)
        for (i in rown) { 
            if (parm %in% names(pred[[i]])) {
                rowcol <- cbind(rep(i, nsess[i]), unlist(labellist[[strata]]))
                tab[rowcol] <- pred[[i]][[parm]][1:nsess[i],fields]
            }
        }
        tab
    }
}
