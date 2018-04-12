###############################################################################
## openCR
## make.table.R
## 2018-02-26 openCR 1.0.0
## 2018-03-27 openCR 1.1.0 select first nsess[i] rows from pred[i]
###############################################################################

make.table <- function (fits, parm = 'phi', fields = 'estimate', ...) {
    if (!inherits(fits, 'openCRlist')) fits <- openCRlist(fits)
    if (is.null(names(fits))) names(fits) <- paste0('fit',1:length(fits))
    nsess <- sapply(fits, function(x) length(x$sessionlabels)) 
    pred <- predict(fits, ...)
    rown <- names(fits)
    coln <- unique(as.character(unlist(lapply(fits, '[[', 'sessionlabels'))))
    tab <- matrix (nrow = length(rown), ncol = length(coln), dimnames=list(model=rown, session = coln))
    tab <- as.table(tab)
    for (i in names(fits)) {
        if (parm %in% names(pred[[i]])) {
            rowcol <- cbind(rep(i, nsess[i]), fits[[i]]$sessionlabels)
            tab[rowcol] <- pred[[i]][[parm]][1:nsess[i],fields]
        }
    }
    tab
}
