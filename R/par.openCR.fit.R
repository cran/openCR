## 2017-05-12
## openCR parallel fit, derived, region.N

par.openCR.fit <- function (arglist, ncores = 1, seed = 123, trace = TRUE,
                          logfile = "logfile.txt", prefix = "") {
    ptm  <- proc.time()
    ## 'inherits' causes R to search in enclosing frames
    if (is.character(arglist))
        arglist <- mget(arglist, inherits = TRUE)
    
    ## force 'trace' to common value across all components of arglist
    arglist <- lapply(arglist, function (x) {x$trace <- trace; x})

    ## check for capthist, mask, dframe mentioned by name
    ## objects are exported to the worker processes as required
    getnames <- function(obj = 'capthist') {
        tmpnames <- sapply(arglist, function(x) if (is.character(x[[obj]])) x[[obj]] else '')
        unique(tmpnames)
    }
    data <- c(getnames('capthist'), getnames('mask'),getnames('dframe'),getnames('details'))
    data <- data[nchar(data)>0]

    ## default details savecall to FALSE across all components of arglist
    arglist <- lapply(arglist, function (x) {
        if (is.null(x$details))
            x$details <- list(savecall = FALSE)
        else if (!('savecall' %in% names(x$details))) {
            x$details[['savecall']] <- FALSE
        }
        x
    })
    
    ## individual fits must use ncores = 1
    if (ncores > 1) {
        ## force 'ncores' to 1 across all components of arglist
        arglist <- lapply(arglist, function (x) {x$ncores <- 1; x})
        clust <- makeCluster(ncores, methods = FALSE, useXDR = .Platform$endian=='big', outfile = logfile)
        clusterSetRNGStream(clust, seed)
        clusterExport(clust, c(data, 'openCR.fit'), environment())
        output <- parLapply(clust, arglist, do.call, what = 'openCR.fit')
        stopCluster(clust)
    }
    else {
        set.seed (seed)
        output <- lapply(arglist, do.call, what = 'openCR.fit')
    }
    
    message('Completed in ', round((proc.time() - ptm)[3]/60,3), ' minutes at ',
        format(Sys.time(), "%H:%M:%S %d %b %Y"))

    if (inherits(output[[1]], 'openCR'))
        output <- openCRlist(output)

    ## apply standard naming convention
    names(output) <- paste0(prefix, names(arglist))

    output
}

# par.derived <- function (secrlist, ncores = 1, ...) {
# 
#     if (!inherits(secrlist, 'secrlist'))
#         stop("requires secrlist input")
#     if (ncores > 1) {
#         clust <- makeCluster(ncores, methods = FALSE, useXDR = .Platform$endian=='big')
#         output <- parLapply(clust, secrlist, derived, ...)
#         stopCluster(clust)
#     }
#     else {
#         output <- lapply(secrlist, derived, ...)
#     }
#     names(output) <- names(secrlist)
#     output
# }
# 
# par.region.N <- function (secrlist, ncores = 1, ...) {
# 
#     if (!inherits(secrlist, 'secrlist'))
#         stop("requires secrlist input")
#     if (ncores > 1) {
#         clust <- makeCluster(ncores, methods = FALSE, useXDR = .Platform$endian=='big')
#         output <- parLapply(clust, secrlist, region.N, ...)
#         stopCluster(clust)
#     }
#     else {
#         output <- lapply(secrlist, region.N, ...)
#     }
#     names(output) <- names(secrlist)
#     output
# }


