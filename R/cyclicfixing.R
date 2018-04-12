## cyclic fixing over real parameters
## cf Schwarz and Arnason 1996 p865
## 2018-02-12 Under development...
## NMOT ESSENTIAL
# cyclic.fit <- function (..., start) {
#     fit0 <- start
#     nreal <- length(parindx)
#     for (pari in 1:nreal) {
#         betai <- unlist(parindx)  
#         betai[parindx[[pari]]] <- NA
#         if (pari>1) details <- list(fixedbeta = fit0$fit$par[betai])
#         else details <- NULL
#         fit0 <- openCR.fit (..., start = fit0, details = details)
#         cat(names(parindx)[pari], " LL = ", logLik(fit0), '\n')
#         
#         ## joint lambda0, sigma?
#     }   
#     
# }