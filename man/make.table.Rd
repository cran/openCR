\name{make.table}
\alias{make.table}

\title{Tabulate Estimates From Multiple Models}

\description{

Session-specific estimates of real parameters (p, phi, etc.) are arranged in a rectangular table.

}

\usage{

make.table(fits, parm = "phi", fields = "estimate", ...)

}

\arguments{
  \item{fits}{openCRlist object}
  \item{parm}{character name of real parameter estimate to tabulate}
  \item{fields}{character column from predict (estimate, SE.estimate, lcl, ucl)}
  \item{...}{arguments passed to \code{\link{predict.openCRlist}}}
}

\details{

The input will usually be from \code{par.openCR.fit}.

}

\value{

A table object.

}

\seealso{

\code{\link{par.openCR.fit}}, \code{\link{openCRlist}}
}

\examples{

arglist <- list(constant = list(capthist=ovenCHp, model=phi~1), 
                session.specific = list(capthist=ovenCHp, model=phi~session))
fits <- par.openCR.fit(arglist, trace = FALSE)
print(make.table(fits), na=".")

}

\keyword{ manip }