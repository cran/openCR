\name{par.openCR.fit}
\alias{par.openCR.fit}
\alias{openCRlist}

\title{Fit Multiple openCR Models}

\description{
  
  This function is a wrappers for \code{\link{openCR.fit}}.

}

\usage{

par.openCR.fit (arglist, ncores = 1, seed = 123, trace = TRUE, logfile = "logfile.txt", 
    prefix = "")

openCRlist (...)

}

\arguments{
  
  \item{arglist}{list of argument lists for \code{secr.fit} or a character
    vector naming such lists}
  
  \item{ncores}{ integer number of cores to be used for parallel processing}

  \item{seed}{integer pseudorandom number seed}

  \item{trace}{logical; if TRUE intermediate output may be logged}
  
  \item{logfile}{character name of file to log progress reports}
  
  \item{prefix}{character prefix for names of output}
  
  \item{\dots}{ openCR objects}
}

\details{

  Any attempt in \code{arglist} to set \code{ncores > 1} for a
  particular secr fit is ignored.

  \code{trace} overrides any settings in \code{arglist}. Reporting of
  intermediate results is unreliable on Windows when \code{ncores > 1}.

  It is convenient to provide the names of the capthist and mask arguments
  in each component of arglist as character values (i.e. in quotes); objects thus
  named are exported from the workspace to each worker process (see Examples).

  \code{openCRlist} forms a special list (class openCRlist) of fitted model (openCR) objects.

}

\value{

For \code{par.openCR.fit} - openCRlist of model fits (see
\code{\link{openCR.fit}}). Names are created by prefixing \code{prefix} to the
names of \code{argslist}. If \code{trace} is TRUE then the total
execution time and finish time are displayed.

}

\seealso{
  \code{\link{openCR.fit}},
  \link{Parallel},
  \code{\link{make.table}}
}

\examples{
\dontrun{

m1 <- list(capthist = ovenCH, model = list(p~1, phi~1)) 
m2 <- list(capthist = ovenCH, model = list(p~session, phi~1))
m3 <- list(capthist = ovenCH, model = list(p~session, phi~session) )
fits <- par.openCR.fit (c('m1','m2','m3'), ncores = 3)
AIC(fits)

}
}
\keyword{ model }

