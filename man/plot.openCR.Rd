\name{plot.openCR}
\alias{plot.openCR}

\title{Plot Estimates}

\description{

Session-specific estimates of the chosen parameter are plotted.

}

\usage{

\method{plot}{openCR} (x, par = "phi", newdata = NULL, add = FALSE, xoffset = 0, ylim = NULL, 
    useintervals = TRUE, CI = TRUE, intermediate.x = TRUE, alpha = 0.05, ...)
    
}

\arguments{
  \item{x}{openCR object from openCR.fit}
  \item{par}{character names of parameter to plot}
  \item{newdata}{dataframe of predictor values for \code{\link{predict}} (optional)}
  \item{add}{logical; if TRUE then points are added to an existing plot}
  \item{xoffset}{numeric offset to be added to all x values}
  \item{ylim}{numeric vector of limits on y-axis}
  \item{useintervals}{logical; if TRUE then x values are spaced according to 
        the intervals attribute}
  \item{CI}{logical; if TRUE then 1-alpha confidence intervals are plotted}
  \item{intermediate.x}{logical; if TRUE then turnover parameters are plotted at 
        the mid point on the x axis of the interval to which they relate}
  \item{alpha}{numeric confidence level default (alpha = 0.05) is 95\% interval}
  \item{\dots}{other arguments passed to \code{\link{points}}}
}

\details{

If ylim is not provided it is set automatically.

Note that the \dots argument is passed only to \code{\link{points}}. If you wish to customize the base plot then do that in advance and use \code{add = TRUE}.

}

\value{

None

}

\seealso{

\code{\link{predict}}

}

\examples{

\dontrun{

fit <- openCR.fit(join(ovenCH), type='CJS', model = phi~session)
plot(fit,'phi', pch = 16, cex=1.3, yl=c(0,1))

}

}

\keyword{ hplot }