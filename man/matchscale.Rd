\name{matchscale}
\alias{matchscale}

\title{Match Kernel}

\description{

Finds scale parameter (move.a) of a movement model that corresponds to desired quantile.

}

\usage{

matchscale(movementmodel, q = 40, p = 0.5, lower = 1e-05, upper = 1e+05, 
   move.b = 1, truncate = Inf)

}

\arguments{
  \item{movementmodel}{character (see \href{../doc/openCR-vignette.pdf}{openCR-vignettes.pdf})}
  \item{q}{desired quantile (distance moved)}
  \item{p}{cumulative probability}
  \item{move.b}{shape parameter of movement kernel}
  \item{lower}{lower bound interval to search }
  \item{upper}{upper bound interval to search }
  \item{truncate}{numeric q value at which distribution truncated}
}

\details{

The \code{truncate} argument must be specified for movementmodel `uniformzi`. For movementmodel `uniform' there is no parameter and the radius of truncation is varied to achieve the requested quantile q corresponding to cumulative probability p.

}

\value{

Numeric value for move.a (scale parameter or zero-inflation in the case of `uniformzi') or truncation radius (`uniform').

}

\seealso{

\code{\link{pkernel}}, \code{\link{make.kernel}}

}

\examples{
matchscale('normal', 40, 0.5)
matchscale('t2D', 40, 0.5, move.b=1)
matchscale('t2D', 40, 0.5, move.b=5)

}
