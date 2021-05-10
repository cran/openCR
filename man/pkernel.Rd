\name{pkernel}
\alias{pkernel}
\alias{dkernel}
\alias{qkernel}

\title{Kernel Distribution Functions}

\description{

Distribution of distance moved for each of the main movement kernels. Theoretical probability density, cumulative distribution function, and quantile function (inverse of the cumulative distribution function).

}

\usage{

pkernel(q, movementmodel = c("normal", "exponential", "t2D"), move.a, move.b, 
    truncate = Inf, lower.tail = TRUE)

dkernel(r, movementmodel = c("normal", "exponential", "t2D"), move.a, move.b, 
    truncate = Inf)

qkernel(p, movementmodel = c("normal", "exponential", "t2D"), move.a, move.b, 
    truncate = Inf, lower.tail = TRUE)

}

\arguments{
  \item{p}{numeric vector of cumulative probabilities (0.5 for median)}
  \item{r}{numeric vector of distance moved}
  \item{q}{numeric vector of quantiles (distance moved)}
  \item{movementmodel}{character (see \href{../doc/openCR-vignette.pdf}{openCR-vignette.pdf} and Note below)}
  \item{move.a}{numeric parameter of movement kernel}
  \item{move.b}{numeric parameter of movement kernel}
  \item{truncate}{numeric q value at which distribution truncated}
  \item{lower.tail}{logical; if TRUE (default), probabilities are P[X <= x] otherwise, P[X > x].}
}

\details{

Formulae are given in openCR-vignette.pdf.

Computation of \code{qkernel} for \code{movementmodel = 'exponential'} uses 
numerical root finding (function \code{\link{uniroot}}).

Truncation (\code{truncate = limit} for finite \code{limit}) adjusts probabilities 
upwards by 1/pkernel(limit,..., truncate = Inf) so that 
pkernel(limit, ..., truncate = limit) equals 1.0.
By default the distribution is not truncated.

}

\value{

For \code{pkernel} --

Vector of cumulative probabilities corresponding to q. The cumulative probability is 1.0 for q > truncate.

For \code{dkernel} --

Vector of probability density at radial distance r (zero for r > truncate).

For \code{qkernel} --

Vector of quantiles (distances moved) corresponding to cumulative probabilities p.

}

\seealso{

\code{\link{make.kernel}}, \code{\link{matchscale}}

}

\examples{
# plot 3 distributions chosen with matchscale to intersect at p = 0.5
q <- 0:100
plot(q, pkernel(q, 'normal', 34), type = 'l', ylab = 'Cumulative probability')
lines(q, pkernel(q, 't2D', move.a = 104, move.b = 5), col = 'darkgreen', lwd = 2)
lines(q, pkernel(q, 't2D', move.a = 40, move.b = 1), col = 'orange', lwd = 2)
points(40, 0.5, pch = 16)
legend(62, 0.36, lty=1, lwd = 2, col = c('black','darkgreen','orange'), 
   legend = c('normal sigma=34', 't2D a=104, b=5', 't2D a=40, b=1'))

# median
abline(v = qkernel(0.5, 'normal', 34))
}
