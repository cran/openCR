\name{make.kernel}
\alias{make.kernel}
\alias{plot.kernel}
\alias{summary.kernel}

\title{Discrete Movement Kernel}

\description{

Functions to create, plot and summarise a discrete representation of a movement kernel.

}

\usage{

make.kernel(movementmodel = c("normal", "exponential", "t2D", "uniform"), 
    kernelradius = 10, spacing, move.a, move.b, clip = FALSE, normalize = TRUE)

\method{plot}{kernel}(x, contour = FALSE, levels = NULL, text = FALSE, ...)

\method{summary}{kernel}(object, ...)

}
  
\arguments{

  \item{movementmodel}{character or function}
  \item{kernelradius}{integer radius of kernel in grid cells}
  \item{spacing}{numeric spacing between cell centres}
  \item{move.a}{numeric parameter of kernel}
  \item{move.b}{numeric parameter of kernel}
  \item{clip}{logical; if TRUE then corner cells are removed}
  \item{normalize}{logical; if TRUE then cell values are divided by their sum}
  
  \item{x}{kernel object from \code{make.kernel}}
  \item{contour}{logical; if TRUE then contour lines are overlaid on any plot}
  \item{levels}{numeric vector of contour levels}
  \item{text}{logical; if TRUE then cell probabilities are overprinted, rounded to 3 d.p.}
  \item{\dots}{other arguments passed to \code{\link{plot.mask}} (not used by summary method)}

  \item{object}{kernel object from \code{make.kernel}}

}

\details{

A kernel object is a type of mask with cell probabilities stored in the covariate `kernelp'. All kernels are truncated at kernelradius x spacing.

Options for \code{movementmodel} are
 \tabular{lll}{
 Kernel \tab Description \tab Parameters \cr
 normal \tab bivariate normal \tab move.a\cr
 exponential \tab negative exponential \tab move.a\cr
 t2D \tab 2-dimensional t-distribution cf Clark et al. (1999) \tab move.a, move.b \cr
 uniform \tab uniform within kernelradius, zero outside \tab (none) \cr
}

Movement models are further described in the \href{../doc/openCR-vignette.pdf}{vignette}.

}

\value{

\code{make.kernel} returns an object of class c('kernel','mask','data.frame').

}

\note{

This implementation in \pkg{openCR} 1.4.0 supercedes the function \code{\link{plotKernel}} that is now deprecated and will be removed in later versions.

}

\references{

Clark, J. S, Silman, M., Kern, R., Macklin, E. and HilleRisLambers, J. (1999) Seed dispersal near and far: patterns across temperate and tropical forests. \emph{Ecology} \bold{80}, 1475--1494.

Nathan, R., Klein, E., Robledo-Arnuncio, J. J. and Revilla, E. (2012) Dispersal kernels: review. In: J. Clobert et al. (eds) \emph{Dispersal Ecology and Evolution}. Oxford University Press. Pp. 187--210.
 
}

\seealso{

\code{\link{mask}}

}

\examples{

k <- make.kernel(movementmodel = 't2D', spacing = 10, move.a = 20, move.b = 1)
summary(k)

}

\keyword{ hplot }
