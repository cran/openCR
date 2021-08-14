\name{make.kernel}
\alias{make.kernel}
\alias{plot.kernel}
\alias{summary.kernel}

\title{Discrete Movement Kernel}

\description{

Functions to create, plot and summarise a discrete representation of a movement kernel.

}

\usage{

make.kernel(movementmodel = c("BVN", "BVE", "BVT","frE", "frG", "frL", 
    "uniform"), kernelradius = 10, spacing, move.a, move.b, 
    sparsekernel = FALSE, clip = FALSE, normalize = TRUE,     
    stat = c('estimate','lcl', 'ucl'), ...)

\method{plot}{kernel}(x, contour = FALSE, levels = NULL, text = FALSE, 
title = NULL, ...)

\method{summary}{kernel}(object, ...)

}
  
\arguments{

  \item{movementmodel}{character or function or openCR object}
  \item{kernelradius}{integer radius of kernel in grid cells}
  \item{spacing}{numeric spacing between cell centres}
  \item{move.a}{numeric parameter of kernel}
  \item{move.b}{numeric parameter of kernel}
  \item{sparsekernel}{logical; if TRUE then only cardinal and intercardinal axes are included}
  \item{clip}{logical; if TRUE then corner cells are removed}
  \item{normalize}{logical; if TRUE then cell values are divided by their sum}
  \item{stat}{character; predicted statistic to use for move.a (openCR object only) }
  
  \item{x}{kernel object from \code{make.kernel}}
  \item{contour}{logical; if TRUE then contour lines are overlaid on any plot}
  \item{levels}{numeric vector of contour levels}
  \item{text}{logical; if TRUE then cell probabilities are overprinted, rounded to 3 d.p.}
  \item{title}{character; if NULL a title is constructed automatically}
  \item{\dots}{other arguments passed to \code{\link{predict.openCR}} (\code{make.model}) or
      \code{\link{plot.mask}} (\code{plot.kernel}) (not used by summary method)}

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
 annular \tab non-zero only at centre and edge cells (after clipping at kernelradius) \tab move.a \cr
 annularR \tab non-zero only at centre and a ring of cells at radius R \tab move.a, move.b \cr
}

The \code{movementmodel} may also be a function or a previously fitted openCR model that includes movement. If a fitted openCR object, parameter values and kernel attributes are derived from that object and other arguments are ignored.

The parameter `move.a' is a scale parameter in metres, except for the `annular' models for which it is the proportion at the centre (probability of not moving). The movement model ``annularR'' uses a variable radius (R = move.b x kernelradius x spacing) and weights each cell according to the length of arc it intersects; ``annularR'' is not currently allowed in \code{openCR.fit}. 

'Sparse' kernels include only those grid cells that lie on 4 axes (N-S, E-W, NW-SE, NE-SW); cell probabilities are adjusted to maintain nearly the same distance distribution as the non-sparse equivalents.

Movement models are further described in the \href{../doc/openCR-vignette.pdf}{vignette}.

\code{plot.kernel()} adds an informative title with font size from the graphical parameter `cex.main'. Set \code{title = ""} to suppress the title.

Useful properties of theoretical (not discretized) kernels may be recovered with \code{\link{matchscale}}, \code{\link{pkernel}}, \code{\link{dkernel}} and \code{\link{qkernel}}.

}

\value{

\code{make.kernel} returns an object of class c('kernel','mask','data.frame').

The kernel object has attributes: 

\tabular{ll}{
Attribute \tab Description \cr
movementmodel \tab saved input \cr
K2 \tab saved kernelradius \cr
move.a \tab saved input \cr
move.b \tab saved input \cr
distribution \tab empirical cumulative distribution function \cr
}

The empirical cumulative distribution is a dataframe with columns for the sorted cell radii `r' and the 
associated cumulative probability `cumprob' (one row per cell).
}

\note{

The plot method for kernels supercedes the function \code{plotKernel} that has been removed.

}

\references{

Clark, J. S, Silman, M., Kern, R., Macklin, E. and HilleRisLambers, J. (1999) Seed dispersal near and far: patterns across temperate and tropical forests. \emph{Ecology} \bold{80}, 1475--1494.

Nathan, R., Klein, E., Robledo-Arnuncio, J. J. and Revilla, E. (2012) Dispersal kernels: review. In: J. Clobert et al. (eds) \emph{Dispersal Ecology and Evolution}. Oxford University Press. Pp. 187--210.
 
}

\seealso{

\code{\link{mask}}, \code{\link{matchscale}}, \code{\link{dkernel}}, \code{\link{pkernel}}, \code{\link{qkernel}}

}

\examples{

k <- make.kernel(movementmodel = 't2D', spacing = 10, move.a = 20, move.b = 1)
summary(k)

}

\keyword{ hplot }
