\name{plotKernel}
\alias{plotKernel}

\title{Plot Movement Kernel}

\description{

Movement between primary sessions is modelled in \pkg{openCR} with a discretized kernel. Each cell of the kernel contains the probability of movement from the central cell. Kernels are `normal' (Gaussian), `exponential' (negative exponential)  or specified with a user-provided function. This function allows you to preview a kernel specification.

}

\usage{

plotKernel(movementmodel = c("normal", "exponential"), kernelradius = 10, spacing, pars, 
    clip = FALSE, plt = TRUE, contour = FALSE, levels = NULL, text = FALSE, ...)

}

\arguments{
  \item{movementmodel}{character or function}
  \item{kernelradius}{integer radius of kernel in grid cells}
  \item{spacing}{numeric spacing between cell centres}
  \item{pars}{numeric vector of 1 or 2 parameter values}
  \item{clip}{logical; if TRUE then corner cells are removed}
  \item{plt}{logical; if TRUE then a plot is produced}
  \item{contour}{logical; if TRUE then contour lines are overlaid on any plot}
  \item{levels}{numeric vector of contour levels}
  \item{text}{logical; if TRUE then cell probabilities are overprinted, roumnded to 3 d.p.}
  \item{\dots}{other arguments passed to plot.mask (e.g. breaks)}
}

\details{

Internally, a mask is generated with kernel probabilities in a covariate, and plotting is done with \code{\link{plot.mask}}.

}

\value{

A dataframe with columns x, y, and kernelp is returned invisibly.

}

\examples{

plotKernel(spacing = 2, k = 10, pars = 10, contour = TRUE, clip = TRUE)

}

\keyword{ hplot }