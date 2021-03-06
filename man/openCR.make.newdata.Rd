\name{openCR.make.newdata}
\alias{openCR.make.newdata}

\title{ Create Default Design Data }

\description{

Internal function used to generate a dataframe containing design data
for the base levels of all predictors in an openCR object.

}

\usage{

openCR.make.newdata(object, all.levels = FALSE)

}

\arguments{

\item{object}{fitted openCR model object}
\item{all.levels}{logical; if TRUE then all covariate factor levels appear in the output}

}

\details{


\code{openCR.make.newdata} is used by \code{predict} in lieu of
user-specified `newdata'. There is seldom any need to call
\code{openCR.make.newdata} directly.

}

\value{

A dataframe with one row for each session, and columns for the
predictors used by \code{object$model}.
  
}

\seealso{
  
\code{\link{openCR.fit}}
  
}

\examples{

## null example (no covariates)
ovenCJS <- openCR.fit(ovenCH)
openCR.make.newdata(ovenCJS)

}

\keyword{ models }
