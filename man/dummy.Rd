% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools.R
\name{dummy}
\alias{dummy}
\title{Transform a factor in a set of dummy variables}
\usage{
dummy(x, ..., keep = FALSE, prefix = NULL, ref = FALSE)
}
\arguments{
\item{x}{a data frame}

\item{...}{series of the data frame, should be factors}

\item{keep}{a boolean, if \code{TRUE}, the original series is kept in
the data frame,}

\item{prefix}{an optional prefix for the names of the computed
dummies,}

\item{ref}{a boolean, if \code{TRUE}, a dummy is created for all the
levels, including the reference level}
}
\value{
a data frame
}
\description{
The normal way to store cathegorical variables in R is to use
factors, each modality being a level of this factor. Sometimes
however, is is more convenient to use a set of dummy variables.
}
\examples{
charitable \%>\% dummy(religion, education)
}
\keyword{misc}
