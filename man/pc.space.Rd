\name{pc.space}
\Rdversion{1.1}
\alias{pc.space}
\alias{plot.pc.space}
\title{Projection of new spectra into a principal component space}
\description{New spectra are projected into the principal component space of a user given spectral sample set in order to check if the new spectra belong to the same population. The \code{prcomp} function from the \code{stats} package is used.}
\usage{
pc.space(base, new)

\method{plot}{pc.space}(x,...)
}
\arguments{
  \item{base}{a numerical spectral matrix containing the population samples.}
  \item{new}{a numerical matrix containing the new samples, which are checked being part of the population of \code{base}.}
  \item{x}{an object of class \code{"pc.space"}.}
  \item{...}{additional arguments.}
}
\details{
A principal component analysis is performed on \code{base} and the spectra from \code{new} are projected/predicted into this score space. The output (the same then using \code{\link{plot.pc.space}}) is plotted that the user can see if the spectra in \code{new} (in red) belong to the population of \code{base}.

\code{plot.pc.space} does the plotting of an object of class \code{"pc.space"}.
}
\value{
\code{pc.space} returns a list with class \code{"pc.space"} containing the following components: 
  \item{prco.base}{a matrix containing the score values of \code{base}.}
  \item{prco.new}{a matrix containing the score values of \code{new}.}
}
\author{Thomas Terhoeven-Urselmans}