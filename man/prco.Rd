\name{prco}
\Rdversion{1.3}
\alias{prco}
\alias{plot.prco}
\alias{summary.prco}
\title{Principal component analysis}
\description{
\code{prco} uses the \code{prcomp} function in the \code{stats} package. In addition the important number of principal components is determined.
}
\usage{
prco(x, sav = "FALSE", save.path = "NULL", output.name = "prco")

\method{plot}{prco}(x,...)

\method{summary}{prco}(object,...)
}
\arguments{
  \item{x}{a numerical, spectral matrix in \code{prco}. An object of class \code{"prco"} in \code{plot.prco}.}
  \item{sav}{a logical value indicating whether the function output (original and compatible spectra) should be saved.}
  \item{save.path}{a character vector giving the path where to save the function output. If \code{"NULL"} (default), the current working directory is taken.}
  \item{output.name}{a character vector giving the name of the function output in case \code{sav} is \code{"TRUE"}.}
  \item{...}{additional arguments.}
 \item{object}{an object of class \code{"prco"}.}
}
\details{
The number of important principal components is selected so that the increase in cumulative explained variance within the next three components is lower than 4 percent.

The \code{plot.prco} function plots two plots: (i) the score values of the important principal components using the \code{pairs} function and (ii) the loading values against the wavebands for the first two principal components.

The \code{summary.prco} function prints the cumulative explained variance of the first ten principal components and the number of important principal components.
}
\value{
\code{prco} returns a list with class \code{"prco"} containing the following components:
  \item{prcomp}{a list of class \code{"prcomp"}.}
  \item{prco}{a numerical matrix containing the score values of the important principal components.}
  \item{i.pc}{numeric indicating the number of important principal components.}
}
\author{Thomas Terhoeven-Urselmans}