\name{trans}
\Rdversion{1.1}
\alias{trans}
\alias{plot.trans}
\title{Spectral transformation}
\description{
\code{trans} transforms spectra either by calculation of (i) derivatives, (ii) continuum removal or (iii) wavelet transform. The used functions are  (i) \code{locpoly} function in \code{"KernSmooth"} package, (ii) \code{chull} and \code{approx} functions in \code{"KernSmooth"} package and (iii) \code{dwt} function in \code{"wavelets"} package , respectively.
}
\usage{
trans(raw, tr = "derivative", order = 1, gap = 21)

\method{plot}{trans}(x,...)
}
\arguments{
  \item{raw}{a matrix containing the raw spectra.}
  \item{tr}{a character naming the transformation method. See details.}
  \item{order}{an integer between 0 and 3 giving the order of derivative. The value 0 performs smoothing based on the \code{gap}.}
  \item{gap}{an integer between 1 and 30 defining the smoothing interval in wavebands.}
 \item{x}{an object of class \code{"trans"}.}
 \item{...}{additional arguments.}
}
\details{
\code{spec.type} uses for spectral transformation (i) the \code{locpoly} function in \code{KernSmooth} package for derivative calculation, (ii) the \code{chull} and \code{approx} functions in \code{"KernSmooth"} package for continuum removal and (iii) the \code{dwt} function in \code{wavelets} package for extraction of wavelet coefficients. Experiences showed for wavelet decomposition that the best ratio of prediction performance and sparse spectral representation is reached when all 128 wavelet coefficients from decomposition level three are taken.

Possible options for \code{tr} are \code{"derivative"}, \code{"continuum removed"} or \code{"wavelet transformed"}.

%%A graph is returned showing the raw and transformed spectra for comparison. Its the same output than from \code{\link{plot.trans}}.
}
\value{
\code{trans} returns a list with class \code{"trans"} containing the following components:
  \item{raw}{a matrix containing the raw spectra.}
  \item{trans}{a matrix containing the transformed spectra.}
  \item{transformation}{a character naming the transformation method.}
}
\author{Thomas Terhoeven-Urselmans}