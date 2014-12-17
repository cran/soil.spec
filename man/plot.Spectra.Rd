\name{plot.Spectra}
\docType{methods}
\alias{plot,Spectra,ANY-method}
\alias{plot.Spectra}
\alias{plot}
\title{Plots spectral signatures}
\description{Plots spectral signatures. Each line represents a single sample.}
\usage{
\S4method{plot}{Spectra,ANY}(x, y, \dots)
}
\arguments{
  \item{x}{object of class \code{"Spectra"}}
  \item{y}{ignored}
  \item{\dots}{other optional arguments}
}
\details{
Wavenumbers are shown on the x-axis at the bottom. The second top x-axis displays the equivalent of wavenumbers interms of wavelengths. 
}
\author{
Andrew Sila and Tomislav Hengl}
\examples{
\dontrun{## Original binary Opus files:
pth = system.file(package = "soil.spec")
lst <- as.list(list.files(path=pth, pattern="*.0$", full.names=TRUE))
xx <- read.opus(lst[[1]])
plot(xx@data)
}
}