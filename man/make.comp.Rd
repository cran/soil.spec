\name{make.comp}
\Rdversion{1.1}
\alias{make.comp}
\title{Making spectra compatible}
\description{
The function \code{make.comp} interpolates spectral absorbance/reflectance values to a new waveband vector by performing cubic spline interpolation using the \code{spline} function from the \code{stats} package.
}
\usage{
make.comp(spec, new, repl.out = "F", sav = "F", path = "", out= "Csm", save.as = "ws")
}
\arguments{
  \item{spec}{a numerical matrix containing the spectra.}
  \item{new}{a numeric vector containing the new wavebands.}
  \item{repl.out}{a logical value indicating whether extrapolated values outside the waveband range of \code{spec} should be replaced by \code{NA} values.}
  \item{sav}{a logical value indicating whether the function output (original and compatible spectra) should be saved.}
  \item{path}{a character vector giving the path where to save the function output. If \code{""} (default), the current working directory is taken.}
  \item{out}{a character vector giving the name of the function output in case \code{sav} is \code{"TRUE"}.}
  \item{save.as}{a character vector indicating the format of the saved output. \code{"ws"} saves the function output named with \code{output.name} as ws. \code{"csv.file"} saves the function output as csv-file.}
}
\details{
Some spectrometer adapt the wavenumber once in a while. Thus, waveband position might be different between spectra measured in different batches. This function makes the spectra compatible. 

Cubic spline interpolation allows to extrapolate values outside the original waveband range. The user can choose if the extrapolated values are replaced by \code{NA} values.

The column names of \code{spec} need to be numeric (no letters e.g. in the first position allowed).
}
\value{
\code{make.comp} returns a list with class \code{"make.comp"} containing the following components:
  \item{original.spectra}{a matrix giving the original spectra.}
  \item{compatible.spectra}{a matrix giving the compatible spectra.}
}
\author{Thomas Terhoeven-Urselmans}