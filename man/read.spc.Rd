\name{read.spc}
\Rdversion{1.3}
\alias{read.spc}
\title{Reads spectral spc-files into R}
\description{
\code{read.spc} reads binary spectral spc-files from a folder into R. The spectra can be made compatible (see details in \code{\link{make.comp}}) either to the first sample wavebands or to the standard wavebands of the ICRAF spectral lab. Information from the scanning method is gathered to check on spectral comparability.
}
\usage{
read.spc(loa = "NULL", save.path = "NULL", sav = "FALSE", output.name = "Spectral matrix", save.as = "workspace", wavenumber = "first.sample")
}
\arguments{
  \item{loa}{a character giving the path name where the spc-files are stored. If \code{"NULL"} (default), the spc-files are stored in the current working directory.}
  \item{save.path}{a character giving the path name where the function output shall be saved.}
  \item{sav}{a logical value indicating whether the function output shall be saved.}
  \item{output.name}{a character giving the function output name, in case \code{sav} is \code{"TRUE"}.}
  \item{save.as}{a character vector indicating the format of the saved output. \code{"workspace"} saves the function output named with \code{output.name} as workspace. \code{"csv.file"} saves the function output as csv-file.}
  \item{wavenumber}{a character giving the way how the samples should be made compatible. If \code{"first.sample"}, all spectra are made compatible to the first read sample. if \code{"ICRAF"}, all spectra are made compatible to the standard wavebands of the ICRAF spectral lab.}
}
\details{
Spectra from the near- and mid-infrared range can be read. In case the spc-files saved in \code{loa} comprise both ranges the function output is given separate for each range.

The function allows to read only spectra in one go, when they have the same material (e.g. soil or plant), were scanned with the same resolution and have the same zero filling. If there are still small differences in the number of wavebands, the spectra are made compatible depending on the argument \code{wavenumber}. In case spectra with different materials shall be read, the user has to decide the material via graphical interface.
}
\value{
\code{read.spc} returns a list with class \code{"read.spc"} containing the following components: 
  \item{spectra}{a numerical matrix containing the read spectra.}
  \item{additional.information}{a data frame containing some scanning method details.}
}
\author{Thomas Terhoeven-Urselmans}