\name{read.spc}
\Rdversion{1.1}
\alias{read.spc}
\title{Reads spectral spc-files into R}
\description{
\code{read.spc} reads binary spectral spc-files from a folder into R. The spectra can be made compatible (see details in \code{make.comp}) either to the first sample wavebands or to the standard wavebands of the ICRAF spectral lab. Information from the scanning method is gathered to check on spectral comparability. The default has been set to ICRAF spectral bands.
}
\usage{
read.spc(loa = "", path = "", sav = "F", out = "Sm", save.as = "workspace", wn = "first")
}
\arguments{
  \item{loa}{a character giving the path name where the spc-files are stored. If \code{""} (default), the spc-files are stored in the current working directory.}
  \item{path}{a character giving the path name where the function output shall be saved.}
  \item{sav}{a logical value indicating whether the function output shall be saved.}
  \item{out}{a character giving the function output name, in case \code{sav} is \code{"T"}.}
  \item{save.as}{a character vector indicating the format of the saved output. \code{"workspace"} saves the function output named with \code{out} as workspace. \code{"csv.file"} saves the function output as csv-file.}
  \item{wn}{a character giving the way how the samples should be made compatible. If \code{"first"}, all spectra are made compatible to the first read sample. if \code{"ICRAF"}, all spectra are made compatible to the standard wavebands of the ICRAF spectral lab.}
}
\details{
Spectra from the near- and mid-infrared range can be read. In case the spc-files saved in \code{loa} comprise both ranges the function output is given separate for each range.

The function allows to read only spectra in one go, when they have the same material (e.g. soil or plant), were scanned with the same resolution and have the same zero filling. If there are still small differences in the number of wavebands, the spectra are made compatible depending on the argument \code{wn}. In case spectra with different materials shall be read, the user has to decide the material via graphical interface.
}
\value{
\code{read.spc} returns a list with class \code{"read.spc"} containing the following components: 
  \item{spectra}{a numerical matrix containing the read spectra.}
  \item{additional.information}{a data frame containing some scanning method details.}
   \item{raw.spectra}{a data frame merging additional.information and spectra; columns 2 to 4 from the additional.information table are excluded to conform to the existing structure of library in existence}
}
\author{Andrew Sila and Thomas Terhoeven-Urselmans}

\examples{
\dontrun{mir<-read.spc(loa="D:/AfSIS/spc files",sav="F",wn="ICRAF")}
\dontrun{raw<-mir$raw.spectra #To extract the object with both additional.information and spectra}
\dontrun{raw[1:4,1:4] #This gives a snapshot of the first four rows and the first columns}
}

