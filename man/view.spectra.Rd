\name{view.spectra}
\alias{view.spectra}
\title{Generate spectral signatures}
\description{
Reads in a data table with infrared data and creates spectral signatures for each sample. }
\usage{
view.spectra(spec)
}
\arguments{
  \item{spec}{
}
}
\details{
Produces spectral signatures for all records found in a table.
Ensure that infrared regions are the columns names of the table from the second column }

\author{
Andrew Sila}
\note{
Ensure the first column contains spectrum ids and the rest contains infrared data.}

\examples{
\dontrun{spec<-read.csv(file.choose())}
\dontrun{view.spectra(spec)}
}
\keyword{spectrum}
