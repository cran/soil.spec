\name{average}
\alias{average}
\Rdversion{1.1}
\title{
Averages data with replicate records}
\description{
Replicate data can be for any type of numeric values. An option to specify if data is spectra or non-spectra is required for this function to determine if to produce spectral signatures of averaged records or not incase of non-spectral data.
}
\usage{
average(spec, type="")
}
\arguments{
  \item{spec}{
Spec is a file with replicate data to be averaged, first column should include duplicated ids to identify rows to be averaged}

  \item{type}{
Specify if file with replicate data has spectral data type or not}
}
\details{
Only numeric columns will be averaged, otherwise \code{NA} values will be returned}

\value{
The original replicate spectra is returned together with averaged table.
A graphical output of the averaged rows is provided and can be saved from the menu if necessary

   \item{rep.spec}{Replicated data before averaging}
  \item{averaged}{Data frame of the average}


}

\author{
Andrew Sila}
\note{
1. If used on different form of data other than IR, no graphics for spectrums are generated}

\note{
2. Ensure that column names have prefix; data tables obtained from \code{read.spc} should first be saved and then read back that way numeric column names will be prefixed automatically. User can control prefix appended to the wavenumbers by naming them explicitly.}
\examples{
\dontrun{spec<-read.csv(file.choose())}
\dontrun{average(spec,type="spectra")}
\dontrun{other<-read.csv(file.choose())}
\dontrun{average(other,type="")}
}

\keyword{Replicate }
\keyword{ Averaged }% __ONLY ONE__ keyword per line
