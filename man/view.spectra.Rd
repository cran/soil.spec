\name{view.spectra}
\alias{view.spectra}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate spectral signatures}
\description{
Reads in a data table with infrared data and creates spectral signatures for each sample. }
\usage{
view.spectra(spec)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{spec}{
%%     ~~Describe \code{spec} here~~
}
}
\details{
Produces spectral signatures for all records found in a table.
Ensure that infrared regions are the columns names of the table from the second column }
%%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%%}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
Andrew Sila}
\note{
Ensure the first column contains spectrum ids and the rest contains infrared data.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
\examples{
#spec<-read.csv(file.choose())
#view.spectra(spec)
}
\keyword{spectrum}
