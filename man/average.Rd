\name{average}
\alias{average}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Average replicate spectral measurements}
\description{
Averages replicate spectral measurement but can also be applied to other data sets as long as the data table are organized in the form; first column with ID to show duplicate row and the remaining columns to be averaged over}
\usage{
average(spec)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{spec}{
Spec is a file with repliate data to be averaged, first column should include duplicated ids to identify rows to be averaged}
}
\details{
Only numeric columns will be averaged, otherwise \code{NA} values will be returned}

\value{
The original replicate spectra is returned together with averaged table.
A graphical output of the averaged rows is provided and can be saved from the menu if necessary
%%  If it is a LIST, use
 \item{rep.spec }{Replicated data before averaging}
 \item{averaged}{Data frame of the average}
 \item{Spectral signatures }{An overlay of all averaged spectraline plot }

}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
Andrew Sila}
\note{
If used on different form of data other than IR, no spectrums are generated and a warning message is given which should be ignored}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Replicate }
\keyword{ Averaged }% __ONLY ONE__ keyword per line
