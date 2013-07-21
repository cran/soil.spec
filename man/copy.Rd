 \name{copy}
\alias{copy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Copying spectral data files}
\description{
A quick tool to determine newly acquired files kept in one computer and copied to another.Files created from a specified period of time and stored in source A are identified and copied to source B.}
\usage{
copy(source, destination, ddmonyyyy)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{source}{
\code{source} This can refer to a pathname found on a local disk or on a remote computer where the data to be copied are stored
}
  \item{destination}{
Same as in \code{source}, but this is where the files meeting the criteria of \code{ddmonyyyy} should be copied to}
  \item{ddmonyyyy}{
Refers to the day, month and year from when files were created; for examples to search for files created  from 30th April 2013 specify as: 30apr2013 }
}
\details{
For the copying to be successful ensure you have right permissions to read from the source path and write permissions to copy files to the specified destination folder}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{Files }{List of files meeting the criterion specified here will be copied to the destination folder}
\item{csvfile}{A csv  file with a list of all the copied files is saved to the destination file indicating the time of creation for each file}
%% ...
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
Andrew Sila}
\note{
To test the speed of writing , try with recent files say those created in the last one to two months}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
\dontrun{new.files<-copy("d:/Measurements/soil","g:/New scans",12mar2013)}
}

\keyword{ source }
\keyword{ destination}% __ONLY ONE__ keyword per line
