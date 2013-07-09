\name{align}
\alias{\S4method{new}{align}}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Combine two tables with columns in unmatching order }
\description{
This tools gives a quick way of combining two tables whose columns are arranged in different order. Decide between the two files which should be treated as the master or the reference and the second to be reorganised to the order of the first before combinining. Second can contain fewer columns than the first}
\usage{
\S4method{new}{align}(file1, file2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file1}{
%%     ~~Describe \code{file1} here~~
}
  \item{file2}{
%%     ~~Describe \code{file2} here~~
}
}
\details{
An object of class \code{file1},\code{file2}, can merged to obtain more records, especially when there new obseravtions or measurements have been generated and should be added to an existing file.}
\value{
If there is at least one common column name, the files merged after aligning the columns in the same order to get a table with more rows
\item{combined}{File1 and file2 put together rowwise}

}
\references{
Unpublished scripts by the author}
\author{
Andrew Sila}
\note{
This is function is included here to help users who would like to spend little time in manipulating data files}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
\examples{
#Read two files;
#file1<-read.csv(file.choose())
#file2<-read.csv(file.choose())

#Call align function
#combined<-align(file1,file2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{reference}
\keyword{master}% __ONLY ONE__ keyword per line
