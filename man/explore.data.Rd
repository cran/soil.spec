\name{explore.data}
\alias{explore.data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Data summary}
\description{
Another way of getting data summary statitics and organized in a user friendly way- tabular form}
\usage{
explore.data(datafile)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datafile}{
%%     ~~Describe \code{datafile} here~~
}
}
\details{
This takes a table and returns the total number of complete cases (excluding missing values), mean,minimum, maximum and standard deviation values for numeric columns in a data table. The output is well formated for a quick and useful exploration view}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{complete.cases}{Number of complete observations}
\item{Average}{Mean of each numeric column}
\item{Min }{Minimum value of each numeric column}
\item{Max }{Maximum value of each numeric column}
\item{Sd}{Standard deviation value of each numeric column}

}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
Andrew Sila}
\note{
If all columns are non-numeric the complete cases will be returned only}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
\examples{

\dontrun{#Start by reading a file with data to be explored}
	\dontrun{soil<-read.csv(file.choose())}
	\dontrun{explore.data(soil)}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{complete.cases}
\keyword{Averages}% __ONLY ONE__ keyword per line
