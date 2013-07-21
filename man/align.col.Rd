\name{align.col}
\alias{align.col}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Combines two files and matches columns}
\description{Combines two files with common data even when columns are not arranged in the same order.
}
\usage{align.col(x,y)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{First file to be used for combination with the second}
  \item{y}{The other file to be used to combine with the first}
}
\details{Files comes with variables in different orders}

\value{
The original replicate spectra is returned together with averaged table.
A graphical output of the averaged rows is provided and can be saved from the menu if necessary
	\item{x.y}{combined file}
}


\author{
Andrew Sila}

\note{
If used on different form of data other than IR, no spectrums are generated and a warning message is given which should be ignored. Also, ensure that column names have prefix; data tables obtained from \code{read.spc} should first be saved and then read back that way numeric column names will be prefixed automatically. User can control prefix appended to the wavenumbers by naming them explicitly.}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Column}
\keyword{ Alignment }% __ONLY ONE__ keyword per line

\examples{
\dontrun{mir<-align.col(x,y)}}