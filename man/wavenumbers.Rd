\name{wavenumbers}
\docType{data}
\alias{wavenumbers}
\title{Standard bands}
\description{Standard band names with upper and lower limits (wavenumbers)}
\usage{
data(wavenumbers)
}
\format{
  \describe{
    \item{BAND}{factor; band names}
    \item{UPPER}{numeric; upper value}
    \item{LOWER}{numeric; lower value}
  }
}
\details{
Standard wavenumbers are required by the \code{\link{read.opus}} function. Each \code{"SpectraModel"} might require different standard bands.
}
\author{Andrew Sila and Tomislav Hengl}
\examples{
data(m.PHIHOX)
## get bands of interest:
w <- attr(m.PHIHOX@model$opt.coef, "dimnames")[[1]]
nc <- as.numeric(sapply(w, function(x){gsub("[^0-9.]", "", x)}))
rn <- rank(nc)
if(any(duplicated(rn))){ warning("Duplicate bands detected") }
UPPER = c((nc[rn][-length(nc)]+nc[rn][-1])/2, nc[rn][length(nc)]+1)
wns <- data.frame(BAND=w, LOWER=c(nc[rn][1]-1, UPPER[-length(nc)]), UPPER=UPPER)
str(wns)
## compare with:
str(wavenumbers)
}
\keyword{data}