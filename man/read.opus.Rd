\name{read.opus-method}
\docType{methods}
\alias{read.opus}
\alias{read.opus,character-method}
\alias{read.opus,list-method}
\title{Reads binary OPUS files containing spectroscopy}
\description{Reads binary \href{http://www.bruker.com/products/infrared-near-infrared-and-raman-spectroscopy/opus-software/downloads.html}{OPUS} files containing  infrared, near infrared and Raman spectroscopy measurements (absorbances) and generates an object of class \code{"SpectraPoints"}.}
\usage{ 
\S4method{read.opus}{character}(file.name, sp = NULL,  
    codes = c("ZFF","RES","SNM","DAT","LWN","FXV","LXV","NPT","MXY","MNY","END","TIM"), 
    plot.spectra = FALSE, print.progress = FALSE, speclib = "ICRAF",
    signif.digit = get("signif.digit", spec.opts), MID, st.wavenumbers = wavenumbers)
\S4method{read.opus}{list}(file.name, \dots)
}
\arguments{
  \item{file.name}{\code{"character"}; file name with extension}
  \item{sp}{\code{"SpatialPoints"}; optional spatial object with coordinates of points with same sample ID's in the coordinates slot row names}
  \item{codes}{\code{"character"}; default OPUS codes}
  \item{plot.spectra}{\code{"logical"}; specifies whether to plot spectral curves every time a new sample is loaded}
  \item{print.progress}{\code{"logical"}; specifies whether to print progress}
  \item{signif.digit}{\code{"integer"}; number of significant digits}
  \item{speclib}{\code{"character"}; specifies which spectral library to base data points}
  \item{MID}{\code{"character"}; metadata ID}  
  \item{st.wavenumbers}{\code{"data.frame"}; standard wavenumbers with band names and upper / lower limits (usually predefined)}
  \item{\dots}{other optional arguments}  
}
\details{This function will read OPUS files from Bruker Optics' Alpha, MPA and Tensor-27 spectrometers. A choice for the type of spectral library to be created from the data tables created is done. To match data points to ICRAF's spectral library set \code{speclib = "ICRAF"}, otherwise to create own spectral library based on OPUS files being converted set \code{speclib="New"}.}
\note{Setting \code{speclib = "New"} may produce spectra with non-overlapping points brought about by slight drifts on the equipement and combinining with different columns can cause problems. If one chooses to create their own library, stable data points observed over time should be used to align all future scans.}
\note{Reading long list of binary files can be time and memory consuming.}
\references{
\itemize{ 
\item \href{http://www.bruker.com/products/infrared-near-infrared-and-raman-spectroscopy/opus-software/downloads.html}{Bruker Coorporation Guide for Infrared Spectroscopy}
}
}
\author{ Andrew Sila and Tomislav Hengl}
\seealso{ \code{\link{SpectraPoints-class}}, \code{\link{predict.SpectraPoints}}, \code{\link{wavenumbers}}
}
\examples{
## Original binary Opus files:
pth = system.file(package = "soil.spec")
lst <- as.list(list.files(path=pth, pattern="*.0$", full.names=TRUE))
file.info(lst[[1]])
xx <- read.opus(lst)
str(xx)
## predict pH
data(m.PHIHOX)
s.xx <- predict(xx, model = m.PHIHOX, prob. = .75)
s.xx
## Note: duplicate samples get unique name by default

## predict all standard soil properties:
nm <- get("attributes", spec.opts)
nm
pr.lst <- NULL
for(k in 1:length(nm)){
  data(list=paste("m.", nm[k], sep=""))
  try( pr.lst[[k]] <- predict(xx, variable=nm[k], 
  model = get(paste("m.", nm[k], sep="")), prob. = .75) )
}
pr <- do.call(cbind, pr.lst)
str(pr)
}
\keyword{methods}

