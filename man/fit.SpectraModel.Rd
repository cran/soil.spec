\name{fit.SpectraModel-method}
\docType{methods}
\alias{fit}
\alias{fit,formula,SpectraPoints,data.frame-method}
\alias{fit,list,SpectraPoints,data.frame-method}
\alias{fit.SpectraModel}
\alias{SpectraModel-class}
\alias{SpectraModel}
\alias{m.ORCCNS}
\alias{m.PHIHOX}
\alias{m.ALUM3S}
\alias{m.ECAM3S}
\alias{m.EXKM3S}
\alias{m.EMGM3S}
\alias{m.ENAM3S}
\alias{m.EXB}
\alias{m.NITCNS}
\alias{m.SNDLDF}
\title{Fits a calibration model for soil spectroscopy}
\description{Fits a calibration model that can be used to convert absorbances to soil properties. Implements the \code{"plsr"} function available from the \href{http://cran.r-project.org/web/packages/pls/}{pls} package.}
\usage{ 
\S4method{fit}{formula,SpectraPoints,data.frame}(formulaString, sampled, reference, 
   idcol = "SAMPLEID", ab.r = get("MIR", envir=spec.opts), 
   CO2.band = get("CO2.band", envir=spec.opts), ncomp, ncomp.max=15, 
   repl = 5, segment.type = "interleaved", prefix = "X", \dots)
\S4method{fit}{list,SpectraPoints,data.frame}(formulaString, sampled, reference, \dots)
}
\arguments{
  \item{formulaString}{\code{"formula"}; formula string}
  \item{sampled}{\code{"SpectraPoints"}; absorbances}   
  \item{reference}{\code{"data.frame"}; reference values estimated in laboratory using at the order of magnitude more precise technique}
  \item{idcol}{\code{"character"}; column with sample IDs}  
  \item{ab.r}{\code{"numeric"}; natural range of wavenumbers of interest}
  \item{CO2.band}{\code{"numeric"}; range of wavenumbers for the CO2 band}
  \item{ncomp}{\code{"integer"}; number of components (if not specified will be determined using cross-validation)}
  \item{ncomp.max}{\code{"integer"}; maximum number of components to be used for model fitting}
  \item{repl}{\code{"integer"}; number of repetitions for \code{"pls::plsr"}}
  \item{segment.type}{\code{"character"}; segment type for \code{"pls::plsr"}}
  \item{prefix}{\code{"character"}; default prefix to be added to the band names}  
  \item{\dots}{other optional arguments}  
}
\details{Returns an object of class \code{"SpectraModel"}, which contains slots: (1) \code{"variable"} (variable name according to the global soil data registry), (2) \code{"Space"} (3D space defined by the 1--3 principal components for absorbances), and (3) \code{"model"} (fitted model; usually of type \code{pls::plsr} and as an output from the partial least squares regression). Fitting models for large data sets can be time and memory consuming.}
\references{
\itemize{ 
\item Mevik, B.-H., Wehrens, R. (2007) \href{http://www.jstatsoft.org/v18/i02/}{The pls Package: Principal Component and Partial Least Squares Regression in R}. Journal of Statistical Software 18(2), 1-24.
\item Africa Soil Information Service (\url{http://africasoils.net})
}
}
\author{ Andrew Sila and Tomislav Hengl }
\seealso{ \code{\link{predict.SpectraPoints}}, \code{\link{spec.opts}}
}
\examples{
\dontrun{
## build a callibration model using AfSIS data:
library(plyr)
library(sp)
library(rgdal)

data(afss)
afss.tbl <- join(afss$sites, afss$horizons)
con <- url("http://gsif.isric.org/lib/exe/fetch.php?media=afss.spec.rda")
load(con)
afss.spec.ab <- Spectra(samples=data.frame(SAMPLEID=afss.spec$SAMPLEID, 
   MID="ICR_SOIL", DateTime=Sys.time()), ab=afss.spec)
sel.r <- !afss.tbl$LONWGS84==100
sp <- data.frame(afss.tbl[sel.r,c("LONWGS84","LATWGS84")])
coordinates(sp) <- ~LONWGS84+LATWGS84
proj4string(sp) <- CRS("+proj=longlat")
attr(sp@coords, "dimnames")[[1]] <- afss.tbl$SAMPLEID[sel.r]
## extend to SprectraPoints:
afss.spec.sp <- SpectraPoints(Spectra=afss.spec.ab, sp=sp)
## select bands of interest:
cutspec <- c(600,2350.8,2379.8,4000)
col.no <- sapply(names(afss.spec.sp@data@ab)[-1], 
   function(x){as.numeric(strsplit(x, "X")[[1]][2])})
sel <- col.no>cutspec[1] & col.no<cutspec[2] | 
       col.no>cutspec[3] & col.no<cutspec[4]
summary(sel)
formulaString <- as.formula(paste("na.omit(log(ORCDRC)) ~ ", 
   paste(names(afss.spec.sp@data@ab[-1])[sel], collapse="+")))
m.ORC <- fit(formulaString, sampled=afss.spec.sp, reference=afss.tbl) 
str(m.ORC)
}
## this model, in fact, is included in the package:
data(m.ORCCNS)
str(m.ORCCNS)
}
\keyword{methods}