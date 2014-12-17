\name{afspec}
\docType{data}
\encoding{latin1}
\alias{afspec}
\title{Sample soil spectroscopy data from Africa Sentinel Sites}
\description{Sample data (113 samples) containing absorbances in the Mid infrared region (MIR) collected at 57 sites (Africa Sentinel Sites).}
\usage{data(afspec)}
\format{
The \code{afspec} data set contains two data frames --- \code{"samples"} and \code{"ab"}. \code{"samples"} table contains the following columns:
  \describe{
  \item{\code{SAMPLEID}}{factor; unique soil sample ID (duplicate IDs are possible as same soil sample can be measured twice)}
  \item{\code{Cluster}}{integer; Africa Sentinel Sites cluster number}
  \item{\code{Plot}}{integer; Africa Sentinel Sites plot number}
  \item{\code{Depth}}{factor; lower depth of the sampled horizon}
  \item{\code{Site}}{factor; site name} 
  \item{\code{Country}}{factor; country name}
  \item{\code{Region}}{factor; region name}
  \item{\code{Longitude}}{numeric; longitude in decimal degrees on the WGS84 datum (\code{Lat} in the LDSF Field Database)}   
  \item{\code{Latitude}}{numeric; latitude in decimal degrees on the WGS84 datum (\code{Lat} in the LDSF Field Database)}
}
\code{"ab"} table contains the following columns:
  \describe{
  \item{\code{SAMPLEID}}{factor; ...}
  \item{\code{m4003.5..m599.8}}{numeric; mid-infrared raw absorbance for denoted wavenumbers (Bruker-LTS instrument)}
  \item{\dots}{other raw absorbance values}
}
}
\author{Africa Soil Information Service (\url{http://africasoils.net})}
\references{
\itemize{
\item Vagen, T., Shepherd, K. D., Walsh, M. G., Winowiecki, L., Desta, L. T., & Tondoh, J. E. (2010) \href{http://worldagroforestry.org/sites/default/files/afsisSoilHealthTechSpecs_v1_smaller.pdf}{AfSIS technical specifications: soil health surveillance}. World Agroforestry Centre, Nairobi, Kenya.
\item World Agroforestry Centre (\url{http://worldagroforestry.org})
}
}
\examples{
library(sp)
library(rgdal)
data(afspec)
sp <- afspec$samples[,c("Longitude","Latitude")]
coordinates(sp) <- ~Longitude+Latitude
proj4string(sp) <- CRS("+proj=longlat +datum=WGS84")
## prepare 'samples' table
samples <- cbind(afspec$samples["SAMPLEID"], 
    MID="AfSIS-MIR", 
    DateTime=Sys.time())
## convert to "SpectraPoints"
afspec.sp <- SpectraPoints(Spectra=Spectra(samples, afspec$ab), sp=sp)
summary(afspec.sp)
get("attributes", spec.opts)
## predict soil pH:
data(m.PHIHOX)
PHIHOX <- predict(afspec.sp, model=m.PHIHOX)
data(m.ORCCNS)
ORCCNS <- predict(afspec.sp, model=m.ORCCNS)
\dontrun{
library(plotKML)
afspec.pnts <- SpatialPointsDataFrame(sp, cbind(PHIHOX, ORCCNS))
data(SAGA_pal)
plotKML(afspec.pnts["PHIHOX"], colour_scale=rev(SAGA_pal[[1]]))
plotKML(afspec.pnts["ORCCNS"], colour_scale=SAGA_pal[[1]])
## plot with global soil legends:
library(GSIF)
data(soil.legends)
brks = c(soil.legends[["PHIHOX"]]$MIN[1]/10, soil.legends[["PHIHOX"]]$MAX/10)
pH.range <- cbind(soil.legends[["PHIHOX"]]$MAX, soil.legends[["PHIHOX"]]$MIN)
levs = paste(signif(rowMeans(pH.range)/10, 2))
afspec.pnts$Col <- cut(afspec.pnts$PHIHOX, breaks=brks, labels=levs)
plotKML(afspec.pnts["Col"], colour_scale=soil.legends[["PHIHOX"]]$COLOR)
}
}
\keyword{datasets}