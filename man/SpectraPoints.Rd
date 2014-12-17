\name{SpectraPoints-class}
\docType{class}
\alias{SpectraPoints-class}
\alias{SpectraPoints}
\alias{Spectra-class}
\alias{Spectra}
\alias{reshape.SpectraPoints}
\alias{summary,SpectraPoints-method}
\alias{validate,SpectraPoints-method}
\alias{validate.SpectraPoints}
\title{A class for spectral absorbance measurements}
\description{A class for spectral absorbance measurements typically in the medium infra-red range (MIR), visible near infra-red (VISMIR) and near-infra-red.}
\section{Slots}{
  \describe{
    \item{\code{metadata}:}{object of class \code{"data.frame"}; metadata table containg relevant information such as \code{"MID"} (metadata ID), \code{"Instrument_name"}, \code{"Instrument_URL"}, \code{"Laboratory_name"}, \code{"Laboratory_contact"} (contact person and his/her e-mail address), \code{"Laboratory_URL"}, \code{"Material_class"} (e.g. soil or vegetation), \code{"Wavenumber_conversion"}, \code{"Wavenlength_unit"}, \code{"Location_error"}}
    \item{\code{data}:}{object of class \code{"Spectra"}; contains absorbance values in wide or long format}
    \item{\code{sp}:}{object of class \code{"SpatialPoints"}; sampling locations}
  }
}
\section{Methods}{
\describe{
	\item{summary}{\code{signature(obj = "SpectraPoints")}: default summary of the object showing number of bands and points}
 	\item{validate}{\code{signature(obj = "SpectraPoints")}: check for overlap in feature space in comparison to the calibration data}
}
}
\details{The class expects by default absorbance, and not reflectance values. For large data sets we advise using the wide data format.
}
\author{ Tomislav Hengl and Andrew Sila }
\seealso{ \code{\link{fit.SpectraModel}}
}
\examples{ 
## generate the class:
library(sp)
library(rgdal)

data(afspec)
## get spatial coordinates:
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
plot(afspec.sp@data)
data(m.PHIHOX)
}
\keyword{classes}
