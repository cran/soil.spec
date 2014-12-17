\name{spec.env}
\alias{spec.env}
\alias{spec.opts}
\title{Soil.spec package specific environmental variables}
\description{Generates internal environment / package specific parameters and settings that can be later on passed to other functions.}
\usage{spec.env(MIR = c(390, 7500), NIRS = c(3900, 12500), 
   NIRP = c(4000,10000), VISNIR1 = c(420, 960), 
   VISNIR2 = c(1020, 1770), 
   VISNIR3 = c(1830, 2480),
   icraf.htsxt = c(3578, 7497.964, 599.76), 
   icraf.alpha = c(2542, 3998.12872, 399.387991), 
   icraf.mpa = c(2307, 12493.2, 3598.69),  
   CO2.band = c(2350.8,2379.8), signif.digit = 5, 
   attributes = c("ORCCNS", "PHIHOX", "ALUM3S", 
   "ECAM3S", "EXKM3S", "EMGM3S", "ENAM3S", "EXB", 
   "NITCNS", "SNDLDF"), 
   mdnames = c("MID", "Instrument_name", "Instrument_URL", 
   "Laboratory_name", "Laboratory_contact", "Laboratory_URL",
   "Material_class", "Wavenumber_conversion", "Wavenlength_unit",
   "Location_error"), show.env = FALSE)
}
\arguments{
  \item{MIR}{\code{"numeric"}, mid-infrared part of spectra}
  \item{NIRS}{\code{"numeric"}, near infrared part of spectra}
  \item{NIRP}{\code{"numeric"}, near infrared part of spectra}
  \item{VISNIR1}{\code{"numeric"}, visible near infrared part of spectra}
  \item{VISNIR2}{\code{"numeric"}, visible near infrared part of spectra}
  \item{VISNIR3}{\code{"numeric"}, visible near infrared part of spectra}
  \item{icraf.htsxt}{\code{"numeric"}, ICRAF MIR part of spectra}
  \item{icraf.alpha}{\code{"numeric"}, ICRAF alpha-MIR part of spectra}
  \item{icraf.mpa}{\code{"numeric"}, ICRAF VISNIR part of spectra}
  \item{CO2.band}{\code{"numeric"}, wavenumbers for the CO2 band}
  \item{signif.digit}{\code{"integer"}, default rounding system for absorbances}
  \item{attributes}{\code{"character"}, standard target variables of interest: \code{"ORCCNS"} (Organic carbon content in soil estimated using the CNS elemental analyzer), \code{"PHIHOX"} (soil pH measured in water), \code{"ALUM3S"} (Aluminium concentration estimated using the Mehlich 3 solution), \code{"ECAM3S"} (exchangeable Calcium determined using using the Mehlich 3 solution), \code{"EXKM3S"} (exchangeable Potassium determined using using the Mehlich 3 solution), \code{"EMGM3S"} (exchangeable Magnesium determined using using the Mehlich 3 solution), \code{"ENAM3S"} (exchangeable Sodium determined using using the Mehlich 3 solution), \code{"EXB"} (sum of exchangeable bases), \code{"NITCNS"} (total nitrogen estimated using estimated using the CNS elemental analyzer), \code{"SNDLDF"} (weight percentage of particles 50-2 mm size determined using the laser diffraction method)}
  \item{mdnames}{\code{"character"}, standard metadata names}
  \item{show.env}{\code{"logical"}, specifies whether the output environment should be printed}
}
\author{Andrew Sila and Tomislav Hengl}
\examples{
# environmental variables:
spec.env()
get("mdnames", envir = spec.opts)
get("attributes", envir = spec.opts)
}
\keyword{options}