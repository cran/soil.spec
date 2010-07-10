\name{isric}
\Rdversion{1.4}
\alias{isric}
\title{Script for publication Terhoeven-Urselmans et al. 2010.}
\description{
The script for the publication: Terhoeven-Urselmans,T., Vagen,T.G., Spaargaren,O. and Shepherd,K.D. 2010. Prediction of soil fertility properties from a globally distributed soil mid-infrared spectral library. Soil Science Society of America Journal. Accepted, doi:10.2136/sssaj2009.0218.
}
\usage{
isric(spec, chem, plo)
}
\arguments{
  \item{spec}{a numerical matrix containing the raw spectra.}
  \item{chem}{a numerical data frame containing the reference values.}
  \item{plo}{a data frame containing two columns. See details}
}
\details{
The first column of \code{plo} named \dQuote{Plot_code} contains the plot code as factor. The second column named \dQuote{Depth} contains the depth-layer as factor. An example is \dQuote{0 to 14} and means that this soil sample was taken from 0 to 14 cm depth.
}
\value{
\code{isric} returns a list with class \code{"isric"} containing the following components:
  \item{Table_1}{a matrix containing the summary of principal component analysis of first derivative spectra. Given as Table 1 in the publication.}
  \item{Table_2_summary}{a table containing the summary of \code{chem}. Given as Table 2 in the publication.}
  \item{Table_2_sd}{a numeric vector containing the standard deviation of \code{chem}. Part of Table 2 in publication}
  \item{Table_3_calibration}{a matrix containing the calibration set statistics. Given as Table 3 in the publication.}
  \item{Table_3_validation}{a matrix containing the validation set statistics. Given as Table 3 in the publication.}
  \item{Table_4_all}{a list containing the correlation between soil property values of the whole data set. Given as upper triangle in Table 4 in the publication.}
  \item{Table_4_calibration}{a list containing the correlation between soil property values of the calibration data set only. Given as lower triangle in Table 4 in the publication.}
  \item{Calibration_predicted}{a data frame containing for the soil property values of the calibration set the (i) untransformed values, (ii) box-cox transformed values, (iii) predicted box-cox transformed values and (iv) back-transformed predicted values.}
  \item{Validation_predicted}{a data frame containing for the soil property values of the validation set the (i) untransformed values, (ii) box-cox transformed values, (iii) predicted box-cox transformed values and (iv) back-transformed predicted values.}
}
\author{Thomas Terhoeven-Urselmans}