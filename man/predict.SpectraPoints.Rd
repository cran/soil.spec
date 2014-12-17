\name{predict.SpectraPoints-method}
\docType{methods}
\alias{predict,SpectraPoints-method}
\alias{predict,data.frame-method}
\alias{predict.SpectraPoints}
\title{Predict soil properties from an object of class \code{"SpectraPoints"}}
\description{Predicts soil properties from an object of class \code{"SpectraPoints"} (or a \code{"data.frame"} containing absorbances and sample IDs) either by using a referent \code{"SpectraModel"} (comes with the package) or by using the user created model. Extends, for example, the \code{"predict.mvr"} function available from the \href{http://cran.r-project.org/web/packages/pls/}{pls} package.}
\usage{ 
\S4method{predict}{SpectraPoints}(object, idcol = "SAMPLEID", model, 
   variable = model@variable, 
   output, validate = FALSE, model.class = "mvr",
   confidence.band = TRUE, prob. = .90, 
   signif.digits = 3, 
   st.wavenumbers=wavenumbers, 
   instr.range = c("ten-mir", "alp-mir", "mpa-nir")[1],
   \dots)
\S4method{predict}{data.frame}(object, idcol = "SAMPLEID", \dots)
}
\arguments{
  \item{object}{object of type \code{"SpectraPoints"}}
  \item{idcol}{\code{"character"}; ID column name} 
  \item{model}{\code{"SpectraModel"}; prediction model}   
  \item{variable}{\code{"character"}; variable name (see 'attributes in \code{spec.env(show.env = TRUE)})}
  \item{output}{\code{"SoilProfileCollection"}; optional output object where the predicted soil properties can be written}
  \item{validate}{\code{"logical"}; runs validation of overlap of points in feature space (see details)}
  \item{model.class}{\code{"character"}; prediction model class} 
  \item{confidence.band}{\code{"logical"}; specifies whether confidence bands should be derived for each prediction}   
  \item{prob.}{\code{"numeric"}; probability for confidence bands}
  \item{signif.digits}{\code{"integer"}; significant digits (all output numbers are typically rounded to three significant digits to save space)}
  \item{st.wavenumbers}{\code{"data.frame"}; standard wavenumbers with band names and upper / lower limits (usually predefined)}
  \item{instr.range}{\code{"character"}; instrument range}   
  \item{\dots}{other optional arguments}  
}
\details{When predicting from a \code{data.frame}, column names must contain wavenumbers and first column must contain sample IDs. Validation of overlap in feature space is done using PCs 1-3 derived from absorbance bands. This method basically compares for each new sample shortest distance (in 3D space; derived using \code{spatstat::nncross.pp3}) to the cloud of referent points and then prints warning if the distance is larger than two times the largest \code{nndist} (distance) for the referent point cloud.}
\references{
\itemize{ 
\item Mevik, B.-H., Wehrens, R. (2007) \href{http://www.jstatsoft.org/v18/i02/}{The pls Package: Principal Component and Partial Least Squares Regression in R}. Journal of Statistical Software 18(2), 1-24.
}
}
\author{ Tomislav Hengl and Andrew Sila}
\seealso{ \code{\link{SpectraPoints-class}}, \code{\link{fit.SpectraModel}}
}
\examples{
data(afspec)
## predict pH using absorbances directly:
data(m.PHIHOX)
x <- predict(afspec$ab, model=m.PHIHOX)
\dontrun{
## how good is this model?
## we look at the calibration data used to produce "m.PHIHOX"
pr <- m.PHIHOX@model$fitted.values
library(plyr)
library(scales)
data(afss)
afss.tbl <- join(afss$sites, afss$horizons)
tbl <- merge(data.frame(predicted=pr, SAMPLEID=attr(pr, "names")), 
   afss.tbl[,c("SAMPLEID","PHIHOX")])
v.r <- range(tbl$PHIHOX)
r.square <- summary(lm(as.formula(paste("PHIHOX ~ predicted - 1")),
   tbl))$adj.r.squared
plot(x=tbl$predicted, y=tbl$PHIHOX, asp=1, xlim=v.r, ylim=v.r, 
   col=alpha("darkgrey", 0.5), xlab="Predicted", ylab="Measured", 
   main=paste("R-square:", formatC(r.square, digits=3), sep=""))
abline(a=0, b=1, lw=2, col="black")
}
}
\keyword{methods}