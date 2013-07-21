\name{regr}
\Rdversion{1.1}
\alias{regr}
\alias{plot.regr}
\alias{pred.regr}
\alias{summary.regr}
\title{Regression functions for soil spectral analysis}
\description{
The next update of this package will fix some problems observed for this method which gives different regression functions (partial least-squares, boosted regression trees, support vector machines) can be chosen for calibrations of one or more constituents. Function settings are optimized for soil spectral analysis, but can be varied. Possible spectral transformations are described in the \code{\link{trans}} function.

}
%This model needs to be brokendown into small pieces manageable in R
%\usage{
%regr(x, y, sw = "NL", sp = "wt", r = "svm", p = "T", pn = 0.3, n, md = "test", drv = 1, bandwidth = 21, val= "CV", ft = "h", lv = 3, dis = "g", nt = 1000, s = 0.01, k = "ra")

%\method{plot}{regr}(x,...)
%\method{summary}{regr}(object,...)
%\method{pred}{regr}(new,model,output.name="test",s="NL")
%}

\arguments{
 \item{x}{a numerical data.frame or matrix containing the raw spectra in \code{regr}. An object of class \code{"regr"} in \code{plot.regr}.}
  \item{y}{a numerical data.frame or matrix containing the constituents.}
  \item{sw}{a character vector giving the path where to save the function output. If \code{"NL"} (default), the current working directory is taken. As well used in \code{pred.regr}.}
  \item{sp}{a character giving the desired spectral transformation. Available are \code{"raw"} (raw spectra), \code{"derivative"} (derivative spectra), \code{"continuum removed"} (continuum removed spectra) and  \code{"wt"} (wavelet coefficients).}
  \item{r}{a character giving the regression method. Available are \code{"pls"} (partial least-squares), \code{"brt"} (boosted regression trees) and \code{"svm"} (support vector machines).}
  \item{p}{a logical indicating whether validation samples should be chosen as a percentage from \code{x} (given in \code{pn}). If \code{"FALSE"} object \code{n} is taken.}
  \item{pn}{a numeric between 0 and 1 giving the percentage of validation samples to choose.}
  \item{n}{an integer giving the number of validation samples when \code{p} is \code{"FALSE"}.}
  \item{md}{a character naming the model output.}
  \item{drv}{an integer between 0 and 3 giving the order of derivative. The value 0 performs smoothing based on \code{bandwidth}.}
  \item{bandwidth}{an integer between 1 and 30 defining the smoothing interval in wavebands.}
  \item{val}{a character defining the type of cross-validation procedure when \code{r} is equal to \code{"pls"}. Available are \code{"none"} (no cross-validation procedure), \code{"CV"} (cross-validation in 10 segments) and \code{"LOO"} (leave-one-out cross-validation).}
  \item{ft}{a character defining the wavelet ft in \code{dwt} function from \code{wavelets} package.}
  \item{lv}{a character defining the lv of wavelet coefficients extraction (1 to 10 possible; 1 yields 512 coefficients, 2 yields 256 coefficients...).}
  \item{dis}{a character giving in the distribution in the \code{gbm.fit} function.}
  \item{nt}{an integer giving the total number of trees to fit in the \code{gbm.fit} function.}
  \item{s}{a character giving in the s parameter in the \code{gbm.fit} function.}
  \item{k}{a character giving in the k used in the \code{svm} function.}
  \item{...}{additional arguments.}
  \item{object}{an object of class \code{"regr"}.}
  \item{new}{a numerical data frame or matrix containing the new spectra.}
  \item{model}{an object of class \code{"regr"}.}
  \item{output.name}{a character naming the prediction output csv-file in \code{pred.regr}.}
}
\details{
Missing values in \code{y} are allowed.

\code{regr} uses the \code{mvr} function in the \code{pls} package for partial least-squares regression, the \code{gbm.fit} function in the \code{gbm} package for boosted regression trees and the \code{svm} function in the \code{e1071} package for support vector machines regression. The number of important PLS latent variables and the svm parameter optimization is done automatically based on experience with soil spectra.

\code{sp} uses for spectral transformation (i) the \code{locpoly} function in \code{KernSmooth} package for derivative calculation, (ii) the \code{chull} and \code{approx} functions in \code{"KernSmooth"} package for continuum removal and (iii) the \code{dwt} function in \code{wavelets} package for extraction of wavelet coefficients. Experiences showed for wavelet decomposition that the best ratio of prediction performance and sparse spectral representation is reached when all 128 wavelet coefficients from decomposition lv three are taken (which is the default).

Settings in the used functions for regression and transformation are chosen based on experience with soil spectra calibrations. It is recommended to take the given default values. Nevertheless, the settings can be adapted to a certain degree. In case you want to use complete functionality use the named functions directly. If \code{r} is \code{"brt"}, the number of samples has to be more than 70.

Column names of \code{x} and \code{new} must contain the wavebands. Wavebands are made automatically compatible if needed (see details in \code{\link{read.spc}).}.

Constituent values are not always normally distributed. This can violate prerequisitives for regression methods. Thus, transformation prior regression can solve this problem. The \code{regr} function uses log, square root and box-cox transformation aside untransformed values and let the user decide graphically which transformation to take for each constituent.

Predictions from \code{pred.regr} are given back with the prediction uncertainty for each individual sample (based on the validation set prediction error). The prediction uncertainty is calculated as the root median square error of prediction (RMedianSEP) using a moving window of in maximum 50 samples with similar predicted values. From the RMedianSEP the confidence interval is calculated.

Predictions are only made if (i) the new spectrum lies within the mahalanobis space of the calibration set, (ii) there is a local neighbor within of 5 and (iii) the predicted value lies within the calibration set range. Otherwise they are set to \code{NA} values. Mahalanobis distance can only be calculated when the number of calibration samples is higher than the number of wavebands/variables.

Calibration statistics contains for each constituents (i) \code{n} the number of samples used in calibration, (ii) \code{r2} the coefficient of determination for the linear regression of measured against predicted values, (iii) \code{a} the slope of the regression line, (iv) \code{bias} the bias, (v) \code{RMSEC} the root means square error of calibration, (vi) \code{RPD} the ratio of constituent standard deviation to RMSEC, (vii) \code{n LV} the number of latent variables used when \code{r} is equal to \code{"pls"}, (viii) \code{n bc out} the number of backtransformed values being \code{NA} values after box-cox transformation and (ix) \code{n trees} the number of trees when \code{r} is equal to \code{"brt"}. Validation statistics contains for each constituents points (i) to (vi). The RMSEC is logically the RMSEP.

The calibration and validation regressions of all constituents are plotted and the statistics printed in the Console.

Nearly each run of \code{regr} yields following warning message: \dQuote{1: in optimize(f = function(lambda)...}. Its related with the box-cox transformation, but does not have any impact or negative side effects.
}
\value{
\code{regr} returns a list with class \code{"regr"} containing the following components excluding the last four ones. \code{pred.regr} returns a list with class \code{"pred.regr"} containing the last four components \code{output.name},  \code{predicted.values},  \code{method}, and  \code{spectral.transformation} (see below):

  \item{model}{a character naming how the model output was named.}
  \item{model}{a list containing the regression output of class \code{"mvr"}, \code{"gbm"} or \code{"svm"}.}
  \item{x.tr}{a matrix containing the transformed spectra.}
  \item{spectral.transformation}{a character naming the spectral transformation.}
  \item{constituents}{a character naming the constituents.}
  \item{constituents.transformation}{a character naming the constituent transformations. Needed for \code{pred.regr}.}
  \item{lambda}{a numeric giving the lambda values in case the box-cox transformation was chosen as constituents transformation. Needed for \code{pred.regr}.}
  \item{method}{a character naming the used regression method.}
  \item{cal.samples}{a list containing the row names of the calibration samples for each soil constituent.}
  \item{val.samples}{a list containing the row names of the validation samples for each soil constituent.}
  \item{cal.statistics}{a matrix containing the calibration statistics for all constituents. See details.}
  \item{cal.mea.pre}{a data frame containing the calibration set measured and predicted values for all constituents.}
  \item{val.statistics}{a matrix containing the validation statistics for all constituents. See details.}
  \item{val.mea.pre}{a data frame containing the validation set measured and predicted values for all constituents.}
  \item{cal.pca}{a list containing objects of the class \code{"prcomp"} for each constituent calibration set. Needed for \code{pred.regr}.}
  \item{mahalanobis}{a list containing numeric vectors having the spectral mahalanobis distance of the constituents calibration sets. Needed for \code{pred.regr}.}
  \item{cal.range}{a list containing numeric vectors having the ranges of the constituents calibration sets. Needed for \code{pred.regr}.}
  \item{rmsep}{a list containing numeric vectors having for each constituent the root median square error of prediction for each validation set sample. See details for further explanation.  Needed for \code{pred.regr}.}
  \item{lm}{a list containing numeric vectors having for each constituent validation set the fitted values calculated by linear regression of measured against predicted values. Needed for \code{pred.regr}.}
  \item{wavebands}{a numeric vector containing the wavebands of \code{x}. Needed for \code{pred.regr}.}
  \item{drv}{an integer giving the order of derivative. Needed for \code{pred.regr}.}
  \item{bandwidth}{an integer defining the smoothing interval in wavebands. Needed for \code{pred.regr}.}
  \item{ft}{a character defining the wavelet ft. Needed for \code{pred.regr}.}
  \item{lv}{a character defining the lv of wavelet coefficients extraction. Needed for \code{pred.regr}.}
  \item{output.name}{a character string giving the name of the saved csv-file from \code{pred.regr}.}
  \item{predicted.values}{a matrix containing the predicted values and its respective confidence interval limits.}
  \item{method}{a character naming the used regression method.}
  \item{spectral.transformation}{a character naming the used spectral transformation method.}
}
\section{Note}{Please note that the usage section was removed due too mant model parameters which exceeds the limit of characters in a single line. The explaination for the concept adopted in the model is included here to provide readers with information about capability of this function. We are committed in revising this usage section to conform with the \code{CRAN} system. The code for this has been tested and works and those who are intrested can contact the author for the code file which can be called locally the same way local functions are loaded in R system }
\author{Thomas Terhoeven-Urselmans}