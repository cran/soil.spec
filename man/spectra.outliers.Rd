\name{spectra.outliers}
\alias{spectra.outliers}
\title{A tool for screening  spectral outliers based on Multidimensional Scaling (MDS) approach}
\description{Reads a data table containing raw spectra which first is transformed using first derivative from Savitzky-Golay algorithm to take off any effects due to measurement differences.    Euclidean distances are computed from the first derivative transformed spectra, other types for metric measure can be used but the function uses Euclidean method only. Multidimensional scaling is then used to represent proximities or Euclidean similarity distances in a more organized structure. Simple visualization are produced on request to show how individual spectral point distribution within the Euclidean space. Mean and standard deviation  for the similarity distances are calculated for determining the most similar points and the most dissimilar ones. There are two types of situations covered here. When checking for outliers within an highly homogenous spectral set where, outliers will be marked for point lying one standard deviation away from the mean distance. This is the case when screening for outliers from a reference standard sample used to keep checking on instrument stability. The other situation is when checking for spectral outliers from a large spectra collection. In this case, the confidence interval for dertermining outliers is increased upto three standard deviations from the mean to allow for any normal sample to sample differences. Therefore, in the function specifying whether spectra comes from standard or non.standard is recommended for the correct confindence interval calculation.
} 
\usage{ spectra.outliers(file.name, spectra="non.standard", 
     plots=TRUE, smethod="euclidean", k=5)
}
\arguments{
  \item{file.name}{\code{"character"}; file name with extension}
  \item{spectra}{\code{"character"}; specify whether spectra is from "non.standard" for actual sample spectra or is a "standard" for reference samples.}
  \item{plots}{\code{"logical"}; specifies whether to produce outlier diagnostic plots or not. Default is set to TRUE}
  \item{smethod}{\code{"character"}; specifies metric method for computing similarity distance. Euclidean method is the one supported }
  \item{k}{\code{"integer"}; number of components to visualize the MDS distance}
}
\value{
\item{dist }{n by k matrix with similarity distances, n=rows of spectral table; k=the number of components used defualut is set to 5}
\item{all.outliers}{Row numbers with spectra marked as outliers}
\item{spec.out}{Data frame with outlier spectra}
\item{spec.normal}{Data frame with remaining spectra after removing outlier spectra}
}
\note{There are two options for getting list of outliers:
	\itemize{
	\item Take all the points marked as outliers falling outside the confidence bands explained above
	\item Get into interactive mode to identify outlying points by clicking all the points appearing as outliers. To stop selection press Ctrl+any keyboard button
}
But whichever the choice, all points outside the confidence band is given in different colour to aid easy identification. User to decide on whether to use the default outliers ir use the interactive mode}
\references{
\itemize{ 
\item Kruskal, J.B. and M. Wish (1978). Multidimensional Scaling. Sage.
\item Cox, T.F., Cox, M.A.A. (2001). Multidimensional Scaling. Chapman and Hall.
}
}
\author{ Andrew Sila and Tomislav Hengl}