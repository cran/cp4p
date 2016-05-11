\name{calibration.plot}
\alias{calibration.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Displaying the "Calibration Plot" of a vector of p-values.
}
\description{
From a proteomics viewpoint, this function displays a graph (the "Calibration Plot") which allows to visually assess the compliance of a differential abundance analysis with FDR control procedure assumptions.

From a statistical viewpoint, this function performs a plot of the cumulative distribution function of 1-p-values. It allows checking whether p-values respect several assumptions of FDR control procedures.
}
\usage{
calibration.plot(p, pi0.method = "pounds", nbins = 20, pz = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{Numeric vector of raw p-values. Raw p-values are assumed without missing values, and between 0 and 1.}
  \item{pi0.method}{
  Numeric value between 0 and 1 corresponding to the proportion of true null hypotheses if known by the user, or the name of an estimation method among \code{"st.boot"}, \code{"st.spline"}, \code{"langaas"}, \code{"jiang"}, \code{"histo"}, \code{"pounds"}, \code{"abh"} or \code{"slim"} (see function \code{\link{estim.pi0}} for details). Default is \code{"pounds"}. If \code{pi0.method="ALL"}, a plot allowing the comparison of the eight estimation methods is displayed.
	}
  \item{nbins}{Number of bins. Parameter used for the \code{"jiang"} and  \code{"histo"} methods in \code{\link{estim.pi0}}. Default is 20.}
  \item{pz}{P-value threshold such as p-values below are associated to false null hypotheses. Used for the \code{"slim"} method in \code{\link{estim.pi0}}. Default is 0.05.}
  %\item{...}{Other arguments to be passed, depending on the chosen estimation method (see \code{\link{estim.pi0}} for details).}
}
\details{
This function provides a graph which displays the cumulative distribution function of 1-p-values as a function of 1-p-values (black curve) as advocated by Schweder and Spjotvoll (1982). 

The blue straight line has a slope equals to the proportion of true null hypotheses (estimated by \code{\link{estim.pi0}}) that is recalled in the caption of the plot. It is close to the black curve for small 1-pvalues if the p-values are independently and uniformly distributed under the null hypothesis.

In addition, two other measures are given in the caption of the graphic. Each has a color that matches that of various areas of the plot and should be carefully consider to assess the well-calibration of p-values (see Giai Gianetto et al. (2016) for details). 

The first measure corresponds to one minus the ratio between the green area and the grey area (referred to as "differentially abundant protein concentration"). The closer to 100\% this measure is, the smaller the false nondiscovery rate is expected. 

The second measure corresponds to the total red area observed on the graph (referred to as "uniformity underestimation"). The smaller this measure is, the more the proportion of true null hypotheses is expected to be not under-estimated (so as to get a conservative p-value adjustment).

Supplementary theoretical justifications on these measures can be found in the tutorial available in the supplementary material of Giai Gianetto et al. (2016).
}
\value{
A list composed of :
\item{pi0}{Numeric value corresponding to the proportion of true null hypotheses (non-differentially abundant proteins or peptides) used for the plot. Numeric vector if \code{pi0.method="ALL"}.}
\item{h1.concentration}{Numeric value corresponding to one minus the ratio between the green area and the grey area. NULL if \code{pi0.method="ALL"}.}
\item{unif.under}{Numeric value corresponding to the total red area observed on the graph (multiplied by 100). NULL if \code{pi0.method="ALL"}.}
}
\references{
Giai Gianetto, Q., Combes, F., Ramus, C., Bruley, C., Couté, Y., Burger, T. (2016). Calibration plot for proteomics: A graphical tool to visually check the assumptions underlying FDR control in quantitative experiments. Proteomics, 16(1), 29-32.

Schweder, T., Spjotvoll, E. (1982). Plots of p-values to evaluate many tests simultaneously. Biometrika, 69(3), 493-502.
}
\author{
Quentin Giai Gianetto <quentin2g@yahoo.fr>
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} 
\seealso{
\code{\link{estim.pi0}}
}
\examples{
#get p-values
data(LFQRatio25)
p=LFQRatio25[,7]

#Plot straight lines whose slopes correspond to different estimates of 
#the proportion of true null hypotheses
r=calibration.plot(p, pi0.method="ALL")
r$pi0
 
#Plot of the graph with the "pounds" method (default)
r=calibration.plot(p)
#Estimate of the proportion of true null hypotheses
r$pi0
#Estimate of the differentially abundant protein concentration 
#(the closer to one, the better)
r$h1.concentration
#Estimate of the "uniformity underestimation" quantity 
#(If null, pi0 is not underestimated.)
r$unif.under

#Plot of the graph using the "slim" method
r=calibration.plot(p, pi0.method="slim")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
