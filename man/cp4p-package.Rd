\name{cp4p-package}
\alias{cp4p-package}
\alias{cp4p}
\docType{package}
\title{Introduction to the CP4P package}
\description{
This package provides tools to check whether a vector of p-values respects the assumptions of classical FDR (false discovery rate) control procedures. 

It is built to be easily used by non-statisticians in the context of quantitative proteomics (yet, it can be applied in other contexts). 

Concretely, it allows estimating the proportion of true null hypotheses (i.e. proportion of non-differentially abundant proteins or peptides in a relative quantification experiment), 
as well as checking whether the p-values are adequately distributed for further FDR  control. 

In addition, the package allows performing an adequately chosen adaptive FDR control procedure to get adjusted p-values.  
}
\details{
\tabular{ll}{
Package: \tab cp4p\cr
Type: \tab Package\cr
Version: \tab 0.3.3\cr
License: \tab GPL-3\cr
Depends:\tab multtest, qvalue, limma, MESS, graphics, stats\cr
}
This package is composed of three functions that take as input a vector of p-values resulting from multiple two-sided hypothesis testing (such as multiple t-tests for equal means for instance). 

First, the function \code{\link{estim.pi0}} allows determining the proportion of true null hypotheses among the set of tests using eight different estimation methods proposed in literature. 

Second, the function \code{\link{calibration.plot}} proposes an intuitive plot of the p-values, so as to visually assess their behavior and well-calibration. 

Third, the function \code{\link{adjust.p}} allows obtaining adjusted p-values in view to perform an adaptive FDR control from a chosen level.

Two proteomic datasets named \code{\link{LFQRatio2}} and \code{\link{LFQRatio25}} allow to use these functions in a concrete framework where the proportion of non-differentially abundant proteins is known.
}
\author{
Quentin Giai Gianetto, Florence Combes, Claire Ramus, Christophe Bruley, Yohann Couté, Thomas Burger

Maintainer: Quentin Giai Gianetto <quentin2g@yahoo.fr>
}
\references{
Giai Gianetto, Q., Combes, F., Ramus, C., Bruley, C., Couté, Y., Burger, T. (2016). Calibration plot for proteomics: A graphical tool to visually check the assumptions underlying FDR control in quantitative experiments. Proteomics, 16(1), 29-32.}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the ~~
%~~ R documentation directory ~~
%\keyword{false discovery rate}
%\keyword{adjusted p-values}
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
%\examples{
%~~ simple examples of the most important functions ~~
%}
