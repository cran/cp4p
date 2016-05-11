\name{estim.pi0}
\alias{estim.pi0}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Estimation of the proportion of true null hypotheses
}

\description{
From a proteomics viewpoint, this function estimates the global proportion of proteins (resp. of peptides) that are non differentially abundant from the tested protein list (resp. from the tested peptide list).
This proportion is later used as a correcting factor to compute the adjusted p-values, that are in turn used to tune a threshold according to a desired false discovery rate.

From a statistical viewpoint, this function allows estimating the proportion of true null hypotheses (pi0) from a vector of raw p-values following eight different estimation methods from the literature.
}

\usage{
estim.pi0(p, pi0.method = "ALL", nbins = 20, pz = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{Numeric vector of raw p-values. Raw p-values are assumed without missing values, and between 0 and 1.}
  \item{pi0.method}{Name of an estimation method for the proportion of true null hypotheses among \code{"st.boot"}, \code{"st.spline"}, \code{"langaas"}, \code{"jiang"}, \code{"histo"}, \code{"pounds"}, \code{"abh"} or \code{"slim"}. Default is \code{"ALL"}: all the eight estimation methods are performed simultaneously.}
  \item{nbins}{Number of bins. Parameter used for the \code{"jiang"} and  \code{"histo"} methods. Default is 20.}
  \item{pz}{P-value threshold such as p-values below are associated to false null hypotheses. Used for the \code{"slim"} method. Wang, Tuominen and Tsai (2011) suggest to take a value between 0.01 and 0.1. Default is 0.05.}
}

\details{
This function allows to estimate the proportion of true null hypotheses following different estimation methods :
\tabular{ll}{
\code{"abh"}\tab the least slope method proposed in Benjamini and Hochberg (2000).\cr
 \tab \cr
\code{"st.spline"}\tab the smoother method described in Storey and Tibshirani (2003). \cr 
 \tab The \code{qvalue} function of R package \code{qvalue} with default tuning is used (Storey (2015)).\cr
  \tab \cr
\code{"st.boot"}\tab the bootstrap method described in Storey et al. (2004).\cr
 \tab The \code{qvalue} function of R package \code{qvalue} with default tuning is used (Storey (2015)).\cr
 \tab \cr
\code{"langaas"}\tab the method described in Langaas, Ferkingstad and Lindqvist (2005) using a convex\cr
  \tab decreasing density estimate for p-values. The \code{convest} function of R package \code{limma}\cr
 \tab  with default tuning is used (Ritchie et al. (2015)).\cr
 \tab \cr
\code{"histo"}\tab the histogram method described in Nettleton, Hwang, Caldo and Wise (2006).\cr
 \tab \cr
\code{"pounds"}\tab the conservative estimate described in Pounds and Cheng (2006).\cr
 \tab \cr
\code{"jiang"}\tab the average estimate method described in Jiang and Doerge (2008).\cr
 \tab \cr
\code{"slim"}\tab the method of Wang, Tuominen and Tsai (2011) using a sliding linear model.\cr
 \tab The default tuning suggested by Wang, Tuominen and Tsai (2011) is used.\cr
 \tab Using their notations, lambda1 is fixed to 0.1, n to 10 and B to 100.
}
To take into account of right censorship on the vector of p-values, each p-value is divided by the maximum p-value present in \code{p}. Accordingly, the p-values of the true null hypotheses are assumed uniformly distributed between 0 and this maximum. This kind of censorship happens in proteomics when a first thresholding is performed on the fold-changes. 

If you want to assume that the p-values are uniformly distributed between 0 and 1, replace \code{p} by \code{c(p,1)} when using \code{estim.pi0}. 
}

\value{
\item{pi0}{Numeric value of the estimated proportion of true null hypotheses from the selected method; Numeric vector if \code{pi0.method="ALL"}.}
}
\references{
Y. Benjamini and Y. Hochberg. On the adaptive control of the false discovery rate in multiple testing with independent statistics. Journal of Educational and Behavioral Statistics, 25(1):60-83, 2000.

H. Jiang and R.W. Doerge. Estimating the proportion of true null hypotheses for multiple comparisons. Cancer informatics, 6:25, 2008.

M. Langaas, B.H. Lindqvist, and E. Ferkingstad. Estimating the proportion of true null hypotheses, with application to dna microarray data. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 67(4):555-572, 2005.

D. Nettleton, J.T.G. Hwang, R.A. Caldo, and R.P. Wise. Estimating the number of true null hypotheses from a histogram of p values. Journal of Agricultural, Biological, and Environmental Statistics, 11(3):337-356, 2006.

S. Pounds and C. Cheng. Robust estimation of the false discovery rate. Bioinformatics, 22(16):1979-1987, 2006.

M.E. Ritchie, B. Phipson, D. Wu, Y. Hu, C.W. Law, W. Shi and G.K. Smyth. “limma powers differential expression analyses for RNA-sequencing and microarray studies.” Nucleic Acids Research, 43(7), pp.e47. 2015.

J.D. Storey, J.E. Taylor, and D. Siegmund. Strong control, conservative point estimation and simultaneous conservative consistency of false discovery rates: a unified approach. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 66(1):187-205, 2004.

J.D. Storey and R. Tibshirani. Statistical significance for genomewide studies. Proceedings of the National Academy of Sciences, 100(16):9440-9445, 2003.

J.D. Storey. qvalue: Q-value estimation for false discovery rate control. R package version 2.0.0, http://qvalue.princeton.edu/, http://github.com/jdstorey/qvalue. 2015.

H.-Q. Wang, L.K. Tuominen, and C.-J. Tsai. SLIM: a sliding linear model for estimating the proportion of true null hypotheses in datasets with dependence structures. Bioinformatics, 27(2):225-231, 2011.
}
\author{
Quentin Giai Gianetto <quentin2g@yahoo.fr>
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{calibration.plot}}, \code{\link{adjust.p}}
}
\examples{
#get p-values
data(LFQRatio2)
p=LFQRatio2[,7]

#estimate the proportion of true null hypotheses with different methods
r=estim.pi0(p)
r$pi0

#estimate the proportion of true null hypotheses with the "abh" method
r=estim.pi0(p, pi0.method="abh")
r$pi0

#compare with one minus the proportion of human proteins 
prop_human=sum(LFQRatio2$Organism=="human")/length(LFQRatio2$Organism)
pi0_true=1-prop_human
pi0_true

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
