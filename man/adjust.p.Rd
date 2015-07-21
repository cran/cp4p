\name{adjust.p}
\alias{adjust.p}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Adjusted p-values for adaptive FDR control}
\description{
This function computes adjusted p-values for adaptive FDR control from a vector of raw (unadjusted) p-values. 
}
\usage{
adjust.p(p, pi0.method = 1, alpha = 0.05, nbins = 20, pz = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{Numeric vector of raw p-values.}
  \item{pi0.method}{Numeric value between 0 and 1 corresponding to the proportion of true null hypotheses (non-differentially abundant proteins or peptides), or the name of an estimation method for this proportion among \code{"st.boot"}, \code{"st.spline"}, \code{"langaas"}, \code{"jiang"}, \code{"histo"}, \code{"pounds"}, \code{"abh"} or \code{"slim"} (see function \code{\link{estim.pi0}} for details). The two-stage Benjamini and Hochberg procedure (Benjamini et al. (2006)) is also available according to an expected FDR given by the \code{alpha} parameter (write \code{pi0.method="bky"}). Default is 1 (classical Benjamini and Hochberg procedure (1995) is performed in this case).}
  \item{alpha}{A nominal type I error rate used for estimating the proportion of true null hypotheses (non-differentially abundant proteins or peptides) in the two-stage Benjamini and Hochberg procedure (used only if \code{pi0.method="bky"}). Default is 0.05.}
  \item{nbins}{Number of bins. Parameter used for the \code{"jiang"} and  \code{"histo"} methods in \code{\link{estim.pi0}}. Default is 20.}
  \item{pz}{P-value threshold such as p-values below are associated to false null hypotheses. Used for the \code{"slim"} method in \code{\link{estim.pi0}}. Default is 0.05.}
  %\item{...}{Arguments to be passed to the chosen estimation method of the proportion of true null hypotheses (non-differentially abundant proteins or peptides). See \code{\link{estim.pi0}} for details.}
}
\details{
The procedure uses an estimation of the proportion of true null hypotheses (non-differentially abundant proteins or peptides), the value or the name of which is precised in input. Next, this estimation is multiplied by the adjusted p-values of the Benjamini and Hochberg procedure (1995) to obtain the final adjusted p-values (see section 3 in Craiu and Sun (2008) for details).

The adjusted p-values of the Benjamini and Hochberg procedure (1995) and of the two-stage Benjamini and Hochberg procedure (Benjamini et al. (2006)) are computed using the R package \code{multtest} (Pollard et al. (2005)).
}
\value{
A list composed of :
  \item{pi0}{The proportion of true null hypotheses (non-differentially abundant proteins or peptides) used to adjust p-values.}
  \item{adjp}{A matrix of raw and adjusted p-values with rows corresponding to each test. First column corresponds to raw p-values and second column to adjusted p-values.}
}
\references{
Y. Benjamini and Y. Hochberg. Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society. Series B (Methodological), 289-300, 1995.

Y. Benjamini, A.M. Krieger, and D.Yekutieli. Adaptive linear step-up procedures that control the false discovery rate. Biometrika, 93(3):491-507, 2006.

R.V. Craiu and L. Sun. Choosing the lesser evil: trade-off between false discovery rate and non-discovery rate. Statistica Sinica, 18:861-879, 2008.

K.S. Pollard, S. Dudoit and M.J. van der Laan. Multiple Testing Procedures: R multtest Package and Applications to Genomics, in Bioinformatics and Computational Biology Solutions Using R and Bioconductor. Springer. 2005.
}
\author{
Quentin Giai Gianetto <quentin2g@yahoo.fr>
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{estim.pi0}}, \code{\link{calibration.plot}}
}
\examples{
#get p-values
data(LFQRatio2)
p=LFQRatio2[,7]

#adjust p-values by estimating the proportion of true null hypotheses 
#with the "pounds" method.
res_pounds=adjust.p(p, pi0.method = "pounds")

#proportion of true null hypotheses with the "pounds" method.
res_pounds$pi0

#plot ajusted p-values in function of raw p-values
plot(res_pounds$adjp)

#adjust p-values by estimating the proportion of true null hypotheses
#using the two-stage Benjamini and Hochberg procedure with a FDR of 0.1.
res_bky=adjust.p(p, pi0.method = "bky", alpha = 0.1)

#proportion of true null hypotheses with the two-stage BH procedure.
res_bky$pi0

#plot adjusted p-values in function of raw p-values
plot(res_bky$adjp)

#compare the two-stage Benjamini and Hochberg procedure 
#with the "pounds" method
plot(res_pounds$adjp[,2],res_bky$adjp[,2])

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
