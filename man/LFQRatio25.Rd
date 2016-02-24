\name{LFQRatio25}
\alias{LFQRatio25}
\docType{data}
\title{Dataset "LFQRatio25"}
\description{This dataset is the final outcome of a quantitative mass spectrometry-based proteomic analysis of two samples containing different concentrations of 48 human proteins (UPS1 standard from Sigma-Aldrich) within a constant yeast background (see Giai Gianetto et al. (2016) for details). It contains the abundance values of the different human and yeast proteins identified and quantified in these two conditions. The conditions A and B represent the measured abundances of proteins when respectively 25fmol and 10fmol of UPS1 human proteins were mixed with the yeast extract before mass spectrometry analyses. Three technical replicates were acquired for each condition.

To identify and quantify proteins, spectra were searched using MaxQuant (version 1.5.1.2) against the Uniprot database, the UPS database and the frequently observed contaminants database. Maximum false discovery rates were set to 0.01 at peptide and protein levels by employing a reverse database strategy.

The abundance values of the dataset were obtained from LFQ values calculated using MaxQuant from MS intensity of unique peptides (see Cox et al. (2014)). The following pre-processing steps were performed to obtain these values using the Perseus toolbox (version 1.5.1.6) available in the MaxQuant environment: log2-transformation of LFQ values, filtering of proteins with at least 3 measured values in one condition and data imputation by replacing missing values with values generated from a normal distribution (see Deeb et al. (2012)).

From a statistical viewpoint, the goal is to find which proteins are differentially abundant between the two conditions among the 1472 quantified proteins. Ideally, the 46 quantified human proteins (out of the original 48 ones) should be concluded as differentially abundant (in such a case the proportion of non-differentially abundant proteins will be pi0=1-46/1472).
}
\usage{data(LFQRatio25)}
\format{A data frame with the 1472 identified proteins in row. 

Columns \code{A.R1}, \code{A.R2} and \code{A.R3} correspond to the (numeric) abundance values of proteins in the three replicates of condition A. 

Columns \code{B.R1}, \code{B.R2} and \code{B.R3} correspond to the (numeric) abundance values of proteins in the three replicates of condition B. 

Column \code{Welch.test.pval} contains the p-values of the Welch t-test between condition A and condition B computed with the Perseus software.

Column \code{Organism} contains categorical values: \code{human} if the protein is identified as human and \code{yeast} otherwise. 

Column \code{Majority.protein.IDs} contains the IDs of proteins.
}
\references{
Cox J., Hein M.Y., Luber C.A., Paron I., Nagaraj N., Mann M. Accurate proteome-wide label-free quantification by delayed normalization and maximal peptide ratio extraction, termed MaxLFQ. Mol Cell Proteomics. 2014 Sep, 13(9):2513-26.

Deeb S.J., D'Souza R.C., Cox J., Schmidt-Supprian M., Mann M. Super-SILAC allows classification of diffuse large B-cell lymphoma subtypes by their protein expression profiles. Mol Cell Proteomics. 2012 May, 11(5):77-89.

Giai Gianetto, Q., Combes, F., Ramus, C., Bruley, C., Couté, Y., Burger, T. (2016). Calibration plot for proteomics: A graphical tool to visually check the assumptions underlying FDR control in quantitative experiments. Proteomics, 16(1), 29-32.
}
\examples{data(LFQRatio25)

#p-values of the Welch t-test between condition A and condition B
p=LFQRatio25[,7]

}
\keyword{datasets}

