\name{LFQRatio2}
\alias{LFQRatio2}
\docType{data}
\title{Dataset "LFQRatio2"}
\description{This dataset is the final outcome of a Nano-LC-MS/MS proteomic analysis consisting to mix a known concentration of UPS1 human proteins within a yeast background. It contains the abundance values of identified proteins in two conditions. The condition A represents the measured abundances of proteins when 25fmol of UPS1 human proteins has been injected in yeast extract while the condition B represents the measured abundances of proteins in the case of 12.5fmol of injected UPS1 human proteins. Three biological replicates have been performed for each condition.

To identify proteins, spectra were searched using MaxQuant (version 1.5.1.2) against the Uniprot database, the UPS database and the frequently observed contaminants database. Maximum false discovery rates were set to 0.01 at peptide and protein levels by employing a reverse database strategy.

The abundance values of the dataset were obtained from LFQ values calculated using MaxQuant from MS intensity of unique peptides (see Cox et al. (2014)). The following pre-processing steps were performed to obtain these values using the Perseus toolbox (version 1.5.1.6) available in the MaxQuant environment: log2-transformation of LFQ values, filtering of proteins with at least 3 measured values in one condition and data imputation by replacing missing values with values generated from a normal distribution (see Deeb et al. (2012)).  

From a statistical viewpoint, the goal is to find which proteins are differentially abundant between the two conditions among the 1481 identified proteins. Ideally, the 38 proteins identified as human should be concluded as differentially abundant (in such a case the proportion of non-differentially abundant proteins will be pi0=1-38/1481).
}
\usage{data(LFQRatio2)}
\format{A data frame with the 1481 identified proteins in row. 

Columns \code{A.R1}, \code{A.R2} and \code{A.R3} correspond to the (numeric) abundance values of proteins in the three biological replicates of condition A. 

Columns \code{B.R1}, \code{B.R2} and \code{B.R3} correspond to the (numeric) abundance values of proteins in the three biological replicates of condition B. 

Column \code{Welch.test.pval} contains the p-values of the Welch t-test between condition A and condition B computed with the Perseus software.

Column \code{Organism} contains categorical values: \code{human} if the protein is identified as human and \code{yeast} otherwise. 

Column \code{Majority.protein.IDs} contains the IDs of proteins.
}
\references{
Cox J., Hein M.Y., Luber C.A., Paron I., Nagaraj N., Mann M. Accurate proteome-wide label-free quantification by delayed normalization and maximal peptide ratio extraction, termed MaxLFQ. Mol Cell Proteomics. 2014 Sep, 13(9):2513-26.

Deeb S.J., D'Souza R.C., Cox J., Schmidt-Supprian M., Mann M. Super-SILAC allows classification of diffuse large B-cell lymphoma subtypes by their protein expression profiles. Mol Cell Proteomics. 2012 May, 11(5):77-89.
}
\examples{data(LFQRatio2)

#p-values of the Welch t-test between condition A and condition B
p=LFQRatio2[,7]

}

\keyword{datasets}

