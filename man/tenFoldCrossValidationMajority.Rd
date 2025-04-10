% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tenFoldCrossValidationMajority.R
\name{tenFoldCrossValidationMajority}
\alias{tenFoldCrossValidationMajority}
\title{10x cross-validation majority classifier}
\usage{
tenFoldCrossValidationMajority(
  countDataRef,
  metaDataRef,
  classColumn,
  higherClassColumn,
  domainColumn,
  sampleColumn,
  nModels = 100,
  maxSamplesPerType = 50,
  nFeatures = 2500,
  nComps = 100,
  maxNeighbors = 25,
  whichSeed = 1,
  outputDir = paste0("./", format(as.Date(Sys.Date(), "\%Y-\%m-\%d"), "\%Y_\%m_\%d")),
  proteinCodingGenes
)
}
\arguments{
\item{countDataRef}{Matrix containing the RNA-transcript per million data. Samples are in the columns,
different RNA-transcripts in the rows.}

\item{metaDataRef}{Metadata file containing the links between the samples and the tumor domain, type and subtype diagnosis.}

\item{classColumn}{Name of column in the metadata file that contains the tumor subtype labels.}

\item{higherClassColumn}{Name of column in the metadata file that contains the tumor type labels.}

\item{domainColumn}{Name of column in the metadata file that contains the tumor domain labels.}

\item{sampleColumn}{Name of column in the metadata file that contains the sample identifiers.}

\item{nModels}{How many models should be created for the Majority classifier?}

\item{maxSamplesPerType}{How many samples should we maximally use per tumor (sub)type?}

\item{nFeatures}{How many of the most variable genes within the dataset should we select for principal component analysis (PCA)?}

\item{nComps}{How many principal components will be selected after PCA?}

\item{maxNeighbors}{What is the maximum number of neighbors to be used for the weighted \emph{k}-nearest neighbor algorithm?}

\item{whichSeed}{For reproducibility, the seed can be specified with this parameter.}

\item{outputDir}{Directory in which you would like to store the R-object containing the results. Default is today's date.}

\item{proteinCodingGenes}{What are the names of the RNA-transcripts that stand for protein-coding genes within our dataset?
Please supply it as a vector. This is needed for ribo-depletion correction model.}
}
\value{
R-object containing the predictions ($classifications), classifications errors ($wrongClassifications),
the probabilities for each classification ($probabilityList),
the metadata file associated to the reference cohort ($metaDataRef),
and metadata for the performed run ($metaDataRun).
}
\description{
Setup to automatically run the 10x stratified cross-validation majority classifier.
}
