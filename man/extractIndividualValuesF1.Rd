% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extractIndividualValuesF1.R
\name{extractIndividualValuesF1}
\alias{extractIndividualValuesF1}
\title{Calculate per tumor (sub)type performance}
\usage{
extractIndividualValuesF1(
  predictionsMM,
  metaDataRef,
  classColumn,
  probabilityThreshold,
  filterOrNot
)
}
\arguments{
\item{predictionsMM}{Predictions for samples by M&M, coming from the function integrateMM.}

\item{metaDataRef}{Metadata file containing the links between the samples and
the tumor (sub)type diagnosis within the reference cohort.}

\item{classColumn}{Name of column in the metadata file that contains the tumor (sub)type labels.}

\item{probabilityThreshold}{What is the probability score threshold you would like to use to call a classification 'confident'?}

\item{filterOrNot}{Do you want to filter the 'confident' classifications only for your calculation?}
}
\value{
Dataframe containing the tumor (sub)type ($tumorType) with its associated precision ($Precision), F1 score ($F1),
sensitivity ($Sensitivity) and recall ($Recall), stratified by population frequency ($nCases).

Note that recall here means the fraction of samples obtaining a 'confident' classification
(probability higher than probabilityThreshold).
}
\description{
Calculate per tumor (sub)type performance
}
