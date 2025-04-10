% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getNonAvailableTiles.R
\name{getNonAvailableTiles}
\alias{getNonAvailableTiles}
\title{Obtain empty tiles for confusion matrix}
\usage{
getNonAvailableTiles(
  predictionsMM,
  probabilityThreshold,
  abbreviations,
  classColumn
)
}
\arguments{
\item{predictionsMM}{Predictions for the tumor (sub)types with their associated probability scores.}

\item{probabilityThreshold}{Probability score that you would like to use as a threshold for calling samples
'confident' or not.}

\item{abbreviations}{Dataframe containing the links between the tumor subtype and their abbreviations ($abbreviationSubtype),
the tumor types and their abbreviations ($abbreviationTumorType), and the domain.}

\item{classColumn}{Specifies the name of the column of the abbreviations dataframe that contains the full tumor type label.}
}
\value{
A dataframe containing the empty tiles for a confusion matrix plot.
}
\description{
Get tiles from the confusion matrix that are never used for a classification.
These tiles are important for the generation of the eventual confusion matrix.
}
