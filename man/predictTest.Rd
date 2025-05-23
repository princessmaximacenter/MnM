% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predictTest.R
\name{predictTest}
\alias{predictTest}
\title{Classify test samples}
\usage{
predictTest(modelList, testData)
}
\arguments{
\item{modelList}{List containing the models generated by the Minority classifier.}

\item{testData}{Test data to be classified. The dataset contains only the variables needed as input for the classifier.}
}
\value{
Dataframe containing the classifications of the tumor subtypes for each sample within the Minority classifier.
Each row is a test sample, with the separate columns containing the classifications from the separate models.
}
\description{
Function to classify new samples using the models generated with the reference cohort.
}
