% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotMeanSDAccuracy.R
\name{plotMeanSDAccuracy}
\alias{plotMeanSDAccuracy}
\title{Plot the results per frequency range of tumor entities}
\usage{
plotMeanSDAccuracy(meanAndSDPlotTrain, meanAndSDPlotTest = NA)
}
\arguments{
\item{meanAndSDPlotTrain}{Dataframe resulting from function 'calculateMeanAndSDAccuracy', for reference cohort data.}

\item{meanAndSDPlotTest}{Dataframe resulting from function 'calculateMeanAndSDAccuracy', for test cohort data.}
}
\value{
Plot showing the accuracy for train and test set for the different frequencies.
}
\description{
Plot the results per frequency range of tumor entities
}
