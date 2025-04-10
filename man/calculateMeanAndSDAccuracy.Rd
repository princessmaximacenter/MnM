% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculateMeanAndSDAccuracy.R
\name{calculateMeanAndSDAccuracy}
\alias{calculateMeanAndSDAccuracy}
\title{Calculate performance of M&M on total train and test set}
\usage{
calculateMeanAndSDAccuracy(
  minorityDir,
  majorityDir,
  metaDataTest = NA,
  subtype = F,
  rounding = F,
  probabilityThreshold
)
}
\arguments{
\item{minorityDir}{Directory in which the minority model(s) are stored.}

\item{majorityDir}{Directory in which the majority model(s) are stored.}

\item{metaDataTest}{Metadata file containing the links between the samples and the tumor type and subtype diagnoses within the test set.}

\item{subtype}{Do you want to obtain the predictions on the tumor subtype classification level (subtype = TRUE)?}

\item{rounding}{Do you want rounded numbers for the performance scores? Default = FALSE.}

\item{probabilityThreshold}{What is the probability score threshold you would like to use to call a classification 'confident'?}
}
\value{
Dataframe containing the average performance of tumor classifications within a certain frequency range (nCases).
Included are the averages for the percentage of correctly and incorrectly classified samples ($meanFractionCorrect and $meanFractionIncorrect),
correctly and incorrectly classified 'confident' samples ($meanFractionCorrectFiltered and $meanFractionIncorrectFiltered),
averaged precision values for the average of all tumor entities within the frequency range ($meanPrecision),
averaged recall values for the average of all tumor entities within the frequency range ($meanRecall),
and averaged F1 scores for the average of all tumor entities within the frequency range ($meanF1).
Furthermore, the standard deviations of all mentioned values are calculated in case multiple runs with different seeds have been performed ($sd...).
Please note that the precision, F1 and recall are calculated for the confident sample classifications only.
The total amount of samples within each frequency range ($meanSamples) is also specified.
Lastly, it's specified whether a cross-validation (Train) or test (Test) type was used.
}
\description{
Calculate performance of M&M on total train and test set
}
