addHigherClassScore <- function(probabilityDFWithHigherClass,
                                higherClassList) {
  probabilityDFWithHigherClass$higherClassScore <- 0
  probabilityDFWithHigherClass$higherClassScore2 <- 0
  probabilityDFWithHigherClass$higherClassScore3 <- 0
  probabilityDFWithHigherClass$higherClassScore4 <- 0

  for (i in seq(1:length(probabilityDFWithHigherClass$predict))) {

    otherModelSubclasses <- higherClassList[[i]]

    totalScore <- probabilityDFWithHigherClass$probability1[i]
    totalScore2 <- 0
    totalScore3 <- 0
    totalScore4 <- 0
    if(is.na(probabilityDFWithHigherClass$predictHigherClass2[i])) {
      totalScore <- totalScore
    } else {
      if(probabilityDFWithHigherClass$predictHigherClass[i] == probabilityDFWithHigherClass$predictHigherClass2[i]) {
        totalScore <- totalScore +
          probabilityDFWithHigherClass$probability2[i]
      } else {
        totalScore2 <- probabilityDFWithHigherClass$probability2[i]
      }
    }

    if(is.na(probabilityDFWithHigherClass$predictHigherClass3[i])) {totalScore <- totalScore} else {
      if (probabilityDFWithHigherClass$predictHigherClass[i] == probabilityDFWithHigherClass$predictHigherClass3[i]) {
        totalScore <- totalScore + probabilityDFWithHigherClass$probability3[i]
      } else if (probabilityDFWithHigherClass$predictHigherClass2[i] == probabilityDFWithHigherClass$predictHigherClass3[i]) {
        totalScore2 <- totalScore2 + probabilityDFWithHigherClass$probability3[i]
      } else {totalScore3 <- probabilityDFWithHigherClass$probability3[i]}

      if(is.na(otherModelSubclasses)[1]) {totalScore <- totalScore} else {

        for (j in seq(1:nrow(otherModelSubclasses))) {
          higherClass <- otherModelSubclasses$higherClass[j]
          subclassScore <- otherModelSubclasses$Score[j]
          if (higherClass ==  probabilityDFWithHigherClass$predictHigherClass[i]) {
            totalScore <- totalScore +
              subclassScore

          } else if (higherClass == probabilityDFWithHigherClass$predictHigherClass2[i] &
                     probabilityDFWithHigherClass$predictHigherClass[i] != probabilityDFWithHigherClass$predictHigherClass2[i]) {
            totalScore2 <- totalScore2 + subclassScore
          } else if (higherClass == probabilityDFWithHigherClass$predictHigherClass3[i] &
                     probabilityDFWithHigherClass$predictHigherClass[i] != probabilityDFWithHigherClass$predictHigherClass3[i] &
                     probabilityDFWithHigherClass$predictHigherClass2[i] != probabilityDFWithHigherClass$predictHigherClass3[i]) {
            totalScore3 <- totalScore3 + subclassScore
          } else {totalScore4 <- totalScore4 + subclassScore}
        }
      }
    }

    probabilityDFWithHigherClass$higherClassScore[i] <- totalScore
    probabilityDFWithHigherClass$higherClassScore2[i] <- totalScore2
    probabilityDFWithHigherClass$higherClassScore3[i] <- totalScore3
    probabilityDFWithHigherClass$higherClassScore4[i] <- totalScore4
  }

  return(probabilityDFWithHigherClass)
}
