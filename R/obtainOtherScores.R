obtainOtherScores <- function(classifierResults, metaDataRef, classColumn, higherClassColumn, nModels = 100, crossValidation = F) {
  predictions <- classifierResults$classifications

  colnames(metaDataRef)[which(colnames(metaDataRef) == classColumn)] <- "class"
  colnames(metaDataRef)[which(colnames(metaDataRef) == higherClassColumn)] <- "higherClass"

  if (crossValidation == T) {
    probabilityList <- classifierResults$probabilityList
    for (i in seq(1:length(probabilityList))) {
      if (i == 1) {
        probabilities <- lapply(probabilityList[[i]], function(x) x / nModels)
      } else {
        probabilities <- c(probabilities, lapply(probabilityList[[i]], function(x) x / nModels))
      }
    }
  } else {
    probabilities <- lapply(classifierResults$probability, function(x) x / nModels)
  }


  combisClassAndHigherClass  <- metaDataRef %>%
    group_by(higherClass, class) %>%
    nest %>% arrange(higherClass) %>%
    select(higherClass, class) %>%
    ungroup


  higherClassList <- list()
  for (i in seq(1:length(probabilities))) {
    allProbs <- sort(probabilities[[i]], decreasing = T)
    combineScores <- NA
    if (length(allProbs) > 3) {
      allProbsLeftover <- allProbs[4:length(allProbs)]
      for (j in seq(1:length(allProbsLeftover))) {
        # Match the higherClass name to the subspecification
        higherClass <- combisClassAndHigherClass %>%
          filter(class == names(allProbsLeftover[j])) %>%
          select(higherClass)
        higherClassScore <- as.numeric(allProbsLeftover[j])
        if (j == 1) {
          higherClassDF <- data.frame(higherClass = higherClass$higherClass,
                                      Score = higherClassScore)
        } else {
          determinedHigherClass <- data.frame(higherClass = higherClass$higherClass,
                                              Score = higherClassScore)
          higherClassDF <-rbind(higherClassDF, determinedHigherClass)
          #if (determinedSubclass) == higherClass
        }
      }
      combineScores <- higherClassDF %>%
        group_by(higherClass) %>%
        nest(.) %>%
        mutate(Score = map(data,  ~ sum(.x$Score))) %>%
        select(higherClass, Score) %>%
        unnest(cols = c(Score))

    }
    higherClassList[[i]] <- combineScores
  }

  return(higherClassList)
}
