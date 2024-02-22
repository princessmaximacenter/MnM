plotPerformanceDomain <- function(precisionDomainTrain,
         precisionDomainTest,
         levelsDomain

) {

  precisionDomainTrain$trainOrTest <- "Train"
  precisionDomainTest$trainOrTest <- "Test"

  precisionDomainsTotal <- rbind(precisionDomainTrain,
                                        precisionDomainTest)
  precisionDomainsTotal$Domain <- factor(precisionDomainsTotal$Domain, levels = levelsDomain)

  precisionDomainsTotal$trainOrTest <- factor(precisionDomainsTotal$trainOrTest, levels = c("Train", "Test"))

  precisionDomainsTotalLonger <- precisionDomainsTotal %>% pivot_longer(cols = c(meanPrecision, meanAccuracy, meanRecall),
                                                                                      names_to = "measurementType")


  precisionDomainsTotalLonger$valuePercent <- paste0(round(precisionDomainsTotalLonger$value *100, digits = 1), "%")


  precisionDomainsTotalLonger <- precisionDomainsTotalLonger %>% pivot_longer(cols = c("sdAccuracy", "sdPrecision", "sdRecall"),
                                                                                            names_to = "whichSD",
                                                                                            values_to = "standardDeviation")

  precisionDomainsTotalLonger$measurementType <- gsub("mean", "", precisionDomainsTotalLonger$measurementType)

  precisionDomainsTotalLonger$whichSD <- gsub("sd", "", precisionDomainsTotalLonger$whichSD)


  precisionDomainsTotalLongerFiltered <- precisionDomainsTotalLonger %>% filter(whichSD == measurementType)


  return(precisionDomainsTotalLongerFiltered)

}
