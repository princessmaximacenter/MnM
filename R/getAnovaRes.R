getAnovaRes <- function(allMaterial, allGenes, classColumn) {
  anovaValues <- data.frame(allGenes = allGenes, p_val = NA, F_val = NA)
  for (i in seq(1:length(allGenes))) {
    ANOVARes <- aov(allMaterial[,allGenes[i]] ~ allMaterial[, classColumn])
    anovaValues[i,"F_val"] <- summary(ANOVARes)[[1]][["F value"]][1]
    anovaValues[i,"p_val"] <- summary(ANOVARes)[[1]][["Pr(>F)"]][1]
  }
  return(anovaValues)
}
