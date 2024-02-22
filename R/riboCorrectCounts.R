#' Correct ribodepletion protocol
#'
#' This functions is used to account for the differences in efficacy of the ribodepletion protocol.
#' It is estimated what is the percentage of protein coding transcripts within the total transcriptome, and subsequently
#' the protein coding fraction is scaled up to 100%, simultaneously removing certain residual ribosomal RNAs from the dataset.
#' @param data Dataframe containing the ribo-depleted count data from an RNA-seq experiment.
#' @param proteinCodingGenes Names of the protein-coding genes within our RNA-seq count data.
#' @param outputDir Directory where the generated model for predicting the protein fraction should be stored.
#' @param saveRiboModels Do you want to save the riboModel?
#'
#' @return list containing the corrected count data ($counts) and the essentials needed to correct new data using the same model ($modelList).

#'
riboCorrectCounts <- function(data,
                              proteinCodingGenes,
                              outputDir,
                              saveRiboModels = T){
  # We look at the protein coding fraction within each patient sample
  # Starting data is the CPM-corrected count data
  proteinCodingFraction <- apply(data[proteinCodingGenes[proteinCodingGenes %in% rownames(data)],],2,sum)/1E6

  # We calculate the mean gene expression
  meanGenes <- apply(data,1,mean)

  # We only select the 5000 genes with the highest expression
  highMeanGenes <- names(meanGenes)[order(meanGenes,decreasing = T)][c(1:5000)]

  # We select only the means of the genes with the highest expression
  meanGenes <- meanGenes[highMeanGenes]

  # Then we calculate the variance of the genes with the highest expression
  varGenes <- apply(data[highMeanGenes,],1,var)

  normalizedData <- apply(data,2,function(x) (x[highMeanGenes]-meanGenes)/varGenes)

  set.seed(1)
  # Do cross-validation for glmnet, generalized linear model, Lasso and Elastic-Net Regularized
  modelCV <- cv.glmnet(x = t(normalizedData),y = 1-proteinCodingFraction,family = "gaussian")


  model <- glmnet(x = t(normalizedData),y = 1-proteinCodingFraction,family = "gaussian",lambda = modelCV$lambda.1se)

  allCoefficients <- as.matrix(coef(model))[which(as.matrix(coef(model)) > 0),]
  relevantCoefficients <- c(allCoefficients[1],allCoefficients[-1][allCoefficients[-1] > 0.01])

  predictProteinCoding <- 1-(relevantCoefficients[1] + t(normalizedData[names(relevantCoefficients)[-1],]) %*% relevantCoefficients[-1])

  # Divide counts for each sample by the percentage of protein coding sequence.
  # The less protein coding reads, the higher the scale-up of eventual reads.
  data <- sapply(c(1:ncol(data)),function(x) smartRound(data[,x]/predictProteinCoding[x],digits =3))
  colnames(data) <- colnames(normalizedData)

  # After scaling-up samples with limited numbers of protein coding reads,
  # the most important ribosomal proteins are removed from the dataset.
  data <- data[!(rownames(data) %in% names(relevantCoefficients)[-1]),]

  #
  modelList <- list("relevantCoefficients"=relevantCoefficients,
                    "meanGenes"=meanGenes,
                    "varGenes"=varGenes)


  riboModelList <- list("counts"=data,"riboModel"=modelList)
  if(saveRiboModels == T) {
    directory <- outputDir
    filename <- paste0(directory, "modelListRiboCounts.rds")
    if (!dir.exists(directory)) {
      dir.create(directory) }
  write_rds(riboModelList, file = filename)
  }
  return(riboModelList)
}
