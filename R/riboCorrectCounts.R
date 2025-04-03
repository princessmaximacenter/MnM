#' Correct ribodepletion protocol
#'
#' This functions is used to account for the differences in efficacy of the ribodepletion protocol.
#' It is estimated what is the percentage of protein coding transcripts within the total transcriptome, and subsequently
#' the protein coding fraction is scaled up to 100%, simultaneously removing certain residual ribosomal RNAs from the dataset.
#' @param data Dataframe containing the ribo-depleted TPM-normalized count data from an RNA-seq experiment.
#' @param proteinCodingGenes What are the names of the RNA-transcripts that stand for protein-coding genes within our dataset?
#' Please supply it as a vector.
#' @param outputDir Directory where the generated model for predicting the protein fraction should be stored.
#' @param saveRiboModels Do you want to save the riboModel? Default is FALSE.
#' @import glmnet
#' @importFrom stats coef
#' @return list containing the corrected count data ($counts) and the essentials needed to correct new data using the same model ($modelList).

riboCorrectCounts <- function(data,
                              proteinCodingGenes,
                              outputDir,
                              saveRiboModels = F){
  # We look at the protein coding fraction within each patient sample
  # Starting data is the CPM-corrected count data
  proteinCodingFraction <- base::apply(data[proteinCodingGenes[proteinCodingGenes %in% base::rownames(data)],],2,base::sum)/1E6

  # We calculate the mean gene expression
  meanGenes <- base::apply(data,1,base::mean)

  # We only select the 5000 genes with the highest expression
  highMeanGenes <- base::names(meanGenes)[base::order(meanGenes,decreasing = T)][c(1:5000)]

  # We select only the means of the genes with the highest expression
  meanGenes <- meanGenes[highMeanGenes]

  # Then we calculate the variance of the genes with the highest expression
  varGenes <- base::apply(data[highMeanGenes,],1,stats::var)

  normalizedData <- base::apply(data,2,function(x) (x[highMeanGenes]-meanGenes)/varGenes)

  base::set.seed(1)
  # Do cross-validation for glmnet, generalized linear model, Lasso and Elastic-Net Regularized
  if (length(proteinCodingFraction) < 30) {
    nfolds <- max(ceiling(length(proteinCodingFraction) / 10), 3)
  } else {
    nfolds <- 10
  }
  modelCV <- glmnet::cv.glmnet(x = base::t(normalizedData),y = 1-proteinCodingFraction,family = "gaussian",
                               nfolds = nfolds)


  model <- glmnet::glmnet(x = base::t(normalizedData),y = 1-proteinCodingFraction,family = "gaussian",lambda = modelCV$lambda.1se)

  allCoefficients <- base::as.matrix(stats::coef(model))[base::which(base::as.matrix(stats::coef(model)) > 0),]
  relevantCoefficients <- c(allCoefficients[1],allCoefficients[-1][allCoefficients[-1] > 0.01])

  predictProteinCoding <- 1-(relevantCoefficients[1] + base::t(normalizedData[base::names(relevantCoefficients)[-1],]) %*% relevantCoefficients[-1])

  # Divide counts for each sample by the percentage of protein coding sequence.
  # The less protein coding reads, the higher the scale-up of eventual reads.
  data <- base::sapply(c(1:base::ncol(data)),function(x) smartRound(data[,x]/predictProteinCoding[x],digits =3))
  base::colnames(data) <- base::colnames(normalizedData)

  # After scaling-up samples with limited numbers of protein coding reads,
  # the most important ribosomal RNAs are removed from the dataset.
  data <- data[!(base::rownames(data) %in% base::names(relevantCoefficients)[-1]),]

  #
  modelList <- base::list("relevantCoefficients"=relevantCoefficients,
                    "meanGenes"=meanGenes,
                    "varGenes"=varGenes)


  riboModelList <- base::list("counts"=data,"riboModel"=modelList)
  if(saveRiboModels == T) {
    directory <- outputDir
    filename <- base::paste0(directory, "modelListRiboCounts.rds")
    if (!base::dir.exists(directory)) {
      base::dir.create(directory) }
    base::saveRDS(riboModelList, file = filename)
  }
  return(riboModelList)
}
