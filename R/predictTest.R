predictTest <- function(modelList, test.data) {
  for (i in seq(1:length(modelList))) {
    model <- modelList[[i]]
    prediction <- predict(model, newdata=test.data)

    if (i == 1) {
      result <- data.frame(fold1 = prediction)
    } else {
      result[, paste0("fold", i)] <- prediction
    }
  }

  return(result)
}
