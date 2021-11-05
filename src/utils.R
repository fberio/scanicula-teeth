# Convert matrix into array and array into matrix
mat2array <- function(x,dimx= 2){
  if (is.vector(x)) x <- t(as.matrix(x))
  if ((dimx > 3) || (dimx < 2)) stop("incorrect dimx parameter")
  if(ncol(x)/dimx != ncol(x)%/%dimx) stop("the number of cols of data is not a multiple of dimx")
  npts <- ncol(x)/dimx
  if (npts - trunc(npts) != 0) stop("dimx and column number are not compatible")
  n <- nrow(x)
  tem <- t(x)
  tem <- array(tem, c(dimx, npts, n))
  tem <- aperm(tem, c(2, 1, 3))
  return(tem)
}

array2mat <- function(arr) {
  bb <- function(x) x <- c(t(x))
  out <- t(apply(arr,3,bb))
  return(out)
}

# Perform PCA
Rpca <- function(mat, lpkcor = TRUE) {
  if (lpkcor) {
    mat <- scale(mat, scale=FALSE) 
    s <- svd(mat, nu=0)
    sdev <- s$d
    val <- s$d * s$d
    vec <- s$u
    scores <- scale(mat, scale=FALSE) %*% s$v
    percent <- 100*val/sum(val)
    cpercent <- cumsum(percent)
  }
  else {
    acp <- prcomp(mat)
    val <- acp$sdev * acp$sdev
    sdev <-acp$sdev
    scores <- acp$x
    vec <- acp$rotation
    percent <- 100*val/sum(val)
    cpercent <- cumsum(percent)
  }
  pca <- list(scores= scores, val= val, sdev= sdev, vec= vec, percent= percent, cpercent= cpercent)
  return(pca)
}

# Functions to compute accuracy, recall, and precision
computeAccuracy <- function(confusionMatrix) {
  return (sum(diag(confusionMatrix))/sum(confusionMatrix))
}

computeRecall <- function(confusionMatrix, groupName){
  return(confusionMatrix[groupName, groupName]/sum(confusionMatrix[groupName,]))
}

computePrecision <- function(confusionMatrix, groupName){
  return(confusionMatrix[groupName, groupName]/sum(confusionMatrix[,groupName]))
}

getTestSets <- function(dataset, nbFolds) {
  stopifnot(nbFolds > 1)
  if (class(dataset) == "data.frame") {
    stopifnot(nrow(dataset) > nbFolds)
    nbItemsPerSet <- round(nrow(dataset) / nbFolds)
    testSets <- list()
    tmp <- dataset
    for (i in 1:(nbFolds-1)) {
      testSets[[i]] <- tmp %>% sample_n(nbItemsPerSet)
      tmp <- suppressMessages(anti_join(tmp, testSets[[i]]))
    }
    testSets[[nbFolds]] <- anti_join(tmp, testSets[[nbFolds - 1]])
  } else {
    if (class(dataset) == "matrix") {
      nbRows <- dim(dataset)[1]
      stopifnot(nbRows > nbFolds)
      x <- 1:nbRows
      testSets <- list()
      tmp <- x
      nbItemsPerSet <- round(dim(dataset)[1] / nbFolds)
      for (i in 1:(nbFolds-1)) {
        indexes <- sample(x=tmp, size=nbItemsPerSet, replace=FALSE)
        testSets[[i]] <- indexes
        tmp <- tmp[-indexes]
      }
      testSets[[nbFolds]] <- tmp
    } else {
      stop("2020 /08/04 Yann : Seuls les data.frame et matrix sont supportés.")
    }
  }
  return(testSets)  
}