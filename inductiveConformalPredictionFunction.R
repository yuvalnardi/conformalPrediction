inductiveConformalPrediction <- function(data, alpha = 0.05, n1 = NULL, 
                                         scoreFun = NULL, scoreFunArg = NULL, conformity = TRUE, 
                                         gridSize = 100, showPlot = TRUE){
  # Performs 2-dimensional Inductive Conformal Prediction using user supplied conformity measure
  #
  # Args:
  #  data:        data matrix, data.frame with 2 columns (variables)
  #  alpha:       significance level (1-alpha is the prediction confidence level)
  #  n1:          split data by the first n1 rows
  #  scoreFun:    the score function 
  #  scoreFunArg: the score function argument (e.g., k for kmeans, epsilon for dbscan, alpha for alpha shape)
  #  conformity:  scoreFun is a "conformity" function if TRUE, otherwise, it is a "non-conformity" function
  #  gridSize:    scalar containing the number of equally spaced points in each direction (x,y)
  #  showPlot:   plot prediction set if TRUE
  #
  # Returns:
  #  the prediction set
  
  ## Error handling
  if (is.null(scoreFun) || is.null(scoreFunArg)){
    stop("Missing conformity/non-conformity function and/or argument !")
  }
  
  ## split data randomly
  n <- nrow(data)
  if (is.null(n1)){
    n1 <- floor(n/2)
  }
  n2 <- n - n1
  sampleIndex <- sample(1:n, n1)
  x1 <- data[sampleIndex, ]
  x2 <- data[- sampleIndex, ]
  
  ## construct scoreFun based on x1
  scoreFunX1 <- do.call(what = scoreFun, args = list(x1, scoreFunArg))
  
  ## evaluate scoreFun (based on x1) on x2
  scoreFunX2 <- NULL
  for (i in 1:n2){
    scoreFunX2 <- c(scoreFunX2, scoreFunX1(x2[i,]))
  }
  if (conformity){
    scoreFunX2 <- sort(scoreFunX2, decreasing = FALSE)
  } else {
    scoreFunX2 <- sort(scoreFunX2, decreasing = TRUE)
  }
  
  ## commpute prediction set threshold (lamnda in {x: g(x) >= lambda})
  if (conformity){
    scoreIndex <- ceiling(alpha*(n2+1)) - 1
  } else {
    scoreIndex <- floor(alpha * (n2+1))
  }
  predThreshold <- scoreFunX2[scoreIndex]
  
  ## evaluate scoreFun (based on x1) on grid
  xPoint <- NULL
  yPoint <- NULL
  for (i in seq(min(x[,1]), max(x[,1]), len = gridSize)){
    for (j in seq(min(x[,2]), max(x[,2]), len = gridSize)){
      scoreIJ <- scoreFunX1( c(i,j) )
      if (conformity){
        if (scoreIJ >= predThreshold){
          xPoint <- c(xPoint, i)
          yPoint <- c(yPoint, j)
        }
      } else {
        if (scoreIJ <= predThreshold){
          xPoint <- c(xPoint, i)
          yPoint <- c(yPoint, j)
        }
      }
    }
  }
  predictionSet <- cbind(xPoint, yPoint)
  
  ## plot data and prediction set
  if (showPlot){
    plot(x, xlab = "X1", ylab = "X2")
    points(predictionSet[,1], predictionSet[,2], col = "red")
  }
  return(predictionSet)
}

################################################################

## example of a conformity measure using K-means (arg = K)
kmeansScoreFun <- function(x, k){
  fit <- kmeans(x, centers = k)
  centers <- fit$centers
  scoreFunction <- function(dataRecord){
    minDistance <- - min(apply(dataRecord-centers, 1, FUN = function(x) {sqrt(sum(x^2))}))
  }
  return(scoreFunction)
}

## example of a conformity measure using 2d KDE (arg = h)
library(ks)
kdeScoreFun <- function(x, h){
  library(MASS)
  
}

################################################################

## example of running "inductiveConformalPrediction"
#P <- inductiveConformalPrediction(x, scoreFun = "kmeansScoreFun", scoreFunArg = 2, conformity = TRUE, gridSize = 100)
#P <- inductiveConformalPrediction(x, scoreFun = "kdeScoreFun", scoreFunArg = 1, conformity = TRUE, gridSize = 100)