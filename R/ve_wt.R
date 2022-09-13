#' Weighted Variance Explained
#' 
#' Calculate weighted or unweighted variance explained.
#' 
#' @details 
#' Unweighted variance explained is calculated by dividing the residual sum of squares between `y_h` and `y` from the total sum of squares of `y`, and subtracting from 1.
#' Weighted variance is obtained by weighting the residual and total sum of squares by the vector `wt` prior to dividing.
#' 
#' @param y_h Vector of predicted values.
#' @param y Vector of true values.
#' @param wt Vector of sample weights.
#' 
#' @examples
#' #generate random data
#' a <- runif(100, 0, 3)
#' 
#' #simulate a linear dependent covariate with noise
#' b <- a + rnorm(100)
#' 
#' #check unweighted VE
#' ve(a, b)
#' 
#' #generate random vector of weights
#' wt = runif(100)
#' 
#' #return weighted VE
#' ve_wt(a, b, wt)
#' 
#' @export
#' 

#calculate unweighted VE
ve <- function(y_h, y){
  SSres = sum((y - y_h)^2)
  SStot = sum((y - mean(y))^2)
  1 - (SSres/SStot)
}

#calculate weighted VE
ve_wt <- function(y_h, y, wt){
  SSres = sum((y - y_h)^2 * wt) / sum(wt)
  SStot = sum((y - mean(y))^2 * wt) / sum(wt)
  1 - ((SSres)/(SStot))
}
