#' Weighted Pearson's correlation
#' 
#' Calculate a weighted version of Pearson's correlation coefficient.
#' 
#' @details 
#' This approach relies on normalizing the weights to obtain "reliability weights", which sum to 1. 
#' These are used to calculate the weighted covariance between predictions `y_h` and true values `y` using their weighted means (also calculated using the normalized weights).
#' The correlation is then calculated using the weighted covariance and the weighted variances for `y_h` and `y`. 
#' 
#' @param y_h Vector of predicted values.
#' @param y Vector of true values.
#' @param wt Vector of sample weights.
#' @param na.rm Logical whether to remove NA values.
#' 
#' @examples
#' #generate random data
#' a <- runif(100, 0, 3)
#' 
#' #simulate a linear dependent covariate with noise
#' b <- a + rnorm(100)
#' 
#' #check unweighted correlation
#' cor(a, b)
#' 
#' #generate random vector of weights
#' wt = runif(100)
#' 
#' #return weighted correlation
#' r_wt(a, b, wt)
#' 
#' #compare to the cov.wt function
#' cov.wt(x = data.frame(a, b), wt = wt, cor = TRUE)$cor[2]
#' 
#' @export
#' 

r_wt <- function(y_h, y, wt, na.rm = FALSE){
  if(na.rm){
    na <- is.na(y_h)|is.na(y)
    y_h <- y_h[!na]
    y <- y[!na]
  }
  
  #normalize weights
  wt <- wt / sum(wt)
  
  #get weighted means
  w_mean <- c(sum(wt*y_h),
              sum(wt*y))
  
  #calculate weighted covariance
  covar <- matrix(
    c(
      y_h - w_mean[1],
      y - w_mean[2]
    ), 
    ncol = 2
  )
  
  covar <- covar[ ,1] * covar[ ,2]
  covar <- sum(covar * wt)
  covar <- covar/(1 - sum(wt^2))
  
  #convert to weighted correlation using the weighted variances
  var_y_h <- sum(wt * (y_h - w_mean[1])^2) / (1 - sum(wt^2))
  var_y <- sum(wt * (y - w_mean[2])^2) / (1 - sum(wt^2))
  
  return(covar/sqrt(var_y_h * var_y))
  
  #useful refs:
  #https://en.wikipedia.org/wiki/Weighted_arithmetic_mean
  #https://stats.stackexchange.com/questions/61225/correct-equation-for-weighted-unbiased-sample-covariance
}
