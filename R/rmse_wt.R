#' Weighted Root Mean Squared Error
#' 
#' Calculate weighted or unweighted root mean squared error.
#' 
#' @details 
#' Weighted RMSE is calculated by taking the root of the weighted MSE, as calualated using `mse_wt`.
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
#' #check unweighted RMSE
#' rmse(a, b)
#' 
#' #generate random vector of weights
#' wt = runif(100)
#' 
#' #return weighted MSE
#' rmse_wt(a, b, wt)
#' 
#' @export
#' 

#calculate unweighted RMSE
rmse <- function(y_h, y){
  sqrt(mse(y_h, y))
}

#calculate weighted RMSE
rmse_wt <- function(y_h, y, wt){
  sqrt(mse_wt(y_h, y, wt))
}
