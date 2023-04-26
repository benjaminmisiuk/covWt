#' Weighted Root Mean Squared Error
#' 
#' Calculate weighted or unweighted root mean squared error.
#' 
#' @details 
#' Weighted RMSE is calculated by taking the root of the weighted MSE, as calualated using `mse_wt`.
#' 
#' @param y Vector of true values.
#' @param y_h Vector of predicted values.
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
#' #check unweighted RMSE
#' rmse(b, a)
#' 
#' #generate random vector of weights
#' wt = runif(100)
#' 
#' #return weighted MSE
#' rmse_wt(b, a, wt)
#' 
#' @export
#' 

#calculate weighted RMSE
rmse_wt <- function(y, y_h, wt, na.rm = FALSE){
  if(na.rm){
    na <- is.na(y_h)|is.na(y)
    y_h <- y_h[!na]
    y <- y[!na]
    wt <- wt[!na]
  }
  
  sqrt(mse_wt(y, y_h, wt))
}

#' @rdname rmse_wt
#' @export
#' 

#calculate unweighted RMSE
rmse <- function(y, y_h, na.rm = FALSE){
  if(na.rm){
    na <- is.na(y_h)|is.na(y)
    y_h <- y_h[!na]
    y <- y[!na]
  }
  
  sqrt(mse(y, y_h))
}

