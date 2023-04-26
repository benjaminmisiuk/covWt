#' Weighted Variance Explained
#' 
#' Calculate weighted or unweighted variance explained.
#' 
#' @details 
#' Unweighted variance explained is calculated by dividing the residual sum of squares between `y` and `y_h` from the total sum of squares of `y`, and subtracting from 1.
#' Weighted variance is obtained by weighting the residual and total sum of squares by the vector `wt` prior to dividing.
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
#' #check unweighted VE
#' ve(b, a)
#' 
#' #generate random vector of weights
#' wt = runif(100)
#' 
#' #return weighted VE
#' ve_wt(b, a, wt)
#' 
#' @export
#' 

#calculate weighted VE
ve_wt <- function(y, y_h, wt, na.rm = FALSE){
  if(na.rm){
    na <- is.na(y_h)|is.na(y)
    y_h <- y_h[!na]
    y <- y[!na]
    wt <- wt[!na]
  }
  
  SSres = sum((y - y_h)^2 * wt) / sum(wt)
  SStot = sum((y - mean(y))^2 * wt) / sum(wt)
  1 - ((SSres)/(SStot))
}

#' @rdname ve_wt
#' @export
#' 

#calculate unweighted VE
ve <- function(y, y_h, na.rm = FALSE){
  if(na.rm){
    na <- is.na(y_h)|is.na(y)
    y_h <- y_h[!na]
    y <- y[!na]
  }
  
  SSres = sum((y - y_h)^2)
  SStot = sum((y - mean(y))^2)
  1 - (SSres/SStot)
}