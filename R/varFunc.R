#' Variogram model functions
#' 
#' Functions to transform distance into semivariance.
#' 
#' @details 
#' Transforms a vector of distances `h` to semivariance according to a semivariogram model function. Note this does not include a nugget component, which should be added to the output of these functions if desired.
#' 
#' @param h Vector of distances.
#' @param a Numerical. Fitted range parameter from a variogram model.
#' 
#' @examples
#' 
#' @export
#' 

varFunc <- function(h, a, fun){
  if(fun == 'Exp') g <- 1 - exp((-h/a))
  if(fun == 'Gau') g <- 1 - exp(-(h/a)^2)
  if(fun == 'Sph'){
    g <- ifelse(
      h > a,
      1,
      ((3*h) / (2*a)) - ((1/2)*(h/a)^3)
    )
  }
  
  return(g)
}