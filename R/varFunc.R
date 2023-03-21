#' Variogram model functions
#' 
#' Functions to transform distance into semivariance.
#' 
#' @details 
#' Transforms a vector of distances `h` to semivariance \eqn{gamma(h)} according to a semivariogram model function. Note this does not include nugget or sill components
#' Use \eqn{c(0)+c(varFunc())}, for example, to plot semivariance, where \eqn{c(0)} is the nugget and $c$ is the partial sill.
#' Currently supported model functions are the exponential, Gaussian, and spherical.
#' 
#' @param h Vector of distances.
#' @param a Numerical. Fitted range parameter from a variogram model.
#' @param fun Character string indicating the variogram model function. Must be one of "Exp", "Gau", or "Sph".
#' 
#' @return A vector of semivariance values.
#' 
#' @examples
#' #set a range parameter (e.g., taken from gstat::fit.variogram())
#' r = 250
#' 
#' #generate some distances to convert to semivariances
#' d <- seq(0, 1000, 50)
#' 
#' varFunc(h=d, a=r, fun = 'Exp')
#' 
#' #look at plots of the different functions
#' plot(varFunc(h=d, a=r, fun = 'Exp'))
#' points(varFunc(h=d, a=r, fun = 'Gau'), col = 'red')
#' points(varFunc(h=d, a=r, fun = 'Sph'), col = 'blue')
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