#' covWt
#' 
#' Calculate covariance weighting for spatial data points based on residual spatial autocorrelation.
#' 
#' @details 
#' Calculates spatial data point weights based on the sum of residual spatial covariances for each point.
#' Provide either a variogram model from [gstat] or the parameters of a model fitted to prediction residuals.
#' Currently supports exponential ("Exp"), Gaussian ("Gau"), and spherical ("Sph") models. 
#' Currently does not support nested variogram models. If you would like to use these, call `covWt` separately for each model and sum them.
#' 
#' @param model variogramModel class, for example as returned from [gstat::fit.variogram()] or [automap::autofitVariogram()].
#' @param m Character. Type of variogram model. Must be one of "Exp", "Gau", or "Sph". Ignored if `model` is provided.
#' @param a Numeric. Range parameter of a variogram model. Ignored if `model` is provided.
#' @param nug Numeric. Nugget parameter of a variogram model. Ignored if `model` is provided.
#' @param psil Numeric. Partial sill of a variogram model. Ignored if `model` is provided.
#' @param dmat Distance matrix. Can be of class "dist" or a symmetric matrix of distances between points with diagonals included.
#' 
#' @return A vector of weights.
#' 
#' @examples
#' 
#' 
#' @export

covWt <- function(x, model, m, a, nug, psil, dmat){
  if(class(dmat) == "dist") dmat <- as.matrix(dmat)
  if(exists("model")){
    if(class(model)[1] != "variogramModel") stop("`model` must be of class 'variogramModel'")
    
    m = as.character(model$model[2])
    a = model$range[2]
    nug = model$psill[1]
    psil = model$psill[2]
  }
  
  #initialize a matrix to store the covariance values
  wt_mat <- matrix(nrow = nrow(dmat), ncol = ncol(dmat))
  
  for(i in 1:nrow(wt_mat)){
    #work with each sample row at a time
    d = dmat[i, ]
    #transform distances from `i` to other points using the selected variogram model
    sv <- varFunc(d, a, m)
    #transform semivariance into covariance and store in the matrix
    wt_mat[i, ] <- nug + psil - (nug + psil * sv) #nug + psil is the sill, nug + psil * auto_mod is semivariance
  }
  
  #get the sum of spatial covariance for each point
  wt = rowSums(wt_mat)
  #if weights are all zero, assign full weight to all points
  if(sum(wt) == 0) wt <- rep(1, length(wt))
  #scale with respect to the minimum
  wt <- min(wt)/wt
  return(wt)
}