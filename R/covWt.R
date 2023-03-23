#' Residual covariance spatial weights
#' 
#' Calculate covariance weighting for spatial data points based on residual spatial autocorrelation.
#' 
#' @details 
#' Calculates spatial data point weights based on the sum of residual spatial covariances for each point. These can be used e.g., for bagging or validation. 
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
#' x <- rgamma(100, 1)
#' y <- rgamma(100, 1)
#' 
#' plot(x, y)
#' 
#' df <- data.frame(x, y)
#' 
#' #get a distance matrix of the observations
#' d <- dist(df)
#' 
#' #we can use parameters from some fitted residual exponential variogram
#' m = "Sph"
#' r = 1
#' nug = 0.1
#' psil = 0.3
#' 
#' wt <- covWt(dmat = d, m = m, a = r, nug = nug, psil = psil)
#' 
#' symbols(x, y, circles = wt, inches = 0.5, xlim = c(0,6), ylim = c(0,6))
#' 
#' #or we could take a residual variogram model from gstat
#' library(gstat)
#' mod <- vgm(psil, "Sph", r, nug)
#' 
#' wt <- covWt(dmat = d, model = mod)
#' 
#' symbols(x, y, circles = wt, inches = 0.5, xlim = c(0,6), ylim = c(0,6))
#' 
#' #we can also use gstat to fit a variogram model
#' data("meuse.all")
#' plot(meuse.all$x, meuse.all$y)
#' 
#' #obtain the residual variogram after modelling the trend (this is not universal kriging)
#' mod <- lm(zinc ~ x + y, data = meuse.all)
#' p <- predict(mod)
#' meuse.all$res <- meuse.all$zinc - p
#' 
#' v <- variogram(object = res~1, locations = ~x+y, data = meuse.all)
#' plot(v)
#' fit <- fit.variogram(v, model = vgm("Exp"))
#' plot(v, fit)
#' 
#' #calculate the distance matrix
#' d = dist(meuse.all[ ,c("x", "y")])
#' 
#' calculate covariance weights using the fitted residual variogram model and plot spatially
#' wt <- covWt(dmat = d, model = fit)
#' 
#' symbols(x = meuse.all$x, y = meuse.all$y, circles = wt, inches = 0.25)
#' 
#' #calculated unweighted and covariance-weighted statistics between the predicted and observed zinc concentrations
#' rmse(p, meuse.all$zinc)
#' rmse_wt(p, meuse.all$zinc, wt)
#' 
#' ve(p, meuse.all$zinc)
#' ve_wt(p, meuse.all$zinc, wt)
#' 
#' @export

covWt <- function(dmat, model = NULL, m, a, nug, psil){
  if(class(dmat) == "dist") dmat <- as.matrix(dmat)
  if(!is.null(model)){
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