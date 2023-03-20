#' Bagger
#' 
#' A generic function for drawing bagging samples for machine learning models.
#' 
#' @details 
#' Draws bootstrap samples (with replacement) from a data frame, allowing for optional specification of probability weights and feature bagging (random selection of predictors). 
#' 
#' @param y Character naming the response variable.
#' @param x Character vector naming the predictor variables.
#' @param data Data frame containing observations of `x` and `y`.
#' @param p Optional vector of probabilities for drawing the observations.
#' @param mtry Integer indicating how many predictors to sample from `x`. Default is all.
#' @param mtry_replace Logical whether to sample the predictors with replacement, is `mtry` is less than the the number of predictors indicated in `x`.
#' 
#' @return A list containing the bootstrapped sample data `$data`, a vector of indices indicating which samples were drawn during bootstrapping `$oob`, 
#' and a vector of which predictor variables were selected `$x`.
#' 
#' @examples
#' data(iris)

#' #define the vector of predictors to include
#' x <- c('Sepal.Length', 'Sepal.Width', 'Petal.Length', 'Petal.Width')
#' 
#' b <- bagger(y = 'Species', x, data = iris)
#' 
#' #compare the number of cases to the original data
#' table(iris$Species)
#' table(b$data$Species)
#' 
#' #down-weight the probability of drawing "setosa" observations
#' iris$prob <- 1
#' iris$prob[iris$Species == 'setosa'] <- 0.5
#' 
#' #re-draw using the probabilities
#' b <- bagger(y = 'Species', x, data = iris, p = iris$prob)
#' 
#' table(iris$Species)
#' table(b$data$Species)
#' 
#' @export


bagger <- function(y, x, data, p = rep(1, nrow(data)), mtry = length(x), mtry_replace = FALSE){
  #bootstrap sample the data frame and record the index of samples in `bs`
  bs <- sample(nrow(data), replace = TRUE, prob = p)
  data_bs <- data[bs, ]
  #record which rows of the data frame were not drawn
  oob <- which(!(1:nrow(data)) %in% unique(bs))
  
  #optionally sample the predictors
  x <- sample(x, mtry, replace = mtry_replace)
  data_mtry <- data_bs[ ,c(y,x)]
  
  return(list(data = data_mtry, oob = oob, x = x))
}