#' Bagger
#' 
#' A generic function for creating a bagging model ensemble. Supports weighted bagging and feature bagging.
#' 
#' @details 
#' Any model that supports the `formula` class may be used, but it might be necessary to specify additional parameters via the `model_params` argument.
#' Models from external libraries must be loaded into the global environment (e.g., using `library()`). Consult individual model documentation
#' for guidance on what goes in `model_params`. Some unique `formula` formats, like in [gam], are not currently supported. You can call and loop `bag` directly 
#' to make a custom bagging model for these.
#' 
#' @param y Character naming the response variable.
#' @param x Character vector naming the predictor variables.
#' @param data Data frame containing observations of `x` and `y`.
#' @param p Optional vector of probabilities for drawing the observations.
#' @param mtry Integer indicating how many predictors to sample from `x`. Default is all.
#' @param mtry_replace Logical whether to sample the predictors with replacement.
#' @param model Character. Name of the model to bag. Must support class `formula`, such as "glm", "earth", "rpart", "loess", "gbm".
#' @param model_params List of parameters to pass to  `model`. See details and examples.
#' @param n Integer. Number of bagging iterations.
#' 
#' @return List containing the models in the ensemble `$models`, the out-of-bag predictions `$oob_preds`, and the number of times each sample is drawn `$n_draws`.
#' 
#' @examples
#' data <- data.frame(iris)
#' 
#' #define the response and a vector of predictors to include
#' y = 'Sepal.Length'
#' x <- c('Species', 'Sepal.Width', 'Petal.Length', 'Petal.Width')
#' 
#' #create a bagging model using multivariate adaptive regression splines
#' library(earth)
#' 
#' B <- bagger(y, x, data = data, model = "earth")
#' 
#' #check the out-of-bag variance explained
#' ve(B$oob_preds, data$Sepal.Length)
#' 
#' #observe how many times each observation was drawn
#' B$n_draws / length(B$models)
#' 
#' #bag again using probability weights
#' #we will calculate them using the oob residuals just as an example
#' r <- abs(data$Sepal.Length - B$oob_preds)
#' wt <- (r - min(r)) / (max(r) - min(r))
#' 
#' B <- bagger(y, x, data = data, model = "earth", p = wt)
#' 
#' #observe the bootstrap draws now
#' B$n_draws / length(B$models)
#' 
#' #you can pass specific model parameters to the bagger
#' B <- bagger(y, x, data = data, model = "glm", model_params = list(family = gaussian(link = 'log')))
#' 
#' #predict with the whole ensemble
#' p <- lapply(B$models, function(x) predict(x, data, type = 'response'))
#' p <- do.call(cbind, p)
#' p <- apply(p, 1, mean, na.rm = TRUE)
#' 
#' @export

bagger <- function(y, x, data, p = rep(1, nrow(data)), mtry = length(x), mtry_replace = FALSE, model, model_params = NULL, n = 500){
  oob_preds <- list()
  mod_bag <- list()
  
  for(i in 1:n){
    #draw a sample
    b <- bag(y, x, data, p, mtry, mtry_replace)
    
    #define the formula
    f <- as.formula(
      paste0(
        y, " ~ ", paste0(b$x, collapse="+")
      )
    )
    
    #run the model
    model_i <- do.call(model, c(list(f, data = quote(b$data)), model_params))
    
    #generate an empty matrix to store the model predictions
    pred <- matrix(nrow = nrow(data))
    
    #predict the out-of-bag samples that were not drawn
    df_oob <- data.frame(data[b$oob, b$x])
    names(df_oob) <- b$x
    pred[b$oob] <- predict(model_i, df_oob)
    
    #record the out-of-bag predictions and the model
    oob_preds[[length(oob_preds) + 1]] <- pred
    mod_bag[[length(mod_bag) + 1]] <- model_i
  }
  
  #take the average of the out-of-bag predictions and count the draws for each observation
  oob_preds <- do.call(cbind, oob_preds)
  n_draws <- apply(oob_preds, 1, function(x) sum(is.na(x)))
  oob_preds <- apply(oob_preds, 1, mean, na.rm=TRUE)
  
  return(list(models = mod_bag, oob_preds = oob_preds, n_draws = n_draws))
}