library(gstat)
library(terra)
library(sp)
library(rgdal)
library(automap)
library(bulkshift)
library(tidyr)
library(spatialEco)
library(randomForest)

#define semivariogram model functions
Gau <- function(h, a) 1 - exp(-(h/a)^2)
Exp <- function(h, a) 1 - exp((-h/a))
Sph <- function(h, a){
  ifelse(
    h > a,
    1,
    ((3*h) / (2*a)) - ((1/2)*(h/a)^3)
  )
}

#setwd()

#read in aligned and masked raster layers from the working directory as a spatRaster
stack <- rast(list.files(full.names = TRUE))
stack

#use eqCov function to get uniform coverage and NAs across layers, if necessary
#stack <- eqCov(stack)

#convert to matrix for faster predictions
mat <- as.matrix(stack)
mat <- mat[complete.cases(mat), ]
gc()

#define parameters for the simulations
n_clust = 5
n = 100
ind_pts = n * 0.2
buff_width = 2000
y_hat = 'smf'
x_hat = c('bath', 'mean_v', 'back', 'x', 'y')

#indicate a directory for writing results
out = 'D:/Documents/R/scratch/'

#if a vector of cluster radii were provided, run them in a loop 
for(l in buff_width){
  
  #also loop through `ind_pts` if a vector was supplied
  for(k in ind_pts){
    
    #initialize results and run selected parameters x100
    results <- list()
    for(j in 1:100){
      cat(j, n, 'points,', k, 'independent,', l, 'm buffer            ', '\r')
      
      #randomly select n_clust locations
      p_clust <- vect(xyFromCell(stack[[y_hat]], bulkshift::bSample(stack[[y_hat]], n_clust)), crs = crs(stack))
      
      #create a spatial buffer of radius l around the cluster locations
      buff <- buffer(p_clust, width = l)
      
      #distribute k points randomly over the raster data
      p1 <- vect(xyFromCell(stack[[y_hat]], bSample(stack[[y_hat]], k)), crs = crs(stack))
      #distribute n-k points within the cluster areas
      p2 <- vect(xyFromCell(mask(stack[[y_hat]], buff), bSample(mask(stack[[y_hat]], buff), n - k)), crs = crs(stack))
      
      #plot(stack[[y_hat]])
      #plot(buff, add = TRUE)
      #plot(p1, add = TRUE)
      #plot(p2, add = TRUE)
      
      #extract environmental data at sample sites
      values(p1) <- terra::extract(stack, p1)[ ,-1]
      values(p2) <- terra::extract(stack, p2)[ ,-1]
      p_c <- rbind(p1, p2)
      
      #convert to data frame for modelling
      df <- as.data.frame(p_c)
      df <- df[complete.cases(df), ]
      
      #define model formula
      form = as.formula(
        paste0(
          y_hat,
          " ~ ", 
          paste0(x_hat, collapse="+")
        )
      )
      
      #unweighted random forest
      rf <- randomForest(form, data = df)
      
      #obtain matrix of out of bag predictions and observed values
      val <- na.omit(cbind(predict(rf), df[ ,y_hat]))
      #predict full matrix over the extent of the raster data
      true_mat <- cbind(predict(rf, mat), mat[ ,y_hat])
      
      #record the "apparent" and "true" unweighted model performance
      result_i <- data.frame(
        y = y_hat,
        n = n,
        n_clust = n_clust,
        a = k,
        b = n - k,
        buff_width = l,
        cor_true = cor(true_mat[ ,1], true_mat[ ,2]),
        mse_true = mse(true_mat[ ,1], true_mat[ ,2]),
        rmse_true = rmse(true_mat[ ,1], true_mat[ ,2]),
        ve_true = ve(true_mat[ ,1], true_mat[ ,2]),
        cor_oob = cor(val[, 1], val[, 2]),
        mse_oob = mse(val[, 1], val[, 2]),
        rmse_oob = rmse(val[, 1], val[, 2]),
        ve_oob = ve(val[, 1], val[, 2])
      )
      rm(true_mat, val)
      
      #calculate out-of-bag residuals
      df$res <- df[ ,y_hat] - predict(rf)
      #convert to sp class
      df_spat <- df; coordinates(df_spat) = ~x+y
      
      #fit a variogram model
      auto_fit <- autofitVariogram(res ~ 1, df_spat, model = c('Exp', 'Gau', 'Sph'))
      #plot(auto_fit)
      
      #extract variogram model parameters
      r = auto_fit$var_model$range[2]
      nug = auto_fit$var_model$psill[1]
      psil = auto_fit$var_model$psill[2]
      
      #calculate distance matrix between sample points
      d_mat <- distance(p_c)
      d_mat <- as.matrix(d_mat)
      
      #initialize a matrix to store the covariance values
      wt_mat <- matrix(nrow = nrow(d_mat), ncol = ncol(d_mat))
      
      for(i in 1:nrow(wt_mat)){
        #work with each sample row at a time
        d = d_mat[i, ]
        
        #transform distances from `i` to other points using the selected variogram model
        if(auto_fit$var_model$model[2] == 'Gau'){
          auto_mod <- Gau(d, r)
        }
        if(auto_fit$var_model$model[2] == 'Exp'){
          auto_mod = Exp(d, r)
        }
        if(auto_fit$var_model$model[2] == 'Sph'){
          auto_mod = Sph(d, r)
        }
        
        #transform semivariance into covariance and store in the matrix
        wt_mat[i, ] <- psil * (1 - auto_mod)
      }
      
      #get the sum of spatial covariance for each point
      wt = rowSums(wt_mat)
      #if weights are all zero, assign full weight to all points
      if(sum(wt) == 0) wt <- rep(1, length(wt))
      #scale with respect to the minimum
      wt <- min(wt)/wt
      
      #assign weights to the data and use for covariance-weighted bagging with random forest
      df <- as.data.frame(p_c)
      df$wt <-  wt
      df <- df[complete.cases(df), ]
      
      rf <- randomForest(form, data = df, weights = df$wt)
      
      #re-calculate the out-of-bag residuals using the weighted model, for validation
      df$res <- df[ ,y_hat] - predict(rf)
      df <- df[complete.cases(df), ]
      
      #convert to sp class
      df_spat <- df; coordinates(df_spat) = ~x+y
      
      #fit a variogram
      auto_fit <- autofitVariogram(res ~ 1, df_spat, model = c('Exp', 'Gau', 'Sph'))
      #plot(auto_fit)
      
      #extract variogram model parameters
      r = auto_fit$var_model$range[2]
      nug = auto_fit$var_model$psill[1]
      psil = auto_fit$var_model$psill[2]
      
      #calculate distance matrix between sample points
      d_mat <- distance(p_c)
      d_mat <- as.matrix(d_mat)
      
      #initialize a matrix to store the covariance values
      wt_mat <- matrix(nrow = nrow(d_mat), ncol = ncol(d_mat))
      
      for(i in 1:nrow(wt_mat)){
        #work with each sample row at a time
        d = d_mat[i, ]
        
        #transform distances from `i` to other points using the selected variogram model
        if(auto_fit$var_model$model[2] == 'Gau'){
          auto_mod <- Gau(d, r)
        }
        if(auto_fit$var_model$model[2] == 'Exp'){
          auto_mod = Exp(d, r)
        }
        if(auto_fit$var_model$model[2] == 'Sph'){
          auto_mod = Sph(d, r)
        }
        
        #transform semivariance into covariance and store in the matrix
        wt_mat[i, ] <- psil * (1 - auto_mod)
      }
      
      #get the sum of spatial covariance for each point
      wt = rowSums(wt_mat)
      #if weights are all zero, assign full weight to all points
      if(sum(wt) == 0) wt <- rep(1, length(wt))
      #scale with respect to the minimum
      wt <- min(wt)/wt
      
      #assign weights to the data and use for validation
      df <- as.data.frame(p_c)
      df$wt <-  wt
      df <- df[complete.cases(df), ]
      
      #obtain matrix of out-of-bag predictions and observed values
      val <- na.omit(cbind(predict(rf), df[ ,y_hat], df$wt))
      #predict full matrix over the extent of the raster data
      true_mat <- cbind(predict(rf, mat), mat[ ,y_hat])
      
      #record the "true" model performance using weighted bagging over the full raster extent 
      result_i$cor_true_wt = cor(true_mat[ ,1], true_mat[ ,2])
      result_i$mse_true_wt = mse(true_mat[ ,1], true_mat[ ,2])
      result_i$rmse_true_wt = rmse(true_mat[ ,1], true_mat[ ,2])
      result_i$ve_true_wt = ve(true_mat[ ,1], true_mat[ ,2])
      
      #record "apparent" model performance using weighted validation with out-of-bag samples
      result_i$cor_val = cov.wt(data.frame(val[ ,1], val[ ,2]), wt = val[ ,3], cor = TRUE)$cor[2]
      result_i$mse_val = mse_wt(val[, 1], val[, 2], val[, 3])
      result_i$rmse_val = rmse_wt(val[, 1], val[, 2], val[, 3])
      result_i$ve_val = ve_wt(val[, 1], val[, 2], val[, 3])

      #record this iteration in the results list
      results[[j]] <- result_i
      
      #clean up the environment
      rm(auto_fit, buff, df_spat, d, nug, psil, r, auto_mod, wt_mat, wt, p_c, p1, p2, val, rf, p_clust, d_mat, result_i, df, true_mat)
      gc()
    }
    
    #after 100x simulation runs for a given set of parameters, output the result to the specified directory
    save(results, file = paste0(out, paste('covwt', y_hat, n, n_clust, l, k, n-k, '.RData', sep = '_')))
  }
}
