library(gstat)
library(terra)
library(sp)
library(rgdal)
library(automap)
library(bulkshift)
library(tidyr)
library(spatialEco)
library(randomForest)
library(covBagging)

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

#aggregated the stack 20x
agg <- aggregate(stack, 20)

#convert to sp grid for conditional Gaussian simulation
grid <- as.data.frame(agg[[y_hat]], xy = TRUE)
grid <- grid[ ,c('x', 'y')]
gridded(grid) = ~x+y
crs(grid) <- crs(agg)

#if a vector of cluster radii were provided, run them in a loop 
for(l in buff_width){
  
  #also loop through `ind_pts` if a vector was supplied
  for(k in ind_pts){
    
    #initialize results and run selected parameters x100 (using `while` because simulations can fail)
    results <- list()
    while(length(results) < 100){
      j = length(results) + 1
      cat(j, n, 'points,', k, 'independent,', l, 'm buffer            ', '\r')
      
      #randomly select n_clust locations
      p_clust <- vect(xyFromCell(stack[[y_hat]], bSample(stack[[y_hat]], n_clust)), crs = crs(stack))
      
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
      
      #calculate out of bag residuals
      df$res <- df[ ,y_hat] - predict(rf)
      #convert to sp class
      df_spat <- df; coordinates(df_spat) = ~x+y
      crs(df_spat) <- crs(agg)
      
      #fit a variogram model
      auto_fit <- autofitVariogram(res ~ 1, df_spat, model = c('Exp', 'Gau', 'Sph'))
      #plot(auto_fit)
      
      #generate KDE for sample points
      suppressWarnings(
        suppressMessages(
          kde <- sp.kde(x = df_spat, newdata = raster::raster(stack[[y_hat]]), scale.factor = 1000000)
        )
      )
      
      #plot(kde)
      
      #extract KDE at sample points
      df$kde <-  raster::extract(kde, df_spat)
      df <- df[complete.cases(df), ]
      
      #bin the density estimates into quantiles
      qt <- quantile(df$kde, seq(0, 1, 0.05))
      df$intvl <- findInterval(df$kde, qt)
      
      #compute the sd of each residual quantile
      lo_df <- data.frame(
        kde_mean = aggregate(df$kde, list(df$intvl), mean)[,2],
        res_sd = aggregate(df$res, list(df$intvl), sd)[,2]
      )
      lo_df <- lo_df[complete.cases(lo_df), ]
      
      #fit loess between sd and KDE mean per quantile bin
      lo <- loess(res_sd ~ kde_mean, data = lo_df, span = 0.5, degree = 0, control = loess.control(surface = "direct"))
      
      #plot(lo_df$kde_mean, lo_df$res_sd); lines(lo_df$kde_mean, lo$fitted)
      rm(lo_df)
      
      #compute muOK to the aggregated grid size
      krg <- krige(
        res ~ 1, 
        remove.duplicates(df_spat), 
        newdata = grid, 
        nmax = 50,
        model = auto_fit$var_model,
        debug.level = 0
      )
      
      muOK <- mean(krg$var1.pred, na.rm=T)
      rm(krg)
      
      #calculate zero-centered unit variance residuals by subtracting muOK and dividing by sd, as predicted by loess model
      df$res_prime <- (df$res - muOK) / predict(lo, newdata = df$kde)
      
      #convert to sp class
      df_spat <- df; coordinates(df_spat) = ~x+y
      crs(df_spat) <- crs(agg)
      
      #fit a variogram model to standardized residuals
      try(
        auto_fit_prime <- autofitVariogram(res_prime ~ 1, df_spat, model = c('Exp', 'Gau', 'Sph'))
      )
      #plot(auto_fit_prime)
      
      #try 500 conditional simulations using standardized residuals with the coarse sp grid and max of 50 local points
      try(
        res_map <- gstat::krige(
          formula = res_prime ~ 1, 
          remove.duplicates(df_spat), 
          newdata = grid,
          model = auto_fit_prime$var_model,
          nsim = 500,
          nmax = 50,
          debug.level = 0
        )
      )
      #spplot(res_map[1:9])
      
      #if simulations did not fail, back-transform and save the result
      if(exists('res_map')){
        
        #convert simulated residual surfaces back to terra
        res_sr <- rast(raster::stack(res_map))
        rm(res_map)
        crs(res_sr) <- crs(agg)
        
        #resample to aggregated grid size
        res_sr <- resample(res_sr, agg)
        
        #back-transform residual prime maps
        kde <- rast(kde); names(kde) <- 'kde_mean'
        kde <- resample(kde, res_sr)
        
        lo_pred <- terra::predict(kde, lo)
        res_sr <- (res_sr * lo_pred) + muOK
        
        #predict the random forest model for the aggregated grid
        agg_pred <- predict(agg, rf)
        #calculate the simulated response surface using the random forest prediction and each simulated residual surface
        yhat_sim <- agg_pred + res_sr
        
        #plot(c(agg[y_hat], yhat_sim[[1]]))
        
        #get the full raster prediction and true response surface
        true_mat <- na.omit(as.matrix(c(predict(agg, rf), agg[y_hat])))
        
        #record the "true" unweighted model performance using the aggregated grids
        result_i$cor_true_agg = cor(true_mat[ ,1], true_mat[ ,2])
        result_i$mse_true_agg = mse(true_mat[ ,1], true_mat[ ,2])
        result_i$rmse_true_agg = rmse(true_mat[ ,1], true_mat[ ,2])
        result_i$ve_true_agg = ve(true_mat[ ,1], true_mat[ ,2])
        
        #calculate validation stats from each conditional simulated response surface
        val500 <- list()
        for(i in 1:nlyr(yhat_sim)){
          
          #for each simulated response surface, get the predictions and simulation as a matrix
          v <- na.omit(as.matrix(c(agg_pred, yhat_sim[[i]])))
          
          #record the validation stats for layer i
          val500[[i]] <- data.frame(
            cor_val = cor(v[ ,1], v[ ,2]),
            mse_val = mse(v[ ,1], v[ ,2]),
            rmse_val = rmse(v[ ,1], v[ ,2]),
            ve_val = ve(v[ ,1], v[ ,2])
          )
        }; rm(v)
        val500 <- do.call(rbind, val500)
        
        #take the mean of the validation statistics over 500 simulations
        result_i$cor_val = mean(val500$cor_val)
        result_i$mse_val = mean(val500$mse_val)
        result_i$rmse_val = mean(val500$rmse_val)
        result_i$ve_val = mean(val500$ve_val)
        
        #record this iteration in the results list
        results[[j]] <- result_i
      }
      
      #clean up the environment
      rm(auto_fit, auto_fit_prime, buff, df_spat, p_c, p1, p2, rf, p_clust, result_i, df, true_mat, res_sr, muOK, kde, lo, lo_pred, qt, agg_pred, val500, yhat_sim)
      gc()
    }
    
    #after 100x simulation runs for a given set of parameters, output the result to the specified directory
    save(results, file = paste0(out, paste('het', y_hat, n, n_clust, l, k, n-k, '.RData', sep = '_')))
  }
}