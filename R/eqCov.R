#' Equal Coverage spatRaster Layers
#' 
#' Generate a spatRaster with multiple layers where coverage and NA values are consistent across layers.
#' 
#' @details 
#' If a cell in any spatRaster layer is NA, all layers in the spatRaster are set to NA at that cell. The result is reduced using `trim` to remove 
#' outer rows or columns with only NAs.
#' 
#' @param x spatRaster containing multiple layers.
#' 
#' @examples
#' #create three spatRasters with uniform values
#' a <- rast(ncols = 10, nrows = 10)
#' b <- rast(ncols = 10, nrows = 10)
#' c <- rast(ncols = 10, nrows = 10)
#' 
#' values(a) <- 1
#' values(b) <- 2
#' values(c) <- 3
#' 
#' #set random NAs in each
#' a[sample(100, 10)] <- NA
#' b[sample(100, 10)] <- NA
#' c[sample(100, 10)] <- NA
#' 
#' #create a multi-layer spatRaster
#' r <- rast(list(a, b, c))
#' plot(r, colNA = 'black')
#' 
#' #set equal coverages
#' r <- eqCov(r)
#' plot(r, colNA = 'black')
#' 
#' #example where an outer row is na
#' a[1:10] <- NA
#' 
#' r <- rast(list(a, b, c))
#' plot(r, colNA = 'black')
#' 
#' r <- eqCov(r)
#' plot(r, colNA = 'black')
#' 
#' @import terra
#' 
#' @export
#' 

eqCov <- function(x){
  #create a copy of x to create the mask
  msk <- rast(x)
  
  #get locations of NA cells for all layers
  for (i in 1:nlyr(msk)) {
    msk[[i]] <- !is.na(x[[i]])
  }
  
  #use `prod` to determine whether any cell value at a given location is NA
  msk<- prod(msk)
  msk[msk == 0] <- NA
  
  #mask and trim x
  x <- mask(x, msk)
  x <- trim(x)
  
  return(x)
}