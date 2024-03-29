% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mse_wt.R
\name{mse_wt}
\alias{mse_wt}
\alias{mse}
\title{Weighted Mean Squared Error}
\usage{
mse_wt(y, y_h, wt, na.rm = FALSE)

mse(y, y_h, na.rm = FALSE)
}
\arguments{
\item{y}{Vector of true values.}

\item{y_h}{Vector of predicted values.}

\item{wt}{Vector of sample weights.}

\item{na.rm}{Logical whether to remove NA values.}
}
\description{
Calculate weighted or unweighted mean squared error.
}
\details{
Weighted MSE is calculated by multiplying the sum of squared errors between \code{y} and \code{y_h} by weights provided in the vector \code{wt},
and dividing by the sum of weights, rather than N.
}
\examples{
#generate random data
a <- runif(100, 0, 3)

#simulate a linear dependent covariate with noise
b <- a + rnorm(100)

#check unweighted MSE
mse(b, a)

#generate random vector of weights
wt = runif(100)

#return weighted MSE
mse_wt(b, a, wt)

}
