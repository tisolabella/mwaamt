## IMPORTS

## library(magrittr) ## For pipe operator %>%


## FUNCTION DEFINITIONS

## singlefit <- function(x, scale, aae) {
##     ## The function for the overall power-law fit
##     return(scale * x^(-aae))
## }

## doublefit <- function(x, scale.1, alpha.1, scale.2, alpha.2) {
##     ## The function for the double power-law fit
##     return(scale.1 * x^(-alpha.1) + scale.2 * x^(-alpha.2))
## }

## ## eq. 1 of the paper
## doublefit <- function(par, alpha.1, x) { 
##     ## The function for the double power-law fit
##     par$A * x^(-alpha.1) + par$B * x^(-par$alpha.2)
## }

## line <- function(x, a, b) {
##     return(a * x + b)
## }

set.resolution <- function(best.parameters, iteration = 2) {
    ## Function to increase the parameter scan resolution
    resolution <- 0.1 / (iteration + 1)
    return.list <- list()
    
    for (best.param in best.parameters) {
        param.list <- c(best.param - 2 * resolution, 
                        best.param - resolution, best.param,
                        best.param + resolution, 
                        best.param + 2 * resolution)
        return.list <- c(return.list, list(param.list))
    }

    return(return.list)
}

## get.chisq <- function(data, expected, sigma, ndf) {
##     ## Calculates the chi squared
##     chi <- sum(((data - expected)^2) / (sigma^2))
    
##     ## Returns the chisquared and the reduced chisq
##     return(c(chi, chi / ndf))
## }

## average <- function(l) {
##     ## Returns the average value of a list
##     s <- sum(l)
##     return(s / length(l))
## }

## weighted.average <- function(l, el) {
##     ## Returns the average of a list weigthed by el
##     s <- sum(l / el^2)
##     ws <- sum(1 / el^2)
##     return(c(s / ws, 1 / sqrt(ws)))
## }

## stddev <- function(l) {
##     ## Returns the standard deviation of a list
##     ll <- l^2
##     avg.sq <- average(ll)
##     sq.avg <- average(l)^2
##     return(sqrt(avg.sq - sq.avg))
## }
