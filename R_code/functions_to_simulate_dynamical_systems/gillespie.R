#' @title Stochastic simulation using the standard Gillespie algorithm
#' @inheritParams deterministic
#' @param tstep equidistant output time step
#' @author Susanne Pieschner based on issb package by Colin Gillespie
#' @return  A matrix whose first column contains the simulation time, the other 
#'  columns contain the process states.
#' @export

gillespie <- function(model, maxtime, tstep=NULL)
{
  sim_time <- 0
  i <- 1
  # initialize time-state matrix
  x <- model$initial
  numCol <- length(model$initial) + 1
  maxSteps <- 100000 
  xmat <- matrix(NA, nrow = maxSteps, ncol = numCol)
  
  stoic <- model$stoic
  get_haz <- model$get_haz
  param <- model$pars
  h <- get_haz(x, param)
  
  while(sim_time <= maxtime && sum(h) > 0){
    xmat[i, ] <- c(sim_time, x)
    # add time to next reaction
    sim_time <- sim_time + rexp(1, sum(h))
    # determine type of next reaction
    j <- sample(length(h), size = 1, prob = h)
    x <- x + stoic[ ,j]
    h <- get_haz(x, param)
    i <- i + 1
  }
  
  if(sum(h) == 0){
    xmat[i, ] <- c(sim_time, x)
    i <- i + 1
  }
  
  xmat[i, ] = c(maxtime, xmat[i-1, 2:numCol])
  
  if(!is.null(tstep)){# return observations with equidistant time step equal to tstep
    time <- seq(from = 0, to = maxtime, by = tstep)
    indizes <- sapply(time, function(x) max(which(xmat[,1]<=x)))
    xmat[1:length(indizes), 2:numCol] <- xmat[indizes, 2:numCol]
    xmat[1:length(indizes), 1] <- time
    xmat <- xmat[1:length(indizes), ]
  }else{
    xmat <- xmat[1:i, ]
  }
  
  colnames(xmat) = c("Time", rownames(stoic))
  return(xmat)
}

