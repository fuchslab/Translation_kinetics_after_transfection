#' @title Simulates the Euler-Maruyama approximation of the SDE solution
#' @inheritParams deterministic
#' @param tstep the step size used for the Euler-Maruyama scheme
#' @author Susanne Pieschner
#' @return  A matrix. The first column contains the simulation time, the other columns contain the species
#' levels
#' @export

diffusion <- function(model, maxtime, tstep)
{
  numSteps <- ceiling(maxtime/tstep) + 1
  x <- model$initial 
  xmat <- matrix(rep(NA, times = numSteps* length(model$initial)) , ncol = length(model$initial))
  xmat[1,] <- x
  stoic <- model$stoic 
  lbound <- model$lbound
  ubound <- model$ubound
  pars <- model$pars
  sqrt_dt <- sqrt(tstep)
  get_haz <- model$get_haz
  # browser()
  
  
  for (i in 2:numSteps){
    # z <- rnorm(ncol(model$stoic), 0, sqrt_dt)
    # #browser()
    # h <- model$get_haz(x, pars, systemSize)
    # x <- x +
    #   stoic_scaled %*% h * tstep +
    #   stoic_scaled %*% diag(sqrt(h)) %*% z
    #
    # #Reflecting barrier
    # #x[x < 0] <- 0#xmat[i-1,][x<0] #???#
    # #browser()
    # x[x < lbound] <- lbound[x < lbound]
    # x[x > ubound] <- ubound[x > ubound]
    # xmat[i,] <- x
    
    #browser()
    h <- model$get_haz(x, pars)
    updated <- FALSE
    j <- 0
    while(!updated){
      z <- rnorm(ncol(model$stoic), 0, sqrt_dt)
      x_prop <- x +
        stoic %*% h * tstep +
        stoic %*% diag(sqrt(h)) %*% z
      j<-j+1
      if(j>1000) browser()
      if(all(c(x_prop >= lbound, x_prop <= ubound))){
        updated <- TRUE
        x <- x_prop
        xmat[i,] <- x
      }
    }
    
    # #browser()
    # h <- model$get_haz(x*systemSize, pars, systemSize)
    # updated <- FALSE
    # j <- 0
    # while(!updated){
    #   z <- rnorm(ncol(model$stoic), 0, sqrt_dt)
    #   x_prop <- x +
    #     stoic_scaled %*% h * tstep +
    #     stoic_scaled %*% diag(sqrt(h)) %*% z
    #
    #   j<-j+1
    #   if(j>1000) browser()
    #   eps <- 10^-( ceiling(log(ubound, base =10)) + 15)
    #   if(all(c(x_prop >= lbound * (1-eps), x_prop <= ubound * (1+eps)))){
    #     if(!all(c(x_prop >= lbound, x_prop <= ubound))){
    #       x_prop[x_prop < lbound] <- lbound[x_prop < lbound]
    #       x_prop[x_prop > ubound] <- ubound[x_prop > ubound]
    #     }
    #     updated <- TRUE
    #     x <- x_prop
    #     xmat[i,] <- x
    #
    #   }
    # }
  }
  
  times <- seq(0, maxtime, by = tstep)
  xmat <- cbind(times, xmat)
  
  if(times[length(times)] < maxtime){
    lastStep <- maxtime - times[length(times)]
    z <- rnorm(ncol(model$stoic), 0, sqrt(lastStep))
    h <- model$get_haz(x, pars, systemSize)
    x <- x +
      stoic %*% h * lastStep +
      stoic %*% diag(sqrt(h)) %*% z
    xmat <- rbind(xmat, c(maxtime, x))
  }
  
  colnames(xmat) <- c("Time", rownames(model$stoic))
  return(xmat)
}


