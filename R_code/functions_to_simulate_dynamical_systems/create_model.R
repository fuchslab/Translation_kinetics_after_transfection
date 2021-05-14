#' @title A function to create a simple model of the translation kinetics after 
#' mRNA transfection
#' @param initial a vector containing the initial conditions of the model.
#' @param pars a vector containing the parameter values.
#' @author Susanne Pieschner
#' @return A list
#' @export

create_model_mRNA1 = function(initial, theta){
  # translation kinetics - model 1
  # This system has two species: mRNA, GFP
  # three reactions: 1) mRNA degradation, 2) translation, 3) GFP degradation
  # three parameters: theta[1] - mRNA degradation rate constant
  #                   theta[2] - translation rate constant
  #                   theta[3] - GFP degradation rate constant
  
  # hazards/ propensities / transition rates
  hazards <- function(x, pars) {
    hazs <- numeric(length(pars))
    hazs[1] <- theta[1] * x[1]
    hazs[2] <- theta[2]     * x[1]
    hazs[3] <- theta[3] * x[2]
    return(hazs)
  }
  
  # stoichiometry matrix
  smat <- matrix(0,nrow=2,ncol=3)
  smat[1,1] <- -1
  smat[2,2] <- 1
  smat[2,3] <- -1
  rownames(smat) <- c("mRNA", "GFP")
  
  # for the linear noise approximation, we need to specify the Jacobian of the 
  # hazard rates:
  # rows x columns = reactions x species
  jac <- function(x, pars)
  {
    fmat <- matrix(0, nrow=3, ncol=2)
    fmat[1,1] <- theta[1]
    fmat[2,1] <- theta[2]
    fmat[3,2] <- theta[3]
    return(fmat)
  }
  
  lbound <- c(0,0)
  ubound <- c(initial[1], Inf)
  
  get_haz <- function(x, p=theta)  hazards(x, p)
  get_jacobian <- function(x, p=theta) jacobian(x, p)
  
  model = list(initial = initial,
               pars = theta,
               stoic = smat,
               get_haz = get_haz,
               get_jacobian = get_jacobian,
               lbound = lbound,
               ubound = ubound)
  return(model)
}
