###
#' This function creates the simulation of Setoguchi et.al. (2008)
#'
#' @description For each simulated dataset, ten covariates (four confounders associated 
#' with both exposure and outcome, three exposure predictors, and three outcome predictors)
#' W_i were generated as standard normal random variables with zero mean and unit variance.
#' Correlations were induced between several of the variables. The binary exposure
#' A has Pr(A=1|W_i) = 1/(1+exp(−beta * f(W_i))). The average exposure probability (in other words,
#' the exposure probability at the average of covariates) was approx. 0.5 and was modeled from W_i
#' according to the scenarios using the formulae provided by Setoguchi et al. 2008 The
#' continuous outcome Y was generated from a linear combination of A and Wi such that 
#' Y = alpha_i*W_i + gamma*A where the effect of exposure, gamma = −0.4.
#'
#' @param n size of the simulation (2000)
#' @param sim simulation scenarios A-G of Setoguchi 2008 (see details)
#' @param seed an RNG seed for repeatability
#' @param cont Continuous (TRUE) simulation of Lee et.al. (2010) or binary (FALSE) simulation of 
#' Setoguchi et.al.(2008).
#' @param plot.it visualize the simulated data. Helpful for debugging.
#' 
#' @details
#' Setoguchi (2008) was interested in simulating a binary response representing occurance of cholorectal 
#' cancer or effect of hormone replacement therapy on fractures. The coefficients in data generation formulae
#' for the simulation (\emph{sim} argument) scenarios are as follows:
#' \itemize{
#' \item A: Model with additivity and linearity
#' \item B: Model with mild non-linearity
#' \item C: Model with moderate non-linearity
#' \item D: Model with mild non-additivity
#' \item E: Model with mild non-additivity and non-linearity
#' \item F: Model with moderate non-additivity
#' \item G: Model with moderate non-additivity and non-linearity
#' }
#' The coefficients used for the formulae are based on claims data, modeling the 
#' use of statins (Setoguchi 2008).
#' 
#' The \emph{cont} argument can switch between Setoguchi binary simulation and the
#' continuous simulation of Lee et.al. (2011).
#' 
#' Setting \emph{plot.it}=TRUE will generate a graphic showing the continuous 
#' Y propensity score, with a horizontal line indicating the threshold cutoff 
#' value used to generate the binary simulation. Points above the line are set
#' to 1, those below are set to 0.
#' 
#' @references 
#' Setuguchi S, Schneeweiss S, Brookhart MA, Glynn RJ, Cook EF. 
#' Evaluation uses of data mining techniques in propensity score estimation: 
#' a simulation study. Pharmacoepidemiology and Drug Safety 2008; 9(4):403-425
#' 
#' Lee B, Lessler J, Stuart E, Improving propensity score weighting using machine learning.
#' Stat Med. 2010 February 10; 29(3):337-346
#' 
#' @export psSimulation
#' 
#' @example
#' # Lee used a continuous response of the y propensity score
#' lee <- psSimulation(n=500, sim="B", cont=TRUE)
#' 
#' # Setoguchi used the binary response based on the propensity score
#' # being lower than a point selected randomly from a uniform 0-1 distribution.
#' setoguchi <- psSimulation(n=2000, sim="B", cont=FALSE)
#' \dontrun{
#' # To generate the simulation of Lee, we need 1000 repititions of each scenario
#' # for datasets of size 500, 1000 and 2000.
#' 
#' # Set n = 500
#' n=500
#' 
#' # Create a list of simulation scenarios (A-G), each list element will contain
#' # a matrix of 1000 simulations (columns) of 500 elements (rows) each.
#' sim500 <- lapply(c("A","B", "C", "D", "E", "F", "G"),
#'    function(scn){
#'      sapply(1:1000, function(ind){psSimulation(n=n, sim=scn)})
#'    })
#'  } 
#'   
#'
psSimulation <- function(n=2000, sim = c("A","B", "C", "D", "E", "F", "G"), cont=TRUE, seed, plot.it=FALSE){
  nvar = 10
  if(!missing(seed)) set.seed(seed)
  ## Covariance matrix, if not provided
  cvar = diag(rep(1, nvar))
  cvar[1,5] = cvar[3,8] = cvar[5,1] = cvar[8,3] = .2
  cvar[2,6] = cvar[4,9] = cvar[6,2] = cvar[9,4] =.9
  if(missing(sim)) sim="A"
  #   cat(sim, "\n")
  sim <- match.arg(sim)
  
  # Beta coefficients (For generating the A parameter)
  bta <- c(.8, -.25, .6, -.4, -.8, -.5, .7)
  bta0 <- 0
  # alpha coefficients (For generating the outcome Y)
  alph <- c( .3, -.36, -.73, -.2, .71, -.19, .26)
  alph0 <- -3.85
  # Gamma coeff (For generating the outcome Y)
  gam <- -.4
  
  ## Create nvar standard normal covariates
  vTerm <- sapply(1:nvar, function(ind){rnorm(n)})
  
  wTerm <- vTerm %*% cvar
  
  ## Binary scale for i=c(1,3,5,6,8,9)
  wTerm[,c(1,3,5,6,8,9)] <- sapply(c(1,3,5,6,8,9), 
                                   function(ind){
                                     xVec <- wTerm[, ind]
                                     pt <- mean(wTerm[, ind])
                                     xVec[which(wTerm[,ind] > pt)] <- 1
                                     xVec[which(wTerm[,ind] <= pt)] <- 0
                                     xVec})
  trueProp <- switch(sim,
                     # Model with additivity and linearity
                     A = 1/(1+ exp(-(bta0 + wTerm[,1:7] %*% bta[1:7]))),
                     # Model with mild non-linearity
                     B = 1/(1+ exp(-(bta0 + wTerm[,1:7] %*% bta[1:7] + wTerm[,2]^2 * bta[2]))),
                     # Model with moderate non-linearity
                     C = 1/(1+ exp(-(bta0 + wTerm[,1:7] %*% bta[1:7]+ 
                                       wTerm[,2]^2 * bta[2]+ wTerm[,4]^2 * bta[4]+ wTerm[,7]^2 * bta[7] ))),
                     # Model with mild non-additivity
                     D = 1/(1+ exp(-(bta0 + wTerm[,1:7] %*% bta[1:7] + wTerm[,2]^2 * bta[2] + .5*bta[1]*wTerm[,1]*wTerm[,3] +
                                       .7*bta[2]*wTerm[,2]*wTerm[,4]+ .5*bta[4]*wTerm[,4]*wTerm[,5]+
                                       .5*bta[5]*wTerm[,5]*wTerm[,6]))),
                     # Model with mild non-additivity and non-linearity
                     E = 1/(1+ exp(-(bta0 + wTerm[,1:7] %*% bta[1:7]+ .5*bta[1]*wTerm[,1]*wTerm[,3] +
                                       .7*bta[2]*wTerm[,2]*wTerm[,4]+ .5*bta[4]*wTerm[,4]*wTerm[,5]+
                                       .5*bta[5]*wTerm[,5]*wTerm[,6]))),
                     # Model with moderate non-additivity
                     F = 1/(1+ exp(-(bta0 + wTerm[,1:7] %*% bta[1:7]+ .5*bta[1]*wTerm[,1]*wTerm[,3] +
                                       .7*bta[2]*wTerm[,2]*wTerm[,4]+ .5*bta[3]*wTerm[,3]*wTerm[,5]+
                                       .7*bta[4]*wTerm[,4]*wTerm[,6]+
                                       .5*bta[5]*wTerm[,5]*wTerm[,7]+ .5*bta[1]*wTerm[,1]*wTerm[,6]+
                                       .7*bta[2]*wTerm[,2]*wTerm[,3]+ .5*bta[3]*wTerm[,3]*wTerm[,4]+ 
                                       .5*bta[4]*wTerm[,4]*wTerm[,5]+ .5*bta[5]*wTerm[,5]*wTerm[,6]
                     ))),
                     # Model with moderate non-additivity and non-linearity
                     G = 1/(1+ exp(-(bta0 + wTerm[,1:7] %*% bta[1:7] + bta[2]*wTerm[,2]*wTerm[,2] +
                                       bta[4]*wTerm[,4]*wTerm[,4] +bta[7]*wTerm[,7]*wTerm[,7] +
                                       .5*bta[1]*wTerm[,1]*wTerm[,3] +
                                       .7*bta[2]*wTerm[,2]*wTerm[,4]+ .5*bta[3]*wTerm[,3]*wTerm[,5]+
                                       .7*bta[4]*wTerm[,4]*wTerm[,6]+
                                       .5*bta[5]*wTerm[,5]*wTerm[,7]+ .5*bta[1]*wTerm[,1]*wTerm[,6]+
                                       .7*bta[2]*wTerm[,2]*wTerm[,3]+ .5*bta[3]*wTerm[,3]*wTerm[,4]+ 
                                       .5*bta[4]*wTerm[,4]*wTerm[,5]+ .5*bta[5]*wTerm[,5]*wTerm[,6]
                     ))),
  )
  
  aTerm <- as.integer(runif(n=1) < trueProp)
  
  yProp <- 1/(1+exp(-(alph0 +wTerm[,c(1:4,8:10)] %*% alph + gam * aTerm)))
  
  yThresh <- runif(n=1)
  if(plot.it){
    plot(yProp, ylim=c(0,1), cex=.5)
    abline(h=yThresh, col=2, lty=2)
  }
  if(cont)
    yTrue <- yProp
  else 
    yTrue<-as.integer(yThresh < yProp)
  
  invisible(list(y= yTrue, w= wTerm, a= aTerm, ps= trueProp))
}
