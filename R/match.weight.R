###
#' match.weight
#'
#' @description matching treatment and control patients by propensity scores matching weights
#'
#' @param treatment indicator variable for treatment/control (1/0)
#' @param propensity a precalculated propensity score on (0,1)
#' 
#' @details
#' 
#' @references 
#'
#' @export match.wieght
#' 
#' @example
#' 
#' \dontrun{}
#'   
#'
match.weight <- function(treatment, 
                         pscore){
  
  # Sanity check the pscore
  if(sum(pscore >=1) + sum(pscore <=0) > 0)
    stop("Mispecification of the propensity score. Values are strictly on (0,1).")
  
  # Sanity check the treatment variable
  if(is.logical(treatment)) treatment = as.numeric(treatment)-1
  if(sum(treatment >1) + sum(treatment <0) > 0)
    stop("Mispecification of the variable. Values are strictly {0,1}.")
  
  weight = min(1-pscore, pscore)/(treatment *pscore +(1-treatment)*(1-pscore))
  
  
}