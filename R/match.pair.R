###
#' match.pair
#'
#' @description matching treatment and control patients pairs based on propensity scores
#'
#' 
#' @param data data.frame containing cases and potential controls. Must 
#' contain the id, group, and match variables.
#' 
#' @param treatment variable name defining cases. treatment= 1 or TRUE if case, 0 or FALSE if control.
#' 
#' @param id unique patient identifier variable name.
#' 
#' @param match matching variable names common to both case and control
#' 
#' @param wts  List of non-negative weights corresponding to each matching variable.  
#' For example wts=10 2 1 corresponding to male, age and birthyr as in the above example.
#' 
#' @param dmaxk  List of non-negative values corresponding to each matching variable.  These numbers are the 
#' largest possible absolute differences compatible with a valid match.  Cases will NOT be matched to a control 
#' if ANY of the INDIVIDUAL matching factor  differences are >DMAXK.  This optional parameter allows one to form 
#' matches of the type male+/-0, age+/-2, birth year+/-5 by specifying DMAXK=0 2 5.
#' 
#' @param dmax   Largest value of Dij considered to be a valid match.  If you want to match exactly on a two-level 
#' factor(such as gender coded as 0 or 1) then assign DMAX to be less than the weight for the factor.  In the 
#' example above, one could use wt=10 for male and dmax=9.  Leave DMAX blank if any Dij is a valid match.  One 
#' would typically NOT use both DMAXK and DMAX.  The only advantage to using both, would be to further restrict
#'  potential matches that meet the DMAXK criteria.
#'  
#' @param dist  Indicates type of distance to calculate.  1=weighted sum(over matching vars) of  absolute 
#' case-control differences(default) 2=weighted Euclidean distance
#' 
#' @param time  Time variable used for risk set matching.  Matches are only valid if the control time > case time. 
#' May need to
#' 
#' @param transf  Indicates whether all matching vars are to be transformed (using the combined case+control data) 
#' prior to computing distances.  0=no(default), 1=standardize to mean 0 and variance 1, 2=use ranks of matching
#'  variables.
#'  
#' @param ncontls  Indicates the number of controls to match to each case.  The default is 1.  With multiple controls
#' per case, the algorithm will first match every case to one control and then again match each case to a second 
#' control, etc.  Controls selected on the first pass will be stronger matches than those selected in later rounds.
#' The output data set contains a variable (cont_n) which indicates on which round the control was selected.   
#' 
#' @param seedca   Seed value used to randomly sort the cases prior to matching. This positive integer will be used 
#' as input to the RANUNI function.  The greedy matching algorithm is order dependent which, among other things means 
#' that cases matched first will be on average more similar to their controls than those matched last(as the number 
#' of control choices will be limited).  If the matching order is related to confounding factors (possibly age or
#' calendar time) then biases may result.  Therefore it is generally considered good practice when using the GREEDY 
#' method to randomly sort both the cases and controls before beginning the matching process.
#' 
#' @param seedco   Seed value used to randomly sort the controls prior to matching using the GREEDY method.  This
#' seed value must also be a positive integer.
#' 
#' @param print= Option to print data for matched cases. Use PRINT=y to print data and PRINT=n or blank to not 
#' print.  Default is y.
#' 
#' @param out=name of SAS data set containing the results of the matching process.  Unmatched cases are not included.
#'   See outnm below.  The default name is __out.  This data set will have the following layout:
#'
#'  Case_id  Cont_id  Cont_n  Dij  Delta_caco MVARS_ca  MVARS_co
#'     1        67      1     5.2  (Differences & actual
#'     1        78      2     6.1   values for matching factors
#'     2        52      1     2.9   for cases & controls)
#'     2        92      2     3.1
#'     .        .       .      .
#'     .        .       .      .
#'
#' @param outnmca=name of SAS data set containing NON-matched cases. Default name is __nmca .
#' 
#' @param outnmco=name of SAS data set containing NON-matched controls. Default name is __nmco .
#'
#' @details
#' 
#' 
#' @references Bergstralh, EJ and Kosanke JL(1995).  Computerized matching of controls.  
#' Section of Biostatistics Technical Report 56.  Mayo Foundation.
#'
#' @export match.pair
#' 
#' @example
#' 
#' \dontrun{}
#'   
#'
match.pairs <- function(data,
                        treatment=,
                        id=,
                        match=,
                        wts=,
                        dmaxk=,
                        dmax=,
                        transf,
                        time=, 
                        dist=,
                        ncontls=,
                        seedca=,
                        seedco=,
                        out=,
                        outnmca=,
                        outnmco=,
                        print=){
  
  
  
}