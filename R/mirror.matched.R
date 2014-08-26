###
#' mirror.matched
#' create a mirrored histogram of propensity score matched and unmatched cases.
#' 
#' @description The propensity score mirrored histogram can be used to verify
#' the matching strategy is balanced across the spectrum of the propensity score
#' (along the x-axis). Colored bins indicate matched sets, uncolored indicates
#' unmatched cases. If the matching set is omitted, or doesn't exist, a mirrored histogram of unmatched control (above zero) vs
#' treatement (below zero) groups is returned.
#' 
#' @param data a data frame consisting of at least a propensity score and treatment column
#' @param ps.name the column name of the propensity score ("propensity")
#' @param treat.name the column name of the treatment effect, coercible into a logical variable ("treatment")
#' @param match.name the column name of the matching information, coercible into a logical variable ("matched")
#' @param colors an optional vector of 2 or 3 colors. By default the matched cases are in green 
#' (light for control, dark for treatment) unmatched are colored white.
#' @param breaks Histogram bin widths associated with the propensity score matching bins. 
#'
#' @return A ggplot histogram object
#' 
#' @examples
#' # create a random set of propensity scores evenly distributed
#' prop <- runif(1000)
#' 
#' # Within each 10% randomly assign some patients to treatment
#' treat<- do.call(c,sapply(1:10, function(ind){rbinom(n=sum(prop>(ind-1)/10 & prop <= ind/10 ), size=1,prob=ind/10)}))
#' match <- pair.matching()
#' 
#' ps.data <- as.data.frame(cbind(prop, treat, match))
#' 
#' # Create and show the graphic
#' mirHist <-mirror.histogram(sta_prob, ps.name="prob", treat.name="statin", match.name="match")
#' show(mirHist)
#' 
#' # Modify the graphics with standard ggplot2 modifiers
#' mirHist+scale_x_continuous(breaks=seq(0,100,20)) +
#' scale_y_continuous(breaks = seq(-25,75, 25), labels=c("25","0","25","50","75")) + theme_bw()
#' 
#' @export mirror.matched
#' @import
#'
#' 
mirror.matched <- function(object, data=NULL, ps.name ="propensity", treat.name="treatment", match.name="matched", colors, breaks){
  
  ## Create vectors of mirrored variables.
  prop.score <- dta[, which(colnames(dta) ==ps.name)]
  treat <- as.logical(dta[, which(colnames(dta) == treat.name)])
  match <- as.logical(dta[, which(colnames(dta) == match.name)])
  
  ## Verify that we have the correct colnames by checking that the lengths of each vector match
  if(length(prop.score) != length(treat) | length(prop.score) != length(match)|
       length(match) != length(treat)){
    stop(paste("Named columns are incorrectly specified. Data colnames,", colnames(dta), 
               "does not include", ps.name, treat.name, "or", "match.name"))
  }
  
  n.obs <- length(prop.score)
  if(n.obs < 1)
    stop("Empty propensity score vector.")
  
  ## Formating the variables
  n.obs <- length(prop.score)
  rng <- range(prop.score)
  if(max(rng) <= 1){
    prop.score <- prop.score *100
    rng <- range(prop.score)
  } 
  
  
  # Create the elements of the histogram
  # control group
  t0 <- as.data.frame(cbind(prop.score[which(!treat)], match[which(!treat)]))
  # treatment group
  t1 <- as.data.frame(cbind(prop.score[which(treat)], match[which(treat)]))
  colnames(t0)<-c("ps", "match")
  colnames(t1)<-c("ps", "match")
  
  
  # Match controls
  m0 <- as.data.frame(t0[which(t0$match==1),])
  
  # matched treatments
  m1 <- as.data.frame(t1[which(t1$match==1),])
  
  
  # Create the histogram breaks dependant on the number of observations.
  if(missing(breaks)) breaks <- seq(0,100,20)
  if(missing(colors)) colors = c("green1", "green4", "white")
  
  mHist <- ggplot() + geom_histogram(aes(x=ps), fill=colors[3], col="black",breaks=breaks, data=as.data.frame(t0))+
    geom_histogram(aes(x=ps), fill=colors[1], col="black",breaks=breaks, data=as.data.frame(m0))+
    geom_histogram(aes(x=ps, y=-..count..), fill=colors[3], col="black",breaks=breaks, data=as.data.frame(t1))+
    geom_histogram(aes(x=ps, y=-..count..), fill=colors[2], col="black",breaks=breaks, data=as.data.frame(m1)) +
    labs(x="Propensity Score", y="Count")
  
  return(mHist)
}