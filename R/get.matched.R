#'  the core function of pair matching
#'
#' @export get.matched
#' 
require(matchit)

get.matched <- function(dat, correct.model = TRUE ) {
  if ( correct.model ) {
    out.ps <- correct.ps.model( dat ) ;
  } else {
    out.ps <- wrong.ps.model( dat ) ; # not used here
  }
  out <- matchit( out.ps$fm$formula, data = dat, method="nearest",
                  caliper=setup$caliper, distance="linear.logit" ) ;
  ## the following code is to get matched data and identify matched pairs
  dat$weights <- out$weights ;
  dat$distance <- out$distance ;
  dat$rownm <- rownames(dat) ;
  match.matrix <- out$match.matrix[ !is.na(out$match.matrix[,1]) &
                                      out$match.matrix[,1] != "-1", 1, drop=FALSE] ;
  trt.id <- as.character( rownames(match.matrix) ) ;
  matched.to <- as.character( match.matrix[,1] )  ;
  pair <- as.character( as.integer(1:length(trt.id)) ) ;
  map <- rbind( cbind( matched.to, pair ), cbind( trt.id, pair ) ) ;
  colnames(map) <- c("rownm","pair") ;
  tmp <- order(map[,"pair"]) ;
  map <- map[tmp, ] ;
  # get the mapping
  mdat <- dat[ dat$rownm %in% c(trt.id, matched.to) , ] ;
  # get matched data set
  mdat2 <- merge(mdat, map, by="rownm", all=TRUE) ;
  tmp <- order(mdat2$pair, mdat2$Z) ;
  mdat2 <- mdat2[tmp, ] ;
  mdat2$ps.lin <- mdat2$distance ;
  after <- mdat2 ;
  after$ps.hat <- inv.logit( after$ps.lin ) ;
  after <- after[ ,c("id","pair")] ;
  # the data set after matching
  before <- dat ; before$ps.lin <- before$distance ;
  before$ps.hat <- inv.logit( before$ps.lin ) ;
  # the data set before matching
  newdat <- merge( x = before, y = after, by = "id", all.x = T, all.y = F ) ;
  newdat$pair <- as.numeric( newdat$pair ) ;
  newdat[ is.na(newdat$pair), "pair" ] <- 0 ;
  newdat$matched <- as.numeric( newdat$pair != 0 ) ;
  newdat <- newdat[ , c("id","Z","Y","X0",
                        unique( c(setup$X.name, setup$V.name,
                                  setup$X4PS.name, setup$X4Y.name) ),
                        "ps.hat","matched","pair")] ;
  # the matching result is in "matched" and "pair"
  # *** NOTE: covariate list is added here !!! ***
  return( newdat ) ;
}
