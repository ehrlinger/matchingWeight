###
#' balance.matched
#'
#' @description tests for balance.
#'
#' @param object a matched object from the match.pair or match.weight functions.
#' 
#' @details
#' 
#' @value 
#' 
#' @references 
#'
#' @export balanced.matched
#' 
#' @example
#' 
#' \dontrun{}
#'   
#'

balanced.matched <- function(object){
  
  n <- nrow(dat) 
  out.ps <- correct.ps.model( dat ) 
  fm.ps <- out.ps$fm 
  X.name <- c( "X0", names(fm.ps$coef)[-1] )  # X0 is a key word
  V.name <- setup$V.name 
  X.mat <- as.matrix( dat[ , X.name, drop=F ] ) 
  V.mat <- as.matrix( dat[ , V.name, drop=F ] ) 
  
  ps.hat <- out.ps$ps.hat 
  Z <- as.numeric( dat[,"Z"] )   # Z is a key word
  W <- pmin( ps.hat, 1-ps.hat )/( Z*ps.hat + (1-Z)*(1-ps.hat) ) 
  # matching weight
  
  tmp1 <- summary(fm.ps)$coef[ ,c(1,2,4)] 
  rownames(tmp1)[1] <- "X0" 
  tmp1 <- cbind( tmp1, CI.lower = tmp1[,"Estimate"] -
                   1.96*tmp1[,"Std. Error"] ) 
  tmp1 <- cbind( tmp1, CI.upper = tmp1[,"Estimate"] +
                   1.96*tmp1[,"Std. Error"] ) 
  ps.model <- tmp1 
  # this is the results of beta from PS model
  
  out1 <- NULL 
  n1 <- sum(Z==1)  n0 <- sum(Z==0) 
  for (j in 1:length(V.name)) {
    this.name <- V.name[j] 
    this.var <- as.numeric( V.mat[ , this.name ] ) 
    Z1.mean <- wtd.mean( x = this.var[Z==1], weights = rep(1, length=n1) ) 
    Z1.var <- wtd.var( x = this.var[Z==1], weights = rep(1, length=n1) ) 
    Z0.mean <- wtd.mean( x = this.var[Z==0], weights = rep(1, length=n0) ) 
    Z0.var <- wtd.var( x = this.var[Z==0], weights = rep(1, length=n0) ) 
    std.diff <- 100*abs(Z1.mean-Z0.mean)/sqrt((Z1.var + Z0.var)/2) 
    # NOTE: std.diff is in absoluate value, i.e., non-negative
    #       Therefore, its mean is strictly positive
    out1 <- rbind( out1, c( Z1.mean, Z1.var, Z0.mean, Z0.var, std.diff ) ) 
  }
  rownames(out1) <- V.name 
  colnames(out1) <- c("Z1.mean","Z1.var","Z0.mean","Z0.var","std.diff.pct") 
  std.diff.before <- out1 
  # standardized difference before applying MW
  
  out1 <- NULL 
  for (j in 1:length(V.name)) {
    this.name <- V.name[j] 
    this.var <- as.numeric( V.mat[ , this.name ] ) 
    Z1.mean <- wtd.mean( x = this.var[Z==1], weights = W[Z==1] ) 
    Z1.var <- wtd.var( x = this.var[Z==1], weights = W[Z==1] ) 
    Z0.mean <- wtd.mean( x = this.var[Z==0], weights = W[Z==0] ) 
    Z0.var <- wtd.var( x = this.var[Z==0], weights = W[Z==0] ) 
    std.diff <- 100*abs(Z1.mean-Z0.mean)/sqrt((Z1.var + Z0.var)/2) 
    # NOTE: std.diff is in absoluate value, i.e., non-negative
    #       Therefore, its mean is strictly positive
    out1 <- rbind( out1, c( Z1.mean, Z1.var, Z0.mean, Z0.var, std.diff ) ) 
  }
  rownames(out1) <- V.name 
  colnames(out1) <- c("Z1.mean","Z1.var","Z0.mean","Z0.var","std.diff.pct") 
  std.diff.after <- out1 
  # standardized difference after applying MW
  
  mu.B1.hat <- as.numeric( t(V.mat) %*% colVec(W*Z) )/sum(W*Z) 
  mu.B0.hat <- as.numeric( t(V.mat) %*% colVec(W*(1-Z)) )/sum(W*(1-Z)) 
  eta.B1.hat <- g.fun( mu.B1.hat, fun.type = setup$V.fun.type ) 
  eta.B0.hat <- g.fun( mu.B0.hat, fun.type = setup$V.fun.type ) 
  beta.hat <- as.numeric( coef(fm.ps) ) 
  theta.hat <- c( eta.B1.hat, eta.B0.hat, beta.hat ) 
  B.hat <- eta.B1.hat - eta.B0.hat 
  n.mu.B1 <- length( mu.B1.hat ) 
  n.mu.B0 <- length( mu.B0.hat )  # NOTE: n.mu.B1 = n.mu.B1 = length(B.hat)
  n.beta <- length( beta.hat ) 
  n.theta <- length(theta.hat) 
  D.mat <- cbind( diag(n.mu.B1), -diag(n.mu.B0),
                  matrix(0, nrow=n.mu.B1, ncol=n.beta) ) 
  # point estimators
  
  calc.ei <- function( X, beta ) {
    # Both X and beta are vectors
    tmp1 <- exp(sum(X*beta)) 
    tmp1/(1+tmp1) 
  }
  
  calc.ei.deriv <- function( X, beta ) {
    # Both X and beta are vectors
    tmp1 <- exp(sum(X*beta)) 
    tmp2 <- tmp1/(1+tmp1) 
    as.numeric( tmp2*(1-tmp2)*X ) 
  }
  
  calc.W.Ze <- function( Z, e ) {
    # Both Z and e are scalar
    if ( e < 0.5 - setup$delta | e > 0.5 + setup$delta ) {
      ans <- min(e, 1-e) 
    } else {
      a <- solve.a() 
      ans <- a[1] + a[2]*e + a[3]*e*e + a[4]*e*e*e 
    }
    ans/( Z*e + (1-Z)*(1-e) ) 
  }
  
  calc.W.deriv.Ze <- function( Z, e ) {
    # Both Z and e are scalar
    tmp0 <- Z*e + (1-Z)*(1-e) 
    if ( e < 0.5 - setup$delta ) {
      ans <- (1-Z)/(tmp0^2) 
    } else if ( e > 0.5 + setup$delta ) {
      ans <- -Z/(tmp0^2) 
    } else {
      a <- solve.a() 
      tmp1 <- a[2] + 2*a[3]*e + 3*a[4]*e*e 
      tmp2 <- a[1] + a[2]*e + a[3]*e*e + a[4]*e*e*e 
      ans <- ( tmp1*tmp0 - (2*Z-1)*tmp2 )/(tmp0^2) 
    }
    ans 
  }
  
  solve.a <- function() {
    delta <- setup$delta 
    tmp1 <- c( 0.5-delta, 1, 0.5-delta, -1 ) 
    tmp2 <- matrix( c( 1, 0.5-delta, (0.5-delta)**2, (0.5-delta)**3 ,
                       0, 1, 2*(0.5-delta), 3*(0.5-delta)**2 ,
                       1, 0.5+delta, (0.5+delta)**2, (0.5+delta)**3 ,
                       0, 1, 2*(0.5+delta), 3*(0.5+delta)**2 ),
                    nrow=4, ncol=4, byrow=T ) 
    a.vec <- as.numeric( solve(tmp2) %*% colVec(tmp1) ) 
    a.vec 
  }
  
  Meat.mat <- Bread.mat <- 0  # the meat and bread of the sandwich variance
  for (i in 1:n) {
    this.V <- as.numeric( V.mat[i,] ) 
    this.X <- as.numeric( X.mat[i,] ) 
    this.Z <- Z[i] 
    
    this.ei <- calc.ei( this.X, beta.hat ) 
    this.ei.deriv <- calc.ei.deriv( this.X, beta.hat ) 
    this.Wi <- calc.W.Ze( this.Z, this.ei ) 
    this.Wi.deriv <- calc.W.deriv.Ze( this.Z, this.ei )*this.ei.deriv 
    
    tmp1 <- this.Wi*this.Z*( this.V - g.inv(eta.B1.hat, setup$V.fun.type) ) 
    tmp2 <- this.Wi*(1-this.Z)*( this.V - g.inv(eta.B0.hat, setup$V.fun.type) ) 
    tmp3 <- (this.Z-this.ei)*this.X 
    this.phi <- c(tmp1, tmp2, tmp3) 
    Meat.mat <- Meat.mat + outer( this.phi, this.phi ) 
    
    Block.11 <- -this.Wi*this.Z*diag( g.inv.deriv(eta.B1.hat, setup$V.fun.type) ) 
    Block.12 <- matrix(0, nrow=n.mu.B1, ncol=n.mu.B0) 
    Block.13 <- this.Z*( colVec( this.V - g.inv(eta.B1.hat, setup$V.fun.type) )
                         %*% rowVec( this.Wi.deriv ) )
    Block.21 <- matrix(0, nrow=n.mu.B0, ncol=n.mu.B1) 
    Block.22 <- -this.Wi*(1-this.Z)*diag( g.inv.deriv(eta.B0.hat, setup$V.fun.type) ) 
    Block.23 <- (1-this.Z)*( colVec( this.V - g.inv(eta.B0.hat, setup$V.fun.type) )
                             %*% rowVec( this.Wi.deriv ) )
    Block.31 <- matrix(0, nrow=n.beta, ncol=n.mu.B1) 
    Block.32 <- matrix(0, nrow=n.beta, ncol=n.mu.B0) 
    Block.33 <- -this.ei*(1-this.ei)*outer(this.X, this.X) 
    Bread.mat <- Bread.mat + rbind( cbind( Block.11, Block.12, Block.13 ) ,
                                    cbind( Block.21, Block.22, Block.23 ) ,
                                    cbind( Block.31, Block.32, Block.33 ) ) 
  }
  B.mat <- Meat.mat/n 
  A.mat <- Bread.mat/n 
  A.mat.inv <- solve(A.mat) 
  Sigma.theta.hat <- ( A.mat.inv %*% B.mat %*% t(A.mat.inv) )/n 
  Sigma.B.hat <- D.mat %*% Sigma.theta.hat %*% t(D.mat) 
  # Sandwich variance estimator
  
  d <- qr(Sigma.B.hat)$rank 
  test.stat <- as.numeric( t(B.hat) %*% solve(Sigma.B.hat) %*% B.hat ) 
  pvalue <- pchisq( test.stat, df = d, lower.tail = FALSE ) 
  # chi-squared test for overall comparison
  
  var.beta.hat <-
    Sigma.theta.hat[(n.mu.B1 + n.mu.B0 + 1):(n.mu.B1 + n.mu.B0 + n.beta),
                    (n.mu.B1 + n.mu.B0 + 1):(n.mu.B1 + n.mu.B0 + n.beta)] 
  sd.beta.hat <- sqrt(diag(var.beta.hat)) 
  beta.hat.lower <- beta.hat - 1.96*sd.beta.hat 
  beta.hat.upper <- beta.hat + 1.96*sd.beta.hat 
  # variance and CI of beta.hat
  
  tmp1 <- sqrt(diag(Sigma.B.hat)) 
  tmp2 <- B.hat/tmp1 
  tmp3 <- 2*pnorm( -abs(tmp2) ) 
  V.test <- cbind( Estimate = B.hat , SD = tmp1 , Z = tmp2 ,
                   pvalue = tmp3 , reject = as.numeric(tmp3 <= 0.05) ) 
  rownames(V.test) <- V.name 
  V.test.reject <- as.numeric( any( V.test[,"reject"] == 1 ) ) 
  # individual test of each of V.name, after g() transformation
  # if any individual Z test rejects, the overall V.test rejects
  
  names(mu.B1.hat) <- paste("mu.B1.", V.name, sep="") 
  names(mu.B0.hat) <- paste("mu.B0.", V.name, sep="") 
  names(eta.B1.hat) <- paste("eta.B1.", V.name, sep="") 
  names(eta.B0.hat) <- paste("eta.B0.", V.name, sep="") 
  names(beta.hat) <- X.name 
  names(theta.hat) <- c( names(eta.B1.hat), names(eta.B0.hat), names(beta.hat) ) 
  names(B.hat) <- c( paste("B.hat.", V.name, sep="") ) 
  rownames(Sigma.theta.hat) <- colnames(Sigma.theta.hat) <- names(theta.hat) 
  rownames(Sigma.B.hat) <- colnames(Sigma.B.hat) <- names(B.hat) 
  rownames(var.beta.hat) <- colnames(var.beta.hat) <- names(beta.hat) 
  names(sd.beta.hat) <- names(beta.hat.lower) <-
    names(beta.hat.upper) <- names(beta.hat) 
  # Add names to point and variance estimators
  
  return( list( X.name = X.name ,
                V.name = V.name ,
                n1 = sum(Z==1) ,
                n0 = sum(Z==0) ,
                eff.n1 = sum( W[Z==1] ) ,
                eff.n0 = sum( W[Z==0] ) ,
                
                std.diff.before = std.diff.before ,
                std.diff.after = std.diff.after ,
                mu.B1.hat = mu.B1.hat ,
                mu.B0.hat = mu.B0.hat ,
                eta.B1.hat = eta.B1.hat ,
                eta.B0.hat = eta.B0.hat ,
                
                beta.hat = beta.hat ,
                var.beta.hat = var.beta.hat ,
                sd.beta.hat = sd.beta.hat ,
                beta.hat.lower = beta.hat.lower ,
                beta.hat.upper = beta.hat.upper ,
                ps.model = ps.model ,
                
                theta.hat = theta.hat ,
                Sigma.theta.hat = Sigma.theta.hat ,
                
                B.hat = B.hat ,
                Sigma.B.hat = Sigma.B.hat ,
                
                V.test = V.test ,
                V.test.reject = V.test.reject ,
                
                test.stat = test.stat,
                d = d,
                pvalue = pvalue,
                reject = as.numeric(pvalue <= 0.05)
  ) ) 
}