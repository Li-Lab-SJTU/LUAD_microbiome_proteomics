## The function of calculating two-part statistics is referred to the paper of Taylor & Pollard (2009). ##
## Refer to: DOI: 10.2202/1544-6115.1425 ##
TwoPartTest <- function(data, group, test="t.test", point.mass=0){
  Index1 <- c(group==1)
  Group1 <- data[Index1]
  Group0 <- data[!Index1]
  n1 <- length(Group1)
  n2 <- length(Group0)
  obs <- c(n1, n2)
  success <- c(sum(Group1!=point.mass), sum(Group0!=point.mass))
  pointmass <- obs-success
  if (sum(success)==0) {
    T2 <- 0
    B2 <- 0
  }
  else if ((success[1]==0)|(success[2]==0)) {
    T2 <- 0
    B2 <- prop.test(pointmass, obs)$statistic
  }
  else if ((success[1]==1)|(success[2]==1)){
    T2 <- 0
    B2 <- prop.test(pointmass, obs)$statistic
  }
  else {
    uniq1 <- length(unique(Group1[Group1!=point.mass]))
    uniq2 <- length(unique(Group0[Group0!=point.mass]))
    if ((uniq1 < 2) & (uniq2 < 2)){
      T2 <- 0
      if (sum(pointmass)==0)
        B2 <- 0
      else
        B2 <- prop.test(pointmass, obs)$statistic
    }
    else if (sum(pointmass)==0){
      B2 <- 0
      if (test=="t.test")
        T2 <- t.test(data~group)$statistic^2
      if (test=="wilcoxon") {
        W <- wilcox.test(data~group, exact=FALSE)$statistic
        mu <- (n1*n2)/2
        sigma <- sqrt((n1*n2*(n1+n2+1))/12)
        T2 <- ((abs(W-mu)-0.5)/sigma)^2
      }
    }
    else {
      B2 <- prop.test(pointmass, obs)$statistic
      contIndex <- data!=point.mass
      cont <- data[contIndex]
      cGroup <- group[contIndex]
      n1c <- sum(cGroup==1)
      n2c <- sum(cGroup==0)
      if (test=="t.test")
        T2 <- t.test(cont~cGroup)$statistic^2
      if (test=="wilcoxon") {
        W <- wilcox.test(cont~cGroup, exact=FALSE)$statistic
        mu <- (n1c*n2c)/2
        sigma <- sqrt((n1c*n2c*(n1c+n2c+1))/12)
        T2 <- ((abs(W-mu)-0.5)/sigma)^2
      }
    }
  }
  X2 <- B2+T2
  if ((T2==0)|(B2==0)) {
    X2pv <- 1-pchisq(X2,1)
  } else {
    X2pv <- 1-pchisq(X2,2)
  }
  ans <- list(statistic=as.numeric(X2), pvalue=as.numeric(X2pv))
  return(ans)
}