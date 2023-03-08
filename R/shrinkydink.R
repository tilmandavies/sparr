

cv.RelRisk <- function(X1,X2=NULL,h=NA,lambda.min=0.1,lambda.max=20,LHmin=TRUE,n.lambda=32){
  if (is.null(X2)){
    levels(X1$marks) <- c("case","control")
    X <- X1
    X2 <- split(X1)$control
    X1 <- split(X1)$case
  } else{
    n1 <- npoints(X1)
    n2 <- npoints(X2)
    X1 <- unmark(X1)
    X2 <- unmark(X2)
    X <- superimpose(X1,X2,check=F)
    marks(X) <- factor(c(rep("case",n1),rep("control",n2)))
  }
  if (is.na(h[1])) h <- sqrt(bw.scott.iso(X1)*bw.scott.iso(X2))
  if (length(h)==1) h <- c(h,h)
  n1 <- npoints(X1)
  n2 <- npoints(X2)
  y <- -(as.numeric(marks(X))-2)
  lambda <- exp(seq(from = log(lambda.min), to = log(lambda.max), length.out = n.lambda))
  f1 <- density.ppp(X,at="points",weights=y,sigma=h[1],positive=T,leaveoneout=F) - y/(2*pi*h[1]^2)
  f2 <- density.ppp(X,at="points",weights=1-y,sigma=h[2],positive=T,leaveoneout=F) - (1-y)/(2*pi*h[2]^2)
  beta.sign <- sign(log(f1)-log(f2)-log(n1/n2))
  CV <- numeric(n.lambda)
  for (i in 1:n.lambda){
    case1 <- (f1/n1-lambda[i]/(n1*2*pi*h[1]^2))/(f2/n2+lambda[i]/(n2*2*pi*h[2]^2))
    case2 <- (f1/n1+lambda[i]/(n1*2*pi*h[1]^2))/(f2/n2-lambda[i]/(n2*2*pi*h[2]^2))
    logRR <- case1*0
    logRR[case1 > 1] <- log(case1[case1 > 1])
    logRR[1/case2 > 1] <- log(case2[1/case2 > 1])
    p <- 1/(n2*exp(-c(logRR))/n1+1)
    CV[i] <- -sum(log(p)*y+log(1-p)*(1-y))
  }
  if(LHmin){
    CVdiff <- diff(CV)
    if(all(CVdiff>=0)) lambda.opt <- lambda[1]
    if(all(CVdiff<=0)) lambda.opt <- lambda[n.lambda]
    if(!all(CVdiff<=0) & !all(CVdiff>=0)) lambda.opt <- lambda[min(which(CVdiff > 0))]
  } else{
    lambda.opt <- lambda[which.min(CV)]
  }
  return(list(lambda=lambda.opt,lambda.grid=lambda,CV=CV))
}


cv.RelRisk.Bithell <- function(X1,X2=NULL,h=NA,lambda.min=0.1,lambda.max=20,LHmin=TRUE,n.lambda=32){
  if (is.null(X2)){
    levels(X1$marks) <- c("case","control")
    X <- X1
    X2 <- split(X1)$control
    X1 <- split(X1)$case
  } else{
    n1 <- npoints(X1)
    n2 <- npoints(X2)
    X1 <- unmark(X1)
    X2 <- unmark(X2)
    X <- superimpose(X1,X2,check=F)
    marks(X) <- factor(c(rep("case",n1),rep("control",n2)))
  }
  if (is.na(h[1])) h <- sqrt(bw.scott.iso(X1)*bw.scott.iso(X2))
  if (length(h)==1) h <- c(h,h)
  n1 <- npoints(X1)
  n2 <- npoints(X2)
  y <- -(as.numeric(marks(X))-2)
  lambda <- exp(seq(from = log(lambda.min), to = log(lambda.max), length.out = n.lambda))
  f1 <- density.ppp(X,at="points",weights=y,sigma=h[1],positive=T,leaveoneout=F) - y/(2*pi*h[1]^2)
  f2 <- density.ppp(X,at="points",weights=1-y,sigma=h[2],positive=T,leaveoneout=F) - (1-y)/(2*pi*h[2]^2)
  CV <- numeric(n.lambda)
  for (i in 1:n.lambda){
    logRR <- log(f1/n1+lambda[i]/(n1*2*pi*h[1]^2)) - log(f2/n2 + lambda[i]/(n1*2*pi*h[1]^2))
    p <- 1/(n2*exp(-c(logRR))/n1+1)
    CV[i] <- -sum(log(p)*y+log(1-p)*(1-y))
  }
  if(LHmin){
    CVdiff <- diff(CV)
    if(all(CVdiff>=0)) lambda.opt <- lambda[1]
    if(all(CVdiff<=0)) lambda.opt <- lambda[n.lambda]
    if(!all(CVdiff<=0) & !all(CVdiff>=0)) lambda.opt <- lambda[min(which(CVdiff > 0))]
  } else{
    lambda.opt <- lambda[which.min(CV)]
  }
  return(list(lambda=lambda.opt,lambda.grid=lambda,CV=CV))
}




