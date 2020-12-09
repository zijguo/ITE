# Calculate the mode
getmode <- function(v) {
  tbl <- table(v)
  if (all(tbl == 1)) {
    # case if all the values are distinct
    median(v)
  } else {
    # normal case if at least one value of v appears multiple times.
    as.numeric(names(which.max(tbl)))
  }
}


#####################################################################################################
################## plugin estimator

###### The following function computes the lasso estimator

# Compute the Lasso estimator:

# - If lambda is given, use glmnet and standard Lasso

# - If lambda is set to the character string "CV", then glmnet with

#   lambda selected by cross-validation is used

# - If lambda is not given or is set to NULL, use square root Lasso
Lasso <- function(X, y, lambda = NULL, intercept = TRUE) {
  p <- ncol(X)
  n <- nrow(X)
  
  htheta <- if (is.null(lambda)) {
    outLas <- cv.glmnet(X, y, family = "gaussian", alpha = 1,
                        intercept = intercept) 
    # Objective : 1/2 * RSS/n + lambda * penalty
    as.vector(coef(outLas, s = outLas$lambda.min))
  } else if (lambda == "CV") {
    outLas <- cv.glmnet(X, y, family = "gaussian", alpha = 1,
                        intercept = intercept)
    # Objective : 1/2 * RSS/n + lambda * penalty
    as.vector(coef(outLas, s = outLas$lambda.1se))
  } else if (lambda == "scalreg") {
    Xc <- if (intercept) {
      cbind(rep(1, n), X)
    } else {
      X
    }
    outLas <- scalreg(Xc, y)
    # return object
    if (intercept) {
      outLas$coefficients
    } else {
      # add a coefficient for the (not estimated) intercept b/c of implementation
      c(0, outLas$coefficients)
    }
  } else {
    outLas <- glmnet(X, y, family = "gaussian", alpha = 1,
                     intercept = intercept)
    # Objective : 1/2 * RSS/n + lambda * penalty
    as.vector(coef(outLas, s = lambda))
  }
  
  if (intercept == TRUE) {
    return(htheta)
  } else {
    return(htheta[2:(p+1)])
  }
}


Initialization.step <- function(X, y, lambda = NULL, intercept = FALSE) {
  n <- nrow(X)
  col.norm <- 1 / sqrt((1 / n) * diag(t(X) %*% X))
  Xnor <- X %*% diag(col.norm)
  
  ### Call Lasso
  htheta <- Lasso(Xnor, y, lambda = lambda, intercept = intercept)
  
  ### Calculate return quantities
  if (intercept == TRUE) {
    Xb <- cbind(rep(1, n), Xnor)
    col.norm <- c(1, col.norm)
  } else {
    Xb <- Xnor
  }
  htheta <- htheta * col.norm
  returnList <- list("lasso.est" = htheta)
  return(returnList)
}



Direction_fixedtuning_lin<-function(X,loading,mu=NULL){
  pp<-ncol(X)
  n<-nrow(X)
  if(is.null(mu)){
    mu<-sqrt(2.01*log(pp)/n)
  }
  loading.norm<-sqrt(sum(loading^2))
  
  if (loading.norm==0){
    H <- cbind(loading, diag(1, pp))
  }else{
    H <- cbind(loading / loading.norm, diag(1, pp))
  }
  v<-Variable(pp+1)
  obj<-1/4*sum((X%*%H%*%v)^2)/n+sum((loading/loading.norm)*(H%*%v))+mu*sum(abs(v))
  prob<-Problem(Minimize(obj))
  result<-solve(prob)
  print("fixed mu")
  print(mu)
  opt.sol<-result$getValue(v)
  cvxr_status<-result$status
  direction<-(-1)/2*(opt.sol[-1]+opt.sol[1]*loading/loading.norm)
  returnList <- list("proj"=direction)
  return(returnList)
}

Direction_searchtuning_lin<-function(X,loading,mu=NULL, resol = 1.5, maxiter = 10){
  pp<-ncol(X)
  n<-nrow(X)
  tryno = 1;
  opt.sol = rep(0,pp+1);
  lamstop = 0;
  cvxr_status = "optimal";
  
  mu = sqrt(2.01*log(pp)/n);
  #mu.initial= mu;
  while (lamstop == 0 && tryno < maxiter){
    ###### This iteration is to find a good tuning parameter
    #print(mu);
    lastv = opt.sol;
    lastresp = cvxr_status;
    loading.norm<-sqrt(sum(loading^2))
    
    if (loading.norm==0){
      H <- cbind(loading, diag(1, pp))
    }else{
      H <- cbind(loading / loading.norm, diag(1, pp))
    }
    v<-Variable(pp+1)
    obj<-1/4*sum((X%*%H%*%v)^2)/n+sum((loading/loading.norm)*(H%*%v))+mu*sum(abs(v))
    prob<-Problem(Minimize(obj))
    result<-solve(prob)
    cvxr_status<-result$status
    if(tryno==1){
      if(cvxr_status=="optimal"){
        incr = 0;
        mu=mu/resol;
        
        opt.sol<-result$getValue(v)
        
        temp.vec<-(-1)/2*(opt.sol[-1]+opt.sol[1]*loading/loading.norm)
        initial.sd<-sqrt(sum((X%*% temp.vec)^2)/(n)^2)*loading.norm
        temp.sd<-initial.sd
      }else{
        incr = 1;
        mu=mu*resol;
      }
    }else{
      if(incr == 1){ ### if the tuning parameter is increased in the last step
        if(cvxr_status=="optimal"){
          
          opt.sol<-result$getValue(v)
          
          lamstop = 1;
        }else{
          mu=mu*resol;
        }
      }else{
        if(cvxr_status=="optimal"&&temp.sd<3*initial.sd){
          mu = mu/resol;
          
          opt.sol <- result$getValue(v)
          
          temp.vec<-(-1)/2*(opt.sol[-1]+opt.sol[1]*loading/loading.norm)
          temp.sd<-sqrt(sum((X%*% temp.vec)^2)/(n)^2)*loading.norm
          #print(temp.sd)
        }else{
          mu=mu*resol;
          opt.sol=lastv;
          lamstop=1;
          tryno=tryno-1
        }
      }
    }
    tryno = tryno + 1;
  }
  direction<-(-1)/2*(opt.sol[-1]+opt.sol[1]*loading/loading.norm)
  step<-tryno-1
  print(step)
  returnList <- list("proj"=direction,
                     "step"=step)
  return(returnList)
}

LF<-function(X,y,loading,intercept=TRUE,init.Lasso=NULL,lambda=NULL,mu=NULL,step=NULL,resol = 1.5,maxiter=10){
  xnew<-loading
  p <- ncol(X);
  n <- nrow(X);
  n_y <- length(y)
  if(n_y!=n)
  {
    print("Check dimensions of X and y")
  }
  else
  {
    data = na.omit(data.frame(y,X))
    X <- as.matrix(data[,-1])
    y <- as.vector(data[,1])
    p <- ncol(X);
    n <- nrow(X);
    col.norm <- 1 / sqrt((1 / n) * diag(t(X) %*% X));
    Xnor <- X %*% diag(col.norm);
    if(is.null(init.Lasso)){
      ####### implement a lasso algorithm to get beta and sigma
      init.Lasso<-Initialization.step(X,y,lambda,intercept)
      htheta<-init.Lasso$lasso.est
    }
    else
    {
      htheta<- init.Lasso
    }
    if (intercept==TRUE){
      Xb <- cbind(rep(1,n),Xnor);
      Xc <- cbind(rep(1,n),X);
      pp <- (p+1);
    } else {
      Xb <- Xnor;
      Xc <- X;
      pp <- p
    }
    sparsity <- sum(abs(htheta) > 0.001)
    sd.est <- sqrt(sum((y - Xb %*% htheta)^2) / max(0.9*n, n - sparsity)) 
### compute the initial estimator
    if(intercept==TRUE){
      loading=rep(0,pp)
      loading[1]=1
      loading[-1]=xnew
    }else{
      loading=xnew
    }
    loading.norm<-sqrt(sum(loading^2))
    lasso.plugin<-sum(loading*htheta)
#####################################################################################################
################## Correction step
    if ((n>=6*p)){
      sigma.hat <- (1/n)*(t(Xc)%*%Xc);
      tmp <- eigen(sigma.hat)
      tmp <- min(tmp$values)/max(tmp$values)
    }else{
      tmp <- 0
    }
    sigma.hat <- (1/n)*(t(Xc)%*%Xc);
    if ((n>=6*p)&&(tmp>=1e-4)){
      direction <- solve(sigma.hat)%*%loading
    }else{
      if(is.null(step)){
        step.vec<-rep(NA,3)
        for(t in 1:3){
          index.sel<-sample(1:n,size=ceiling(0.5*min(n,p)), replace=FALSE)
          Direction.Est.temp<-Direction_searchtuning_lin(Xc[index.sel,],loading,mu=NULL, resol, maxiter)
          step.vec[t]<-Direction.Est.temp$step
        }
        step<-getmode(step.vec)
      }
      print(paste("step is", step))
      Direction.Est<-Direction_fixedtuning_lin(Xc,loading,mu=sqrt(2.01*log(pp)/n)*resol^{-(step-1)})
      while(is.na(Direction.Est)&&(step>0)){
        step<-step-1
        Direction.Est <- Direction_fixedtuning_lin(Xc, loading, mu = sqrt(2.01 * log(pp) / n) * resol^{-(step - 1)})
      }
      direction<-Direction.Est$proj
    }
    correction = t(Xc%*%direction)%*%(y - Xc%*%htheta)/n;
    debias.est=lasso.plugin+correction*loading.norm
    se<-sd.est*sqrt(sum((Xc%*%direction)^2)/(n)^2)*loading.norm
    returnList <- list("prop.est" = debias.est,
                       "se" = se,
                       "proj"=direction,
                       "plug.in"=lasso.plugin
    )
    return(returnList)
  }
}


ITE<-function(X1,y1,X2,y2,loading,intercept=TRUE,init.Lasso1=NULL,init.Lasso2=NULL,lambda1=NULL,lambda2=NULL,mu1=NULL,mu2=NULL,step1=NULL,step2=NULL,resol = 1.5,maxiter=10){
  Est1<-LF(X1,y1,loading,intercept=intercept,init.Lasso=init.Lasso1,lambda=lambda1,mu=mu1,step=step1,resol = 1.5,maxiter=10)
  Est2<-LF(X2,y2,loading,intercept=intercept,init.Lasso=init.Lasso2,lambda=lambda2,mu=mu2,step=step2,resol = 1.5,maxiter=10)
  diff.point<-Est1$prop.est - Est2$prop.est
  diff.sd<-sqrt((Est1$se)^2 + (Est2$se)^2)
  alpha=0.05
  diff.CI<-c(diff.point-qnorm(1-alpha/2)*diff.sd,diff.point+qnorm(1-alpha/2)*diff.sd)
  diff.decision<-(diff.point>qnorm(1-alpha)*diff.sd)
  returnList <- list("prop.est" = diff.point,
                     "se" = diff.sd,
                     "CI"=diff.CI,
                     "decision"=diff.decision
  )
  return(returnList)
}
