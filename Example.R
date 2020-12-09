###### loading the package
library(MASS)
library(CVXR)
library(AER)
library(Matrix);
library(glmnet);
###### read the source code
source('~/Dropbox/Projects-Collaboration/Treatment Selection/JRSSB-Third-Revision/Submission-Code-Data/ITE_Linear.R', encoding = 'UTF-8')
###### set up the model generation parameters
n1 = 91
p = 171
n2 = 92
A1gen <- function(rho,p){
A1=matrix(0,p,p)
for(i in 1:p){
 for(j in 1:p){
   A1[i,j]<-rho^(abs(i-j))
}
}
return(A1)
}
mu <- rep(0,p)
mu[1:5] <- c(1:5)/5
rho = 0.5
Cov <- (A1gen(rho,p))/2
Cov2<- (A1gen(0.5,p))
beta1 <- rep(0,p)
beta1[1:10] <- c(1:10)/5
beta2 <- rep(0,p)
beta2[1:5] <- c(1:5)/10
####### generate the loading
loading <- MASS::mvrnorm(1,rep(0,p),Cov2)/2
for(i in 11:p){
  loading[i]<-loading[i]/10
}
true.ITE<-sum((beta1-beta2)*loading)
##### generate the data
X1 <- MASS::mvrnorm(n1,mu,Cov)
X2 <- MASS::mvrnorm(n2,mu,Cov)
y1 = X1%*%beta1 + rnorm(n1)
y2 = X2%*%beta2 + rnorm(n2)
######## implement our proposed algorithm
########## Here are the inputs
#### X1 denotes the design for the first sample
#### y1 denotes the outcome for the first sample 
#### X2 denotes the design for the second sample 
#### y2 denotes the outcome for the second sample
#### loading is xnew (a p-dimension vector denoting the covariate observation for the new observation)
Est <- ITE(X1 = X1, y1 = y1, X2 = X2, y2 = y2,loading = loading, intercept = TRUE)
######### Here are the outputs 
###### the point estimator of ITE
Est$prop.est
###### the standard error of ITE
Est$se
###### the 95% interval for ITE
Est$CI
##### the test of H_0: the first treatment is worse than the second treatment 
##### This is only a one-side test and hence the CI is more informative
Est$decision


