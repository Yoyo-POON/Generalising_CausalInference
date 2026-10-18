generate_weights<-function(Smod,Amod, data){
  S1data<-subset(data, S==1)
  w_reg<-glm(Smod, family="binomial", data=data)
  ps<- predict(w_reg,newdata=data, type="response")
  w_reg2<-glm(Amod, family="binomial", data=S1data)
  pa<- predict(w_reg2,newdata=data, type="response")
  w = (data$A*data$S)/(ps*pa) + ((1-data$A) *data$S) /(ps*(1-pa))
  data$w<-w
  list<-list(w=w,dat=data, Smod=w_reg, Amod=w_reg2, ps=ps, pa=pa)
  return(list)
}



OM_est<-function(data, outcome = "c"){
  
  S1data_A1<-subset(data, S==1 & A==1)
  
  if ( outcome == "c"){
    OR1mod<-glm(formula=Y~X1+X2+X3,family="gaussian", data=S1data_A1)
  } else {
    OR1mod<-glm(formula=Y~X1+X2+X3,family="binomial", data=S1data_A1)
  }
  
  p1<- predict(OR1mod,newdata=data, type="response")
  
  S1data_A0<-subset(data, S==1 & A==0)
  
  if ( outcome == "c"){
    OR0mod<-glm(formula=Y~X1+X2+X3,family="gaussian", data=S1data_A0)
  } else {
    OR0mod<-glm(formula=Y~X1+X2+X3,family="binomial", data=S1data_A0)
  }
  
  p0<- predict(OR0mod,newdata=data, type="response")
  diff<-mean(p1)-mean(p0)
  list<-list(mu1=mean(p1), mu0=mean(p0), diff= diff,
             p1=p1, p0=p0)
  return(list)
}

IPW1_est<-function(data){
  A<-data$A
  S<-data$S
  w<-data$w
  Y<-data$Y
  mu1<- mean(A*S*w*Y)
  mu0<-mean((1-A)*S*w*Y)
  diff=mu1-mu0
  estimates<-data.frame(estimator="IPW1",mu1=mu1, mu0=mu0, diff=diff)
  return(estimates)
}

IPW2_est<-function(data){
  S1data<-subset(data, S==1)
  A<-data$A
  S<-data$S
  Y<-data$Y
  w<-data$w
  sum1_IPW<-sum(S*A*w*Y)
  sum0_IPW<-sum(S*(1-A)*w*Y)
  norm1<-(sum(S*A*w))^-1
  norm0<-(sum(S*(1-A)*w))^-1
  mu1<-norm1*sum1_IPW
  mu0<-norm0*sum0_IPW
  diff<-mu1-mu0
  estimates<-data.frame(estimator="IPW2",mu1=mu1, mu0=mu0, diff=diff)
  return(estimates)
}

AIPW1_est<-function(data){
  A<-data$A
  S<-data$S
  Y<-data$Y
  p1<-data$p1
  p0<-data$p0
  w<-data$w
  mu1<-mean(S*A*w*(Y-p1) + p1)
  mu0<-mean(S*(1-A)*w*(Y-p0) + p0)
  diff <-mu1- mu0
  estimates<-data.frame(estimator="AIPW1",mu1=mu1, mu0=mu0, diff=diff)
  return(estimates)
}

AIPW2_est<-function(data){
  A<-data$A
  S<-data$S
  Y<-data$Y
  p1<-data$p1
  p0<-data$p0
  w<-data$w
  sum1_AIPW2<-sum(S*A*w*(Y-p1))
  sum0_AIPW2<-sum(S*(1-A)*w*(Y-p0))
  norm1<-(sum(S*A*w))^-1
  norm0<-(sum(S*(1-A)*w))^-1
  mu1<-norm1*sum1_AIPW2 + mean(p1)
  mu0<-norm0*sum0_AIPW2 + mean(p0)
  diff<-mu1- mu0
  estimates<-data.frame(estimator="AIPW2",mu1=mu1, mu0=mu0, diff=diff)
  return(estimates)
}


AIPW3_est<-function(data, outcome){
  S1data_A1<-subset(data, S==1 & A==1)
  
  if ( outcome == "c"){
    OR1mod<-glm(formula=Y~X1+X2+X3,family="gaussian", data=S1data_A1, weights=w)
  } else {
    OR1mod<-glm(formula=Y~X1+X2+X3,family="binomial", data=S1data_A1, weights=w)
  }
  
  p1<- predict(OR1mod,newdata=data, type="response")
  
  
  S1data_A0<-subset(data, S==1 & A==0)
  
  if ( outcome == "c"){
    OR0mod<-glm(formula=Y~X1+X2+X3,family="gaussian", data=S1data_A0, weights=w)
  } else {
    OR0mod<-glm(formula=Y~X1+X2+X3,family="binomial", data=S1data_A0, weights=w)
  }
  
  p0<- predict(OR0mod,newdata=data, type="response")
  diff<-mean(p1)-mean(p0)
  estimates<-data.frame(estimator="AIPW3",mu1=mean(p1), mu0=mean(p0), diff=diff)
  return(estimates)
}



mu_true <- function(theta0, theta1, outcome = "c"){
  
  #####
  #  the true mean that from data generating algorithm
  #####
  
  set.seed(1234)
  
  
  X1 <- rnorm(1e6)
  X2 <- rnorm(1e6)
  X3 <- rnorm(1e6)
  X  <- cbind(1, X1, X2, X3)
  # X  <- cbind(1, X1, X2, X3, X1*X2, X1*X3, X2*X3) # for model misspecification cases
  
  if (outcome == "c") {
    mu1 <- mean(X %*%theta1)
    mu0 <- mean(X %*%theta0)
  } else {
    mu1_in <- 1/(1+exp(-(X %*%theta1))) # Y^(a=1) ~ Bernoulli --> mean = prob of Y^(a=1) = E[Y|X, S=1, A=1]
    mu0_in <- 1/(1+exp(-(X %*%theta0))) # mean = prob of Y^(a=0) = E[Y|X, S=1, A=0]
    
    mu1 <- mean(mu1_in) # mean = E[Y] = mean(E[Y|X, S=1, A=1])
    mu0 <- mean(mu0_in) # mean = E[Y] = mean(E[Y|X, S=1, A=0])
  }
  diff<-mu1- mu0
  
  # remove storage for X to release memory
  rm(X)
  gc()
  
  
  return(list(mu1 = mu1, mu0 = mu0, diff = diff))
}



trial_est<-function(data){
  #####
  # The raw mean that from simulation data, without any model/ weightings
  # Trial-only: S = 1 only
  #####
  mu1 <- mean(data[data$S == 1 & data$A == 1, "Y"])
  mu0 <- mean(data[data$S == 1 & data$A == 0, "Y"])
  diff <-mu1- mu0
  estimates<-data.frame(estimator="Trial-only",mu1=mu1, mu0=mu0, diff=diff)
  return(estimates)
}
