
# EXAMPLE CODE FOR NON-DESTRUCTIVE ONE-SHOT DEVICES UNDER GENERAL LIFETIMES


source("dpd estimation all.R")
source("waldtest.R")

N = 200
stressvec = c(30,40) 
ITvec= c(6,10,14,18,20,24,28,32,36,40,44,48,52) 
tauvec = c(18,52) 

model.distribution = "gamma"
theta = c(5.3, -0.05, 1.5) 
theta.cont = c(5.7, -0.05, 1.5) 

beta.list = c(0.2,0.4,0.6,0.8,1)

n <- simulate.sample(theta, theta.cont, N, tauvec, stressvec, ITvec, model.distribution, seed = sample(100,1))


# MDPDE
estimators = estimate.function(N, tauvec, stressvec, ITvec, n, initial=c(5, -0.02, 1), model.distribution = model.distribution)
estimators

#WALD TEST FOR TESTING H0: theta_0 = -0.05 VS H1: theta_0 != -0.05
m.function = function(theta){return(c(theta[2]+0.05))}
M=c(0,1,0) 
signif=0.05
reject_point = qchisq(1-signif, 1)

for(ind in 1:length( c(0,beta.list) ) ){
  W = Wtest(estimators[[ind]],M,N,tauvec,stressvec,ITvec, c(0,beta.list)[ind], model.distribution) 
  print( paste( "For beta = ", c(0,beta.list)[ind], "the null is rejected", W>reject_point) )
  
}
