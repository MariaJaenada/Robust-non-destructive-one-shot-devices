
# Wald-type test for the SSALT model under general lifetime distributions
# Computes main matrices involved and asymptotic distribution 

library("hypergeo")

density.distribution.gamma <-function(lambda, tauvec, stressvec, t){
  
  nu = lambda[length(lambda)]
  
  k = length(tauvec) #number of stress level/number of time of stress change
  if(length(lambda)!= (k+1)){stop("number of stress change = number stress levels = number of dist. parameters")}
  
  #derivative of the probabilities
  tauvec1 = c(0,tauvec[1:(k-1)]); stressvec1 = c(0,stressvec[1:(k-1)]) #shifted positions to substract
  
  hi1 = c(0)
  for (i in 1:(length(stressvec)-1)){
    #calculate h_i
    suma = 0
    for(k in 1:i){suma = suma + (-1/lambda[i+2-k]+1/lambda[i+1-k])*tauvec[i+1-k]}
    hi1[i+1] = lambda[i+1]*suma
  }
  #create position vector of the times of stress change
  stress.position = 1+sum(t >= tauvec) #> or >=
  
  if(stress.position > length(tauvec)){stress.position = stress.position-1} #greatest position
    
  ttras = (t+hi1[stress.position])
  if (t <= tauvec[length(tauvec)]){ 
    fsol = dgamma(ttras, shape = nu, scale = lambda[stress.position]) #(nu/lambda[stress.position])*ttras^(nu-1)*exp(-ttras^nu)
  }else{
    stop("lifetime out of experimental limits")
  }
  
  return(fsol)
}


cumulative.distribution.gamma <-function(lambda, tauvec, stressvec, t){
  
  nu = lambda[length(lambda)]
  
  k = length(tauvec) #number of stress level/number of time of stress change
  if(length(lambda)!= (k+1)){stop("number of stress change = number stress levels = number of dist. parameters")}
  
  #derivative of the probabilities
  tauvec1 = c(0,tauvec[1:(k-1)]); stressvec1 = c(0,stressvec[1:(k-1)]) #shifted positions to substract
  
  hi1 = c(0)
  for (i in 1:(length(stressvec)-1)){
    #calculate h_i
    suma = 0
    for(k in 1:i){suma = suma + (-1/lambda[i+2-k]+1/lambda[i+1-k])*tauvec[i+1-k]}
    hi1[i+1] = lambda[i+1]*suma
  }
  #create position vector of the times of stress change
  stress.position = 1+sum(t >= tauvec) #> or >=
  
  if(stress.position > length(tauvec)){stress.position = stress.position-1} #greatest position
  
  ttras = (t+hi1[stress.position])
  if (t <= tauvec[length(tauvec)]){ 
    fsol = pgamma(ttras, shape = nu, scale = lambda[stress.position]) 
  }else{
    stop("lifetime out of experimental limits")
  }
  
  return(fsol)
}

density.distribution.weibull <-function(lambda, tauvec, stressvec, t){
  
  nu = lambda[length(lambda)]
  
  k = length(tauvec) #number of stress level/number of time of stress change
  if(length(lambda)!= (k+1)){stop("number of stress change = number stress levels = number of dist. parameters")}
  
  hi1 = shifted.times(tauvec, lambda)
  
  #create position vector of the times of stress change
  stress.position = 1+sum(t >= tauvec) #> or >=
  
  if(stress.position > length(tauvec)){stress.position = stress.position-1} #greatest position
  
  ttras = (t+hi1[stress.position])/lambda[stress.position]
  if (t <= tauvec[length(tauvec)]){ 
    fsol = (nu/lambda[stress.position])*ttras^(nu-1)*exp(-ttras^nu)
  }else{
    stop("lifetime out of experimental limits")
  }
  
  return(fsol)
}

density.distribution.exponential <-function(lambda, tauvec, stressvec, t){
  
  nu = lambda[length(lambda)]
  
  k = length(tauvec) #number of stress level/number of time of stress change
  if(length(lambda)!= (k+1)){stop("number of stress change = number stress levels = number of dist. parameters")}
  
  hi1 = shifted.times(tauvec, lambda)
  
  #create position vector of the times of stress change
  stress.position = 1+sum(t >= tauvec) #> or >=
  
  if(stress.position > length(tauvec)){stress.position = stress.position-1} #greatest position
  
  ttras = (t+hi1[stress.position])
  if (t <= tauvec[length(tauvec)]){ 
    fsol = exp(-ttras/lambda[stress.position])/(lambda[stress.position])
  }else{
    stop("lifetime out of experimental limits")
  }
  
  return(fsol)
}

Wmatrix.gamma <- function(estimate, tauvec, stressvec, ITvec, distribution){
  
  lambda = c(exp(estimate[1] + estimate[2]*stressvec), estimate[3])
  k = length(tauvec) #number of stress level/number of time of stress change
  
  if(length(lambda)!= (k+1)){stop("number of stress change = number stress levels = number of dist. parameters")}
  
  #derivative of the probabilities
  #shifted positions
  tauvec1 = c(0,tauvec[1:(k-1)])
  stressvec1 = c(0,stressvec[1:(k-1)]) 
  
  hi1 =  hi1ast = c(0)

  for (i in 1:(length(stressvec)-1)){
    #calculate h_i
    suma = 0
    for(j in 1:i){suma = suma + (-1/lambda[i+2-j]+1/lambda[i+1-j])*tauvec[i+1-j]}
    hi1[i+1] = lambda[i+1]*suma
    #calculate h_i^ast
    suma = 0
    for(j in 1:i){suma = suma + (stressvec[i+2-j]/lambda[i+2-j]-stressvec[i+1-j]/lambda[i+1-j])*tauvec[i+1-j]}
    hi1ast[i+1] = lambda[i+1]*suma + hi1[i+1]*stressvec[i+1]
  }
  
  #create position vector of the times of stress change
  stress.position = rep(1,length(ITvec))
  for(i in 1:k){stress.position[ITvec >= tauvec[i]] = i+1}
  stress.position[length(ITvec)] = i
  
  #derivative matrix
  W = matrix(NA, nrow = length(ITvec)+1, ncol = length(estimate)) #length(theta) = 3
  density.distribution = match.fun(paste0("density.distribution.", distribution))
  model.distribution = match.fun(paste0("cumulative.distribution.", distribution))
  
  zj1_1 = zj1_2 = zj1_3 = 0 #initial values of -1 wj = zj-zj1
  for(j in 1:length(ITvec)){
    densitytj = density.distribution(lambda, tauvec, stressvec,ITvec[j])
    distributionj =  model.distribution(lambda, tauvec, stressvec, ITvec[j])
    i = stress.position[j]
    zj_1 = densitytj*(-ITvec[j]-hi1[i])
    zj_2 = densitytj*((-ITvec[j]-hi1[i])*stressvec[i]+hi1ast[i])
    alpha = lambda[length(lambda)]
    zj_3 = distributionj*(log((ITvec[j]+hi1[i])/lambda[i]) - digamma(alpha) ) - ((ITvec[j]+hi1[i])/lambda[i])^alpha*genhypergeo(c(alpha, alpha), c(alpha+1, alpha+1), -((ITvec[j]+hi1[i])/lambda[i]))/(alpha^2*gamma(alpha))
    W[j,] = c(zj_1-zj1_1, zj_2-zj1_2, zj_3-zj1_3)
    
    #update
    zj1_1 = zj_1; zj1_2 = zj_2; zj1_3 = zj_3
  }
  
  W[j+1,]=c(-zj_1, -zj_2, -zj_3)
  return(W)
}


Wmatrix.lognormal <- function(estimate, tauvec, stressvec, ITvec, distribution){
  
  lambda = c(estimate[1] + estimate[2]*stressvec, estimate[3])
  k = length(tauvec) #number of stress level/number of time of stress change
  
  if(length(lambda)!= (k+1)){stop("number of stress change = number stress levels = number of dist. parameters")}
  
  #derivative of the probabilities
  tauvec1 = c(0,tauvec[1:(k-1)]); stressvec1 = c(0,stressvec[1:(k-1)]) #shifted positions to substract
  
  hi1 = shifted.times(tauvec, exp(lambda))
  hi1ast = c(0)
  for (i in 1:(length(stressvec)-1)){
    #calculate h_i^ast
    suma = 0
    for(j in 1:i){suma = suma + (stressvec[i+2-j]/exp(lambda[i+2-j])-stressvec[i+1-j]/exp(lambda[i+1-j]))*tauvec[i+1-j]}
    hi1ast[i+1] = exp(lambda[i+1])*suma + hi1[i+1]*stressvec[i+1]
  }
  
  #create position vector of the times of stress change
  stress.position = rep(1,length(ITvec))
  for(i in 1:k){stress.position[ITvec >= tauvec[i]] = i+1}
  stress.position[length(ITvec)] = i
  
  #derivative matrix
  W = matrix(NA, nrow = length(ITvec)+1, ncol = length(estimate)) #length(theta) = 3
  
  zj1_1 = zj1_2 = zj1_3 = 0 #initial values of -1 wj = zj-zj1
  for(j in 1:length(ITvec)){
    sigma = lambda[length(lambda)]
    i = stress.position[j]
    
    densitytj = (-1/sigma)*dnorm((log(ITvec[j]+hi1[i])-lambda[i])/sigma)
    
    zj_1 = densitytj
    zj_2 = densitytj*(stressvec[i]-hi1ast[i]/(ITvec[j]+hi1[i]))
    zj_3 = densitytj*(log(ITvec[j]+hi1[i])-lambda[i])/sigma
    W[j,] = c(zj_1-zj1_1, zj_2-zj1_2, zj_3-zj1_3)
    
    #update
    zj1_1 = zj_1; zj1_2 = zj_2; zj1_3 = zj_3
  }
  
  W[j+1,]=c(-zj_1, -zj_2, -zj_3)
  return(W)
}

Dbeta <- function(estimate, tauvec, stressvec, ITvec, exponent, distribution){
  
  if (model.distribution == "lognormal"){
    lambda = c(estimate[1]+ estimate[2]*stressvec, estimate[3]) #for logormal mui is linearly related
  }else{
    lambda = c(exp(estimate[1]+ estimate[2]*stressvec), estimate[3])
  }
  
  pi = theoretical.probability(lambda, tauvec, stressvec, ITvec, distribution)
  return(diag(pi^exponent))
}

matrixJ <- function(estimate, tauvec, stressvec, ITvec, beta,distribution){
  
  Wmatrix = match.fun(paste0("Wmatrix.", distribution))
  W <- Wmatrix(estimate, tauvec, stressvec, ITvec, distribution)
  D1 <- Dbeta(estimate, tauvec, stressvec, ITvec, beta-1, distribution)
  return(t(W)%*%D1%*%W)
}

matrixK <- function(estimate, tauvec, stressvec, ITvec, beta, distribution){
  Wmatrix = match.fun(paste0("Wmatrix.", distribution))
  W <- Wmatrix(estimate, tauvec, stressvec, ITvec, distribution)
  D2 <- Dbeta(estimate, tauvec, stressvec, ITvec, 2*beta-1, distribution)
  
  if (model.distribution == "lognormal"){
    lambda = c(estimate[1]+ estimate[2]*stressvec, estimate[3]) #for logormal mui is linearly related
  }else{
    lambda = c(exp(estimate[1]+estimate[2]*stressvec), estimate[3])
  }
  
  pi = theoretical.probability(lambda, tauvec, stressvec, ITvec, distribution)
  PI = (pi^beta)%*%(t(pi)^beta)
  return(t(W)%*%(D2-PI)%*%W)
}

Wtest <- function(estimate,M, N, tauvec, stressvec, ITvec, beta, distribution){
  
  J = matrixJ(estimate, tauvec, stressvec, ITvec, beta, distribution)
  K = matrixK(estimate, tauvec, stressvec, ITvec, beta, distribution) #beta determines the estimator to which I am calculating the asymptotic matrix
  
  m.theta = m.function(estimate) #dimension 1,2,3
  
  W= N*t(m.theta)%*%solve(t(M)%*%solve(J)%*%K%*%solve(J)%*%M)%*%m.theta
  return(W)
}

asymptotic <- function(estimate, theta, N, tauvec, stressvec, ITvec, beta, distribution){
  J = matrixJ(estimate, tauvec, stressvec, ITvec, beta, distribution)
  K = matrixK(estimate, tauvec, stressvec, ITvec, beta, distribution) #beta determines the estimator to which I am calculating the asymptotic matrix
  
  varcov = solve(J)%*%K%*%solve(J)
  asymp0 = (sqrt(N)/sqrt(varcov[1,1]))*(estimate[1]-theta[1])
  asymp1 = (sqrt(N)/sqrt(varcov[2,2]))*(estimate[2]-theta[2])
  asymp2 = (sqrt(N)/sqrt(varcov[3,3]))*(estimate[3]-theta[3])
  
  return(c(asymp0,asymp1, asymp2))
}

asymptotic.variances <- function(estimate,  N, tauvec, stressvec, ITvec, beta, distribution){
  J = matrixJ(estimate, tauvec, stressvec, ITvec, beta, distribution)
  K = matrixK(estimate, tauvec, stressvec, ITvec, beta, distribution) #beta determines the estimator to which I am calculating the asymptotic matrix
  
  varcov = solve(J)%*%K%*%solve(J)
  
  return(sqrt(diag(varcov))/sqrt(N))
}

asymptotic.IC.accuracy <- function(theta, theta.cont,N, tauvec, stressvec, ITvec,distribution,  B=100){
  
  beta.list = c(0.2,0.4,0.6,0.8,1)
  
  count = 0
  betas = c(0,beta.list)
  asympIC = matrix(0, nrow =3, ncol=length(betas))
  zalpha = qnorm(0.975)
  for(b in 1:B){
    n = simulate.sample(theta,  theta.cont, N, tauvec, stressvec, ITvec, distribution, seed = b)
    estimators = estimate.function(N, tauvec, stressvec, ITvec, n, model.distribution = distribution)
    
    if(sum(sapply(estimators, function(x) sum(is.na(x))))<1){ #none of the estimators is nan
      
      for(ind in 1:length(betas)){
        asympIC[,ind] = asympIC[,ind]+ 2*zalpha*asymptotic.variances(estimators[[ind]], N, tauvec, stressvec, ITvec, betas[ind], distribution)
        }
      count = count+1
    }
    asympIC = asympIC/count
  }
  print(count)
  return(asympIC)
}

asymptotic.distribution <- function(theta, theta.cont,N, tauvec, stressvec, ITvec,distribution,  B=100){
  
  beta.list = c(0.2,0.4,0.6,0.8,1)

  count = 0
  betas = c(0,beta.list)
  asymp = matrix(NA, nrow =B, ncol=length(betas))
  
  for(b in 1:B){
    n = simulate.sample(theta,  theta.cont, N, tauvec, stressvec, ITvec, distribution, seed = b)
    estimate.function = match.fun(paste0("estimate.function.", distribution))
    estimators = estimate.function(N, tauvec, stressvec, ITvec, n)
    
    if(sum(sapply(estimators, function(x) sum(is.na(x))))<1){ #none of the estimators is nan
      
      for(ind in 1:length(betas)){
        asymp[b, ind] = asymptotic(estimators[[ind]], theta, N, tauvec, stressvec, ITvec, betas[ind], distribution)
      }
      count = count+1
    }
  }
  print(count)
  
  for(ind in 1:length(betas)){
    print("media y desviación"); print(mean(asymp[, ind]));  print(sd(asymp[, ind]))
  }
  
  #plot results
  for (j in 1:length(betas)){
    x11()
    hist(asymp[,j],freq=FALSE)
  }
  #return(asymp)
}
