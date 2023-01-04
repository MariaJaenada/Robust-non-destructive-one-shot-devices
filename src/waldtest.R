
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
  
  
  # hi1 = c(0,  (-1+lambda[2]/lambda[1])*tauvec[1] )
  # hi1ast = c(0,  hi1[2]*stressvec[2] + lambda[2]*(stressvec[2]/lambda[2] - stressvec[1]/lambda[1])*tauvec[1] )
  
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
      #densitytj*(ITvec[j]+hi1[i])/lambda[length(lambda)]*log((ITvec[j]+hi1[i])/lambda[i])
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
        #asymp0[b,ind] = 2*zalpha*asymp[1]
        #asymp1[b,ind] = 2*zalpha*asymp[2]
        #asymp2[b,ind] = 2*zalpha*asymp[3]
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

# empirical level or power 
signif= 0.05

level_estimate <- function(theta,  theta.cont, M, N, tauvec, stressvec, ITvec,distribution, signif, B=50){
  
  beta.list = c(0.2,0.4,0.6,0.8,1)
  reject_point = qchisq(1-signif, 1) #qchisq(1-signif,1)
  
  count = 0
  betas = c(0,beta.list)
  level = rep(0,length(betas)+1)
  
  for(b in 1:B){
    n = simulate.sample(theta,  theta.cont, N, tauvec, stressvec, ITvec, distribution, seed = b)
    #print(b);print(n);print(sum(n[-length(n)]))
    
    #OPTIMUM
    # pilot=c(0,0,0)
    # for(p in 1:length(estimators)){ #no se como usar funciones superiores sapply
    #   pilot = pilot+ estimators[[p]]/length(estimators)
    # }
    # optimum = IJW(n = n, pilot = pilot) ##
    p=c(0,0,0)
    for(beta.init in c(0,0.4,0.8)){
      pilot = tryCatch(optimr(par= c(4, -0.02, 1), DPDloss,  N = N, tauvec = tauvec, stressvec = stressvec, ITvec = ITvec, n = n,
                          beta = beta.init, model.distribution = distribution, method ="BFGS"),
                   error=function(sol){sol$code=3;return(NA)})$par
      p = p + pilot/3
    }
    
    optimum = list(0.1, p) #IJW(n = n, pilot = pilot, model.distribution = model.distribution) 
    #list(0.1, pilot) #
    #estimate.function = match.fun(paste0("estimate.function.", distribution))
    estimators = estimate.function(N, tauvec, stressvec, ITvec, n, model.distribution = distribution)
    if(sum(sapply(estimators, function(x) sum(is.na(x))))<1){ #none of the estimators is nan
      
      for(ind in 1:length(betas)){
        W = Wtest(estimators[[ind]],M,N,tauvec,stressvec,ITvec, betas[ind], distribution) #each Wtest is constructed by approx the matrices with its own estimate
        #print(W)
        level[ind] = level[ind]+(W>reject_point)
      }
      #optimum
      W = Wtest(optimum[[2]],M,N,tauvec,stressvec,ITvec, optimum[[1]],distribution) 
      level[ind+1] = level[ind+1] + (W>reject_point)
      count = count+1
    }
  }
  
  level= level/count
  
  return(level)
}


pvalue_estimate <- function(theta,  theta.cont, M, N, tauvec, stressvec, ITvec,distribution, signif, B=50){
  
  beta.list = c(0.2,0.4,0.6,0.8,1)
  reject_point = qchisq(1-signif, 1) #qchisq(1-signif,1)
  
  count = 0
  betas = c(0,beta.list)
  pvalue = rep(0,length(betas)+1)
  
  for(b in 1:B){
    n = simulate.sample(theta,  theta.cont, N, tauvec, stressvec, ITvec, distribution, seed = b)
    #print(b);print(n);print(sum(n[-length(n)]))
    
    #OPTIMUM
    # pilot=c(0,0,0)
    # for(p in 1:length(estimators)){ #no se como usar funciones superiores sapply
    #   pilot = pilot+ estimators[[p]]/length(estimators)
    # }
    # optimum = IJW(n = n, pilot = pilot) ##
    p=c(0,0,0)
    for(beta.init in c(0,0.4,0.8)){
      pilot = tryCatch(optimr(par= c(4, -0.02, 1), DPDloss,  N = N, tauvec = tauvec, stressvec = stressvec, ITvec = ITvec, n = n,
                          beta = beta.init, model.distribution = distribution, method ="BFGS"),
                   error=function(sol){sol$code=3;return(NA)})$par
      p = p + pilot/3
    }
    
    optimum = list(0.1, p) #IJW(n = n, pilot = pilot, model.distribution = model.distribution) 
    # #
    #estimate.function = match.fun(paste0("estimate.function.", distribution))
    estimators = estimate.function(N, tauvec, stressvec, ITvec, n, model.distribution = distribution)
    if(sum(sapply(estimators, function(x) sum(is.na(x))))<1){ #none of the estimators is nan
      
      for(ind in 1:length(betas)){
        W = Wtest(estimators[[ind]],M,N,tauvec,stressvec,ITvec, betas[ind], distribution) #each Wtest is constructed by approx the matrices with its own estimate
        #print(W)
        pvalue[ind] = pvalue[ind]+(1-pchisq(W,1))#(W>reject_point)
      }
      #optimum
      W = Wtest(optimum[[2]],M,N,tauvec,stressvec,ITvec, optimum[[1]],distribution) 
      pvalue[ind+1] = pvalue[ind+1] + (1-pchisq(W,1))
      count = count+1
    }
  }
  
  pvalue= pvalue/count
  
  return(pvalue)
}

beta.error_estimate <- function(theta, elle, theta.cont, M, N, tauvec, stressvec, ITvec,distribution, signif, B=50){
  
  #type II error: Probability of acepting under the alternative hypothesis
  #theta: true parameter BELONGING TO THE NULL HYPOTHESIS: closest element
  #elle: vector direction of the violating null restriction. For the moment elle=number violating only second component
  #theta.cont: Contamination of the true theta = theta+1/sqrt(N)*elle GIVEN IN TERMS OF CONTAMINATION OF THE CLOSEST
  
  beta.list = c(0.2,0.4,0.6,0.8,1)
  reject_point = qchisq(1-signif, 1)
  velle = c(0,1,0)*elle
  ncp = drop(M%*%velle)^2
  
  count = 0
  betas = c(0,beta.list)
  typeII = rep(0,length(betas)+1)
  
  theta.a = theta+1/sqrt(N)*velle
  theta.a.cont = theta.cont + 1/sqrt(N)*velle
  
  for(b in 1:B){
    n = simulate.sample(theta.a,  theta.a.cont, N, tauvec, stressvec, ITvec, distribution, seed = b)
    #print(b);print(n);print(sum(n[-length(n)]))
    
    #OPTIMUM
    # pilot=c(0,0,0)
    # for(p in 1:length(estimators)){ #no se como usar funciones superiores sapply
    #   pilot = pilot+ estimators[[p]]/length(estimators)
    # }
    # optimum = IJW(n = n, pilot = pilot) ##
    p=c(0,0,0)
    for(beta.init in c(0,0.4,0.8)){
      pilot = tryCatch(optimr(par= c(4, -0.02, 1), DPDloss,  N = N, tauvec = tauvec, stressvec = stressvec, ITvec = ITvec, n = n,
                          beta = beta.init, model.distribution = distribution, method ="BFGS"),
                   error=function(sol){sol$code=3;return(NA)})$par
      p = p + pilot/3
    }
    
    optimum = list(0.1, p) #IJW(n = n, pilot = pilot, model.distribution = model.distribution) 
    
    estimators = estimate.function(N, tauvec, stressvec, ITvec, n, initial=pilot, model.distribution = distribution)
    if(sum(sapply(estimators, function(x) sum(is.na(x))))<1){ #none of the estimators is nan
      
      for(ind in 1:length(betas)){
        W = Wtest(estimators[[ind]],M,N,tauvec,stressvec,ITvec, betas[ind], distribution) #each Wtest is constructed by approx the matrices with its own estimate
        
        typeII[ind] = typeII[ind]+ pchisq(W, df = 1, ncp = ncp) # TOO MUCH DIFFERENT RESULTS pchisq(W, df = 1, ncp = ncp)
        #print(pchisq(W, df = 1, ncp = ncp))
      }
      #optimum
      W = Wtest(optimum[[2]],M,N,tauvec,stressvec,ITvec, optimum[[1]],distribution) 
      typeII[ind+1] = typeII[ind+1] + pchisq(W, df = 1, ncp = ncp)
      count = count+1
    }
  }
  
  typeII= typeII/count
  
  return(typeII)
}

power_estimate <- function(theta, elle.grid, theta.cont, M, N, signif, B=50){
  
  #type II error: Probability of acepting under the alternative hypothesis
  #theta: true parameter BELONGING TO THE NULL HYPOTHESIS: closest element
  #elle: vector direction of the violating null restriction. For the moment elle=number violating only second component
  #theta.cont: Contamination of the true theta = theta+1/sqrt(N)*elle GIVEN IN TERMS OF CONTAMINATION OF THE CLOSEST
  
  reject_point = qchisq(1-signif, 1)
  typeII = c()
  for (l in 1:length(elle.grid)){
    
    velle = c(0,1,0)*elle.grid[l]
    ncp = drop(M%*%velle)^2
    
    #theta.a = theta+1/sqrt(N)*velle
    #theta.a.cont = theta.cont + 1/sqrt(N)*velle
    
    typeII[l] = pchisq(reject_point, df=1, ncp=ncp) 
  }
  return(typeII)
}
  

chi.distance <- function(theta.null,  theta.real, M, N, tauvec, stressvec, ITvec, signif, B=50){
  
  beta.list = c(0.2,0.4,0.6,0.8,1)
  distance = 0
  
  for(b in 1:B){
    n = simulate.sample(theta.null,  theta.null, N, tauvec, stressvec, ITvec, seed = b)
    #print(b);print(n);print(sum(n[-length(n)]))
    n1 = simulate.sample(theta.real,  theta.real, N, tauvec, stressvec, ITvec, seed = b)
    distance = distance+mean((n-n1)^2)
  }
  
  distance = distance/B
  
  return(distance)
}

# OBTAIN THE LEVEL FOR GRAPHICS 
level.contamination <-function(theta, m, d, N, tauvec, stressvec, ITvec, signif, reduction.param.vec, B=50){

  beta.list = c(0.2,0.4,0.6,0.8,1)
  level.mat = matrix(0, nrow =  length(reduction.param.vec), ncol = length(beta.list)+1)

  for(red.param in 1:length(reduction.param.vec)){
    level.mat[red.param,] = level_estimate(theta, m, d, N, tauvec, stressvec, ITvec, signif, reduction.param.vec[red.param], B)
  }
  return(level.mat)
}

level.samplesize <-function(theta, m, d, Nvec, tauvec, stressvec, ITvec, signif, reduction.param, B=50){

  beta.list = c(0.2,0.4,0.6,0.8,1)
  level.mat = matrix(0, nrow =  length(Nvec), ncol = length(beta.list)+1)

  for(samplesize in 1:length(Nvec)){
    level.mat[samplesize,] = level_estimate(theta, m, d, Nvec[samplesize], tauvec, stressvec, ITvec, signif, reduction.param, B)
  }
  return(level.mat)
}

# OBTAIN THE GRAPHICS
graphics.contamination <-function(theta, theta.cont.vec, M, N, tauvec, stressvec, ITvec, distribution, signif, B=50, aylabel = "level", axlabel = expression(tilde(a)["0"])){
  
  beta.list = c(0.2,0.4,0.6,0.8,1)
  level.mat = matrix(0, nrow =  length(theta.cont.vec), ncol = length(beta.list)+2)
  
  for(cont in 1:length(theta.cont.vec)){
    level.mat[cont,] = pvalue_estimate(theta,theta.cont.vec[[cont]], M, N, tauvec, stressvec, ITvec, distribution, signif, B)
  }

  
  #plot results
  #if(aylabel == "level"){position = "topleft"}else{position = "bottomleft"}
  position = "topright"
  x11()
  #colores=c(1,3,6,4,2,7,5,8,9,10,11,12)
  colores= c("#0D0D0D","#3B3B3B", "#1874CD","#87CEFA", "#008B00", "#FF7F00", "#EE7942", "#CD2990", "#EE7AE9", "#CD3333")
  formas=c(17,3,15,5,2,7,1,8,9,10,11,12)

  epsilons = epsilons
  plot(epsilons, level.mat[,1],col=colores[1],"o",lty=1, pch=formas[1],lwd=1.6, cex=1.05, ylim=c(min(level.mat),max(level.mat)), xlab = axlabel, ylab = aylabel)
  for (j in 2:(length(beta.list)+1)){
    lines(epsilons ,level.mat[,j],col=colores[j],"o",lty=j, pch=formas[j],lwd=1.6,cex=1.05)
  }
  legend(position, title= expression(beta), inset=c(0,0.0), c("MLE",as.character(c(0.2,0.4,0.6,0.8,1))), col = colores, lty = c(1,2:6),pch= formas, cex =0.75 )
  #title("")
  return(level.mat)
}

graphics.samplesize <-function(theta, theta.cont, M, Nvec, tauvec, stressvec, ITvec, distribution, signif, B=50, aylabel = "level"){
  
  beta.list = c(0.2,0.4,0.6,0.8,1)
  level.mat = matrix(0, nrow =  length(Nvec), ncol = length(beta.list)+2)
  
  for(samplesize in 1:length(Nvec)){
    level.mat[samplesize,] = level_estimate(theta, theta.cont, M, Nvec[samplesize], tauvec, stressvec, ITvec, distribution, signif,  B)
  }
  # 
  # #plot results
  # if(aylabel == "level"){position = "topleft"}else{position = "bottomleft"}
  # #x11()
  # #colores=c(1,3,6,4,2,7,5,8,9,10,11,12)
  # colores= c("#0D0D0D","#3B3B3B", "#1874CD","#87CEFA", "#008B00", "#FF7F00", "#EE7942", "#CD2990", "#EE7AE9", "#CD3333")
  # formas=c(17,3,15,5,2,7,1,8,9,10,11,12)
  # 
  # plot(Nvec, level.mat[,1],col=colores[1],"o",lty=1, pch=formas[1],lwd=1.6, cex=1.05, ylim=c(min(level.mat),max(level.mat)), xlab = expression(N), ylab = aylabel)
  # for (j in 2:length(beta.list)+2){
  #   lines(Nvec ,level.mat[,j],col=colores[j],"o",lty=j, pch=formas[j],lwd=1.6,cex=1.05)
  # }
  # legend(position, title= expression(beta), inset=c(0,0.0), c("MLE",as.character(c(0.2,0.4,0.6,0.8,1)), "Opt"), col = colores, lty = c(1,2:6),pch= formas, cex =0.75 )
  # #title("")
  return(level.mat)
}

