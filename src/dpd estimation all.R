# MULTINOMIAL REPRESENTATION OF NON-DESTRUCTIVE ONE-SHOT DEVICES UNDER GENERAL LIFETIMES
# Parameters
# ---------------------------------------
# model.distribution or distribution = {exponential, Weibull, gamma, lognormal}
# N : Number of devices
# n : sample
# tauvec : vector of times of stress change
# ITvec : vector of inspection times
# theta : parameter vector of general lifetime scale = exp(theta[1]+theta[2]x), common shape = theta[3]
# lambda: lifetime distribution parameters (different for each distribution) 
# --------------------------------------

library(optimx)

survival.units <- function(N,n){
  
  # Compute survival units on each step of the experiment
  
  survival = c(N)
  for(step in 2:length(n)){survival[step] = survival[step-1]-n[step-1]}
  return(survival)
}

shifted.times <- function(tauvec, lambda){
  
  #compute vector of shited times for the CE model
  
  k = length(tauvec)
  tauvec1 = c(0,tauvec[1:(k-1)]) #shifted positions
  #stressvec1 = c(0,stressvec[1:(k-1)]) #shifted stress levels
  
  #calculate h_i vector
  hi1 = c(0)
  for (i in 1:(length(tauvec)-1)){
    suma = 0
    for(k in 1:i){suma = suma + (-1/lambda[i+2-k]+1/lambda[i+1-k])*tauvec[i+1-k]}
    hi1[i+1] = lambda[i+1]*suma
  }
  return(hi1)
}

model.distribution.exponential <-function(lambda, tauvec, stressvec, t){
  
  #compute lifetime for the step-stress model under exponential distributions
  
  alpha = lambda[length(lambda)] #useless but maintained to generalize the code
  k = length(tauvec) #number of stress level/number of time of stress change
  
  if(length(lambda)!= (k+1)){
    stop("number of stress change = number stress levels = number of dist. parameters")
  }
  
  hi1 = shifted.times(tauvec, lambda)
  
  #create position vector of the times of stress change
  stress.position = sum(t > tauvec)+1
  ttras = (t+hi1[stress.position])
  
  if (t <= tauvec[length(tauvec)]){ 
    DFsol = 1-exp(-ttras/lambda[stress.position]) 
    
  }else{stop("lifetime out of experimental limits")}
  
  return(DFsol)
}

model.distribution.weibull <-function(lambda, tauvec, stressvec, t){
  
  alpha = lambda[length(lambda)]
  k = length(tauvec) #number of stress level/number of time of stress change
  
  if(length(lambda)!= (k+1)){
    stop("number of stress change = number stress levels = number of dist. parameters")
  }
  
  hi1 = shifted.times(tauvec, lambda)
  
  #create position vector of the times of stress change
  stress.position = sum(t > tauvec)+1
  ttras = (t+hi1[stress.position])
  
  if (t <= tauvec[length(tauvec)]){ 
    DFsol = 1-exp(-(ttras/lambda[stress.position])^alpha)
  }else{stop("lifetime out of experimental limits")}
  
  return(DFsol)
}

model.distribution.gamma <-function(lambda, tauvec, stressvec, t){
  
  alpha = lambda[length(lambda)]
  k = length(tauvec) #number of stress level/number of time of stress change
  
  if(length(lambda)!= (k+1)){
    stop("number of stress change = number stress levels = number of dist. parameters")
  }
  
  hi1 = shifted.times(tauvec, lambda)
  
  #create position vector of the times of stress change
  stress.position = sum(t > tauvec)+1
  ttras = (t+hi1[stress.position])
  
  if (t <= tauvec[length(tauvec)]){ 
    DFsol = pgamma(ttras, shape = alpha, scale = lambda[stress.position])
  }else{stop("lifetime out of experimental limits")}
  
  return(DFsol)
}

model.distribution.lognormal <-function(lambda, tauvec, stressvec, t){
  
  sigma = lambda[length(lambda)]
  k = length(tauvec) #number of stress level/number of time of stress change
  
  if(length(lambda)!= (k+1)){
    stop("number of stress change = number stress levels = number of dist. parameters")
  }
  
  hi1 = shifted.times(tauvec, exp(lambda))
  
  #create position vector of the times of stress change
  stress.position = sum(t > tauvec)+1
  ttras = (t+hi1[stress.position])
  
  if (t <= tauvec[length(tauvec)]){ 
    DFsol = pnorm((log(ttras)-lambda[stress.position])/sigma, 0,1)
  }else{stop("lifetime out of experimental limits")}
  
  return(DFsol)
}
  
theoretical.probability <- function(lambda, tauvec, stressvec, ITvec, 
                                    model.distribution = "gamma"){
  
  f <- paste0("model.distribution.", model.distribution)
  model.distribution = match.fun(f)
  
  th = c(model.distribution(lambda, tauvec, stressvec, ITvec[1]))
  for(step in 2:length(ITvec)){
    th[step] = model.distribution(lambda, tauvec, stressvec, ITvec[step])-model.distribution(lambda, tauvec, stressvec, ITvec[step-1])
  }
  th[step+1] = 1- model.distribution(lambda, tauvec, stressvec, ITvec[step])
  return(th)
}

simulate.sample <- function(theta, theta.cont, N, tauvec, stressvec, ITvec,  
                            model.distribution, seed){
  
  #theta.cont = c(theta[1]*reduction.param, theta[2], theta[3]) 
  if (model.distribution == "lognormal"){
    lambda = c(theta[1]+ theta[2]*stressvec, theta[3]) #for logormal mui is linearly related
    lambda.cont =c(theta.cont[1]+theta.cont[2]*stressvec,theta.cont[3])
  }else{
    lambda = c(exp(theta[1]+ theta[2]*stressvec), theta[3])
    lambda.cont =c(exp(theta.cont[1]+theta.cont[2]*stressvec),theta.cont[3])
  }
  #lambda.cont = c(reduction.param*lambda[1], lambda[2], lambda[3])
  
  set.seed(seed)
  
  th = theoretical.probability(lambda, tauvec, stressvec, ITvec, model.distribution)
  th[th<1e-8]=0
  
  #contaminate j interval
  j = 3
  
  f <- paste0("model.distribution.", model.distribution)
  model.distribution = match.fun(f)
  
  th[j] = model.distribution(lambda, tauvec, stressvec, ITvec[j]) - model.distribution(lambda.cont, tauvec, stressvec, ITvec[j-1])
  th[j] = max(th[j],0)
  
  #th[length(th)] = max(0,1- sum(th[1:(length(th)-1)])) #correct so the probability vector sums up to 1
  th = th/sum(th)
  n = drop(rmultinom(1, N, prob = th)) 
  
  return(n)
}

loglikelihood <- function(theta, N, tauvec, stressvec, ITvec, n, model.distribution){
  
  #theta.cont = c(theta[1]*reduction.param, theta[2], theta[3]) 
  if (model.distribution == "lognormal"){
    lambda = c(theta[1]+ theta[2]*stressvec, theta[3]) #for logormal mui is linearly related
  }else{
    lambda = c(exp(theta[1]+ theta[2]*stressvec), theta[3])
  }
  th = theoretical.probability(lambda, tauvec, stressvec, ITvec, model.distribution)
  like = 0
  for(step in 1:length(n)){ like = like + n[step]*log(th[step]) }
  
  return(-like)
}

dbeta <- function(hat.p, th.p, beta){
  return(th.p^(beta+1)-((beta+1)/beta)*(hat.p*(th.p^beta)))
}

DPDloss <- function(theta, N, tauvec, stressvec, ITvec, n, beta, model.distribution){
  
  #theta.cont = c(theta[1]*reduction.param, theta[2], theta[3]) 
  if (model.distribution == "lognormal"){
    lambda = c(theta[1]+ theta[2]*stressvec, theta[3]) #for logormal mui is linearly related
  }else{
    lambda = c(exp(theta[1]+ theta[2]*stressvec), theta[3])
  }
  
  th = theoretical.probability(lambda, tauvec, stressvec, ITvec, model.distribution)
  if (beta == 0){DPD=loglikelihood(theta, N, tauvec, stressvec, ITvec, n, model.distribution)
  }else{
    DPD = 0
    for(step in 1:length(n)){ DPD = DPD + dbeta(n[step]/N,th[step], beta)}
  }
  return(DPD)
}

weighted.phi <- function(theta, N, tauvec, stressvec,ITvec, n, beta, model.distribution){
  
  #theta.cont = c(theta[1]*reduction.param, theta[2], theta[3]) 
  if (model.distribution == "lognormal"){
    lambda = c(theta[1]+ theta[2]*stressvec, theta[3]) #for logormal mui is linearly related
  }else{
    lambda = c(exp(theta[1]+ theta[2]*stressvec), theta[3])
  }
  
  th = theoretical.probability(lambda, tauvec, stressvec, ITvec, model.distribution)
  # survival = survival.units(N,n)
  p = n/N
  # p[p==0] = 0.0001 ; p[p==1] = 0.999 #avoiding log(0)
  
  phi=0
  for(step in 1:length(n)){
    phiterm =p[step]*((p[step]/th[step])^beta -1)
    phi = phi + phiterm
    
  }
  return((1/(beta*(beta+1)))*phi)
}

estimate.function <-function(N, tauvec, stressvec, ITvec, n, initial = c(4, -0.02, 1), beta.list = c(0.2,0.4,0.6,0.8,1), model.distribution){
  
  estimates = list()
  
  #initial point depend of the distribution
  if(model.distribution == "gamma"){

    #initial estimates are found from exponential distribution
    init =  tryCatch(optimr(par= initial, loglikelihood,  N = N, tauvec = tauvec, stressvec = stressvec, ITvec = ITvec, n = n, model.distribution = "exponential", method ="BFGS"),
                    #lower = c(-Inf,-Inf,0), upper = c(Inf, 0, Inf), method ="L-BFGS-B"),
                    error=function(sol){sol$code=3;return(NA)})
    if(!is.na(init)){initial = init$par}else{initial = initial}

  }
  # else if(model.distribution == "lognormal"){
  #   
  #   last = length(n)
  #   
  #   #Define mui and sigma mean and sd of log(t+hi) ~ N(mui,sigma)
  #   time.position.left=1
  #   hi1 = shifted.times(t)
  #   for(stress.position in c(1:length(tauvec)) ){
  #       time.position.right = sum(ITvec <= tauvec[stress.position])
  #       ti = ITvec[time.position.left:time.position.right]
  #       time.position.left = time.position.right
  #       
  #       
  #       #ttras = (t+hi1[stress.position])
  #   }
  #   logtimes = (log(ITvec)+c(0,log(ITvec)[-1]))/2 #interval id
  #   
  #   #divide data by intervals
  #   #first compute mu using mean of logtimes (without taking account stress level)
  #   mu = sum(logtimes*n[-last])/sum(n[-last])
  #   
  #   sigma = sum((logtimes-mean(logtimes))^2*n[-last])/sum(n[-last])
  #   logdata = log(n)
  #   mus = n/
  #   initial = c(2,2) 
  # }
  
  #MLE
  MLE =  tryCatch(optimr(par= initial, loglikelihood,  N = N, tauvec = tauvec, stressvec = stressvec, ITvec = ITvec, n = n, model.distribution =  model.distribution, method ="BFGS"),
                         #lower = c(-Inf,-Inf,0), upper = c(Inf, 0, Inf), method ="L-BFGS-B"),
                  error=function(sol){sol$code=3;return(NA)})
  if(!is.na(MLE$par)){estimates[["MLE"]] = MLE$par}else{estimates[["MLE"]] = NA}
  MLE
  #initial = MLE$par
  for(beta in beta.list){
    #initial = (initial + estimates[[length(estimates)]])/2
    DPDE =  tryCatch(optimr(par= initial,DPDloss,  N = N, tauvec = tauvec, stressvec = stressvec, ITvec = ITvec, n = n, beta= beta, model.distribution = model.distribution, method ="BFGS"),
                            #lower = c(-Inf,-Inf,0), upper = c(Inf, 0, Inf), method ="L-BFGS-B"),
                     error=function(sol){sol$code=3;return(NA)})
    #print(DPDE)
    if(!is.na(DPDE$par)){estimates[[paste0("DPD", beta)]] = DPDE$par}else{estimates[[paste0("DPD", beta)]] = NA}
  }
  
  # for(beta in c(-0.8,-0.5,0.6,1)){
  #   
  #   minphi =tryCatch(optimr(par= initial, weighted.phi,  N = N, tauvec = tauvec, stressvec = stressvec, ITvec = ITvec, n = n, beta= beta),
  #                      error=function(sol){sol$code=3;return(NA)})
  #   if(!is.na(minphi)){estimates[[paste0("phi", beta)]] = minphi$par}else{estimates[[paste0("phi", beta)]] = NA}
  #   
  # }
  
  return(estimates)
}

RMSE.function <- function(theta, theta.hat){return(mean(((theta-theta.hat)/theta)^2) )}

simulate <-function(theta, theta.cont.vec,  N, tauvec, stressvec, ITvec, model.distribution, B=50){
  
  beta.list = c(0.2,0.4,0.6,0.8,1)
  nbetas = length(beta.list)+2
  res.mat = matrix(0, nrow =  length(theta.cont.vec), ncol = nbetas)
  
  for(cont in 1:length(theta.cont.vec)){
    count = 0
    #print(cont)
    for(b in 1:B){
      
      n = simulate.sample(theta, theta.cont.vec[[cont]], N, tauvec, stressvec, ITvec,  seed = b, model.distribution = model.distribution)
    
      p=c(0,0,0)
      for(beta.init in c(0,0.4,0.8)){
        pilot = tryCatch(optimr(par= c(4, -0.02, 1), DPDloss,  N = N, tauvec = tauvec, stressvec = stressvec, ITvec = ITvec, n = n,
                            beta = beta.init, model.distribution = model.distribution, method ="BFGS"),
                     error=function(sol){sol$code=3;return(NA)})$par
        p = p + pilot/3
      }

      optimum = list("beta"=0.1, "estimator" = p) #IJW(n = n, pilot = pilot, model.distribution = model.distribution)
      #print(optimum)
         #
      estimators = estimate.function(N, tauvec, stressvec, ITvec, n, initial=c(5, -0.02, 1), model.distribution = model.distribution)
      if(sum(sapply(estimators, function(x) sum(is.na(x))))<1){ #none of the estimators is nan
        res.mat[cont,-nbetas] = res.mat[cont,-nbetas]+ unlist(lapply(estimators, RMSE.function, theta = theta))
        count = count+1
      }
      res.mat[cont,nbetas] = res.mat[cont,nbetas]+ RMSE.function(theta = theta, optimum[[2]])
    }
    res.mat[cont,] = res.mat[cont,]/count
    #print(count)
  }
  return(res.mat)
}
