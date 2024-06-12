

# needed packages

# packages needed
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
list.of.packages <- c("boot", "circular", "Renext", "pracma", "fitdistrplus","readr","survival",
                      "parallel","doParallel","compiler", "tidyverse")

ipak(list.of.packages)







########################################################################################## usefull functions Author: Aurélien Nicosia email:
########################################################################################## nicosia.aurelien@gmail.com Date: 28/07/2015


# packages needed
require("boot")
library("boot")
require("circular")
library("circular")
require("Renext")
library(Renext)
library(pracma)



#' Conditional likelihood function of the Step Selection Function 
#' 
#' Function that compute the conditional likelihood function of the Step Selection Function
#' at time step \code{id}
#' 
#'  
#' \describe{
#'   Compute the conditional likelihood function of the Step Selection Function
#'    }
#' 
#' 
#' @param beta  vector of parameters
#' @param X data at time step \code{id}
#' @param id time step           
#'  
#' @return \item{1/S}{ the coonditional likelihood } 
#' 
#' @references 
#' 
#' @author Aurélien Nicosia 


Lt <- function(beta, X, id) {
  S = 0
  for (i in (1:length(X[, 1]))) S = S + exp(beta %*% (X[i, 
                                                        ] - X[id, ]))
  return(1/S)
}



#' Conditional likelihood function of the Step Selection Function
#' 
#' Function that compute the intermediate function Qduring the E step in the EM algorithm.
#' 
#'  
#' \describe{
#'   Compute the intermediate function (Q, in the EM algorithm)
#'    }
#' 
#' 
#' @param beta  vector of parameters
#' @param E_smooth matrix T*K of smooth probabilities
#' @param df.matrix data frame as list(Case, strata, x, Times)         
#' @param k state          
#'  
#' @return \item{L}{ Intermediate function } 
#' 
#' @references 
#' 
#' @author Aurélien Nicosia 
#' 
LL <- function(beta, E_smooth, df.marix, k) {
  
  x <- as.matrix(df.marix$x)
  y <- as.vector(df.marix$y)
  strata <- as.factor(df.marix$strata)
  steps <- length(levels(strata))
  
  L <- 0
  
  for (t in (1:(steps))) {
    lev = levels(strata)[t]
    stratet <- which(strata == lev)
    id = which(y[stratet] == 1)
    S = 0
    X = x[stratet, ]
    for (i in (1:length(X[, 1]))) S = S + exp(beta %*% (X[i, 
                                                          ] - X[id, ]))
    LV = 1/S
    if (is.finite(log(S))) {
      L <- L + E_smooth[t, k] * log(S)
    }
  }
  L
}

############################# Filtering-Smoothing algorithm ###############


#' Filtering-Smoothing algorithm
#' 
#' Function that compute filtering and smoothing probabilities
#' 
#'  
#' \describe{
#'   During the M step of the EM algorithm, we need to filter and smooth
#'   the states.
#'    
#' 
#' @param param  list of parameters of the model 
#'  param <- list(beta =...,P = ...)
#' 
#' @param df.matrix data frame as list(Case, strata, x, Times)         
#'  
#' @return \item{E_smooth}{ Matrix T*K of the smooth probabilities \code{P(S_{tk}=1|F_{T})}} 
#' @return \item{E_smooth0}{ Vector K of the smooth probabilities \code{P(S_{0,k}=1|F_{T})}} 
#' 
#' @return \item{E_filter}{ Matrix T*K of the filter probabilities \code{P(S_{tk}=1|F_{t-1})} } 
#' @return \item{LogLike}{ Value of the log-likelihood function w.r.t the parameters \code{param} } 
#' @references 
#' 
#' @author Aurélien Nicosia 
FilterSmooth <- function(param, df.marix) {
  
  # initialisation:
  pi0 <- param$pi0
  P <- param$P
  beta <- as.matrix(param$beta)
  K <- length(beta[, 1])
  
  x <- as.matrix(df.marix$x)
  strate = df.marix$strata
  p <- length(x[1, ])
  
  y <- as.vector(df.marix$y)
  strata <- as.factor(df.marix$strata)
  lev = levels(strata)
  steps <- length(levels(strata))
  
  
  
  E_smooth0 <- rep(0, K)
  E_Filter <- matrix(0, steps, K)
  LogLike <- 0
  i <- 1
  P_trans <- as.vector(pi0 %*% P)
  levt = lev[1]
  strate1 <- which(strata == levt)
  id = which(y[strate1] == 1)
  
  f <- Lt(beta, X = x[strate1, ], id = id)
  f1 <- as.vector(exp(log(f) - max(log(f))))
  denom <- (P_trans) %*% f1
  E_Filter[1, ] <- as.vector((P_trans * f1))/denom
  
  if (is.na(log((P_trans) %*% (f))) == FALSE) {
    LogLike <- log((P_trans) %*% (f))
  }
  
  
  # 1. filtering Algorithm
  
  for (i in (2:steps)) {
    
    P_trans <- as.vector(E_Filter[i - 1, ] %*% P)
    levi = lev[i]
    stratei <- which(strata == levi)
    id = which(y[stratei] == 1)
    f <- Lt(beta, x[stratei, ], id)
    f1 <- as.vector(exp(log(f) - max(log(f))))
    denom <- (P_trans) %*% f1
    E_Filter[i, ] <- as.vector(P_trans * f1)/denom
    if (is.na(log((P_trans) %*% (f))) == FALSE) {
      LogLike <- LogLike + log((P_trans) %*% (f))
    }
  }
  
  # 2. Smoothing Algorithm
  E_smooth <- matrix(0, steps, K)
  # initialisation
  E_smooth[steps, ] <- E_Filter[steps, ]
  
  for (t in (steps - 1):1) {
    for (l in (1:K)) {
      for (k in (1:K)) {
        
        denom <- (E_Filter[t, ]) %*% P[, k]
        num <- (P[l, k] * E_Filter[t, l] * E_smooth[t + 
                                                      1, k])
        if (!is.na(num/denom)) {
          E_smooth[t, l] <- E_smooth[t, l] + num/denom
        }
      }
      
    }
  }
  for (l in (1:K)) {
    for (k in (1:K)) {
      
      denom0 <- pi0 %*% P[, k]
      num0 <- P[l, k] * pi0[l] * E_smooth[1, k]
      if (!is.na(num0/denom0)) {
        E_smooth0[l] <- E_smooth0[l] + num0/denom0
      }
    }
  }
  out <- list(E_Filter = E_Filter, E_smooth = E_smooth, E_smooth0 = E_smooth0, 
              LogLike = LogLike)
  class(out)
  out
  
}


############################# Log-Likelihood of the observed data ###############


#' Observed Log-likelihood function with a Markovian hidden process. 
#' 
#' Function that compute the observed log-likelihood given the assumption that the underlying hidden
#' process is Markovian.
#' 
#'  
#' \describe{
#'   
#'    
#' 
#' @param theta  vector of all the parameters of the model in this order: P and beta
#' 
#' @param df.matrix data frame as list(Case, strata, x, Times)         
#'
#' @param pi0 vector of initial distribution of the hidden process
#' 
#'  
#' @return \item{LogLike}{ Log-likelihood value } 
#' 
#' @details \code{LogLike} is the value of the observed log-likelihood function at parameter \code{theta} 
#' with the Markovian assumption of the hidden structure
#' @references 
#' 
#' @author Aurélien Nicosia 

logL <- function(theta, df.marix, pi0) {
  
  x <- as.matrix(df.marix$x)
  p = length(x[1, ])
  y <- as.vector(df.marix$y)
  strata <- as.factor(df.marix$strata)
  steps <- length(levels(strata))
  K <- length(pi0)
  beta <- NULL
  P <- NULL
  for (l in (1:K)) {
    P <- rbind(P, c(theta[(1 + (l - 1) * (K - 1)):(1 + (l - 
                                                          1) * (K - 1) + K - 2)], 1 - sum(theta[(1 + (l - 1) * 
                                                                                                   (K - 1)):(1 + (l - 1) * (K - 1) + K - 2)])))
    beta <- rbind(beta, theta[((K * (K - 1) + 1 + (l - 1) * 
                                  (p))):((K * (K - 1) + 1) + (l - 1) * (p) + (p - 1))])
  }
  E_smooth0 <- rep(0, K)
  E_Filter <- matrix(0, steps, K)
  LogLike <- 0
  i <- 1
  P_trans <- as.vector(pi0 %*% P)
  
  strate1 <- which(strata == levels(strata)[1])
  f <- Lt(beta, x[strate1, ], y[strate1])
  f1 <- as.vector(exp(log(f) - max(log(f))))
  denom <- (P_trans) %*% f1
  E_Filter[1, ] <- as.vector((P_trans * f1))/denom
  
  if (is.na(log((P_trans) %*% (f))) == FALSE) {
    LogLike <- log((P_trans) %*% (f))
  }
  
  
  # 1. filtering Algorithm
  
  for (i in (2:steps)) {
    
    P_trans <- as.vector(E_Filter[i - 1, ] %*% P)
    strate <- which(strata == levels(strata)[i])
    f <- Lt(beta, x[strate, ], y[strate])
    f1 <- as.vector(exp(log(f) - max(log(f))))
    denom <- (P_trans) %*% f1
    E_Filter[i, ] <- as.vector(P_trans * f1)/denom
    if (is.na(log((P_trans) %*% (f))) == FALSE) {
      LogLike <- LogLike + log((P_trans) %*% (f))
    }
  }
  
  
  return(LogLike)
}


############################# EM-algorithm ###############



#' Fit the HMM-SSF using initial parameter
#' 
#' Function that fit the HMM-SSF on data
#' using initial parameters given \code{param}
#' 
#'  
#' \describe{
#'   The function produce an estimation of the parameters of the HMM-SSF
#'    
#' @param formula formula that correspond to the conditional logistic model. The formula has to be in the same form
#' as the \code{coxph} function, i.e., of the form \code{case.status~explanatory_variable+strata(matched.set)  
#' 
#' @param dat data frame containing the variables in the model.
#' 
#' @param param  list of all the parameters of the model in this order: P and beta.
#' 
#' @param df.matrix data frame as list(Case, strata, x, Times)         
#' 
#' @param EMMax numbers of EM algorithm's maximum iteration
#' 
#' @param precision of the convergence, i.e. EM converge if the minimum
#' between two consecutives estimations is less than 10^(-\code{precision})  
#' 
#' @param EMMin numbers of EM algorithm's minimum iteration 
#' 
#' @return \item{beta}{ Matrix p*K of coefficient beta associated to the SSF model. } 
#' 
#' @return \item{LB}{ Value of the maximized log-likelihood function } 
#' @return \item{s}{ Numbers of EM algirithm's iteration} 
#' @return \item{P}{ estimated transition matrix} 
#' @return \item{pi0}{ initial probability distribution} 
#' @return \item{smooth}{ Matrix T*K of the smooth probabilities \code{P(S_{tk}=1|F_{T})}}
#' @return \item{thetat}{ history of estimated parameters in the EM algorithm} 
#' @references 
#' 
#' @details 
#' \code{param} should be list(pi0 =, P =, beta =
#' @author Aurélien Nicosia 


EstimSSFMultiState <- function(formula, dat, param, df.marix, EMMax, 
                               precision, EMMin) {
  
  
  
  # initialisation:
  pi0 <- param$pi0
  P <- param$P
  beta <- as.matrix(param$beta)
  K <- length(pi0)
  
  x <- as.matrix(df.marix$x)
  p <- length(x[1, ])
  
  y <- as.vector(df.marix$y)
  strata <- as.factor(df.marix$strata)
  steps <- length(levels(strata))
  
  thetat <- E_smooth <- E_Filter <- E_smooth0 <- paramop <- param1 <- param2 <- NULL
  LB <- 2
  LBa <- 1
  s <- 1
  
  for (i in (1:K)) {
    param1 <- c(param1, P[i, -K])
    param2 <- c(param2, (beta[i, ]))
  }
  paramopold <- c(param1, param2)
  paramopnew <- paramopold + 1
  
  # EM algorithm
  
  while (norm(as.matrix(paramopnew - paramopold), type = "M") > 
         10^(-precision[1]) && s < EMMin[1]) {
    
    paramopold <- paramopnew
    
    ########### E-step
    
    if (s > EMMax[1]) {
      break
    }
    
    # Posterior expectation
    FS <- FilterSmooth(param, df.marix)
    E_smooth0 <- FS$E_smooth0
    E_Filter <- FS$E_Filter
    E_smooth <- FS$E_smooth
    LBa <- FS$LogLike
    Pa <- param$P
    
    ############ M-step
    
    num <- matrix(0, K, K)
    P <- matrix(0, K, K)
    for (h in (1:K)) {
      denom <- sum(E_smooth[-steps, h]) + E_smooth0[h]
      
      for (k in (1:K)) {
        for (t in (2:steps)) {
          v <- E_Filter[t - 1, ] %*% Pa
          num[h, k] <- num[h, k] + E_smooth[t, k] * Pa[h, 
                                                       k] * E_Filter[t - 1, h]/v[k]
        }
        
        v <- pi0 %*% Pa
        num[h, k] <- num[h, k] + E_smooth[1, k] * Pa[h, 
                                                     k] * pi0[h]/v[k]
        P[h, k] <- num[h, k]/denom
        
      }
    }
    
    for (i in (1:K)) {
      W = NULL
      for (t in (1:length(E_smooth[, 1]))) {
        W = c(W, c(rep(E_smooth[t, i], sum(strata == 
                                             strata[1]))))
      }
      dat$W = W
      op <- coxph(formula, data = dat, weights = W)
      beta[i, ] <- op$coefficients
    }
    
    
    param <- list(pi0 = pi0, beta = beta, P = P)
    
    FSL <- FilterSmooth(param, df.marix)
    LB <- FSL$LogLike
    
    # check for spurious maximum
    
    a <- eigen(t(P))
    ind <- which.max(abs(a$values))
    pstatio <- as.matrix(abs(a$vectors[, ind])/norm(as.matrix(a$vectors[, 
                                                                        ind]), "o"))
    
    # look for a spurious maxima
    if (min(pstatio) < 1e-04 | max(abs(beta)) > 1000 | min(P) < 
        10^(-5)) {
      break
    }
    s <- s + 1
    
    
    paramop <- param1 <- param2 <- NULL
    
    for (i in (1:K)) {
      param1 <- c(param1, P[i, -K])
      param2 <- c(param2, t(beta[i, ]))
    }
    paramop <- c(param1, param2)
    paramopnew <- paramop
    thetat <- rbind(thetat, paramop)
  }
  
  FSL <- FilterSmooth(param, df.marix)
  LB <- FSL$LogLike
  out <- list(beta = beta, LB = LB, s = s, P = P, pi0 = pi0, 
              thetat = thetat, smooth = FSL$E_smooth)
  class(out)
  return(out)
  
  
}


############################# Global maximum of the likelihood function ###############



#' Fit the multi-state SSF model
#' 
#' Function that find the global maximum of the likelihood of
#'  the multi-state SSF model on data
#' 
#'  
#' \describe{
#'   The function produce an estimation of the parameters of the multi-state SSF model.
#'   The method use several starting point of the EM algorithm to rich the global maxima of the 
#'   likelihood function. 
#'    
#'    
#' 
#' 
#' @param formula formula that correspond to the conditional logistic model. The formula has to be in the same form
#' as the \code{coxph} function, i.e., of the form \code{case.status~explanatory_variable+strata(matched.set)  
#' 
#' @param dat data frame containing the variables in the model.
#' 
#' @param df.matrix data frame as list(Case, strata, x, Times)     
#' 
#' @param K the number of hidden states to consider.   
#' 
#' @param Initparam  list of all the initial parameters of the model in this order: Pinit and betainit to consider. 
#' NULL by default, in this case the function check for multiple initial conditions.
#' 
#' @param nb_init number of short EM algorithm.
#' 
#' @param EMMax vectors of size 2 of the numbers of EM algortihm's maximum iteration  (short, long)
#' 
#' @param precision vector of size 2 of the convergence of the precision of short and 
#'        long run algorithm, i.e. EM converge if the minimum
#'        between two consecutives estimations is less than 10^(-\code{precision})     
#' 
#' @param EMMin vectors of size 2 of the numbers of EM algortihm's minmum iteration  (short, long)
#' 
#' @param precision vector of size 2 of the convergence of the precision of short and 
#'        long run algorithm, i.e. EM converge if the minimum
#'        between two consecutives estimations is less than 10^(-\code{precision})                                   
#'  
#'
#' @return \item{Markov}{ List of estimation with a Markovian hidden structure } 
#' 
#' @return \item{EM.itermax}{ Number of EM algorithm iterations} 
#' 
#' @return \item{V}{ Estimated covariance matrix of the estimated coefficients} 
#' 
#' @return \item{theta}{ Vector of estimated parameters} 
#' 
#' @return \item{smooth}{ Matrix T*K of the smooth probabilities \code{P(S_{tk}=1|F_{T})}}
#' 
#' @references 
#' 
#' @details The list Markov is structured as:
#' \item{likelihood.Markov} Likelihood function.
#' \item{beta.Markov} matrix of the estimated parameters associated to Step Selection Function.
#' 
#' Concerning the hidden structure, depending on the nature of the hidden process, the output is
#' 
#' \item{P.Markov} estimated transition matrix with standard errors. Note that for each state k, only k-1 
#' probabilites are displayed, i.e. P(S_t=j|S_{t-1}=k), for j=1...K-1.
#' 
#' @author Aurélien Nicosia 


GlobalMaximaSSF <- function(formula, df.marix, dat, K, Initparam = NULL, 
                            nb_init = 30, EMMax = c(50, 2000), precision = c(2, 8), EMMin = c(10, 
                                                                                              1000)) {
  
  pi0 <- Initparam$pi0
  betainit <- Initparam$betainit
  Pinit <- Initparam$Pinit
  x <- as.matrix(df.marix$x)
  y <- as.vector(df.marix$y)
  strata <- as.factor(df.marix$strata)
  steps <- length(levels(strata))
  p <- length(x[1, ])
  # 1 etat
  
  
  if (K > 1) {
    
    # initialisation
    
    LBprec <- -exp(100)
    condition <- "FALSE"
    thetat1 <- thetat <- NULL
    
    if (is.null(pi0)) {
      w <- runif(K)
      pi0 <- w/sum(w)
      
    }
    
    j <- 1
    print("short-run EM algortim")
    while (j < nb_init) {
      
      if (is.null(betainit)) {
        beta <- NULL
        for (i in (1:K)) {
          conc <- runif(1, min = 2, max = 30)
          beta <- rbind(beta, c(conc, conc * runif(p, 
                                                   min = -1, max = 1)))
          
          
        }
      }
      if (!is.null(betainit)) {
        beta <- betainit + t(replicate(K, runif(p, -0.1, 
                                                0.1)))
      }
      if (is.null(Pinit)) {
        P <- NULL
        for (i in (1:K)) {
          vect <- runif(K)
          P <- rbind(P, vect/sum(vect))
        }
      }
      param <- list(beta = beta, P = P, pi0 = pi0)
      # Estimation
      S <- EstimSSFMultiState(formula, dat, param, df.marix, 
                              EMMax[1], precision[1], EMMin[1])
      out = S
      
      a <- eigen(t(P))
      LB <- out$LB
      ind <- which.max(abs(a$values))
      pstatio <- as.matrix(abs(a$vectors[, ind])/norm(as.matrix(a$vectors[, 
                                                                          ind]), "o"))
      
      if (min(pstatio) > 1e-04 && max(abs(out$beta)) < 
          1000 && min(P) > 10^(-5) && LB > LBprec) {
        condition <- "TRUE"
        betaop <- out$beta
        Pop <- out$P
        
        LBop <- out$LB
        thetat <- out$thetat  
        LBprec <- LBop
      }
      print(paste0(round(j/nb_init * 100), "% of the short-run EM algorithm done"))
      j <- j + 1
      
    }
    print(paste0(round(j/nb_init * 100), "% of the short-run EM algorithm done"))
    thetat1 <- thetat
    
    
    if (condition == "FALSE") {
      LB <- -exp(100)
      s <- 0
      print("Not enough initial parameters")
      break
    }
    if (condition == "TRUE") {
      print("long-run EM algorithm...")
      param <- list(beta = betaop, P = Pop, pi0 = pi0)
      S <- EstimSSFMultiState(formula, dat, param, df.marix, 
                              EMMax[2], precision[2], EMMin[2])
      
      
      out = S
      betaop <- out$beta
      Pop <- out$P
      s <- out$out
      LB <- out$LB
      thetat <- rbind(thetat1, out$thetat)
      theta1 <- theta2 <- NULL
      for (i in (1:K)) {
        theta1 <- c(theta1, Pop[i, -K])
        theta2 <- c(theta2, (betaop[i, ]))
      }
      theta <- c(theta1, theta2)
    }
    thetat <- thetat[-1, ]
    
    if (length(thetat) == 0) {
      print("need more initial condition")
      break
    }
    
    V3 <- optim(theta, fn = logL, gr = NULL, df.marix = df.marix, pi0 = param$pi0, 
                method = "L-BFGS-B", control = list(fnscale = -1), 
                hessian = TRUE, lower = theta - 10^-4 * rep(1, length(theta)), 
                upper = theta + 10^-4 * rep(1, length(theta)))
    if (det(V3$hessian) != 0) {
      V <- -(solve(V3$hessian))
    }
    
    if (det(V3$hessian) == 0) {
      V <- -pinv(V3$hessian)
    }
    
    parambetafinal <- NULL
    p = length(x[1, ])
    for (i in (1:K)) {
      parambeta <- cbind(betaop[i, ], sqrt(diag(V[(K * 
                                                     (K - 1) + 1 + (i - 1) * (p)):(K * (K - 1) + p + 
                                                                                     (i - 1) * (p)), (K * (K - 1) + 1 + (i - 1) * 
                                                                                                        (p)):(K * (K - 1) + p + (i - 1) * (p))])))
      colnames(parambeta) <- c(paste0("estimate ", "(k=", 
                                      i, ")"), paste0("s.e. ", "(k=", i, ")"))
      parambetafinal <- cbind(parambetafinal, parambeta)
      
    }
    
    paramPfinal <- NULL
    for (i in (1:K)) {
      paramP <- cbind(as.matrix(Pop[i, -K]), as.matrix(sqrt(diag(as.matrix(V[(1 + 
                                                                                (i - 1) * (K - 1)):((i) * (K - 1)), (1 + (i - 
                                                                                                                            1) * (K - 1)):((i) * (K - 1))])))))
      colnames(paramP) <- c(paste0("estimate ", "(P", i, 
                                   ".)"), paste0("s.e. ", "(P", i, ".)"))
      paramPfinal <- cbind(paramPfinal, paramP)
      
    }
    rownames(paramPfinal) <- paste0("P.", 1:(K - 1))
    
    Markov <- list(P.Markov = paramPfinal, beta.Markov = parambetafinal, 
                   likelihood.Markov = LB)
    out <- list(Markov = Markov, 
                EM.itermax = s, V = V, theta = theta, smooth =S$smooth )
    
  }
  class(out)
  out
  
}


