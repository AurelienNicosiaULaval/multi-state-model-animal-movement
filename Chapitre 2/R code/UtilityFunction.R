
#' Euclidean distance between two points
#' 
#' Calculate the euclidean distance between two points.
#'  
#' @param x0, y0  coordinates of one point.
#' @param x1, y1  coordinates of the other point.
#' 
#' @return The euclidean distance between the two points given as input. 
#' 
#' @author Sophie Baillargeon 
getDistance <- function(x0, y0, x1, y1) {
  as.vector(sqrt((x1 - x0)^2 + (y1 - y0)^2))
}


#' Angle between two points
#' 
#' Calculate the angle between a vector pointing North whose origin is the point (x0, y0)
#' and a vector whose origin is also point (x0, y0) and which contains point (x1, y1).
#' This angle consider a clockwise rotation from the vector pointing north.
#' Therefore, it is a positive number whose value is between 0 and 2*pi (units=radians).
#'  
#' @param x0, y0  coordinates of the point of origin.
#' @param x1, y1  coordinates of the other point.
#' 
#' @return The angle between the two points given as input, in radians. 
#' 
#' @author  Sophie Baillargeon
getAngle <- function(x0, y0, x1, y1) {
  at <- atan2(y1 - y0, x1 - x0)
  as.vector(at)
}






#' Weigthed von Mises distribution in the state \code{k}
#' 
#' Function that compute the weigthed von Mises distribution in the state \code{k}
#' with respect to the smooth probabilities \code{E.smooth}.
#' 
#'  
#' \describe{
#'   During the M step of the EM algorithm, we have to maximize the weigthed
#'    (by the smooth probabilities) von Mises distribution. 
#'    }
#' 
#' @param kappa  matrix of parameters associated to the von Mises distributions
#' 
#' @param E.smooth  matrix T*K of smooth probabilities
#' @param data data frame of the directions (y), distances (d), explanatory angles variables (x)
#'            and explanatory real variables (z)
#' @param k state           
#'  
#' @return \item{L}{ the weighted log-likelihood function related to the direction } 
#' 
#' @references 
#' 
#' @author Aurélien Nicosia 
LL.vonMises <- function(kappa, E.smooth, data, k) {
  x <- data$x
  y <- data$y
  z <- data$z
  steps <- length(y)
  y0 <- y[1]
  L <- 0
  t <- 1
  
  conc <- sqrt((kappa %*% (c(sin(y0), z[t, ] * sin(x[t, ]))))^2 + 
    (kappa %*% (c(cos(y0), z[t, ] * cos(x[t, ]))))^2)
  l.t0 <- E.smooth[t, k] * (-log(2 * pi * besselI(conc, 0, 
    expon.scaled = FALSE)) + kappa %*% (c(cos(y[t] - y0), 
    z[t, ] * cos(y[t] - x[t, ]))))
  if (!is.na(l.t0)) 
    L <- L + l.t0
  for (t in (2:(steps))) {
    conc <- sqrt((kappa %*% (c(sin(y[t - 1]), z[t, ] * sin(x[t, 
      ]))))^2 + (kappa %*% (c(cos(y[t - 1]), z[t, ] * cos(x[t, 
      ]))))^2)
    l.tk <- E.smooth[t, k] * (-log(2 * pi * besselI(conc, 
      0, expon.scaled = FALSE)) + kappa %*% (c(cos(y[t] - 
      y[t - 1]), z[t, ] * cos(y[t] - x[t, ]))))
    if (!is.na(l.tk)) 
      L <- L + l.tk
  }
  
  return(L)
}



#' Weigthed distribution of distance 
#' 
#' Function that compute the weigthed distribution of the distances 
#' 
#' 
#'  
#' \describe{
#'   During the M step of the EM algorithm, we have to maximize the weigthed
#'    (by the smooth probabilities) distribution of distances. 
#'    }
#' 
#' @param param  parameters c(log(shape),log(scale)) of the distances
#' 
#' @param d  vector of observed distances
#' @param weight vector of weight of each observation
#' @param dist type of distribution on the step length. The choices are 'gamma' (Default) or 'weibull'
#'  
#' @return \item{L}{ the weighted distribution related to the distance } 
#' 
#' @references 
#' 
#' @author Aurélien Nicosia 
l.dist = function(param, d, weight, dist = "gamma") {
  shape = exp(param[1])
  scale = exp(param[2])
  if (dist == "gamma") 
    L = sum(weight * dgamma(d, shape = shape, scale = scale, 
      log = TRUE))
  if (dist == "weibull") 
    L = sum(weight * dweibull(d, shape = shape, scale = scale, 
      log = TRUE))
  return(L)
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
#'  param <- list(pi0 = ..., kappa = ...,thetalin = ...,P = ...)
#' 
#' @param data data frame of the directions (y), distances (d), explanatory angles variables (x)
#'            and explanatory real variables (z)
#' @param type of the model: 'angle-dist' for bivariate model on direction-step length
#'  and 'angle' for univariate model on direction (Default).           
#'  
#' @param dist type of distribution on the step length. The choices are 'gamma' (Default) or 'weibull'
#'  
#' @return \item{E.smooth}{ Matrix T*K of the smooth probabilities \code{P(S_{tk}=1|F_{T})}} 
#' 
#' @return \item{E.filter}{ Matrix T*K of the filter probabilities \code{P(S_{tk}=1|F_{t-1})} } 
#' @return \item{LogLike}{ Value of the log-likelihood function w.r.t the parameters \code{param} } 
#' @references 
#' 
#' @author Aurélien Nicosia 
FilterSmooth <- function(param, data, type = "angle", dist = "gamma") {
  
  # initialization:
  pi0 <- param$pi0
  P <- param$P
  kappa <- param$kappa
  thetalin <- param$thetalin
  K <- length(pi0)
  p <- length(kappa[1, ]) - 1
  
  # information on the data set
  x <- as.matrix(data$x)
  y <- as.matrix(data$y)
  z <- as.matrix(data$z)
  d <- as.matrix(data$d)
  steps <- length(y)
  y0 <- y[1]
  E.smooth0 <- rep(0, K)
  E.filter <- matrix(0, steps, K)
  LogLike <- 0
  
  
  i <- 1
  P.trans <- as.vector(pi0 %*% P)
  l <- as.vector(sqrt((kappa %*% (c(sin(y0), z[i, ] * sin(x[i, 
    ]))))^2 + (kappa %*% (c(cos(y0), z[i, ] * cos(x[i, ]))))^2))
  f <- as.vector((1/(2 * pi * (besselI(l, 0, expon.scaled = FALSE)))) * 
    exp(kappa %*% (c(cos(y[i] - y0), z[i, ] * cos(y[i] - 
      x[i, ])))))
  g <- rep(1, length(f))
  if (type == "angle-dist" & dist == "gamma") 
    g <- dgamma(d[1], shape = thetalin[, 1], scale = thetalin[, 
      2])
  if (type == "angle-dist" & dist == "weibull") 
    g <- dweibull(d[1], shape = thetalin[, 1], scale = thetalin[, 
      2])
  f.renorm <- as.vector(exp(log(f * g) - max(log(f * g))))
  denom <- (P.trans) %*% f.renorm
  E.filter[1, ] <- as.vector((P.trans * f.renorm))/denom
  
  if (!is.na(log((P.trans) %*% (f * g)))) 
    LogLike <- log((P.trans) %*% (f * g))
  
  
  # 1. filtering Algorithm
  for (i in (2:steps)) {
    P.trans <- as.vector(E.filter[i - 1, ] %*% P)
    l <- as.vector(sqrt((kappa %*% (c(sin(y[i - 1]), z[i, 
      ] * sin(x[i, ]))))^2 + (kappa %*% (c(cos(y[i - 1]), 
      z[i, ] * cos(x[i, ]))))^2))
    f <- as.vector((1/(2 * pi * (besselI(l, 0, expon.scaled = FALSE)))) * 
      exp(kappa %*% (c(cos(y[i] - y[i - 1]), z[i, ] * cos(y[i] - 
        x[i, ])))))
    g <- rep(1, length(f))
    if (type == "angle-dist" & dist == "gamma") 
      g <- dgamma(d[i], shape = thetalin[, 1], scale = thetalin[, 
        2])
    if (type == "angle-dist" & dist == "weibull") 
      g <- dweibull(d[i], shape = thetalin[, 1], scale = thetalin[, 
        2])
    f.renorm <- as.vector(exp(log(f * g) - max(log(f * g))))
    denom <- (P.trans) %*% f.renorm
    E.filter[i, ] <- as.vector(P.trans * f.renorm)/denom
    if (!is.na(log((P.trans) %*% (f * g)))) 
      LogLike <- LogLike + log((P.trans) %*% (f * g))
  }
  
  # 2. Smoothing Algorithm
  E.smooth <- matrix(0, steps, K)
  E.smooth[steps, ] <- E.filter[steps, ]
  
  for (t in (steps - 1):1) {
    for (l in (1:K)) {
      for (k in (1:K)) {
        denom <- (E.filter[t, ]) %*% P[, k]
        num <- (P[l, k] * E.filter[t, l] * E.smooth[t + 
          1, k])
        if (!is.na(num/denom)) 
          E.smooth[t, l] <- E.smooth[t, l] + num/denom
      }
    }
  }
  for (l in (1:K)) {
    for (k in (1:K)) {
      denom0 <- pi0 %*% P[, k]
      num0 <- P[l, k] * pi0[l] * E.smooth[1, k]
      if (!is.na(num0/denom0)) 
        E.smooth0[l] <- E.smooth0[l] + num0/denom0
    }
  }
  out <- list(E.smooth0 = E.smooth0, E.filter = E.filter, E.smooth = E.smooth, 
    LogLike = LogLike)
  class(out) = "Filter-Smooth"
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
#' @param theta  vector of all the parameters of the model in this order: P, kappa and thetalin
#' 
#' @param data data frame of the directions (y), distances (d), explanatory angles variables (x)
#'            and explanatory real variables (z)
#' @param pi0 vector of initial distribution of the hidden process
#' 
#' @param nb.target number of targets                         
#'            
#' @param type of the model: 'angle-dist' for bivariate model on direction-step length
#'  and 'angle' for univariate model on direction (Default).           
#'  
#' @param dist type of distribution on the step length. The choices are 'gamma' (Default) or 'weibull'
#'  
#' @return \item{LogLike}{ Log-likelihood value } 
#' 
#' @details \code{LogLike} is the value of the observed log-likelihood function at parameter \code{theta} 
#' with the Markovian assumption of the hidden structure
#' @references 
#' 
#' @author Aurélien Nicosia 

logL <- function(theta, data, pi0, nb.target, type, dist) {
  x <- data$x
  y <- data$y
  z <- data$z
  d <- data$d
  steps <- length(y)
  y0 <- y[1]
  p <- nb.target
  K <- length(pi0)
  kappa <- NULL
  P <- NULL
  for (l in (1:K)) {
    P <- rbind(P, c(theta[(1 + (l - 1) * (K - 1)):(1 + (l - 
      1) * (K - 1) + K - 2)], 1 - sum(theta[(1 + (l - 1) * 
      (K - 1)):(1 + (l - 1) * (K - 1) + K - 2)])))
    kappa <- rbind(kappa, theta[((K * (K - 1) + 1 + (l - 
      1) * (p + 1))):((K * (K - 1) + 1) + (l - 1) * (p + 
      1) + (p))])
  }
  thetalin <- theta[c((length(theta) - K * 2 + 1):length(theta))]
  E.filter <- matrix(0, steps, K)
  LogLike <- 0
  x <- as.matrix(x)
  z <- as.matrix(z)
  cond <- "OK"
  i <- 1
  P.trans <- as.vector(pi0 %*% P)
  l <- as.vector(sqrt((kappa %*% (c(sin(y0), z[i, ] * sin(x[i, 
    ]))))^2 + (kappa %*% (c(cos(y0), z[i, ] * cos(x[i, ]))))^2))
  f <- as.vector((1/(2 * pi * (besselI(l, 0, expon.scaled = FALSE)))) * 
    exp(kappa %*% (c(cos(y[i] - y0), z[i, ] * cos(y[i] - 
      x[i, ])))))
  g <- rep(1, length(f))
  if (type == "angle-dist") {
    if (dist == "gamma") 
      g <- dgamma(d[1], shape = thetalin[c(1, 3)], scale = thetalin[c(2, 
        4)])
    if (dist == "weibull") 
      g <- dweibull(d[1], shape = thetalin[c(1, 3)], scale = thetalin[c(2, 
        4)])
  }
  f.renorm <- as.vector(exp(log(f * g) - max(log(f * g))))
  denom <- (P.trans) %*% f.renorm
  E.filter[1, ] <- as.vector((P.trans * f.renorm))/denom
  
  if (is.na(log((P.trans) %*% (f * g))) == FALSE) 
    LogLike <- log((P.trans) %*% (f * g))
  # 1. Algorithme filtering
  for (i in (2:steps)) {
    P.trans <- as.vector(E.filter[i - 1, ] %*% P)
    l <- as.vector(sqrt((kappa %*% (c(sin(y[i - 1]), z[i, 
      ] * sin(x[i, ]))))^2 + (kappa %*% (c(cos(y[i - 1]), 
      z[i, ] * cos(x[i, ]))))^2))
    f <- as.vector((1/(2 * pi * (besselI(l, 0, expon.scaled = FALSE)))) * 
      exp(kappa %*% (c(cos(y[i] - y[i - 1]), z[i, ] * cos(y[i] - 
        x[i, ])))))
    if (type == "angle-dist") {
      if (dist == "gamma") 
        g <- dgamma(d[i], shape = thetalin[c(1, 3)], 
          scale = thetalin[c(2, 4)])
      if (dist == "weibull") 
        g <- dweibull(d[i], shape = thetalin[c(1, 3)], 
          scale = thetalin[c(2, 4)])
    }
    f.renorm <- as.vector(exp(log(f * g) - max(log(f * g))))
    denom <- (P.trans) %*% f.renorm
    E.filter[i, ] <- as.vector(P.trans * f.renorm)/denom
    if (is.na(log((P.trans) %*% (f * g))) == FALSE) 
      LogLike <- LogLike + log((P.trans) %*% (f * g))
  }
  
  return(LogLike)
}



#' Transition matrix os semi_Markov process
#' 
#' Function that compute the transition matrix of the approximated semi-Markov process
#' 
#'
#' 
#' @param m  vector of the size of the Dwell times
#' 
#' @param p  vector of probabilities of the Dwell times.
#' @return \item{Gamma}{ Transition matrix } 
#' 
#' @references 
#' 
#' @author Roland Langrock 


gen.Gamma <- function(m, p) {
  Gamma <- diag(sum(m)) * 0
  ## state aggregate 1
  if (m[1] == 1) {
    Gamma[1, 1] <- 1 - dnbinom(0, size = p[1], prob = p[3])
    Gamma[1, 2] <- 1 - Gamma[1, 1]
  }
  if (m[1] > 1) {
    Gamma[1, m[1] + 1] <- dnbinom(0, size = p[1], prob = p[3])
    Gamma[1, 2] <- 1 - Gamma[1, m[1] + 1]
    for (i in 2:(m[1] - 1)) {
      cc <- rep(1, sum(m))
      for (k in 1:(i - 1)) {
        cc[k] <- Gamma[k, k + 1]
      }
      dd <- prod(cc)
      if (dd > 1e-12) 
        Gamma[i, m[1] + 1] <- dnbinom(i - 1, size = p[1], 
          prob = p[3])/dd
      if (dd < 1e-12) 
        Gamma[i, m[1] + 1] <- 1
      Gamma[i, i + 1] <- 1 - Gamma[i, m[1] + 1]
    }
    cc <- rep(1, sum(m))
    for (k in 1:(m[1] - 1)) {
      cc[k] <- Gamma[k, k + 1]
    }
    dd <- prod(cc)
    if (dd > 1e-12) 
      Gamma[m[1], m[1] + 1] <- dnbinom(m[1] - 1, size = p[1], 
        prob = p[3])/dd
    if (dd < 1e-12) 
      Gamma[m[1], m[1] + 1] <- 1
    Gamma[m[1], m[1]] <- 1 - Gamma[m[1], m[1] + 1]
  }
  ## state aggregate 2
  if (m[2] == 1) {
    Gamma[2, 2] <- 1 - dnbinom(0, size = p[2], prob = p[4])
    Gamma[2, 1] <- 1 - Gamma[2, 2]
  }
  if (m[2] > 1) {
    Gamma[m[1] + 1, 1] <- dnbinom(0, size = p[2], prob = p[4])
    Gamma[m[1] + 1, m[1] + 2] <- 1 - Gamma[m[1] + 1, 1]
    for (i in 2:(m[2] - 1)) {
      cc <- rep(1, sum(m))
      for (k in 1:(i - 1)) {
        cc[k] <- Gamma[m[1] + k, m[1] + k + 1]
      }
      dd <- prod(cc)
      if (dd > 1e-12) 
        Gamma[m[1] + i, 1] <- dnbinom(i - 1, size = p[2], 
          prob = p[4])/dd
      if (dd < 1e-12) 
        Gamma[m[1] + i, 1] <- 1
      Gamma[m[1] + i, m[1] + i + 1] <- 1 - Gamma[m[1] + 
        i, 1]
    }
    cc <- rep(1, sum(m))
    for (k in 1:(m[2] - 1)) {
      cc[k] <- Gamma[m[1] + k, m[1] + k + 1]
    }
    dd <- prod(cc)
    if (dd > 1e-12) 
      Gamma[m[1] + m[2], 1] <- dnbinom(m[2] - 1, size = p[2], 
        prob = p[4])/dd
    if (dd < 1e-12) 
      Gamma[m[1] + m[2], 1] <- 1
    Gamma[m[1] + m[2], m[1] + m[2]] <- 1 - Gamma[m[1] + m[2], 
      1]
  }
  Gamma
}



#' Observed Log-likelihood function with a semi-Markovian hidden process. 
#' 
#' Function that compute the observed log-likelihood given the assumption that the underlying hidden
#' process is semi-Markovian.
#' 
#'  
#' \describe{
#'   
#'    
#' 
#' @param theta  vector of all the parameters of the model in this order: P, kappa and thetalin
#' 
#' @param data data frame of the directions (y), distances (d), explanatory angles variables (x)
#'            and explanatory real variables (z)
#' @param pi0 vector of initial distribution of the hidden process
#' 
#' @param nb.target number of targets  
#' 
#' @param m vector or size K of the accurency of the approximation                       
#'            
#' @param type of the model: 'angle-dist' for bivariate model on direction-step length
#'  and 'angle' for univariate model on direction (Default).           
#'  
#' @param dist type of distribution on the step length. The choices are 'gamma' (Default) or 'weibull'
#'  
#' @return \item{LogLike}{ Log-likelihood value } 
#' 
#' @details \code{LogLike} is the value of the observed log-likelihood function at parameter \code{theta} 
#' with the Markovian assumption of the hidden structure
#' @references 
#' 
#' @author Aurélien Nicosia 


logLSemi <- function(theta, y0, data, pi0, nb.target, m = c(30, 
  30), type = "angle-dist", dist = "gamma") {
  if (m[1] == 1) 
    theta[1] = 1
  if (m[2] == 1) 
    theta[3] = 1
  
  
  theta[c(1, 3)] <- exp(theta[c(1, 3)])
  theta[c(2, 4)] <- inv.logit(theta[c(2, 4)])
  n <- length(theta)
  theta[c((n - 3):n)] <- exp(theta[c((n - 3):n)])
  p <- nb.target
  K <- length(pi0)
  x <- data$x
  y <- data$y
  z <- data$z
  d <- data$d
  steps <- length(y)
  p1 <- c(theta[1], theta[2])
  p2 <- c(theta[3], theta[4])
  kappa <- rbind(theta[5:(5 + (p))], theta[(5 + p + 1):(length(theta) - 
    4)])
  thetalin <- theta[c((length(theta) - 3):length(theta))]
  
  P <- gen.Gamma(m, cbind(c(p1[1], p2[1]), c(p1[2], p2[2])))
  kappacopy <- NULL
  thetalincopy <- NULL
  for (i in (1:m[1])) {
    kappacopy <- rbind(kappacopy, kappa[1, ])
    thetalincopy <- rbind(thetalincopy, thetalin[1:2])
  }
  for (i in (1:m[2])) {
    kappacopy <- rbind(kappacopy, kappa[2, ])
    thetalincopy <- rbind(thetalincopy, thetalin[3:4])
  }
  
  kappa <- kappacopy
  thetalin <- thetalincopy
  param = list(pi0 = pi0, kappa = kappacopy, thetalin = thetalincopy, 
    P = P)
  data <- list(x = x, y = y, z = z, d = d)
  type = "angle-dist"
  
  E.filter <- matrix(0, steps, K)
  LogLike <- 0
  x <- as.matrix(x)
  z <- as.matrix(z)
  cond <- "OK"
  # initialization
  i <- 1
  P.trans <- as.vector(pi0 %*% P)
  l <- as.vector(sqrt((kappa %*% (c(sin(y0), z[i, ] * sin(x[i, 
    ]))))^2 + (kappa %*% (c(cos(y0), z[i, ] * cos(x[i, ]))))^2))
  f <- as.vector((1/(2 * pi * (besselI(l, 0, expon.scaled = FALSE)))) * 
    exp(kappa %*% (c(cos(y[i] - y0), z[i, ] * cos(y[i] - 
      x[i, ])))))
  if (type == "angle-dist") {
    if (dist == "gamma") 
      g <- dgamma(d[1], shape = thetalin[, 1], scale = thetalin[, 
        2])
    if (dist == "weibull") 
      g <- dweibull(d[1], shape = thetalin[, 1], scale = thetalin[, 
        2])
  }
  f.renorm <- as.vector(exp(log(f * g) - max(log(f * g))))
  denom <- (P.trans) %*% f.renorm
  E.filter[1, ] <- as.vector((P.trans * f.renorm))/denom
  
  if (is.na(log((P.trans) %*% (f * g))) == FALSE) 
    LogLike <- log((P.trans) %*% (f * g))
  # 1. Filtering algorithm
  
  for (i in (2:steps)) {
    P.trans <- as.vector(E.filter[i - 1, ] %*% P)
    l <- as.vector(sqrt((kappa %*% (c(sin(y[i - 1]), z[i, 
      ] * sin(x[i, ]))))^2 + (kappa %*% (c(cos(y[i - 1]), 
      z[i, ] * cos(x[i, ]))))^2))
    f <- as.vector((1/(2 * pi * (besselI(l, 0, expon.scaled = FALSE)))) * 
      exp(kappa %*% (c(cos(y[i] - y[i - 1]), z[i, ] * cos(y[i] - 
        x[i, ])))))
    if (type == "angle-dist") {
      if (dist == "gamma") 
        g <- dgamma(d[i], shape = thetalin[, 1], scale = thetalin[, 
          2])
      if (dist == "weibull") 
        g <- dweibull(d[i], shape = thetalin[, 1], scale = thetalin[, 
          2])
    }
    f.renorm <- as.vector(exp(log(f * g) - max(log(f * g))))
    denom <- (P.trans) %*% f.renorm
    E.filter[i, ] <- as.vector(P.trans * f.renorm)/denom
    if (is.na(log((P.trans) %*% (f * g))) == FALSE) 
      LogLike <- LogLike + log((P.trans) %*% (f * g))
  }
  return(LogLike)
}



############################# EM-algorithm ###############


#' Fit the general hidden random walf model using initial parameter
#' 
#' Function that fit the General hidden random walk model on data
#' using initial parameters given \code{param}
#' 
#'  
#' \describe{
#'   The function produce an estimation of the parameters of the General Hidden Random Walk model
#'    
#'    
#' 
#' @param param  list of all the parameters of the model in this order: P, kappa and thetalin.
#' 
#' @param data data frame of the directions (y), distances (d), explanatory angles variables (x)
#'            and explanatory real variables (z)
#' @param EMMax numbers of EM algortihm's maximum iteration
#' 
#' @param EMMin numbers of EM algortihm's minimum iteration 
#' 
#' @param precision of the convergence, i.e. EM converge if the minimum
#' between two consecutives estimations is less than 10^(-\code{precision})                                   
#'  
#' @param type of the model: 'angle-dist' for bivariate model on direction-step length
#'  and 'angle' for univariate model on direction (Default).  
#'  
#' @param dist type of distribution on the step length. The choices are 'gamma' (Default) or 'weibull'
#'  
#' @return \item{beta}{ Matrix p*K of normalized coefficient beta associated to the targets. } 
#' 
#' @return \item{LL}{ Value of the maximized log-likelihood function } 
#' @return \item{s}{ Numbers of EM algirithm's iteration} 
#' @return \item{P}{ estimated transition matrix} 
#' @return \item{pi0}{ initial probability distribution} 
#' @return \item{kappa}{ matrix (p+1)*K of the parameters associated to the direction} 
#' @return \item{thetalin}{ matrix 2*K of the parameters associated to the distance} 
#' @return \item{thetat}{ history of estimated parameters in the EM algorithm} 
#' @references 
#' 
#' @details \code{p} reprensents the number of targets. 
#' The vector beta_k in the state \code{k} represents kappa[k, -1]/kappa[k, 1].
#' \code{param} should be list(pi0 =, P =, kappa =,  thetalin=, )
#' 
#' @author Aurélien Nicosia 



fitGHRW <- function(param, data, EMMax = 50, precision = 2, EMMin = 10, 
  type = "angle-dist", dist = "gamma") {
  
  # initialization:
  pi0 <- param$pi0
  P <- param$P
  kappa <- as.matrix(param$kappa)
  thetalin <- param$thetalin
  K <- length(pi0)
  p <- length(kappa[1, ]) - 1
  
  x <- as.matrix(data$x)
  y <- as.vector(data$y)
  z <- data$z
  d <- data$d
  steps <- length(y)
  y0 <- y[1]
  
  thetat <- E.smooth <- E.filter <- E.smooth0 <- paramop <- param1 <- param2 <- NULL
  LL <- 2
  LL.prev <- 1
  s <- 1
  
  paramop <- param1 <- param2 <- param3 <- NULL
  
  for (i in seq(along = pi0)) {
    param1 <- c(param1, P[i, -K])
    param2 <- c(param2, t(kappa[i, ]))
    param3 <- c(param3, t(thetalin[i, ]))
  }
  paramopold <- c(param1, param2, param3)
  paramopnew <- paramopold + 1
  
  # Em algorithm
  
  while (norm(as.matrix(paramopnew - paramopold), type = "M") > 
    10^(-precision) && s < EMMin) {
    
    paramopold <- paramopnew
    ########### E-step
    if (s > EMMax) {
      break
    }
    
    FS <- FilterSmooth(param, data, type, dist)
    E.smooth0 <- FS$E.smooth0
    E.filter <- FS$E.filter
    E.smooth <- FS$E.smooth
    LL.prev <- FS$LogLike
    Pa <- param$P
    
    ############ M-step
    
    num <- matrix(0, K, K)
    P <- matrix(0, K, K)
    for (h in seq(along = pi0)) {
      denom <- sum(E.smooth[-steps, h]) + E.smooth0[h]
      for (k in seq(along = pi0)) {
        for (t in c(2:length(y[]))) {
          v <- E.filter[t - 1, ] %*% Pa
          num[h, k] <- num[h, k] + E.smooth[t, k] * Pa[h, 
          k] * E.filter[t - 1, h]/v[k]
        }
        
        v <- pi0 %*% Pa
        num[h, k] <- num[h, k] + E.smooth[1, k] * Pa[h, 
          k] * pi0[h]/v[k]
        P[h, k] <- num[h, k]/denom
      }
    }
    
    for (i in seq(along = pi0)) {
      op <- optim(param$kappa[i, ], LL.vonMises, E.smooth = E.smooth, 
        data = data, k = i, control = list(fnscale = -1), 
        method = "L-BFGS-B")
      kappa[i, ] <- op$par
      if (type == "angle-dist") 
        op2 <- optim(log(param$thetalin[i, ]), l.dist, 
          d = d, weight = E.smooth[, i], dist = dist, 
          control = list(fnscale = -1), method = "L-BFGS-B")
      thetalin[i, ] <- exp(op2$par)
    }
    
    
    ind <- order(-kappa[, 1])
    kappa <- kappa[ind, ]
    thetalin <- thetalin[ind, ]
    P <- P[ind, ind]
    
    param <- list(pi0 = pi0, kappa = kappa, thetalin = thetalin, 
      P = P)
    FSL <- FilterSmooth(param, data, type, dist)
    LL <- FSL$LogLike
    
    # check for spurious maximum
    beta <- kappa[, 2:(p + 1)]/kappa[, 1]
    a <- eigen(t(P))
    ind <- which.max(abs(a$values))
    pstatio <- as.matrix(abs(a$vectors[, ind])/norm(as.matrix(a$vectors[, 
      ind]), "o"))
    critkappa <- max(abs(kappa[, 1]))
    # look for a spurious maxima
    if (min(pstatio) < 1e-04 | critkappa > 100 | max(abs(beta)) > 
      100 | min(P) < 10^(-5) | min(thetalin) < 0.01 | max(thetalin) > 
      20) 
      break
    
    s <- s + 1
    
    paramop <- param1 <- param2 <- param3 <- NULL
    
    for (i in seq(along = pi0)) {
      param1 <- c(param1, P[i, -K])
      param2 <- c(param2, t(kappa[i, ]))
      param3 <- c(param3, t(thetalin[i, ]))
    }
    paramop <- c(param1, param2, param3)
    paramopnew <- paramop
    thetat <- rbind(thetat, paramop)
  }
  
  FSL <- FilterSmooth(param, data, type, dist)
  LL <- FSL$LogLike
  out <- list(beta = beta, LL = LL, s = s, P = P, pi0 = pi0, 
    kappa = kappa, thetalin = thetalin, thetat = thetat)
  class(out) = "fitGHRW"
  out
  
}









############################# Global maximum of the likelihood function ###############


#' Fit the general hidden random walk model
#' 
#' Function that find the global maximum of the likelihood of
#'  the General hidden random walk model on data
#' 
#'  
#' \describe{
#'   The function produce an estimation of the parameters of the General Hidden Random Walk model.
#'   The method use several starting point of the EM algorithm to rich the global maxima of the 
#'   likelihood function. 
#'    
#'    
#' 
#' 
#' @param data data frame of the directions (y), distances (d), explanatory angles variables (x)
#'            and explanatory real variables (z)
#' @param EMMax vectors of size 2 of the numbers of EM algortihm's maximum iteration  (short, long)
#' 
#' @param EMMin vectors of size 2 of the numbers of EM algortihm's minmum iteration  (short, long)
#' 
#' @param precision vector of size 2 of the convergence of the precision of short and 
#'        long run algorithm, i.e. EM converge if the minimum
#'        between two consecutives estimations is less than 10^(-\code{precision})                                   
#'  
#' @param type of the model: 'angle-dist' for bivariate model on direction-step length
#'  and 'angle' for univariate model on direction (Default).  
#'  
#' @param dist type of distribution on the step length. The choices are 'gamma' (Default) or 'weibull'
#'  
#' @param semi indicates TRUE for the semi-Markov estimation. Default is FALSE.
#'
#' @return \item{Markov}{ List of estimation with a Markovian hidden structure } 
#' 
#' @return \item{SemiMarkov}{ List of estimation with a semi-Markovian hidden structure  } 
#' @return \item{EM.itermax}{ Number of EM algorithm iterations} 

#' @references 
#' 
#' @details The list Markov and SemiMarkov output have the same structure, as:
#' \item{fit} of the model with \code{likelihood}, \code{AIC}, \code{BIC}
#' \item{kappa} matrix of the estimated parameters associated to the direction with standard errors.
#' \item{lambda} matrix of the estimated parameters associated to the distance with standard errors.
#' 
#' Concerning the hidden structure, depending on the nature of the hidden process, the output is
#' 
#' \item{P.Markov} estimated transition matrix with standard errors. Note that for each state k, only k-1 
#' probabilites are displayed, i.e. P(S_t=j|S_{t-1}=k), for j=1...K-1.
#' \item{Dwell.SemiMarkov} estimated parameters of the Dwell time in each state of the semi_Markov
#' process
#' 
#' @author Aurélien Nicosia 

GHRandomWalk <- function(data, K, paraminitiaux = NULL, nb_init = 50, 
  EMMax = c(50, 2000), precision = c(2, 8), EMMin = c(10, 1000), 
  type = "angle-dist", dist = "gamma", semi = FALSE) {
  
  pi0 <- paraminitiaux$pi0
  kappainit <- paraminitiaux$kappainit
  thetalininit <- paraminitiaux$thetalininit
  L = length(thetalininit[1, ])
  Pinit <- paraminitiaux$Pinit
  
  x <- as.matrix(data$x)
  y <- as.vector(data$y)
  z <- data$z
  d <- data$d
  steps <- length(y)
  y0 <- y[1]
  
  # 1 etat
  if (K == 1) 
    print("please use the consensus function")
  if (K > 1) {
    # initialization
    LLprec <- -exp(100)
    condition <- FALSE
    thetat1 <- thetat <- NULL
    if (is.null(pi0)) {
      w <- runif(K)
      pi0 <- w/sum(w)
    }
    # Estimation for nb_init initial parameters
    print("short-run EM algorithm")
    j <- 1
    while (j < nb_init) {
      if (is.null(kappainit)) {
        kappa <- NULL
        thetalin <- NULL
        for (i in seq(along = pi0)) {
          conc <- runif(1, min = 2, max = 30)
          kappa <- rbind(kappa, c(conc, conc * runif(p, 
          min = -1, max = 1)))
        }
      }
      if (type == "angle-dist") {
        if (is.null(thetalininit)) 
          thetalin <- matrix(runif(K * L, min = 1, max = 10), 
          K, L)
      }
      if (!is.null(kappainit)) 
        kappa <- kappainit + t(replicate(K, runif(p + 
          1, -1, 1)))
      
      if (!is.null(thetalininit)) {
        if (type == "angle-dist") 
          thetalin <- thetalininit + matrix(runif(K * 
          L, min = 1, max = 5), K, L)
      }
      if (is.null(Pinit)) {
        P <- NULL
        for (i in seq(along = pi0)) {
          vect <- runif(K)
          P <- rbind(P, vect/sum(vect))
        }
      }
      param <- list(kappa = kappa, thetalin = thetalin, 
        P = P, pi0 = pi0)
      
      # Estimation
      S <- fitGHRW(param, data, EMMax[1], precision[1], 
        EMMin[1], type, dist)
      P <- S$P
      a <- eigen(t(P))
      LL <- S$LL
      ind <- which.max(abs(a$values))
      pstatio <- as.matrix(abs(a$vectors[, ind])/norm(as.matrix(a$vectors[, 
        ind]), "o"))
      critkappa <- min(abs(S$kappa[, 1]))
      
      # check for False maxima
      if (min(pstatio) > 1e-04 && critkappa < 100 && max(abs(S$beta)) < 
        100 && min(P) > 10^(-5) && min(thetalin) > 0.01 && 
        LL > LLprec) {
        condition <- TRUE
        kappaop <- S$kappa
        Pop <- S$P
        thetalinop <- S$thetalin
        LLop <- S$LL
        thetat <- S$thetat
        LLprec <- LLop
      }
      print(paste0(round(j/nb_init * 100), "% of the short-run EM algorithm done"))
      j <- j + 1
      
    }
    thetat1 <- thetat
    
    # Long-run EM algorithm
    
    if (condition == "FALSE") {
      LL <- -exp(100)
      s <- 0
      print("Not enough initial parameters! Please try again.")
      break
    }
    if (condition == "TRUE") {
      print("long-run EM algorithm...")
      param <- list(kappa = kappaop, thetalin = thetalinop, 
        P = Pop, pi0 = pi0)
      S <- fitGHRW(param, data, EMMax[2], precision[2], 
        EMMin[2], type, dist)
      kappaop <- S$kappa
      Pop <- S$P
      thetalinop <- S$thetalin
      beta <- S$beta
      s <- S$s
      LL <- S$LL
      param.Markov <- list(kappa = kappaop, thetalin = thetalinop, 
        P = Pop, pi0 = pi0)
      
      thetat <- rbind(thetat1, S$thetat)
      theta1 <- theta2 <- theta3 <- NULL
      for (i in seq(along = pi0)) {
        theta1 <- c(theta1, Pop[i, -K])
        theta2 <- c(theta2, (kappaop[i, ]))
        theta3 <- c(theta3, thetalinop[i, ])
      }
      theta <- c(theta1, theta2, theta3)
      theta.Markov = theta
      nb.target <- length(kappaop[1, ]) - 1
    }
    thetat <- thetat[-1, ]
    
    if (length(thetat) == 0) {
      print("need more initial condition")
      break
    }
    V3 <- optim(theta, fn = logL, gr = NULL, data = data, 
      pi0 = param$pi0, nb.target = nb.target, type = type, 
      dist = dist, method = "L-BFGS-B", control = list(fnscale = -1), 
      hessian = TRUE, lower = theta - 10^-4 * rep(1, length(theta)), 
      upper = theta + 10^-4 * rep(1, length(theta)))
    if (det(V3$hessian) != 0) 
      V <- -(solve(V3$hessian))
    if (det(V3$hessian) == 0) 
      V <- -pinv(V3$hessian)
    
    paramkappafinal <- NULL
    for (i in seq(along = pi0)) {
      se = sqrt(diag(V[(K * (K - 1) + 1 + (i - 1) * (p + 
        1)):(K * (K - 1) + 1 + p + (i - 1) * (p + 1)), 
        (K * (K - 1) + 1 + (i - 1) * (p + 1)):(K * (K - 
          1) + 1 + p + (i - 1) * (p + 1))]))
      paramkappa <- cbind(kappaop[i, ], se, round(1 - pnorm(abs(kappaop[i, 
        ]/se)), 4))
      colnames(paramkappa) <- c(paste0("estimate ", "(k=", 
        i, ")"), paste0("s.e. ", "(k=", i, ")"), paste0("p.value ", 
        "(k=", i, ")"))
      paramkappafinal <- cbind(paramkappafinal, paramkappa)
    }
    paramPfinal <- NULL
    for (i in seq(along = pi0)) {
      paramP <- cbind(as.matrix(Pop[i, -K]), as.matrix(sqrt(diag(as.matrix(V[(1 + 
        (i - 1) * (K - 1)):(1 + (i - 1) * (K - 1)), (1 + 
        (i - 1) * (K - 1)):(1 + (i - 1) * (K - 1))])))))
      colnames(paramP) <- c(paste0("estimate ", "(P", i, 
        ".)"), paste0("s.e. ", "(P", i, ".)"))
      paramPfinal <- cbind(paramPfinal, paramP)
    }
    rownames(paramPfinal) <- paste0("P.", 1:(K - 1))
    
    paramdistfinal <- NULL
    V.dist = diag(as.matrix(V[c((length(theta) - 2 * K + 
      1):length(theta)), c((length(theta) - 2 * K + 1):length(theta))]))
    for (i in seq(along = pi0)) {
      paramdist <- cbind(as.matrix(thetalinop[i, ]), as.matrix(sqrt(V.dist[c(1 * 
        (i == 1) + 3 * (i == 2), 2 * (i == 1) + 4 * (i == 
        2))])))
      colnames(paramdist) <- c(paste0("estimate ", "(k=", 
        i, ")"), paste0("s.e. ", "(k=", i, ")"))
      paramdistfinal <- cbind(paramdistfinal, paramdist)
      
    }
    
    # semi-Markov result for a 2-states
    if (K > 2 & semi) 
      cat("Semi-Markovian approximation with more than 2 state is too much demanding... please select K=2")
    if (K == 2 & semi) {
      p1op <- c(log(1), logit(Pop[1, 2]))
      p2op <- c(log(1), logit(Pop[2, 1]))
      theta <- c(p1op, p2op, kappaop[1, ], kappaop[2, ], 
        log(thetalinop[1, ]), log(thetalinop[2, ]))
      nb.target <- length(kappaop[1, ]) - 1
      m <- c(30, 30)
      w <- runif(sum(m))
      pi0 <- w/sum(w)
      minL = -exp(16)
      for (l in (1:nb_init)) {
        thetaop = theta + runif(length(theta), -1, 1)
        V3 <- optim(thetaop, fn = logLSemi, gr = NULL, 
          y0 = y0, data = data, pi0 = pi0, nb.target = nb.target, 
          m = m, dist = dist, method = "L-BFGS-B", control = list(fnscale = -1, 
          maxit = 10^(8)), hessian = TRUE)
        if (V3$value > minL) {
          minL = V3$value
          V2 = V3
        }
      }
      V <- sqrt(diag(-(solve(V2$hessian))))
      
      p1op <- c(exp(V2$par[1]), inv.logit(V2$par[2]))
      sp1op <- c(p1op[1] * V[1], p1op[2] * (1 - p1op[2]) * 
        V[2])
      p2op <- c(exp(V2$par[3]), inv.logit(V2$par[4]))
      sp2op <- c(p2op[1] * V[3], p2op[2] * (1 - p2op[2]) * 
        V[4])
      
      kappaop <- V2$par[c(5:((5 + 2 * nb.target + 1)))]
      skappaop <- V[c(5:((5 + 2 * nb.target + 1)))]
      n <- length(V)
      thetalinop <- exp(V2$par[c((n - 3):n)])
      sthetalinop <- thetalinop * V[c((n - 3):n)]
      paramkappafinalS <- NULL
      paramkappafinalS <- cbind(kappaop[1:(nb.target + 
        1)], skappaop[1:(nb.target + 1)], kappaop[(nb.target + 
        2):(2 * (nb.target + 1))], skappaop[(nb.target + 
        2):(2 * (nb.target + 1))])
      
      colnames(paramkappafinalS) <- c("estimate (k=1)", 
        "s.e.", "estimate (k=2)", "s.e.")
      
      Dwell <- rbind(c(p1op[1], sp1op[1], p1op[2], sp1op[2]), 
        c(p2op[1], sp2op[1], p2op[2], sp2op[2]))
      colnames(Dwell) <- c("size", "s.e", "prob", "s.e")
      rownames(Dwell) <- c("p1", "p2")
      
      paramdistfinalS <- NULL
      paramdistfinalS <- cbind(thetalinop[c(1, 2)], sthetalinop[c(1, 
        2)], thetalinop[c(3, 4)], sthetalinop[c(3, 4)])
      colnames(paramdistfinalS) <- c("lambda (k=1)", "s.e", 
        "lambda (k=2)", "s.e")
      
      P <- gen.Gamma(m, cbind(c(p1op[1], p2op[1]), c(p1op[2], 
        p2op[2])))
      Pa <- P
      kappaop <- rbind(t(paramkappafinalS[, 1]), t(paramkappafinalS[, 
        3]))
      
      kappacopy <- NULL
      thetalincopy <- NULL
      for (i in (1:m[1])) {
        kappacopy <- rbind(kappacopy, kappaop[1, ])
        thetalincopy <- rbind(thetalincopy, thetalinop[1:2])
      }
      for (i in (1:m[2])) {
        kappacopy <- rbind(kappacopy, kappaop[2, ])
        thetalincopy <- rbind(thetalincopy, thetalinop[3:4])
      }
      param = list(pi0 = pi0, kappa = kappacopy, thetalin = thetalincopy, 
        P = P)
      data <- list(x = x, y = y, z = z, d = d)
      FS <- FilterSmooth(param, data, type, dist)
      
      E.smooth0 <- FS$E.smooth0
      E.filter <- FS$E.filter
      E.smooth <- FS$E.smooth
      LL.prev <- FS$LogLike
    }
    
    
    fit.Markov = list(likelihood.Markov = LL, AIC.Markov = -2 * 
      LL + 2 * length(theta.Markov), BIC.Markov = -2 * 
      LL + log(length(y)) * length(theta.Markov))
    if (semi) 
      fit.SemiMarkov = list(likelihood.SemiMarkov = LL.prev, 
        AIC.SemiMarkov = -2 * LL.prev + 2 * length(theta), 
        BIC.SemiMarkov = -2 * LL.prev + log(length(y)) * 
          length(theta))
    
    Markov <- list(param = param.Markov, P.Markov = paramPfinal, 
      kappa.Markov = paramkappafinal, lambda.Markov = paramdistfinal, 
      fit.Markov = fit.Markov)
    if (semi) 
      SemiMarkov <- list(Dwell.SemiMarkov = Dwell, kappa.SemiMarkov = paramkappafinalS, 
        lambda.SemiMarkov = paramdistfinalS, fit.SemiMarkov = fit.SemiMarkov)
    out <- list(Markov = Markov, param = param, EM.itermax = s)
    if (semi) 
      out <- list(Markov = Markov, param = param, SemiMarkov = SemiMarkov, 
        EM.itermax = s)
  }
  class(out) <- "GHRandomWalk"
  out
}


#' @param x An object, produced by the \code{\link{consensus}} function, to print.
#' @param \dots Further arguments to be passed to \code{print.default}. 
print.GlobalMaxima <- function(x, ...) {
  cat("\nMarkov specification:")
  cat("\nMax log-likelihood and classical criteria:", x$fit.Markov, 
    "\n")
  cat("\nParameters of directions :\n")
  print.default(x$fit.Markov$kappa.Markov, print.gap = 2, quote = FALSE, 
    right = TRUE, ...)
  cat("\n")
  cat("\nParameters of distances :\n")
  print.default(x$fit.Markov$lambda.Markov, print.gap = 2, 
    quote = FALSE, right = TRUE, ...)
  cat("\n")
  cat("\nTransition matrix :\n")
  print.default(x$fit.Markov$P.Markov, print.gap = 2, quote = FALSE, 
    right = TRUE, ...)
  cat("\n")
  invisible(x)
}


####################### consensus 1 state
####################### ################################################'



#' Consensus Model
#' 
#' Function that fit a consensus model for angular variables and its \code{print} method.
#' 
#' The control argument is a list that can supply any of the following components: 
#' \describe{
#'   \item{\code{pginit}}{ The approximate number of points on the grid of possible initial beta values tried when 
#'                         \code{initbeta} is not given. The default is 1000, which runs quickly. A large value of 
#'                         \code{pginit} makes the function slower.  }
#'   \item{\code{maxiter}}{ The maximum number of iterations. The default is 1000.  }
#'   \item{\code{mindiff}}{ The minimum difference between two max cosine to be reached. It defines the convergence criteria:
#'                          if the difference between the max cosine for the updated parameters values and the max
#'                          cosine for the parameters values at the previous iteration is below \code{mindiff}, convergence is
#'                          reached. The default is 0.000001.   }
#' }
#' 
#' @param formula  A formula with the dependent angle on the left of the ~ operator and terms specifying
#'                 the explanatory variables on the right. These terms must be written \code{x:z}, where
#'                 \code{x} is an explanatory angle which relative importance migth depend on the
#'                 positive variable \code{z}. It is not mandatory to specify a \code{z} variable for each
#'                 explanatory angle. For \code{model='simplified'}, the first explanatory angle listed is 
#'                 the reference direction (if a \code{z} variable was specified for this angle, it is ignored).
#' @param data  An optional data frame, list or environment 
#'              containing the variables in the model formula. If not found in data, the variables are taken from 
#'              \code{environment(formula)}, typically the environment from which \code{angular} is called. 
#' @param model  A character string, either \code{'complete'} for the complete model with an intercept (the default) 
#'               or \code{'simplified'} for the simplified model without an intercept.
#' @param initparam  A numerical vector, initial values for the parameters. The default is to use the best initial
#'                   values found among some values tried on a grid of possible values for the parameters.
#' @param control A list of control parameters. See Details.
#' 
#' @return \item{MaxLL}{ the maximum value of the log likeliohood } 
#' @return \item{parameters}{ the parameter estimates and their standard errors (obtained from two definitions) } 
#' @return \item{varcov1}{ the estimated variance covariance matrix for the parameter estimates (obtained from the first definition) }
#' @return \item{varcov2}{ the estimated variance covariance matrix for the parameter estimates (obtained from the second definition) }
#' @return \item{parambeta}{ the beta parameter estimates and their standard errors (obtained by linearization) } 
#' @return \item{varcovbeta1}{ the estimated variance covariance matrix for the beta parameter estimates (obtained by linearization) }
#' @return \item{varcovbeta2}{ the estimated variance covariance matrix for the beta parameter estimates (Sandwhich form) }
#' @return \item{autocorr}{ the autocorrelation of the residuals \eqn{\sin(y_i-\mu_i)}{sin(yi-mui)}  }
#' @return \item{iter.detail}{ the iteration details }
#' @return \item{converge}{ an indicator of convergence }
#' @return \item{call}{ the function call }
#' 
#' @author Sophie Baillargeon, Louis-Paul Rivest and Aurélien Nicosia
consensus <- function(formula, data, model = "simplified", weights = NULL, 
  initbeta = NULL, control = list()) {
  
  call <- mfcall <- match.call()
  model <- model[1]
  
  
  
  ### information from formula and data
  mfargs <- match(c("formula", "data"), names(mfcall), 0L)
  mfcall <- mfcall[c(1L, mfargs)]
  mfcall[[1L]] <- as.name("model.frame")
  mf <- eval(mfcall, parent.frame())
  # useful objects
  nobs <- nrow(mf)
  nomterms <- attr(attr(mf, "terms"), "term.labels")
  nterms <- length(nomterms)
  p <- if ("simplified" == model) 
    nterms - 1 else nterms
  nparam <- if ("simplified" == model) 
    p + 1 else p + 2
  # paramname <- paste0('kappa', 0:p)
  paramname <- nomterms
  if ("complete" == model) 
    paramname <- c(paramname, paste0("beta", p + 1))
  # first column = response variable
  y <- as.vector(mf[, 1])
  # explanatory variables
  noms <- strsplit(nomterms, split = ":")
  noms <- do.call(rbind, noms)
  if ("simplified" == model) {
    x0 <- mf[, noms[1, 1]]
    noms <- noms[-1, , drop = FALSE]
  }
  matx <- as.matrix(mf[, noms[, 1], drop = FALSE])
  if (ncol(noms) == 1) {
    matz <- matrix(1, ncol = ncol(matx), nrow = nrow(matx))  # all z are 1
  } else {
    matz <- as.matrix(mf[, noms[, 2], drop = FALSE])
    matz[, noms[, 2] == noms[, 1]] <- 1  # unspecified z are 1
  }
  weight = rep(1, nobs) * (is.null(weights)) + (!is.null(weights)) * 
    weights
  
  ### log-likelihood function
  LL <- function(param) {
    angleref <- if ("simplified" == model) 
      x0 else rep(param[p + 2], nobs)
    # length of the vector
    sinmu <- param[1] * sin(angleref) + (matz * sin(matx)) %*% 
      param[2:(p + 1)]
    cosmu <- param[1] * cos(angleref) + (matz * cos(matx)) %*% 
      param[2:(p + 1)]
    long <- as.vector(sqrt(sinmu^2 + cosmu^2))
    # predicted value form the model
    mui <- as.vector(atan2(sinmu, cosmu))
    # log likelihood
    term1 <- param[1] * cos(y - angleref) + (matz * cos(y - 
      matx)) %*% param[2:(p + 1)]
    # LL <- sum(weight*term1) - sum(weight*log(besselI(long, 0,
    # expon.scaled = FALSE)))
    LL <- sum(term1) - sum(log(besselI(long, 0, expon.scaled = FALSE)))
    
    
    list(LL = LL, long = long, mui = mui)
  }
  # Function that update parameter of the log-likelihood
  # function
  paramUpdate <- function(paramk, long, mui) {
    angleref <- if ("simplified" == model) 
      x0 else rep(paramk[p + 2], nobs)
    matx0 <- cbind(angleref, matx)
    matz0 <- cbind(rep(1, nobs), matz)
    # score vector
    Along <- as.vector(besselI(long, 1, expon.scaled = FALSE)/besselI(long, 
      0, expon.scaled = FALSE))
    matu <- matz0 * (cos(y - matx0) - cos(matx0 - mui) * 
      Along)
    if ("complete" == model) {
      
      matu <- cbind(matu, paramk[1] * sin(y - angleref) - 
        sin(mui - angleref) * Along)
    }
    vecs <- colSums(matu)
    names(vecs) <- paramname
    # Fisher information matrix
    Xc <- matz0 * cos(matx0 - mui)
    Xs <- matz0 * sin(matx0 - mui)
    if ("complete" == model) {
      
      Xc <- cbind(Xc, paramk[1] * sin(mui - paramk[p + 
        2]))
      Xs <- cbind(Xs, paramk[1] * cos(mui - paramk[p + 
        2]))
    }
    Dc <- diag(1 - Along/long - Along^2, nrow = nobs, ncol = nobs)
    Ds <- diag(Along/long, nrow = nobs, ncol = nobs)
    matI <- t(Xc) %*% Dc %*% Xc + t(Xs) %*% Ds %*% Xs
    colnames(matI) <- rownames(matI) <- paramname
    # update of parameters
    dparam <- as.vector(solve(matI, vecs))
    paramk1 <- paramk + dparam
    list(paramk1 = paramk1, dparam = dparam, matu = matu, 
      matI = matI)
  }
  # initial values of parameters beta we try 10 000 differents
  # possible values
  if (is.null(initbeta)) {
    pginit <- if (is.null(control$pginit)) 
      1000 else control$pginit
    pg <- round(pginit^(1/nparam))
    possparam <- rep(list(seq(-1, 1, length.out = pg + 2)[-c(1, 
      pg + 2)]), p + 1)
    if ("complete" == model) 
      possparam[[nparam]] <- seq(0, 2 * pi, length.out = pg + 
        2)[-c(1, pg + 2)]
    possVal <- cbind(expand.grid(possparam), NA)
    colnames(possVal) <- c(paramname, "LL")
    maxLL <- function(param) LL(param = param)$LL
    possVal[, nparam + 1] <- apply(possVal[, 1:nparam], 1, 
      maxLL)
    paramk <- unlist(possVal[which.max(possVal[, nparam + 
      1]), 1:nparam])
  } else {
    if (length(initbeta) != nparam) 
      stop("for the requested model, 'initparam' must be of length ", 
        nparam)
    paramk <- initbeta
  }
  # log likelihood value with respect to initial parameter
  calcul <- LL(param = paramk)
  maxLLk <- calcul$LL
  long <- calcul$long
  mui <- calcul$mui
  # Initialization
  iter <- iter.sh <- 0
  maxiter <- if (is.null(control$maxiter)) 
    1000 else control$maxiter
  mindiff <- if (is.null(control$mindiff)) 
    1e-06 else control$mindiff
  conv <- FALSE
  # Initialisation of the matrix of informations during the
  # iterations
  iter.detail <- matrix(NA, nrow = maxiter + 1, ncol = nparam + 
    3)
  colnames(iter.detail) <- c(paramname, "maxLL", "iter", "nitersh")
  iter.detail[1, ] <- c(paramk, maxLLk, iter, iter.sh)
  # start of the fit
  while (!conv && iter <= maxiter) {
    # update of the parameters
    maj <- paramUpdate(paramk = paramk, long = long, mui = mui)
    paramk1 <- maj$paramk1
    dparam <- maj$dparam
    # computation of the log-likelihood
    calcul <- LL(param = paramk1)
    maxLLk1 <- calcul$LL
    long <- calcul$long
    mui <- calcul$mui
    # if the criteria as decrease then do step halving
    iter.sh <- 0
    while (maxLLk1 < maxLLk) {
      iter.sh <- iter.sh + 1
      paramk1 <- paramk + dparam/(2^iter.sh)
      calcul <- LL(param = paramk1)
      maxLLk1 <- calcul$LL
      long <- calcul$long
      mui <- calcul$mui
      if (iter.sh >= maxiter) 
        break
    }
    # Does the criteria increase more than mindiff?
    if (maxLLk1 < maxLLk) {
      conv <- FALSE
      warning("the algorithm did no converge, it failed to maximize the log likelihood")
      break
    } else {
      conv <- if (maxLLk1 - maxLLk > mindiff) 
        FALSE else TRUE
      paramk <- paramk1
      maxLLk <- maxLLk1
      iter <- iter + 1
      iter.detail[iter + 1, ] <- c(paramk, maxLLk, iter, 
        iter.sh)
    }
  }
  if (iter > maxiter + 1) {
    warning("the algorithm did not converge, the maximum number of iterations was reached")
  } else {
    iter.detail <- iter.detail[1:(iter + 1), , drop = FALSE]
  }
  
  ### Computation of standard errors and Fisher information
  ### matrix for the final parameters.
  if (maxLLk == maxLLk1) {
    maj <- paramUpdate(paramk = paramk, long = long, mui = mui)
  }
  matu <- maj$matu
  matI <- maj$matI
  # parametric estimation of the covariance matrix
  v1 <- solve(matI)
  # non parametric estimation
  mid <- matrix(0, ncol = nparam, nrow = nparam)
  for (i in 1:nobs) {
    mid <- mid + t(matu[i, , drop = FALSE]) %*% matu[i, , 
      drop = FALSE]
  }
  v2 <- v1 %*% mid %*% v1
  
  ### Results for the betas
  paramb <- paramk[2:(p + 1)]/paramk[1]
  matDeriv <- rbind(-paramk[2:(p + 1)]/paramk[1]^2, diag(1/paramk[1], 
    nrow = p, ncol = p))
  vb <- t(matDeriv) %*% v1[1:(p + 1), 1:(p + 1)] %*% matDeriv
  vb2 <- t(matDeriv) %*% v2[1:(p + 1), 1:(p + 1)] %*% matDeriv
  # names(paramb) <- colnames(vb) <- rownames(vb)
  # <-colnames(vb2) <- rownames(vb2)<- paste0('beta', 1:p)
  names(paramb) <- colnames(vb) <- rownames(vb) <- colnames(vb2) <- rownames(vb2) <- paramname[-1]
  
  ### Output
  zvalue <- abs(paramk)/sqrt(diag(v1))
  p <- round(2 * pnorm(abs(paramk)/sqrt(diag(v1)), lower.tail = FALSE), 
    5)
  parameters <- cbind(paramk, sqrt(diag(v1)), zvalue, p)
  # colnames(parameters) <- c('estimate', paste('stderr', 1:2,
  # sep=''))
  colnames(parameters) <- c("estimate", "stderr", "z value", 
    "P(|z|>.)")
  rownames(parameters) <- paramname
  parambeta <- cbind(paramb, sqrt(diag(vb)), sqrt(diag(vb2)))
  colnames(parambeta) <- c("estimate", paste("stderr", 1:2, 
    sep = ""))
  res <- sin(y - mui)
  autocorr <- acf(res, plot = FALSE)
  out <- list(MaxLL = maxLLk, parameters = parameters, varcov1 = v1, 
    varcov2 = v2, parambeta = parambeta, varcovbeta1 = vb, 
    varcovbeta2 = vb2, autocorr = autocorr, matx = matx, 
    matz = matz, y = y, long = long, mui = mui, iter.detail = iter.detail, 
    call = call)
  class(out) <- "consensus"
  out
}


#' @param x An object, produced by the \code{\link{consensus}} function, to print.
#' @param \dots Further arguments to be passed to \code{print.default}. 
print.consensus <- function(x, ...) {
  cat("\nMaximum log-likelihood :", x$MaxLL, "\n")
  cat("\nParameters:\n")
  print.default(x$parameters, print.gap = 2, quote = FALSE, 
    right = TRUE, ...)
  cat("\n")
  # cat('\nBeta Parameters:\n') print.default(x$parambeta,
  # print.gap = 2, quote = FALSE, right=TRUE, ...) cat('\n')
  invisible(x)
}





##### plot smooth proba on the trajectory (K=2 states)
##### #######################



#' Plot smooth probabilities on trajectory
#' 
#' Function that smooth eacgh step of the trajectory by the smooth probabilites
#' This plot can be used to caracterize spatially the behavior of the animal.
#' 
#'  
#'    
#' 
#' @param param  list of all the parameters of the model in this order: P, kappa and thetalin.
#' @param data data frame of the directions (y), distances (d), explanatory angles variables (x)
#'            and explanatory real variables (z)
#' @param long vectors of longtitude coordinates of the trajectory (x axis)
#' 
#' @param lat vectors of latitude coordinates of the trajectory (y axis)
#' 
#' @param type of the model: 'angle-dist' for bivariate model on direction-step length
#'  and 'angle' for univariate model on direction (Default).  
#'  
#' @param dist type of distribution on the step length. The choices are 'gamma' (Default) or 'weibull'
#' 
#' @references 
#' 
#' 
#' @author Aurélien Nicosia 

PlotEstim <- function(data, param, long, lat, type = "angle-dist", 
  dist = "gamma") {
  
  
  K <- length(param$pi0)
  p <- length(param$kappa[1, ]) - 1
  steps <- length(y)
  FS <- FilterSmooth(param, data, type, dist)
  E.smooth <- FS$E.smooth
  explo <- 1
  position <- cbind(long, lat)
  
  Sx1 <- max(position[, 1]) - min(position[, 1])
  Sy1 <- max(position[, 2]) - min(position[, 2])
  x1 <- seq(min(position[, 1]) - 0.1 * Sx1, max(position[, 
    1]) + 0.1 * Sx1, length.out = 1024)
  y1 <- seq(min(position[, 2]) - 0.1 * Sy1, max(position[, 
    2]) + 0.1 * Sy1, length.out = 1024)
  Map <- matrix(1, length(x1), length(y1))
  image(x1, y1, (Map), col = "white")
  
  n <- 5
  colorter <- rgb(c(seq(0, 1, 1/n), rep(1, n)), c(abs(seq(0, 
    1 - 1/n, 1/n)), 1, abs(seq(1 - 1/n, 0, -1/n))), c(rep(1, 
    n), abs(seq(1, 0, -1/n))))
  
  for (i in (2:(length(position[, 1]) - 1))) {
    classe <- floor(E.smooth[i, explo] * 10) + 1
    segments(position[i - 1, 1], position[i - 1, 2], position[i, 
      1], position[i, 2], lwd = 1, col = colorter[classe])
  }
}



