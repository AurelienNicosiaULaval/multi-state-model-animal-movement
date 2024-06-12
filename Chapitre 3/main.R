############################## Run a two state HMM-SSF and general random walk model
############################## Author: Aurélien nicosia



# Please load the appropriate workspace

source("UtilityFunction_HMMSSF.R")


# HMM-SSF ------------


dat <- read_csv("1051complet.csv")

# initialization parameter of the EM algorithm with a one
# state model

fit <- coxph(Surv(Times, Cote) ~ cos.persis + dist.neg + dist.log + 
  eau + Dec + meadow + road + strata(idstrate), data = dat, 
  x = TRUE)
fit

# data set for the fit

df.marix <- list(y = dat$Cote, strata = dat$idstrate, x = as.matrix(fit$x), 
  Times = dat$Times)
K <- 2
Initparam <- list(betainit = rbind(fit$coefficients, fit$coefficients))

# Fit
SSF <- GlobalMaximaSSF(Surv(Times, Cote) ~ cos.persis + dist.neg + 
  dist.log + eau + Dec + meadow + road + strata(idstrate), 
  df.marix, dat, K, Initparam, nb_init = 4, precision = c(2, 
    4), EMMin = c(3, 50))

SSF



# General Random Walk Model ---------

source("UtilityFunction_RandomWalkModel.R")

jeu <- data.frame(read.table(file = "caribou.txt", header = TRUE));

# initialization parameter of the EM algorithm with a one state model
n <- length(jeu$y);
jeu$yprec <- c(0, jeu$y[ - n]);
mc <- consensus(formula = y ~ yprec + xcut + xcenter, data = jeu);
mc

kappac <- cbind(mc$parameters[, 1], mc$parameters[, 3]);   
dmom <- ini.mixexp2(jeu$d);
mel.d <- mledist(jeu$d, distr = "gamma");
pi0 <- c(0.7, 0.3);
kappainit <- rbind(c(kappac[1, 1], kappac[2:3, 1]),
                   c(kappac[1, 1], kappac[2:3, 1]));
thetalininit <- rbind(c(mel.d$estimate[1], 1/mel.d$estimate[2]),
                      c(mel.d$estimate[1], 1/mel.d$estimate[2]));
colnames(thetalininit) <- c("shape", "scale");
paraminitiaux <- list(pi0 = pi0, kappainit = kappainit,
                      thetalininit = thetalininit, Pinit = NULL);

# data set for the fit

target <- cbind(jeu$xcut, jeu$xcenter);
p <- length(target[1, ]);
z <- matrix(1,n,p);
jeu <- list(x = as.matrix(target), y=as.vector(jeu$y),
            z = z, d = as.vector(jeu$d));


# FIT
K <- 2
GM <- GHRandomWalk(data = jeu, K, paraminitiaux, nb_init = 5, dist = "gamma");
GM


# Plot of the behavior
PlotEstim(jeu, param = GM$param)