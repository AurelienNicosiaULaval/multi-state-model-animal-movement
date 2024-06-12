##############################
#
#  Run a two state general random walk model
#
# Author: Aurélien nicosia
#
###############################

# needed packages

# packages needed
install.packages(c("boot","circular","Renext","pracma","fitdistrplus"))
library(boot)
library(circular)
library(Renext)
library(pracma)
library(fitdistrplus)

# Please load the appropriate workspace 

source("UtilityFunction.R")
jeu <- data.frame(read.table(file = "caribou.txt", header = TRUE));

# initialization parameter of the EM algorithm with a one state model
n <- length(jeu$y);
jeu$yprec <- c(0, jeu$y[ - n]);
mc <- consensus(formula = y ~ yprec + xcut + xcenter, data = jeu);


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

  
steps <- length(jeu$y);
position <- matrix(0, steps, 2);
position[1, ] <- c(0, 0);
for (i in (2:steps))
{
  position[i, ] = position[i - 1, ] + d[i] * 1000 * cbind(cos(y[i]), sin(y[i]));
}

# Plot of the behavior
PlotEstim(jeu, param = GM$param, long = position[, 1], lat = position[, 2])
    


