# generate a dataset, using code from Giganti & Shepherd (2020)
set.seed(455)
library(mvtnorm)
# Fixed values
obs <- 4000
ID <- rep(1:obs)
X1X2cov <- (-0.25)
X1var <- 1
X2var <- 1
Astar.var <- 1
A.var <- 1
beta.star <- c(1, 0, 0)
gamma.star <- c(-3, 0, 0, 0.5)
beta <- c(0, -1, 0.5, 0.9, 0.5)
gamma <- c(-5.5, -2, 1, 0.0, 5.0, 0.5)

# Generate values
# Covariates X and X2
X1X2 <- rmvnorm(obs, mean = c(0, 0),
                sigma = matrix(c(X1var, X1X2cov, X1X2cov, X2var), 2, 2)
                )
X1 <- X1X2[,1]
X2 <- X1X2[,2]
Xmat <- cbind(1,X1,X2)
rownames(Xmat) <- ID
colnames(Xmat) <- c("Intercept","X1","X2")

####################################
#Unvalidated data
####################################
#A*
A.star <- rnorm(obs, Xmat %*% beta.star, Astar.var)

#D*
LP.D.star <- Xmat %*% gamma.star[1:ncol(Xmat)] +
  A.star * gamma.star[length(gamma.star)]
p.D.star <- exp(LP.D.star) / (1 + exp(LP.D.star))

#Indicator of event (D*) at time t
D.star <- rbinom(obs, 1, p.D.star)

####################################
#Validated data
####################################
#A
LP.A <- Xmat %*% beta[1:ncol(Xmat)] + A.star * beta[ncol(Xmat) + 1] +
  D.star * beta[ncol(Xmat) + 2]
A <- rnorm(obs, LP.A, A.var)

#D
LP.D <- Xmat %*% gamma[1:ncol(Xmat)] + A.star * gamma[ncol(Xmat) + 1] +
  D.star * gamma[ncol(Xmat) + 2] + A * gamma[ncol(Xmat) + 3]
p.D <- exp(LP.D) / (1 + exp(LP.D))

#Indicator of event (D) at time t
D <- rbinom(obs, 1, p.D)

subsampleN <- 1000
inthresh <- 2

#These values are fixed for a given simulation, but may vary across simulations
#Specify the audit/subsample size
n.sample <- subsampleN

#Specify the inclusion threshold
inclusionthreshold <- inthresh

#Randomly sample records to represent audited cohort
ChooseSubset <- sample(1:obs, n.sample, replace = F)
sampled <- as.numeric((1:obs) %in% ChooseSubset)

# Original dataset
dat <- dat.original <- data.frame(ID,X1,X2,A.star,D.star,A,D)
dat$Intercept <- 1
dat$A <- ifelse(sampled == 0, NA, dat$A)
dat$D <- ifelse(sampled == 0, NA, dat$D)
dat$D <- factor(dat$D, levels = c(0, 1))

giganti_data <- dat
usethis::use_data(giganti_data, overwrite = TRUE)
