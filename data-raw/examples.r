library(R2jags)
library(lattice)
library(MASS)

myburnin <- 100
mysample <- 100

expit <- function(arg) {
  exp(arg)/(1+exp(arg))
}

#generate data
n <- 1000
xzcorr <- 0.25

x1 <- rnorm(n)
x2xb <- x1
x2pr <- expit(x2xb)
x2 <- 1*(runif(n)<x2pr)
x3xb <- 0.5*x1-0.5*x2
x3lambda <- exp(x3xb)
x3 <- rpois(n, x3lambda)

resvar <- 10

y <- x1+x2+x3+rnorm(n,sd=resvar^0.5)

x1misspr <- expit(y-5)
x1[(runif(n)<x1misspr)] <- NA
x2[(runif(n)<0.25)] <- NA
x3[(runif(n)<0.25)] <- NA
mydata <- data.frame(y,x1,x2,x3)

cca <- lm(y~x1+x2+x3, mydata)


bayesmiss(mydata, smformula="y~x1+x2+x3",method=c("","norm","logit","pois"),order=c(0,1,2,3),nChains=5)



#generate data
n <- 1000
xzcorr <- 0.25

x1 <- rnorm(n)
x2xb <- x1
x2pr <- expit(x2xb)
x2 <- 1*(runif(n)<x2pr)
x3 <- ceiling(3*runif(n))

resvar <- 10

y <- x1+x2+rnorm(n,sd=resvar^0.5)

x1misspr <- expit(y-5)
x1[(runif(n)<x1misspr)] <- NA
x2[(runif(n)<0.25)] <- NA
x3[(runif(n)<x1misspr)] <- NA
mydata <- data.frame(y,x1,x2,x3)

cca <- lm(y~x1+x2+x3, mydata)


bayesmiss(mydata, smoutcome="y",method=c("norm","norm","logit","ologit"),order=c(4,1,2,3))

library(rjags,coda)
jags.data <- as.list(mydata)
jags.data$n <- n
jags.data$beta_x1_mean <- rep(0, 1)
jags.data$beta_x1_prec <- 0.0001*diag(1)
jags.data$tau_x1_alpha <- 0.5
jags.data$tau_x1_beta <- 0.5
jags.data$beta_x2_mean <- rep(0, 2)
jags.data$beta_x2_prec <- 0.0001*diag(2)
jags.data$beta_y_mean <- rep(0, 3)
jags.data$beta_y_prec <- 0.0001*diag(3)
jags.data$tau_y_alpha <- 0.5
jags.data$tau_y_beta <- 0.5
jags.params <- c("beta_y","tau_y")
jagsmodel <- jags.model(data=jags.data, file="bayesmissmod.bug",n.chains=5)
burnin <- coda.samples(jagsmodel, variable.names=jags.params, n.iter=200)
mainsample <- coda.samples(jagsmodel, variable.names=jags.params, n.iter=200)
summary(mainsample)
gelman.diag(mainsample)


#generate data
n <- 10000

x1 <- rnorm(n)
x2xb <- x1
x2pr <- expit(x2xb)
x2 <- 1*(runif(n)<x2pr)
x3tilde <- x1+x2+rlogis(n)
x3 <- rep(1,n)
x3[x3tilde>0] <- 2
x3[x3tilde>1] <- 3

resvar <- 10

y <- x1+x2+x3+rnorm(n,sd=resvar^0.5)

x1misspr <- expit(y-5)
#x1[(runif(n)<x1misspr)] <- NA
#x2[(runif(n)<0.25)] <- NA
x3[(runif(n)<x1misspr)] <- NA
mydata <- data.frame(y,x1,x2,x3)

cca <- lm(y~x1+x2+x3, mydata)


bayesmiss(mydata, smformula="y~x1+x2+x3",method=c("","","","ologit"),order=c(0,0,0,1),nChains=5)


