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


bayesmiss(mydata, omformula="y~x1+x2+x3",method=c("","norm","logit","pois"),order=c(0,1,2,3),nChains=5)

#generate data
n <- 1000
xzcorr <- 0.25

x1 <- rnorm(n)
x2xb <- x1
x2pr <- expit(x2xb)
x2 <- 1*(runif(n)<x2pr)
#x3 <- ceiling(3*runif(n))
x3 <- 1*(runif(n)<0.5)

resvar <- 10

y <- x1+x2+x3+rnorm(n,sd=resvar^0.5)

x1misspr <- expit(y-5)
x1[(runif(n)<x1misspr)] <- NA
x2[(runif(n)<0.25)] <- NA
x3[(runif(n)<0.25)] <- NA
mydata <- data.frame(y,x1,x2,x3)

cca <- lm(y~x1+x2+x3, mydata)


bayesmiss(mydata, omformula="y~x1+x2+x3",method=c("","norm","logit","logit"),order=c(0,1,2,3),nChains=5)



#generate data
n <- 1000
x1 <- 1*(runif(n)<0.5)
x2 <- 1*(runif(n)<0.5)
x3 <- 1*(runif(n)<0.5)
x4 <- 1*(runif(n)<0.5)
resvar <- 10
y <- x1+x2+x3+x4+rnorm(n,sd=resvar^0.5)
#make covariates missing
x1[(runif(n)<0.25)] <- NA
x2[(runif(n)<0.25)] <- NA
x3[(runif(n)<0.25)] <- NA
mydata <- data.frame(y,x1,x2,x3)

#cca <- lm(y~x1+x2+x3, mydata)


bayesmiss(mydata, omformula="y~x1+x2+x3",method=c("","logit","logit","logit"),order=c(0,1,2,3),nChains=5)

