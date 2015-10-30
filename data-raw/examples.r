
myburnin <- 100
mysample <- 100

expit <- function(arg) {
  exp(arg)/(1+exp(arg))
}

#generate data
n <- 5000
xzcorr <- 0.25
numCovs <- 4
SigmaMat <- array(xzcorr, dim=c(numCovs,numCovs))
diag(SigmaMat) <- 1
cov <- mvrnorm(n, mu=rep(0,numCovs), Sigma=SigmaMat)
x1 <- cov[,1]
x2 <- cov[,2]
x3xb <- x1+x2
x3pr <- expit(x3xb)
x3 <- 1*(runif(n)<x3pr)
z <- cov[,4]
resvar <- 1

y <- x1+x2+x3+z+rnorm(n,sd=resvar^0.5)

x1misspr <- expit(y-5)
x1[(runif(n)<x1misspr)] <- NA
x2[(runif(n)<0.1)] <- NA
x3[(runif(n)<0.1)] <- NA
mydata <- data.frame(y,z,x1,x2,x3)

cca <- lm(y~x1+x2+x3+z, mydata)


bayesmiss(mydata, omformula="y~x1+x2+x3+z",method=c("","","norm","norm","logit"),order=c(0,0,1,2,3),nChains=5)

