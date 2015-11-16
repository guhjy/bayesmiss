#' Bayesian regression with missing covariates
#'
#' \code{bayesmiss} generates JAGS model code and an R script to perform
#' Bayesian regression, allowing for missingness in covariates.
#'
#' \code{bayesmiss} faciliates running Bayesian regression models in which there
#' are missing values in some of the covariates. The function generates two files,
#' one a JAGS model file, and one a R script file.
#'
#' The JAGS model defines the
#' user's substantive model, and models as required to handle the missing covariates.
#'
#' The R script file generated contains commands to generate the required JAGS data
#' and parameter objects, which are passed when calling \code{jags}.
#'
#' The \code{order} argument is used to specify the order in which the joint
#' distribution of the partially observed covariates is to be modelled. Models
#' for the partially observed covariates automatically condition on substantive model
#' covariates which are fully observed.
#'
#' The R script file generated contains commands to generate the required JAGS data
#' and parameter objects, which are passed when calling \code{jags}. Note that
#' \code{bayesmiss} doesn't actually run this code - the user, having ensured that
#' the \code{R2jags} package is installed, must run this code to fit the model.
#' Note that the MCMC options in the call to \code{jags} are just suggested default.
#' It is up to the user to ensure, via the usual diagnostics for MCMC, that a
#' sufficient number of iterations have been run to ensure convergence of the chains.
#'
#' If it is desired to add interactions or non-linear covariate effects, the easiest
#' approach is to first run \code{bayesmiss} omitting these terms, and then modify
#' the JAGS code file and R code specifying the priors as needed.
#'
#' @param originaldata the data frame upon which the analysis is to be performed.
#'
#' @param smtype a string specifying the substantive model type
#'
#' @param smformula an expression for the substantive model type
#'
#' @param method a vector of strings specifying the imputation method for each of
#' the columns in \code{originaldata}. Currently the possible values are "norm"
#' (normal linear model), "logit" (logistic regression), "mlogit" (multinomial
#' logistic regression, and "ologit" (ordinal logistic regression". The elements
#' corresponding to fully observed variables should be the empty string "".
#' Variables imputed using mlogit or ologit should be stored as factors, coded
#' 1:K, where K is the number of categories.
#'
#' @param order a vector specifying the order in which the joint model for missing
#' covariates should be constructed. e.g. c(0,0,1,2) specifies that the first two
#' variables are fully observed and that the third and fourth are partially observed,
#' with the covariate model constructed by modelling the marginal distribution of
#' the third variable and the conditional distribution of the fourth given the third.
#'
#' @return \code{bayesmiss} generates two files in the current working directory:
#' \code{bayesmissmod.bug} is a JAGS model file for the constructed model, and \code{bayesmissRscript.r}
#' is an R file containing R code for generating the required JAGS parameters and data
#' objects, and a call to the \code{jags} function of the \code{R2jags} package for
#' fitting the model.
#'
bayesmiss <- function(originaldata,smtype,smformula,method,order) {
  n <- dim(originaldata)[1]
  #create matrix of response indicators
  r <- 1*(is.na(originaldata)==0)

  #outcome model
  outcomeCol <- which(colnames(originaldata)==as.formula(smformula)[[2]])
  outcomename <- as.character(as.formula(smformula)[[2]])

  omcovnames <- attr(terms(as.formula(smformula)), "term.labels")
  if (smtype=="lm") {
    omDist <- paste("      ",outcomename,"[i] ~ dnorm(mu[i], outcome_tau)",sep="")
    omLinPred <- "      mu[i] <- beta[1]"
  }
  else if (smtype=="poisson") {
    omDist <- paste("      ",outcomename,"[i] ~ dpois(mu[i])",sep="")
    omLinPred <- "      log(mu[i]) <- beta[1]"
  }
  else stop("You must supply a valid value as the smtype argument.")
  for (i in 1:length(omcovnames)) {
    omLinPred <- paste(omLinPred, " + beta[",i+1,"]*", omcovnames[i], "[i]", sep="")
  }

  #start constructing model code string
  modelCode <- c(
    "model {",
    "   for (i in 1:n) {",
    omDist,
    omLinPred)
  if (smtype=="lm") {
    priorCode <- c("   beta ~ dmnorm(betamean,betaprec)","   outcome_tau ~ dgamma(tau_alpha, tau_beta)")
  }
  else if (smtype=="poisson") {
    priorCode <- "   beta ~ dmnorm(betamean,betaprec)"
  }

  fullObsVars <- which((colSums(r)==n) & (colnames(originaldata) %in% omcovnames))

  #start adding to R script file
  rScriptText <- "library(rjags,coda)"
  rScriptText <- c(rScriptText,paste("jags.data <- as.list(",deparse(substitute(originaldata)),")", sep=""),
  "jags.data$n <- n",
  paste("jags.data$betamean <- rep(0, ",length(omcovnames)+1,")",sep=""),
  paste("jags.data$betaprec <- 0.0001*diag(",length(omcovnames)+1,")",sep=""))
  if (smtype=="lm") {
    rScriptText <- c(rScriptText,"jags.data$tau_alpha <- 0.5","jags.data$tau_beta <- 0.5")
  }

  #models for partially observed covariates
  numPartialVars <- max(order)
  prevMissVars <- vector(mode="integer",0)

  if (numPartialVars>0) {
    for (var in 1:numPartialVars) {
      targetCol <- which(order==var)
      missName <- colnames(originaldata)[targetCol]
      missCovNames <- colnames(originaldata)[c(fullObsVars,prevMissVars)]

      if (method[targetCol]=="norm") {
        missDist <- paste("      ",missName,"[i] ~ dnorm(mu_",missName,"[i], tau_",missName,")", sep="")
        missLinPred <- paste("      ","mu_",missName,"[i] <- gamma_",missName,"[1]", sep="")
      }
      else if (method[targetCol]=="logit") {
        missDist <- paste("      ",missName,"[i] ~ dbern(mu_",missName,"[i])", sep="")
        missLinPred <- paste("      ","logit(mu_",missName,"[i]) <- gamma_",missName,"[1]", sep="")
      }
      else if (method[targetCol]=="pois") {
        missDist <- paste("      ",missName,"[i] ~ dpois(mu_",missName,"[i])", sep="")
        missLinPred <- paste("      ","log(mu_",missName,"[i]) <- gamma_",missName,"[1]", sep="")
      }
      else if (method[targetCol]=="mlogit") {
        missDist <- paste("      ",missName,"[i] ~ dcat(pi_",missName,"[i,])", sep="")
        numCats <- nlevels(originaldata[,targetCol])

        missLinPred <- paste("      ","pi_",missName,"[i,1] <- 1", sep="")
        for (catNum in 2:numCats) {
          missLinPred <- paste(missLinPred, " - pi_",missName,"[i,",catNum,"]", sep="")
        }
        for (catNum in 2:numCats) {
          missLinPred <- c(missLinPred, paste("      ","pi_",missName,"[i,",catNum,"] <- exp(xb_",missName,"[i,",catNum-1,"])/denom_",missName,"[i]", sep=""))
        }
        denomExpr <- paste("      ","denom_",missName,"[i] <- 1", sep="")
        for (catNum in 2:numCats) {
          denomExpr <- paste(denomExpr,"+exp(xb_",missName,"[i,",catNum-1,"])", sep="")
        }
        missLinPred <- c(missLinPred,denomExpr)
        for (catNum in 2:numCats) {
          xbExpr <- paste("      ","xb_",missName,"[i,",catNum-1,"] <- gamma_",missName,"_",catNum-1,"[1]", sep="")
          if (length(missCovNames)>0) {
            for (i in 1:length(missCovNames)) {
              xbExpr <- paste(xbExpr, " + gamma_",missName,"_",catNum-1,"[",i+1,"]*", missCovNames[i], "[i]", sep="")
            }
          }
          missLinPred <- c(missLinPred, xbExpr)
        }
      }
      else if (method[targetCol]=="ologit") {
        missDist <- paste("      ",missName,"[i] ~ dcat(pi_",missName,"[i,])", sep="")
        numCats <- nlevels(originaldata[,targetCol])

        missLinPred <- paste("      ","xb_",missName,"[i] <- 0", sep="")
        if (length(missCovNames)>0) {
          for (i in 1:length(missCovNames)) {
            missLinPred <- paste(missLinPred, " + gamma_",missName,"[",i,"]*", missCovNames[i], "[i]", sep="")
          }
        }
        piExpr <- paste("      ","pi_",missName,"[i,1] <- 1", sep="")
        for (catNum in 2:numCats) {
          piExpr <- paste(piExpr, " - pi_",missName,"[i,",catNum,"]", sep="")
        }
        for (catNum in 2:(numCats-1)) {
          piExpr <- c(piExpr, paste("      ","pi_",missName,"[i,",catNum,"] <- 1/(1+exp(-k_",missName,"_",catNum," + xb_",
                                    missName,"[i])) - 1/(1+exp(-k_",missName,"_",catNum-1," + xb_",missName,"[i]))", sep=""))
        }
        piExpr <- c(piExpr, paste("      ","pi_",missName,"[i,",numCats,"] <- 1 - 1/(1+exp(-k_",missName,"_",numCats-1," + xb_",missName,"[i]))", sep=""))
        missLinPred <- c(missLinPred, piExpr)
      }
      else stop(paste("Method ",method[targetCol]," not recognised.",sep=""))


      if ((method[targetCol]=="norm") | (method[targetCol]=="logit") | (method[targetCol]=="pois")) {
        if (length(missCovNames)>0) {
          for (i in 1:length(missCovNames)) {
            missLinPred <- paste(missLinPred, " + gamma_",missName,"[",i+1,"]*", missCovNames[i], "[i]", sep="")
          }
        }
      }

      #append to modelCode
      modelCode <- c(modelCode,"",missDist,missLinPred)

      #append to priorCode
      if (method[targetCol]=="mlogit") {
        priorCode <- c(priorCode,"")
        for (catNum in 2:numCats) {
          priorCode <- c(priorCode,paste("   gamma_",missName,"_",catNum-1," ~ dmnorm(gamma_",missName,"_",catNum-1,"_mean,gamma_",missName,"_",catNum-1,"_prec)", sep=""))
        }
      }
      else {
        priorCode <- c(priorCode,"",paste("   gamma_",missName," ~ dmnorm(gamma_",missName,"_mean,gamma_",missName,"_prec)", sep=""))
        if (method[targetCol]=="norm") {
          priorCode <- c(priorCode, paste("   tau_",missName," ~ dgamma(tau_alpha, tau_beta)", sep=""))
        }
        else if (method[targetCol]=="ologit") {
          #priors for cutpoints in ordinal logistic regression
          priorCode <- c(priorCode, paste("   k_",missName,"_1 ~ dnorm(k_",missName,"_1_mean, k_",missName,"_1_prec)", sep=""))
          for (catNum in 2:(numCats-1)) {
            priorCode <- c(priorCode, paste("   k_",missName,"_",catNum," <- k_",missName,"_1 + kinc_",missName,"_",catNum-1,  sep=""))
            priorCode <- c(priorCode, paste("   kinc_",missName,"_",catNum-1," ~ dlnorm(kinc_",missName,"_",catNum-1,"_mean, kinc_",missName,"_",catNum-1,"_prec)",sep=""))
          }
        }
      }

      #append priors to JAGS data list
      if (method[targetCol]=="mlogit") {
        for (catNum in 2:numCats) {
          rScriptText <- c(rScriptText, paste("jags.data$gamma_",missName,"_",catNum-1,"_mean <- rep(0, ",length(missCovNames)+1,")", sep=""),
                           paste("jags.data$gamma_",missName,"_",catNum-1,"_prec <- 0.0001*diag(",length(missCovNames)+1,")", sep=""))
        }
      }
      else if (method[targetCol]=="ologit") {
        rScriptText <- c(rScriptText, paste("jags.data$gamma_",missName,"_mean <- rep(0, ",length(missCovNames),")", sep=""),
                         paste("jags.data$gamma_",missName,"_prec <- 0.0001*diag(",length(missCovNames),")", sep=""))
        rScriptText <- c(rScriptText, paste("jags.data$k_",missName,"_1_mean <- 0", sep=""))
        rScriptText <- c(rScriptText, paste("jags.data$k_",missName,"_1_prec <- 0.0001", sep=""))
        for (catNum in 2:(numCats-1)) {
          rScriptText <- c(rScriptText, paste("jags.data$kinc_",missName,"_",catNum-1,"_mean <- 0", sep=""))
          rScriptText <- c(rScriptText, paste("jags.data$kinc_",missName,"_",catNum-1,"_prec <- 0.0001", sep=""))
        }
      }
      else {
        rScriptText <- c(rScriptText, paste("jags.data$gamma_",missName,"_mean <- rep(0, ",length(missCovNames)+1,")", sep=""),
          paste("jags.data$gamma_",missName,"_prec <- 0.0001*diag(",length(missCovNames)+1,")", sep=""))
      }

      prevMissVars <- c(prevMissVars,targetCol)
    }
  }

  #write model to file
  modFileName <- "bayesmissmod.bug"
  fileConn <- file(modFileName )
  modelCode <- c(modelCode, "   }","")
  priorCode <- c(priorCode, "}")
  writeLines(c(modelCode, priorCode), fileConn)
  close(fileConn)

  print(paste("Your JAGS model has been written to: ", modFileName, sep=""))

  #run bayes analysis
  rScriptFileName <- "bayesmissRScript.r"
  fileConn <- file(rScriptFileName)
  rScriptText <- c(rScriptText, "jags.params <- c(\"beta\", \"outcome_tau\")")

  rScriptText <- c(rScriptText, paste("jagsmodel <- jags.model(data=jags.data, file=\"",modFileName,"\",n.chains=5)", sep=""))
  rScriptText <- c(rScriptText, "burnin <- coda.samples(jagsmodel, variable.names=c(\"beta\"), n.iter=200)")
  rScriptText <- c(rScriptText, "mainsample <- coda.samples(jagsmodel, variable.names=c(\"beta\"), n.iter=200)")
  rScriptText <- c(rScriptText, "summary(mainsample)")
  rScriptText <- c(rScriptText, "gelman.diag(mainsample)")
  writeLines(rScriptText, fileConn)
  close(fileConn)

  print(paste("Your R script has been written to: ", rScriptFileName, sep=""))

}
