bayesmiss <- function(originaldata,omformula,method,order,nIter=200,nChains=5) {
  n <- dim(originaldata)[1]
  #create matrix of response indicators
  r <- 1*(is.na(originaldata)==0)

  #outcome model
  outcomeCol <- which(colnames(originaldata)==as.formula(omformula)[[2]])
  outcomename <- as.character(as.formula(omformula)[[2]])
  outcomeModDist <- paste("      ",outcomename,"[i] ~ dnorm(mu[i], outcome_tau)",sep="")
  omcovnames <- attr(terms(as.formula(omformula)), "term.labels")
  omDist <- paste("      ",outcomename,"[i] ~ dnorm(mu[i], outcome_tau)",sep="")
  omLinPred <- "      mu[i] <- beta[1]"
  for (i in 1:length(omcovnames)) {
    omLinPred <- paste(omLinPred, " + beta[",i+1,"]*", omcovnames[i], "[i]", sep="")
  }

  #start constructing model code string
  modelCode <- c(
    "model {",
    "   for (i in 1:n) {",
    omDist,
    omLinPred)
  priorCode <- c("beta ~ dmnorm(betamean,betaprec)","outcome_tau ~ dgamma(tau_alpha, tau_beta)")

  fullObsVars <- which((colSums(r)==n) & (colnames(originaldata) %in% omcovnames))
  partialVars <- which(method=="norm")

  #add lines, if needed, for passively imputed covariates
  #passiveVars <- which((method!="") & (method!="norm"))
  #if (length(passiveVars)>0) {
  #  modelCode <- c(modelCode, paste("      ",colnames(originaldata)[passiveVars],"[i] <- ",method[passiveVars],sep=""))
  #}

  #start adding to R script file
  rScriptText <- c(paste("jags.data <- as.list(",deparse(substitute(originaldata)),")", sep=""),
  "jags.data$n <- n",
  paste("jags.data$betamean <- rep(0, ",length(omcovnames)+1,")",sep=""),
  paste("jags.data$betaprec <- 0.0001*diag(",length(omcovnames)+1,")",sep=""),
  "jags.data$tau_alpha <- 0.5",
  "jags.data$tau_beta <- 0.5")

  #models for partially observed covariates
  numPartialVars <- max(order)
  prevMissVars <- vector(mode="integer",0)

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
    for (i in 1:length(missCovNames)) {
      missLinPred <- paste(missLinPred, " + gamma_",missName,"[",i+1,"]*", missCovNames[i], "[i]", sep="")
    }

    #append to modelCode
    modelCode <- c(modelCode,missDist,missLinPred)

    #append to priorCode
    priorCode <- c(priorCode,paste("gamma_",missName," ~ dmnorm(gamma_",missName,"_mean,gamma_",missName,"_prec)", sep=""))
    if (method[targetCol]=="norm") {
      priorCode <- c(priorCode, paste("tau_",missName," ~ dgamma(tau_alpha, tau_beta)", sep=""))
    }

    #append priors to JAGS data list
    rScriptText <- c(rScriptText, paste("jags.data$gamma_",missName,"_mean <- rep(0, ",length(missCovNames)+1,")", sep=""),
      paste("jags.data$gamma_",missName,"_prec <- 0.0001*diag(",length(missCovNames)+1,")", sep=""))

    prevMissVars <- c(prevMissVars,targetCol)
  }


  #write model to file
  modFileName <- "bayesMissMod.bug"
  fileConn <- file(modFileName )
  modelCode <- c(modelCode, "   }","")
  priorCode <- c(priorCode, "}")
  writeLines(c(modelCode, priorCode), fileConn)
  close(fileConn)

  print(paste("Your JAGS model has been written to: ", modFileName, sep=""))

  #run bayes analysis
  rScriptFileName <- "bayesMissRScript.r"
  fileConn <- file(rScriptFileName)
  rScriptText <- c(rScriptText, "jags.params <- c(\"beta\", \"outcome_tau\")")
  rScriptText <- c(rScriptText, paste("modelFit <- jags(data=jags.data, parameters.to.save=jags.params, n.iter=",nIter,
                   ", n.thin=1, model.file=\"",modFileName,"\", n.chains=",nChains,")", sep=""))
  rScriptText <- c(rScriptText, "modelFit")
  writeLines(rScriptText, fileConn)
  close(fileConn)

  print(paste("Your R script has been written to: ", rScriptFileName, sep=""))

  #run analysis
  #source(rScriptFileName, local=TRUE)
}
