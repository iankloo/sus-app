#
# This is a Plumber API. You can run the API by clicking
# the 'Run API' button above.
#
# Find out more about building APIs with Plumber here:
#
#    https://www.rplumber.io/
#

library(plumber)
library(data.table)
library(rstan)
library(resample)

source('R/functions.R')

#' function to make sure requests aren't blocked - known issue with plumber
#' @filter cors
cors <- function(req, res) {
  res$setHeader("Access-Control-Allow-Origin", "*")
  if (req$REQUEST_METHOD == "OPTIONS") {
    res$setHeader("Access-Control-Allow-Methods","*")
    res$setHeader("Access-Control-Allow-Headers", req$HTTP_ACCESS_CONTROL_REQUEST_HEADERS)
    res$status <- 200 
    return(list())
  } else {
    plumber::forward()
  }
  
}

#* Process SUS data
#* @post /sus
function(req) {
  dat <- jsonlite::fromJSON(req$postBody)
  dat <- data.table(dat)
  
  #get sus scores
  df_sub <- dat[, c(2:ncol(dat)), with = FALSE]
  cols <- colnames(df_sub)
  df_sub <- data.table(apply(df_sub, 2, as.numeric))
  colnames(df_sub) <- cols
  odds <- seq(1, ncol(df_sub) - 1, by = 2)
  evens <- seq(2, ncol(df_sub), by = 2)
  sus_scores <- apply(df_sub, 1, sus_converter, odds, evens)
  
  
  sub <- list(N=length(sus_scores),
              J=1,#always set to 1 if using single sus result
              y=sus_scores,
              g=rep(1,length(sus_scores)))#always rep 1 if single sus result
  
  new.fit <- stan(file = 'R/stanmod.stan', data = sub,refresh=0)
  #mod <- readRDS('stanmod.rds')
  #new.fit <- stan(model_code = 'mod', data = sub,refresh=0)
  
  #samps <- as.vector(extract(new.fit)$mu)
  #lower<-summary(new.fit, pars = c("mu"), probs = c(0.05, 0.95))$summary[[4]]
  #upper<-summary(new.fit, pars = c("mu"), probs = c(0.05, 0.95))$summary[[5]]
  
  #fit <- stan(file = 'BayesCode.stan', data = bayes.dat,refresh=0)

  bayes.est<-rstan::extract(new.fit,pars="mu[1]")$`mu[1]`
  sig.est<-rstan::extract(new.fit,pars="sigma")$sigma
  alpha=(0-bayes.est)/sig.est
  beta=(100-bayes.est)/sig.est
  ex.val<-bayes.est + 
    sig.est*(dnorm(alpha)-dnorm(beta))/
    (pnorm(beta)-pnorm(alpha))

  bayes.ci<-quantile(ex.val,probs=c(.025,.975))

  
  ci <- matrix(c(bayes.ci[1], bayes.ci[2]), nrow = 1)
  bayes <- list(sus_scores = sus_scores, ci = ci[1,], replicates = ex.val, mean = mean(ex.val))
  
  
  
  N.bootstrap <- 2000
  ci.hw.increase <- 0
  b <- bootstrap(sus_scores, mean(sus_scores), R = N.bootstrap)  
  ci <- CI.bca(b, probs = c(0.025-ci.hw.increase, 0.975+ci.hw.increase), expand = TRUE)
  boot <- list(sus_scores = sus_scores, ci = ci[1,], replicates = b$replicates, mean = mean(sus_scores))

  #out <- list(sus_scores = sus_scores, ci = ci[1,], replicates = samps)
  out <- list(bayes, boot)
  return(out)
}
