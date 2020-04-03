# simulation study to see if the ratio estimator is working properly

# Set up parallelization engines
doParallel = FALSE
if(doParallel) {
  library(doMC)  # for parallel model fitting
  #library(foreach)

  # see http://viktoriawagner.weebly.com/blog/five-steps-to-parallel-computing-in-r
  detectCores() 
  cl <- makeCluster(2)
  # Need to export some libraries to the cluster
  # see http://stackoverflow.com/questions/18981932/logging-with-plyr-in-parallel-true
  clusterEvalQ(cl, library(unmarked))
  registerDoMC(3) 
} 


library(SightabilityModel)

source("SightabilityModel/R/Sight.Est.Ratio.R")
source("SightabilityModel/R/Wong.Est.Ratio.R")


library(ggplot2)
library(plyr)

set.seed(12343343)
# First a single stratum of with 50 blocks, sampling 10 blocks, sightability of .5 (binomial)

N.blocks <- 50
prob.sight <- 0.5

population <- plyr::ldply(1:N.blocks, function(block){
  # number of groups is uniform between 10 and 20
  n.groups <- trunc(runif(1, min=10, max=21))
  # number of bulls is uniform between 5 and 10
  n.bulls  <- trunc(runif(n.groups, min=5, max=10))
  n.cows   <- rgeom(n.bulls, .2)
  # sightability of 0.5 using binomial
  n.bulls.obs <- rbinom(n.groups, n.bulls, prob.sight)
  n.cows.obs  <- rbinom(n.groups, n.cows , prob.sight)
  
  # total animals observed
  n.total  <- n.bulls + n.cows
  n.total.obs <- n.bulls.obs + n.cows.obs
  stratum =1 
  data.frame(stratum=stratum, block=block, group=1:n.groups,
             n.bulls, n.cows, n.bulls.obs, n.cows.obs,
             n.total, n.total.obs)
})

# actual values
head(population)
nrow(population)

length(unique(population$block))
sum(population$n.total)
sum(population$n.bulls)
sum(population$n.cows)




# now for the simulation study
n.sim <- 1000

n.sightability <- 100
n.sample.blocks<- N.blocks*.2

sim.res <- plyr::llply(1:n.sim, function(sim,population, n.sightability, n.sample.blocks){
   stratum.info <- plyr::ddply(population, "stratum", plyr::summarize,
                            Nh=length(unique(block)))
   # sightability trials. Random sample of groups
   sight.trials <- population[sample(1:nrow(population), n.sightability, replace=TRUE),]
   fit <- glm( cbind(sight.trials$n.total.obs, sight.trials$n.total-n.total.obs) ~ 1, data=sight.trials, family=binomial(link=logit))
   
   beta <- coef(fit)
   beta.cov <- vcov(fit)
   
   # random sample of blocks
   block.sample <- sample(1:length(unique(population$block)), n.sample.blocks)
   obs.data <- population[ population$block %in% block.sample,]
   obs.data$total       <- obs.data$n.total.obs
   obs.data$numerator   <- obs.data$n.bulls.obs
   obs.data$denominator <- obs.data$n.cows.obs
   obs.data$subunit     <- obs.data$block
   
   # number of blocks
   sample.stratum.info <- plyr::ddply(obs.data, "stratum", plyr::summarize,
                                      nh=length(unique(block)))
   sample.stratum.info <- merge(sample.stratum.info, stratum.info)
   # make sure that totals are estimated properly
   #   browser()
   total.est <- Sight.Est( observed ~1, odat=obs.data, 
                           sampinfo=sample.stratum.info,
                           bet=beta,
                           varbet = beta.cov)
   
   # ratio of bulls to cows
   ratio.est <- Sight.Est.Ratio(observed ~1, odat=obs.data, 
                           sampinfo=sample.stratum.info,
                           bet=beta,
                           varbet = beta.cov)
   #browser()
   list(total.est=total.est, ratio.est=ratio.est)
}, population=population, n.sightability=n.sightability, n.sample.blocks=n.sample.blocks, .parallel=doParallel)

if(doParallel) stopCluster(cl) # stop parallel processing


# extract estimates from simulation study

sim.est <- plyr::ldply(sim.res, function(x){
  #browser()
  res <- data.frame(total.hat      =x$total.est$est["tau.hat"],
                    total.hat.se    =sqrt(x$total.est$est["VarTot"]),
                    bulls.hat      =x$ratio.est$numerator$est["tau.hat"],
                    bulls.hat.se   =sqrt(x$ratio.est$numerator$est["VarTot"]),
                    cows.hat       =x$ratio.est$denominator$est["tau.hat"],
                    cows.hat.se    =sqrt(x$ratio.est$denominator$est["VarTot"]),
                    ratio.hat      =x$ratio.est$est["ratio.hat"],
                    ratio.hat.se   =sqrt(x$ratio.est$est["VarRatio"]))
  row.names(res) <- NULL
  res
})

head(sim.est)


compare.sim <- function(df,actual, decdigit=0){
   # given the estimate and estiamted se compare to the true value
   # and compare the mean variance to actual variance
   # and compare the mean se to the actual se
   mean.est <- mean(df[,1])
   per.bias     <- (mean.est - actual)/actual*100
   mean.var <- mean(df[,2]^2)
   actual.var<- var(df[,1])
   var.ratio   <- mean.var/actual.var
   actual.se <- sd(df[,1])
   mean.se   <- mean(df[,2])
   se.ratio <- mean.se/actual.se
   cat("Mean of estimates ",mean.est, " vs. actual ", actual, " with %bias ", per.bias, "\n")
   cat("Mean of var est   ",mean.var, " vs. actual var", actual.var, "with ratio of ", var.ratio, "\n")
   cat("Mean of se est    ",mean.se , " vs actual se  ", actual.se , "with ratio of ", se.ratio, "\n")
}

compare.sim( sim.est[,c("total.hat","total.hat.se")], actual=sum(population$n.total))
compare.sim( sim.est[,c("bulls.hat","bulls.hat.se")], actual=sum(population$n.bulls))
hist(sim.est$bulls.hat)
hist(sim.est$bulls.hat.se^2)

compare.sim( sim.est[,c("cows.hat" ,"cows.hat.se" )], actual=sum(population$n.cows))
compare.sim( sim.est[,c("ratio.hat" ,"ratio.hat.se" )], actual=sum(population$n.bulls)/sum(population$n.cows))

# look at the components and total covariance value of numerator-hat and denominator-hat
cov.components.est <- plyr::ldply(sim.res, function(x){
  #browser()
  res <- x$ratio.est$cov.nhat.dhat.components
  row.names(res) <- NULL
  res
})



head(cov.components.est)

# check the individual components of the covariance term
# Term 1
mean(cov.components.est$E2.4.6.t1)
mean(sim.est$bulls.hat*sim.est$cows.hat)

# Term 2
mean(cov.components.est$E2.4.6.t2)
sum(population$n.bulls*population$n.cows)

# term3
mean(cov.components.est$E2.4.6.t3)
sum(plyr::daply(population, "block", function(x){
   temp <- outer(x$n.bulls, x$n.cows, FUN="*")
   diag(temp)<- 0
   sum(temp)
})
)

# term 4
mean(cov.components.est$E2.4.6.t4)
temp <- 0
for(i in 1:N.blocks){
  for(j in 1:N.blocks){
     if(i != j){
        bulls <- population$n.bulls[ population$block==i]
        cows  <- population$n.cows [ population$n.cows==j]
        temp <- temp + sum(outer(bulls, cows, FUN="*"))
     }
  }
}
temp    

# sum of 2nd, 3rd, 4th component
mean( cov.components.est$E2.4.6.t2+cov.components.est$E2.4.6.t3+cov.components.est$E2.4.6.t4)
sum(population$n.bulls)*sum(population$n.cows)

# check the covariance
mean(cov.components.est$cov.nhat.dhat)
cov(sim.est$bulls.hat, sim.est$cows.hat)

