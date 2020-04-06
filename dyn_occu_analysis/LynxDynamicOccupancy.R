####################################
##
##  R script to run dynamic occupancy model from:
##  Hostetter NJ, D Ryan, D Grosshuesch, T Catton, S Malick-Wahls, TA Smith, and B Gardner. 
##  Quantifying spatiotemporal occupancy dynamics and multi-year core use areas at a species range boundary. 
##  Diversity and Distributions.
##
##  Updated: 20-Jan-2020
##  
####################################

# Required packages
library(jagsUI)

# load formatted data
load(file = "~/Hostetter_etal_2020_DiversityAndDistributions/data/data.Rdata")

# data to run model
jags.data <- data[["jags.data"]]




##
## JAGS model
##

# Specify model in BUGS language
sink("Dynocc.jags")
cat("
model {

## priors

# intercepts
pbeta0.psi~dbeta(1,1)
pbeta0.gamma~dbeta(1,1) 
pbeta0.phi~dbeta(1,1)
pbeta0.p~dbeta(1,1)

beta0.psi<-logit(pbeta0.psi)
beta0.gamma<-logit(pbeta0.gamma)
beta0.phi<-logit(pbeta0.phi)
beta0.p<-logit(pbeta0.p)


#survival and colonization vary by year
for(t in 1:(nyear-2)){
 beta.phi.yr[t]~dnorm(0,0.01) 
 beta.gamma.yr[t]~dnorm(0,0.01) 
}
 beta.phi.yr[nyear-1]<- -1*(sum(beta.phi.yr[1:(nyear-2)]))
 beta.gamma.yr[nyear-1]<- -1*(sum(beta.gamma.yr[1:(nyear-2)]))

#habitat covariates for initial occu, persistence, and colonization
for(ss in 1:4){ 
 beta.psi[ss]~dnorm(0,0.01)   # covariates on initial occu
 beta.gamma[ss]~dnorm(0,0.01) # covariates on colonization
 beta.phi[ss]~dnorm(0,0.01)   # covariates on persistence
}

#continuous covariates on detection (effort, doy, doy^2)
for(ss in 1:3){  
 beta.p[ss]~dnorm(0,0.01)   
}

#snow conditions on detection
for(ss in 1:2){
 beta.p.snow[ss]~dnorm(0,0.01) 
}
 beta.p.snow[3]<- -1*(sum(beta.p.snow[1:2]))

#road conditions on detection
 beta.p.road[1]~dnorm(0,0.01) 
 beta.p.road[2]<- -1*(beta.p.road[1])


# Ecological process across all sites
for (s in 1:nsite){
   #Occasion 1: initial occu.
   z[s,1] ~ dbern(psi[s,1])
   psi[s,1]<-1/(1+exp(-(beta0.psi + beta.psi[1]*Spcteverg[s] + beta.psi[2]*SAvgSnowDepth[s] + 
                        beta.psi[3]*SLidarHigh[s] + beta.psi[4]*SLidarMid[s]  ))) 

   #Occasion >1: dynamic processes of persistence and colonization
   for (t in 2:nyear){
      z[s,t] ~ dbern(z[s,t-1]*phi[s,t-1] + (1-z[s,t-1])*gamma[s,t-1])
 
      # persistence
      phi[s,t-1]<-1/(1+exp(-(beta0.phi + beta.phi.yr[t-1] + beta.phi[1]*Spcteverg[s] + beta.phi[2]*SAvgSnowDepth[s] + 
                        beta.phi[3]*SLidarHigh[s] + beta.phi[4]*SLidarMid[s] ))) 

      # colonization
      gamma[s,t-1]<-1/(1+exp(-(beta0.gamma + beta.gamma.yr[t-1] + beta.gamma[1]*Spcteverg[s] + beta.gamma[2]*SAvgSnowDepth[s] + 
                        beta.gamma[3]*SLidarHigh[s] + beta.gamma[4]*SLidarMid[s]  ))) 

      # Monitor site- and occasion-specific occupancy  probability
      psi[s,t] <- psi[s,t-1]*phi[s,t-1] + (1-psi[s,t-1])*gamma[s,t-1]
 
      # monitor turnover rate as probability that a site changes occupancy status between years (MacKenzie et al. ed 2, pg 362)
      turnover[s,t-1] <- (psi[s,t-1]*(1 - phi[s,t-1]) ) +  ( (1-psi[s,t-1])*gamma[s,t-1] )

   } #t

  # monitor numbers of years a site was occupied 
  yrsOcc[s] <- sum(z[s,1:nyear]) 

 } #s

# Observation process: only loop over sites and occasions with surveys (speeds things up) 
for (s in 1:nSites){
 for(t in 1:nyrsSamp[s]){                      
  for(k in 1:nSurveysSamp[s,yearSite[s,t]]){                 # where nSurveysSamp defines the number of occasion surveyed at site s, in year yearSite[s,t]
    p[s,k,yearSite[s,t]]<-(1/(1+exp(-(beta0.p + beta.p.snow[SnowCond[s,k,yearSite[s,t]]] +    # Intercept and snow conditions
                           beta.p[1]*effort[s,k,yearSite[s,t]] +                              # effort
                           beta.p[2]*SDaySurv[s,k,yearSite[s,t]] +                            # day-of-year
                           beta.p[3]*SDaySurv2[s,k,yearSite[s,t]]+                            # day-of-year squared
                           beta.p.road[RoadCondMat[s,yearSite[s,t]]] ))))*                    # road conditions (good/poor)
                           z[surveyedCells[s],yearSite[s,t]]                                  # occupancy state

    y[s,k,yearSite[s,t]] ~ dbern(p[s,k,yearSite[s,t]])
  
  } #k
 } #t
} #s

#monitor numbers of occupied sites across entire study area
for(t in 1:nyear){
 n.occ[t] <- sum(z[1:nsite,t]) 
 pct.occ[t] <- n.occ[t]/nsite
}


}
",fill = TRUE)
sink()


##
## END JAGS model
##



# MCMC settings
nc <- 6; nAdapt=5000; nb <- 5000; ni <- 50000+nb; nt <- 5

# Initial values
zik <- data[["zik"]] # inital values for z to keep things from locking up
inits <- function(){ list(z = zik, 
                         pbeta0.psi=runif(1,.1,.25),beta.psi=runif(4,-.25,.25),
                         pbeta0.phi=runif(1,.75,1),beta.phi=runif(4,-.25,.25),beta.phi.yr=c(runif(3,-.25,.25),NA),
                         pbeta0.gamma=runif(1,0.1,.25),beta.gamma=runif(4,-.25,.25),beta.gamma.yr=c(runif(3,-.25,.25),NA),
                         pbeta0.p=runif(1,.1,.25),beta.p.snow=c(1,0,NA),beta.p.road=c(.5,NA), beta.p=c(runif(1,.5,1),runif(2,0,.5))
                     )}

# Parameters monitored
params <- c("beta0.psi","beta.psi", 
            "beta0.phi","beta.phi","beta.phi.yr",
            "beta0.gamma","beta.gamma","beta.gamma.yr",
            "beta0.p","beta.p","beta.p.snow","beta.p.road",  
            "n.occ","psi", "turnover","yrsOcc") 

# Call JAGS from R 
out <- jags(jags.data, inits, params, "Dynocc.jags", 
            n.adapt=nAdapt, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=T)


# Summarize posteriors
# e.g., 
round(out$summary[c("beta0.psi",paste("beta.psi[",1:4,"]", sep=""),
            "beta0.phi",paste("beta.phi[",1:4,"]", sep=""),paste("beta.phi.yr[",1:4,"]", sep=""),
            "beta0.gamma",paste("beta.gamma[",1:4,"]", sep=""),paste("beta.gamma.yr[",1:4,"]", sep=""),
            "beta0.p",paste("beta.p[",1:3,"]", sep=""),paste("beta.p.snow[",1:3,"]", sep=""),paste("beta.p.road[",1:2,"]", sep="")),
           c("50%", "2.5%", "97.5%","Rhat", "n.eff") ],3)
