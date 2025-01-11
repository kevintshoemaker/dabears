# Bear project script ---------------------
# This script runs a CMR model for northern NV black bear population in JAGS
#    and summarizes model results: abundance, natural mortality, road mortality, management removals. 

rm(list=ls())

# load packages ------------------------------

library(jagsUI)   # need to install jags first
library(lubridate)
library(Hmisc)
library(IPMbook)
library(R2ucare)
library(coda)

# read in data  ---------------------------

load("beardata_bundle_v6.RData")

# run ucare tests on live captures and release...

test3sr(as.matrix(caphist_liverel),rep(1,times=nrow(caphist_liverel)))   # slight transience detected (could be mortalities?)
test2ct(as.matrix(caphist_liverel),rep(1,times=nrow(caphist_liverel)))   # weak 'trap dependence' detected (should we model this?)
test3sm(as.matrix(caphist_liverel),rep(1,times=nrow(caphist_liverel)))   # no overdispersion detected
test2cl(as.matrix(caphist_liverel),rep(1,times=nrow(caphist_liverel)))   # no overdispersion detected... interesting!


# load functions ------------------------------

# thisvar=allvars[3]
visualizeRelation <- function(thisvar,thisproc,bottom=T){
  thisvar2 <- allvars2[allvars==thisvar]
  
  allproc2 <- c("Detection prob", "Nat. survival (phi)",
                "Road mort prob", "Mgmt Removal prob")
  ylab <- allproc2[which(allproc==thisproc)]
  
  factor <- ifelse(thisvar%in%c("snotel","freezedate"),F,T)
  
  if(thisvar=="freezedate")  xlab = "Timing of last frost (ordinal date)"
  if(thisvar=="snotel") xlab = "Snowpack (tenths of an inch)"
  if(thisvar=="is.male"){
    xlab = "Sex"; xticks = c("Male","Female")
  } 
  if(thisvar=="is.juv"){
    xlab = "Age"; xticks = c("Juvenile","Adult")
  } 
  
  if(factor){
    xvec <- c(0,1) 
  }else if(thisvar=="freezedate"){
    xvec <- seq(min(data.for.jags[[thisvar]]),max(data.for.jags[[thisvar]]),length=100)
    realmean <- mean(firstfreeze1$jday)   # NOTE: fix this to grab the real means and standard deviations
    realsd <- sd(firstfreeze1$jday)
    xvecr <- xvec*realsd + realmean
    # month.name[month(ymd("2000-01-01")+days(floor(xvec*realsd + realmean)))]
  }else{
    xvec <- seq(min(data.for.jags[[thisvar]]),max(data.for.jags[[thisvar]]),length=100)
    realmean <- mean(SWE[-1])   # NOTE: fix this to grab the real means and standard deviations
    realsd <- sd(SWE[-1])
    xvecr <- xvec*realsd + realmean
  }
  
  thisint <- qlogis(mod$sims.list[[sprintf("%s0",thisproc)]])
  thisslop <- mod$sims.list[[sprintf("%s.%s.eff",thisproc,thisvar2)]]
  # if(thisproc=="p"){
  #   thisyreff <- mod$sims.list[["p.raneff.yr"]][,5]
  # }else{
  #   thisyreff <- numeric(length(thisint))
  # }
  
  predmat <- sapply(xvec,function(t) plogis(thisint + thisslop*t  )  )  # + thisyreff
  dim(predmat)
  
  qs <- t(apply(predmat,2,function(t) quantile(t,c(0.025,0.5,0.975)) ))
  
  if(factor){
    bars <- qs[,2][c(2,1)] 
    #names(bars) <- rev(as.character(as.logical(xvec)))
    names(bars) <- xticks
    x <- barplot(bars,ylim=c(0,range(qs)[2]*1.3),ylab=ylab,xlab=xlab)
    errbar(x,qs[,2][c(2,1)],qs[,3][c(2,1)],qs[,1][c(2,1)],add=T)   #xlim=c(-0.5,1.5)
    #axis(1,at=c(0,1),labels=c("no","yes") )
  }else if(thisvar=="freezedate"){
    if(bottom==F) xlab=""
    plot(xvec,qs[,2],type="l",lwd=2,col=gray(0.2),ylim=c(min(qs),max(qs)),
         xlab=xlab,ylab=ylab,xaxt="n")
    lines(xvec,qs[,1],lty=2)
    lines(xvec,qs[,3],lty=2)
    thissq <- seq(min(xvec),max(xvec),length=5)
    thissq2 <- thissq*realsd + realmean
    if(bottom==T) axis(1,at=thissq,labels=month.name[month(ymd("2000-01-01")+days(floor(thissq2*realsd + realmean)))] )
    if(bottom==F) axis(1,at=thissq,labels=rep("",times=length(thissq)) )
    rug(data.for.jags[[thisvar]])
  }else{
    if(bottom==F) xlab=""
    plot(xvec,qs[,2],type="l",lwd=2,col=gray(0.2),ylim=c(min(qs),max(qs)),
         xlab=xlab,ylab=ylab,xaxt="n")
    lines(xvec,qs[,1],lty=2)
    lines(xvec,qs[,3],lty=2)
    thissq <- seq(min(xvec),max(xvec),length=5)
    thissq2 <- thissq*realsd + realmean
    if(bottom==T) axis(1,at=thissq,labels=round(thissq2) )
    if(bottom==F) axis(1,at=thissq,labels=rep("",times=length(thissq)) )
    rug(data.for.jags[[thisvar]])
  }
  
}

# make sure individuals are marked as having been observed on their first year of capture...

sapply(1:ninds, function(t) caphist_ms[t,max(1,firsts[t]-1)] )  # okay,this is true...
firsts

# JAGS code ----------------------------------

fn = "cmr_jags_bear_v6.txt"

cat("

model{

### capture probability model  (prob of bear detected in urban setting)

for(yr in 1:(nyears)){
  p0[yr] ~ dunif(0,1)  # mean detection prob within surveyed areas
  logit.p0[yr] <- log(p0[yr]/(1-p0[yr]))    # convert to logit (intercept for logit linear model with standardized covariates)
}

p.male.eff ~ dnorm(0,1)          # effect of sex on mean capture probability (slightly regularized prior)
p.juv.eff ~ dnorm(0,1)           # effect of age on mean capture probability (slightly regularized prior)
p.snot.eff ~ dnorm(0,1)          # effect of snotel on mean capture probability (slightly regularized prior)
p.freeze.eff ~ dnorm(0,1)        # effect of freeze date on mean capture probability (slightly regularized prior)
p.prev.eff ~ dnorm(0,1)          # effect of being previously detected (with slightly regularized prior)

for(i in 1:ninds){  # loop through juveniles for state space model
  logit(p[i,firsts[i]]) <- logit.p0[firsts[i]] + p.male.eff*is.male[i] + p.prev.eff + 
          p.snot.eff*snotel[firsts[i]] + p.freeze.eff*freezedate[firsts[i]] + 
          p.juv.eff*is.juv[i,firsts[i]]  
  for(yr in (firsts[i]+1):nyears){
    isprev[i,yr] <- ifelse(y[i,yr-1]==2,1,0)  # my first use of ifelse in JAGS, hope it works! 
    logit(p[i,yr]) <- logit.p0[yr] + p.male.eff*is.male[i] + p.prev.eff*isprev[i,yr] + 
          p.snot.eff*snotel[yr] + p.freeze.eff*freezedate[yr] + p.juv.eff*is.juv[i,yr]  
  }
}

### natural survival process

  # survival priors
phi0 ~ dunif(0.3,.96)                 # mean/intercept survival
logit.phi0 <- log(phi0/(1-phi0))
phi.male.eff ~ dnorm(0,1)          # effect of sex on mean survival  (slightly regularized prior)
phi.juv.eff ~ dnorm(0,1)           # effect of age on mean survival  (slightly regularized prior)
phi.urb.eff ~ dnorm(0,1)           # effect of being urban on mean survival  (slightly regularized prior)
phi.snot.eff ~ dnorm(0,1)          # effect of max snow depth on survival (slightly regularized prior)
phi.freeze.eff ~ dnorm(0,1)        # effect of freeze date on survival

# TODO: add random year effect?

      #### survival model (logit-linear)
for(i in 1:ninds){  # loop through individuals
  for(yr in firsts[i]:(nyears-1)){
    logit(phi[i,yr]) <- logit.phi0 + phi.male.eff*is.male[i] +            # this should thought be interpreted as natural survival/mort                          
                                         phi.juv.eff*is.juv[i,yr] + phi.snot.eff*snotel[yr] + phi.freeze.eff*freezedate[yr]
     
  }
}

# Prob of management removal: prob of being captured and removed each year  (use similar structure for prob of road mortality)

pmr0 ~ dunif(0,.4)                  # mean prob of management removal
logit.pmr0 <- log(pmr0/(1-pmr0))   # convert to logit (intercept for logit linear model with standardized covariates)
pmr.male.eff ~ dnorm(0,1)          # effect of sex on management removal probability  (slightly regularized prior)
pmr.juv.eff ~ dnorm(0,1)           # effect of age on management removal probability  (slightly regularized prior)
pmr.snot.eff ~ dnorm(0,1)          # effect of snotel on management removal probability (slightly regularized prior)
pmr.freeze.eff ~ dnorm(0,1)        # effect of freeze date on pmr

for(i in 1:ninds){
  for(yr in (firsts[i]+1):nyears){
    logit(pmr[i,yr])  <- logit.pmr0 + pmr.male.eff*is.male[i] + pmr.snot.eff*snotel[yr] + pmr.freeze.eff*freezedate[yr] +  
                             pmr.juv.eff*is.juv[i,yr]      # probability of management remove conditional on being alive
  }
}

### Prob of road mortality

rmort0 ~ dunif(0,0.2)               # mean prob of road mortality
logit.rmort0 <- log(rmort0/(1-rmort0))   # convert to logit (intercept for logit linear model with standardized covariates)
rmort.male.eff ~ dnorm(0,1)          # effect of sex on road mortality probability  (slightly regularized prior)
rmort.juv.eff ~ dnorm(0,1)           # effect of age on road mortality probability  (slightly regularized prior)
rmort.snot.eff ~ dnorm(0,1)          # effect of snotel on road mortality probability (slightly regularized prior)
rmort.freeze.eff ~ dnorm(0,1)        # effect of freeze date on road mort prob

for(i in 1:ninds){
  for(yr in (firsts[i]+1):nyears){
    logit(rmort[i,yr])  <- logit.rmort0 + rmort.male.eff*is.male[i] + rmort.snot.eff*snotel[yr] +  rmort.freeze.eff*freezedate[yr] +
                             rmort.juv.eff*is.juv[i,yr] #+ [OTHER EFFECTS??]     # probability of management remove conditional on being alive
  }
}

p1[1] <- 1    # alive
p1[2] <- 0    # removed by management
p1[3] <- 0    # dead on road
p1[4] <- 0    #

# LIKELIHOOD (multi-state)

for(i in 1:ninds){
    
  z[i,firsts[i]] ~ dcat(p1)   # DATA NODE, is coded as 0 in year of known death and all years thereafter

  for(yr in (firsts[i]+1):lasts[i]){
  
    p.z[1,i,yr,1] <- 1-sum(p.z[1,i,yr,2:4]) #  # transition from alive last year to alive this year  phi[i,y-1] * (1-rmort[i,y]) * (1- (p[i,y] * pmr[i,y]) )  
    p.z[1,i,yr,2] <- p[i,yr] * pmr[i,yr]   # transition from alive last year to captured and removed this year
    p.z[1,i,yr,3] <- rmort[i,yr]  #*(1-p[i,yr] * pmr[i,yr])    # transition from alive last year to dead on road this year
    p.z[1,i,yr,4] <- (1-phi[i,yr-1])   # *(1- rmort[i,yr]*(1-p[i,yr] * pmr[i,yr]))   # probability of natural mortality
    p.z[2,i,yr,1] <- 0   # removed this year to alive next year
    p.z[2,i,yr,2] <- 0
    p.z[2,i,yr,3] <- 0
    p.z[2,i,yr,4] <- 1
    p.z[3,i,yr,1] <- 0
    p.z[3,i,yr,2] <- 0
    p.z[3,i,yr,3] <- 0
    p.z[3,i,yr,4] <- 1
    p.z[4,i,yr,1] <- 0
    p.z[4,i,yr,2] <- 0
    p.z[4,i,yr,3] <- 0
    p.z[4,i,yr,4] <- 1
    
    z[i,yr] ~ dcat(p.z[z[i,yr-1],i,yr,] )   # this year's state
    
    p.y[1,i,yr,1] <- 1-p[i,yr]                     # prob of being not seen if alive
    p.y[1,i,yr,2] <- p[i,yr]                       # prob of being seen and released if alive
    p.y[1,i,yr,3] <- 0                            # prob of being seen and removed if alive
    p.y[1,i,yr,4] <- 0                            # prob of being seen dead on road if alive (assume all DOR are reported)
    p.y[2,i,yr,1] <- 0
    p.y[2,i,yr,2] <- 0
    p.y[2,i,yr,3] <- 1                            # 100% chance of been seen when removed from population  (management removal)
    p.y[2,i,yr,4] <- 0
    p.y[3,i,yr,1] <- 0
    p.y[3,i,yr,2] <- 0
    p.y[3,i,yr,3] <- 0
    p.y[3,i,yr,4] <- 1                            # 100% chance of being seen DOR if dead on road
    p.y[4,i,yr,1] <- 1                            # 100% chance of being not seen if long dead
    p.y[4,i,yr,2] <- 0
    p.y[4,i,yr,3] <- 0
    p.y[4,i,yr,4] <- 0
    
    y[i,yr] ~ dcat(p.y[z[i,yr],i,yr,])
    
  }
}


### Compute abundance using H-T estimator

for(i in 1:ninds){
  for(yr in 1:(firsts[i]-1)){
    invpcap[i,yr] <- 0
  }
  # invpcap[i,firsts[i]] <- 1/p[i,firsts[i]]
  for(yr in (firsts[i]):nyears){
    invpcap[i,yr] <- 1/p[i,yr] # prob captured and released given alive this year
  }
}

for(yr in 1:nyears){
  N[yr] <- inprod(caphist_liverel[1:ninds,yr],invpcap[1:ninds,yr])
}

}
    
",file=fn)



# PACKAGE DATA FOR JAGS -----------------------

firstfreeze <- read.table("RenoFreezeDates.txt", skip = 3, sep = ",", header = T)

firstfreeze$lastdate <- mdy(firstfreeze$Last)

head(firstfreeze)
firstfreeze$jday <- yday(firstfreeze$lastdate)
firstfreeze$year <- year(firstfreeze$lastdate)
yearswewant <- 1998:2022
firstfreeze1 <- firstfreeze[which(firstfreeze$year%in%yearswewant),]

# write.csv(firstfreeze1, "freezedates_for_kelley.csv",row.names = F)

caphist_ms[is.na(caphist_ms)] = 1
data.for.jags = list(
  ninds = ninds,
  nyears = nyears,
  y = caphist_ms,
  caphist_liverel = caphist_liverel,
  is.male = is.male,
  is.juv = is.juv,
  snotel = as.numeric(scale(SWE)), 
  freezedate = as.vector(scale(firstfreeze1$jday)),
  z = status,
  firsts = firsts,
  lasts = lasts
)

# SET PARAMS TO MONITOR -----------------
params.to.monitor = c(
  "p0",
  "p.male.eff",
  "p.snot.eff",
  "p.freeze.eff",
  "p.juv.eff",
  "p.prev.eff",
  "phi0",
  "phi.male.eff",
  "phi.juv.eff",
  "phi.snot.eff",
  "phi.freeze.eff",
  "N",
  "pmr0",
  "pmr.male.eff",
  "pmr.juv.eff",
  "pmr.snot.eff",
  "pmr.freeze.eff",
  "rmort0",
  "rmort.male.eff",
  "rmort.juv.eff",
  "rmort.snot.eff",
  "rmort.freeze.eff"
)


# INITIALIZATION FOR JAGS ------------------

inits.for.jags <- function(){
  list(
    p0 = runif(nyears,0.3,0.4),
    p.male.eff = rnorm(1,0,0.05),
    p.snot.eff = rnorm(1,0,0.05),
    p.freeze.eff = rnorm(1,0,0.05),
    p.juv.eff = rnorm(1,0,0.05),
    p.prev.eff = rnorm(1,0,0.05),
    phi0 = runif(1,0.6,0.7),
    phi.male.eff = rnorm(1,0,0.1),
    phi.juv.eff = rnorm(1,0,0.1),
    phi.snot.eff = rnorm(1,0,0.1),
    phi.freeze.eff = rnorm(1,0,0.1)
  )
}

# inits.for.jags()


# RUN JAGS  ---------------------

# ?jags

debug=F

if(debug){
  mod <- jags(data.for.jags, inits.for.jags, params.to.monitor, fn,
              n.chains=1, n.adapt=100, n.iter=500, n.burnin=100, n.thin=1,
              parallel=FALSE, n.cores=1, verbose=TRUE)
}else{
  # real model run
  
  mod <- jags(data.for.jags, inits.for.jags, params.to.monitor, fn,
              n.chains=3, n.adapt=1000, n.iter=10000, n.burnin=5000, n.thin=5,
              parallel=TRUE, n.cores=3, verbose=TRUE)
  save(mod,file=sprintf("JagsResults_v6_%s.RData",Sys.Date() ))
}

# load("JagsResults_v6_2024-02-11.RData")


  # rm(mod, mod_freeze, mod_snow)

  # load("JagsResults_2023-02-01.RData")
  # mod_freeze <- mod
  # load("JagsResults_2023-05-19_snow.RData")   
  # mod_snow <- mod

# visualize results and assess convergence -----------------------

# mod<-mod_freeze
# mod$DIC   # DIC for freeze model is 9176.544
# mod<-mod_snow
# mod$DIC   # DIC for freeze and snow model is 9172.044 

# try model with snow depth for 

samples <- mod$samples

mod$DIC    # 3327 for model with firstyears included, both snowpack and final frost
           # 2822 for model with firstyears removed, both snowpack and final frost

# mod$samples
# mod
# "p0",
# "p.male.eff",
# "p.juv.eff",
# "urb0",
# "urb.snow.eff",
# "phi0",
# "phi.male.eff",
# "phi.juv.eff",
# "phi.bcs.eff",
# "phi.urb.eff",
# "alive"


plot(samples[,"N[10]"])

plot(samples[,"p0[9]"])
plot(samples[,"p.male.eff"])   # males much less likely to be captured
plot(samples[,"p.juv.eff"])    # slightly negative after removing first years
plot(samples[,"p.snot.eff"])  # nothing- could be removed
plot(samples[,"p.freeze.eff"])  # could be removed
plot(samples[,"p.prev.eff"])     # strong effect of previous capture

# plot(samples[,"p.bcs.eff"])


# plot(samples[,"urb0"])

# plot(samples[,"urb.snow.eff"])

plot(samples[,"phi0"])
plot(samples[,"phi.male.eff"])  # near zero, not great convergence, could be removed
plot(samples[,"phi.juv.eff"])   # lower survival for juveniles
plot(samples[,"phi.snot.eff"])     # not great convergence- slightly higher survival but corsses zero...  could be removed?
plot(samples[,"phi.freeze.eff"])   # not great convergence- slightly negative

# plot(samples[,"phi.bcs.eff"])

# plot(samples[,"phi.urb.eff"])
# plot(samples[,"alive"])
#plot(samples[,"alive[1,5]"])

# plot(samples[,"p.year.prec"])

# plot(samples[,"p.raneff.yr[19]"])

plot(samples[,"pmr0"])
plot(samples[,"pmr.male.eff"])
plot(samples[,"pmr.juv.eff"])  # juveniles more likely to be removed
plot(samples[,"pmr.snot.eff"])  # positive, overlaps with zero 
plot(samples[,"pmr.freeze.eff"]) # very strong positive relationship

plot(samples[,"rmort0"])     
plot(samples[,"rmort.male.eff"])  # males more likely to succumb to road mort
plot(samples[,"rmort.juv.eff"])   # juveniles more likely to be killed on roads
plot(samples[,"rmort.snot.eff"])    # no effect of snowpack
plot(samples[,"rmort.freeze.eff"])   # slight positive effect of freeze


#gelman.diag(samples)

# coda::gelman.diag(samples[,"urb0"])
coda::gelman.diag(samples[,"phi0"])
# gelman.diag(samples[,"urb0"])


# abundance over time


nyears

names(mod$sims.list)

library(Hmisc)
library(coda)

abundmat <- array(NA, dim=c(nyears,3))
y=1
for(y in 1:nyears){
  abundmat[y,1] <- median(mod$sims.list$N[,y])
  temp <- HPDinterval(as.mcmc(mod$sims.list$N[,y]))
  abundmat[y,2] <- temp[,"lower"]
  abundmat[y,3] <- temp[,"upper"]
}

# lines with ci
png("figs/abund_trajectory_HT2.png",5,4,units = "in",res = 300)
plot(1:nyears,abundmat[,1],type="l",lwd=3,ylim=c(0,800),xaxt="n",ylab="Abundance",xlab="Year" )
# lines(1:nyears,abundmat[,2],lty=2)
# lines(1:nyears,abundmat[,3],lty=2)
Hmisc::errbar(1:nyears,abundmat[,1],abundmat[,3],abundmat[,2],add=T)
axis(1,at=1:nyears,labels=1998:2022)
dev.off()

colnames(abundmat) <- c("med","lb","ub")
write.csv(abundmat, "newabund2.csv",row.names = F)   # write to CSV for Kelley

## Here we probably want to start coding and visualizing the different effects 
## the effects plots 
## we want to look at the plots 
## add some tyoe of table to look at 
# ?errbar

# Univariate effects plots ---------------------------

names(mod$sims.list)

allvars <- c("is.male","is.juv","snotel","freezedate")    # list of covariates
allvars2 <- c("male","juv","snot","freeze")
allproc <- c("p","phi","rmort","pmr")        # list of processes

thisvar <- allvars[1]
thisproc <- allproc[1]
factor <- T

# visualize all relationships
v=1;p=1
for(v in 1:length(allvars)){
  for(p in 1:length(allproc)){
    thisfile <- sprintf("figs/univariateEffPlot_%s_%s.svg",allproc[p],allvars[v])
    svg(thisfile,4,4)
     visualizeRelation(allvars[v],allproc[p])
    dev.off()
  }
}


# make figure with all frost and SWE relationships
{
  thisfile <- "figs/allEffects2.png"
  png(thisfile,5,6,units = "in",res = 300)

  layout(matrix(c(1:6),ncol=2,byrow=T),
         heights=c(10,7,10))
  
  par(mai=c(0.1,0.7,0.8,0.1))
  visualizeRelation(allvars[4],allproc[2],bottom=F)
  mtext("Freeze Date",line=2)
  visualizeRelation(allvars[3],allproc[2],bottom=F)
  mtext("SWE",line=2)
  
  par(mai=c(0.1,0.7,0.1,0.1))
  # visualizeRelation(allvars[4],allproc[2],bottom=F)
  # visualizeRelation(allvars[3],allproc[2],bottom=F)
  
  visualizeRelation(allvars[4],allproc[3],bottom=F)
  visualizeRelation(allvars[3],allproc[3],bottom=F)
  
  par(mai=c(0.8,0.7,0.1,0.1))
  visualizeRelation(allvars[4],allproc[4],bottom=T)
  visualizeRelation(allvars[3],allproc[4],bottom=T)
  dev.off()
}


# Make table of all parameters --------------------

params.to.monitor

paramsTable <- data.frame(
  parm = params.to.monitor
)
paramsTable$median <- 0
paramsTable$LB <- 0
paramsTable$UB <- 0

paramsTable <- paramsTable[!paramsTable$parm%in%c("N","alive","p.raneff.yr"),]

library(HDInterval)
r=1
for(r in 1:nrow(paramsTable)){
  # grab posterior
  thisparam <- paramsTable$parm[r]
  thispost <- mod$sims.list[[thisparam]]
  paramsTable$median[r] <- median(thispost)
  thishpd <- HDInterval::hdi(as.vector(thispost))
  paramsTable$LB[r] <- thishpd["lower"]
  paramsTable$UB[r] <- thishpd["upper"]
}

write.csv(paramsTable,file="paramsTable2.csv",row.names=F)


paramsTable3 <- data.frame(parm=params.to.monitor)
paramsTable3$text <- character(nrow(paramsTable3))
paramsTable3 <- paramsTable3[!paramsTable3$parm%in%c("N","alive","p.raneff.yr"),]

r=1
for(r in 1:nrow(paramsTable)){
  # grab posterior
  thisparam <- paramsTable3$parm[r]
  this <- format(c(paramsTable$median[r],paramsTable$LB[r],paramsTable$UB[r]),digits=2,nsmall=1,trim=T)
  paramsTable3$text[r] <- sprintf("%s (%s - %s)",this[1],this[2],this[3])
}

write.csv(paramsTable3,file="paramsTable3.csv",row.names=F)

## plot management removals as a function of snotel data

# allyears
# data.for.jags$snotel
# 
# snotel <- data.for.jags$snotel
# names(snotel) <- allyears
# snotel
# 
# removed <- apply(data.for.jags$removed,2,sum,na.rm=T) 
# names(removed) <- allyears
# removed
# 
# hbc <- apply(data.for.jags$HBC,2,sum,na.rm=T) 
# names(hbc) <- allyears
# hbc
# 
# 
# plot(snotel,removed,type="p")
# 
# plot(snotel,hbc,type="p")








### how can we start to plot and save figures 
# save.image()

## summaries for paper ----------


quantile(mod$sims.list$p0,c(0.05,0.5,0.95) )


