rm(list=ls())

setwd("~/Dropbox/GitHub/Nevada-Black-Bear-Abundance")
load(paste0(getwd(),"/Results/Data_Analysis_Workspace.RData"))


################################################################################
### Plot
################################################################################

##
## N
##

par(mfrow=c(1,1))

# ## Males
# plot(years,unique(estM),type='l', ylim=c(0,400),col=4)
# lines(years,unique(lclM),col=4,lty=2)
# lines(years,unique(uclM),col=4,lty=2)
# 
# ## Females
# lines(years,unique(estF),type='l',col=2)
# lines(years,unique(lclF),col=2,lty=2)
# lines(years,unique(uclF),col=2,lty=2)
# 
# ## 2022 estimate
# tail(unique(estF)+unique(estM),1)
# tail(unique(lclF)+unique(lclM),1)
# tail(unique(uclF)+unique(uclM),1)

plot(years,unique(estF)+unique(estM),
     type='l',col=1,ylim=c(0,800),
     ylab="Population Size",
     xlab="Year"
)
lines(years,unique(lclF)+unique(lclM),lty=1,col=1)
lines(years,unique(uclF)+unique(uclM),lty=1,col=1)
legend(1998,600,legend=c("Lackey et al. 2013",
                         "Sedinger 2014",
                         "Sedinger 2018",
                         "Sultaire et al. 2023",
                         "This study"),pch=16,col=c(2:5,1))
points(2008,262,pch=16,col=2)
segments(2008,262,2008,262+31,col=2)
segments(2008,262,2008,262-31,col=2)

points(2014,445,pch=16,col=3)
segments(2014,445,2014,445+31,col=3)
segments(2014,445,2014,445-31,col=3)

points(2018,431,pch=16,col=4)
segments(2018,431,2018,431+33,col=4)
segments(2018,431,2018,431-33,col=4)

points(2020,418,pch=16,col=5)
segments(2020,418,2020,239,col=5)
segments(2020,418,2020,740,col=5)

points(tail(years,1),mean,pch=16,col=1)
segments(tail(years,1),mean,tail(years,1),lcl,col=1)
segments(tail(years,1),mean,tail(years,1),ucl,col=1)

lines(years,apply(ch.m,2,sum),type='l',lty=3)


################################################################################
### Plot 2
################################################################################

##
## phi
##

par(mfrow=c(1,1))

# ## Males
# plot(years,unique(estM),type='l', ylim=c(0,400),col=4)
# lines(years,unique(lclM),col=4,lty=2)
# lines(years,unique(uclM),col=4,lty=2)
# 
# ## Females
# lines(years,unique(estF),type='l',col=2)
# lines(years,unique(lclF),col=2,lty=2)
# lines(years,unique(uclF),col=2,lty=2)
# 
# ## 2022 estimate
# tail(unique(estF)+unique(estM),1)
# tail(unique(lclF)+unique(lclM),1)
# tail(unique(uclF)+unique(uclM),1)

plot(years,unique(estF.phi),
     type='l',col=2,ylim=c(0.2,1),
     ylab="Survival rate",
     xlab="Year",lwd=2
)
lines(years,unique(lclF.phi),lty=3,col=2)
lines(years,unique(uclF.phi),lty=3,col=2)

lines(years,unique(estM.phi),lty=1,col=4,lwd=2)
lines(years,unique(lclM.phi),lty=3,col=4)
lines(years,unique(uclM.phi),lty=3,col=4)

legend(2003,.4,legend=c("Female","Male"),col=c(2,4),lwd=2)


################################################################################
### Plot
################################################################################

##
## lambda
##

par(mfrow=c(1,1))

# ## Males
# plot(years,unique(estM),type='l', ylim=c(0,400),col=4)
# lines(years,unique(lclM),col=4,lty=2)
# lines(years,unique(uclM),col=4,lty=2)
# 
# ## Females
# lines(years,unique(estF),type='l',col=2)
# lines(years,unique(lclF),col=2,lty=2)
# lines(years,unique(uclF),col=2,lty=2)
# 
# ## 2022 estimate
# tail(unique(estF)+unique(estM),1)
# tail(unique(lclF)+unique(lclM),1)
# tail(unique(uclF)+unique(uclM),1)

plot(years,lambda,
     type='l',col=1,ylim=c(0,2),
     ylab=expression(lambda),
     xlab="Year",lwd=2
)
abline(h=1,lty=3)









