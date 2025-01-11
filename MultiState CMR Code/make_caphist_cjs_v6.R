
# read in and process bear data
# last modified: 2/11/24


# Clear workspace ----------------------------

rm(list=ls())


# Load packages ---------------------------

library(lubridate)
library(readxl)


# Load data -------------------------------

bears<-read_excel("Bear data_Master2022.xlsx", sheet = 1,skip=1)  # read in master bear data
temp <- read_excel("Bear data_Master2022.xlsx", sheet = 1,skip=0)   #add col names
names(bears) <- names(temp)
# bears[1:5,1:5]

# Process data ---------------------------

bears$DATE <- ymd(bears$DATE)
bears$YEAR <- year(bears$DATE)
bears$MONTH <- month(bears$DATE)
bears$YDAY <- yday(bears$DATE)
bears$SEX <- toupper(bears$SEX)
bears$SEX[bears$SEX%in%c("?","UNK")] <- NA
sort(unique(bears$SEX))

bears$EST. <- toupper(bears$EST.)
bears$AGE <- bears$EST.
bears$AGE[grepl("10\\+",bears$EST.)] <- '10'
bears$AGE[grepl("15\\+",bears$EST.)] <- '15'
bears$AGE[grepl("MOS",bears$EST.)] <- '0'
bears$AGE[grepl("MO",bears$EST.)] <- '0'
bears$AGE[grepl("M",bears$EST.)] <- '0'
bears$AGE[grepl("0.5",bears$EST.)] <- '0'
bears$AGE[grepl("WKS",bears$EST.)] <- '0'
bears$AGE[grepl("WEEKS",bears$EST.)] <- '0'
bears$AGE[grepl("UNK",bears$EST.)] <- NA
bears$AGE[grepl("ADULT",bears$EST.)] <- '10'

bears$AGE[bears$EST.=="?"] <- NA

sort(unique(bears$AGE))
sort(unique(bears$EST.))

temp <- subset(bears,CORRECTED=="101")
temp$AGE

bears$AGE <- as.numeric(bears$AGE)
hist(bears$AGE)

temp <- subset(bears,!is.na(`BEAR #`))
length(which(temp$AGE==0))/length(which(!is.na(temp$AGE)))  

sort(unique(bears$AGE))

bears$BCS[bears$BCS=="?"] <- NA
bears$BCS[bears$BCS=="NA"] <- NA
bears$BCS[bears$BCS=="UNK"] <- NA
bears$BCS[bears$BCS=="G"] <- '3'
sort(unique(bears$BCS))

bears$BCS <- as.numeric(bears$BCS)
hist(bears$BCS)

bears$`MORTALITY TYPE`[bears$`MORTALITY TYPE`=="UNK"] <- NA
bears$`MORTALITY TYPE`[bears$`MORTALITY TYPE`=="ILLEGAL"] <- "ILL"
allmorts <- sort(unique(bears$`MORTALITY TYPE`))
bears[which(bears$`BEAR #`=="P10"),] <- NA


hist(bears$AGE)


# Make initial list of bear IDs  

allids <- sort(unique(bears$CORRECTED))
allids <- setdiff(allids ,c("NA","NOTCH",'UNK',"BB006","TBD"))
# View(bears[which(bears$CORRECTED=="RED8"),])
allids <- allids[!grepl("CDFW",allids)]

allids <- allids[!is.na(allids)]

ninds <- length(allids)    # 699 bears!
ninds

# bears2 <- subset(bears,CORRECTED%in%allids)
# nrow(bears2)   # 1465 observations

# remove cubs from master data -----------------------

bears <- subset(bears,!AGE==0)   # 1650 bears

# Make final list of bear IDs  

allids <- sort(unique(bears$CORRECTED))
allids <- setdiff(allids ,c("NA","NOTCH",'UNK',"BB006","TBD"))
# View(bears[which(bears$CORRECTED=="RED8"),])
allids <- allids[!grepl("CDFW",allids)]

allids <- allids[!is.na(allids)]

ninds <- length(allids)    # 594 bears!
ninds

bears2 <- subset(bears,CORRECTED%in%allids)
nrow(bears2)   # 1220 observations

# Get ready to make caphist -------------------------------


ninds <- length(allids)    # 594 bears
ninds

### define time frame for sessions and occasions
### make caphist array with dimension ninds, nyears 

realyears <- 1998:2022
nyears <- length(realyears)

# Make capture history arrays (and ancillary data structures) -------------------------------

caphist_ms <- array(1,dim=c(ninds,nyears))         # dim: nind, nyears
caphist_liverel <- array(0,dim=c(ninds,nyears))    # live releases only (for H-T abundance estimation)
status <- array(NA,dim=c(ninds,nyears))

sexarray <- numeric(ninds)
sexarray[] <- NA
agearray <- array(NA,dim=c(ninds,nyears))     # cubs, yearlings, juveniles adult?  Or represent as real age in years??
bcsarray <- array(1,dim=c(ninds,nyears))         # dim: nind, nyears

count_as_cap <- c("MGMT","PDEP","ACC","HBC")  # only these capture types count as a capture for the model
dontcount_as_cap <- setdiff(allmorts,count_as_cap)

torm <- c()  # list of individuals to remove
firsts <- numeric(ninds)
lasts <- numeric(ninds)

# STATES   1: alive, 2: management removal, 3: dead on road, 4: dead
# OBSERVATIONS:   1:not seen,  2: seen and released alive, 3: management removal, 4: seen dead on road

i=1;y=1
for(i in 1:ninds){
  temp <- subset(bears,CORRECTED==allids[i])
  tempsex <- temp$SEX[temp$SEX%in%c("M","F")]
  if(length(tempsex)>0) sexarray[i] <- tempsex[1]
  firstyr <- min(temp$YEAR,na.rm=T)
  firsts[i] <- match(firstyr,realyears)
  lastyr <- max(temp$YEAR,na.rm=T)
  lastndx <- match(lastyr,realyears)
  firstage <- min(temp$AGE)
  agearray[i,firsts[i]:nyears] <- firstage+0:(nyears-firsts[i])
  for(y in firsts[i]:lastndx){
    thisyear <- realyears[y]
    temp2 <- subset(temp,YEAR==thisyear)
    if(nrow(temp2)>0){  # if any observations this year
      dead <- ifelse(any(temp2$`MORTALITY TYPE`%in%allmorts),1,0)  # is it dead?
      if(dead==1){  # if dead
        if(thisyear==firstyr){       # if it's the first year this bear is in the study, remove it
          torm = c(torm,allids[i])   # remove if died in first year
        }else{  # if after the first year...
          mort_type <- names(which.max(table(na.omit(temp2$`MORTALITY TYPE`))))
          if(mort_type=="MGMT"){
            status[i,firsts[i]:max(firsts[i],y-1)] <- 1
            status[i,y] <- 2
            caphist_ms[i,y] <- 3
            lasts[i] <- min(nyears,y+1)
          }else if(mort_type=="HBC"){
            status[i,firsts[i]:max(firsts[i],y-1)] <- 1
            status[i,y] <- 3
            caphist_ms[i,y] <- 4
            lasts[i] <- min(nyears,y+1)
          }else{    # if some other mortality type
            status[i,firsts[i]:max(firsts[i],y-1)] <- 1
            status[i,y] <- NA
            caphist_ms[i,y:nyears] <- NA
            lasts[i] <- y-1   # remove from the likelihood if died from some other cause
          } # end going through causes of mortality
        } # end first-year check
      }else{  # if alive
        status[i,firsts[i]:y] <- 1
        caphist_ms[i,y] <- 2
        caphist_liverel[i,y] <- 1
      }
      
      bcsarray[i,y] <- mean(temp2$BCS,na.rm=T)
    } # end if observed this year
  }
}

lasts[lasts==0] <- nyears

length(torm)   # 85 individuals removed at this stage...
tokeep <- setdiff(allids,torm)
tokeep_ndx <- match(tokeep,allids)
agearray <- agearray[tokeep_ndx,]
allids <- tokeep
ninds <- length(allids)
bcsarray <- bcsarray[tokeep_ndx,] 
caphist_liverel <- caphist_liverel[tokeep_ndx,]
caphist_ms <- caphist_ms[tokeep_ndx,]
status <- status[tokeep_ndx,]
sexarray <- sexarray[tokeep_ndx]
firsts <- firsts[tokeep_ndx]
lasts <- lasts[tokeep_ndx]

rownames(caphist_ms) = allids
rownames(caphist_liverel) = allids
colnames(caphist_ms) = realyears
colnames(caphist_liverel) = realyears


names(sexarray) = allids
rownames(agearray) = allids
colnames(agearray) = realyears

agearray2 <- agearray
agearray2[agearray>=5] <- 5

# Make additional data structures for JAGS --------------------------------

#### this makes nyears and nocc 
allyears <- realyears
nyears <- length(allyears)
ninds <- length(allids)

is.male =ifelse(sexarray=="M",1,0)

# agearray2

juv.cutoff <- 4
is.juv <- ifelse(agearray2<=juv.cutoff,1,0)

 # read in snow data
# rm(precip)
dat<-read.csv("sierra_data.csv" )

ndx <- match(allyears,dat$Year)
SWE<-dat$average.1[ndx]
names(SWE)<-allyears

# snowtel <- precip
# par(mfrow=c(2,1))
# plot(x= c(allyears), y = snowtel[,as.character(allyears)],type="l")
# plot(x= c(allyears), y = SWE[,as.character(allyears)],type="l")
# cor(snowtel[,as.character(allyears)],SWE[,as.character(allyears)])

## choose to use SWE

## First freeze data -------------

firstfreeze <- read.table("RenoFreezeDates.txt", skip = 3, sep = ",", header = T)
firstfreeze$lastdate <- mdy(firstfreeze$Last)

head(firstfreeze)
firstfreeze$jday <- yday(firstfreeze$lastdate)
firstfreeze$year <- year(firstfreeze$lastdate)

firstfreeze1 <- firstfreeze[match(allyears,firstfreeze$year),]


# SAVE ALL KEY DATA --------------
save(
  ninds,
  nyears,
  is.male,
  is.juv,
  SWE,
  caphist_liverel,
  caphist_ms,
  status,
  firsts,
  lasts,
  agearray2,
  allids,
  allyears,
  firstfreeze1,
  file="beardata_bundle_v6.RData")


### KTS Ended here -------------


rm(list = ls())
load("beardata_bundle_cjs.RData")
mgmt2d <- array(0, dim=c(ninds,nyears))
hbc2d <- array(0, dim=c(ninds,nyears))

i=1;y=1
for(i in 1:ninds){
  for(y in 1:nyears){
    test1 <- any(mgmtarray[i,y,]==1)
    if(!is.na(test1)) if(test1) mgmt2d[i,y] <- 1  
    
    test1 <- any(hbcarray[i,y,]==1)
    if(!is.na(test1)) if(test1) hbc2d[i,y] <- 1  
  }
}

rownames(mgmt2d) <- rownames(caphist2d)
rownames(hbc2d) <- rownames(caphist2d)

### what are we going to do with the bears that we dont have in the model 

# does the total number of bears removed corr with the snotel or maybe last freeze
## Do we have a signal 
## how can we plot the total number of HBC, MNGT and look at a simple regression
## Josh can try to do that 
## the above can not be incorp to the current model but can still have something that is overall useful 


################### 
## here I will want to try to load in the new Data from Heather that has the total number of deaths by year 
beardeaths <- read_excel("Unknown Mortalities.xlsx")
beardeaths <- beardeaths[-1,]

dat<-read.csv("sierra_data.csv" )


# head(dat)
# precip<-matrix(data=dat$Prec_avg,nrow=1,nrow(dat))
# colnames(precip)<-as.character(dat$Year[1:26])

# rownames(precip)<-"AvgPrecip"

SWE<-matrix(data = dat$average.1[1:26],nrow=1,26)
colnames(SWE)<-as.character(dat$Year[1:26])
rownames(SWE)<-"AveSWE"
## choose to use SWE
# SAVE ALL KEY DATA
firstfreeze <- read.table("RenoFreezeDates.txt", skip = 3, sep = ",", header = T)

firstfreeze$lastdate <- mdy(firstfreeze$Last)

head(firstfreeze)
firstfreeze$jday <- yday(firstfreeze$lastdate)
firstfreeze$year <- year(firstfreeze$lastdate)

yearswewant <- 1998:2022
firstfreeze1 <- firstfreeze[which(firstfreeze$year%in%yearswewant),]
### this is a plot to visualize the last freeze and looking at the SWE
plot(firstfreeze1$jday, SWE[,-1])
cor(firstfreeze1$jday, SWE[,-1])

#### now I want to look at the total number of deaths and see if we have any type of relation.... 
sumofHBC <- c()
# sumofallmort <- c()
sumofMGMT <- c()
# sumofbears <- list()
year <- 2022
# beardeaths[which(beardeaths$Year==year),]
nrow(beardeaths)   ## 343 inds that were not in final model 

#### Here I need to find when and where HBC and MGMT add into the death count 
year <- 1998
counter=1
for(year in 1998:2022){
  sumofHBC[[paste0("y",year)]] <- sum(length(which(beardeaths$Year==year&beardeaths$`MORTALITY TYPE`=="HBC"))+
                                      apply(hbc2d,2,sum,na.rm=T)[counter]  )
  # if(year<2022){
  #   sumofallmort[[paste0("y",year)]] <- sum(length(which(beardeaths$Year==year))+length(which(lasts==(year-1997))))  
  # }else{
  #   sumofallmort[[paste0("y",year)]] <- sum(length(which(beardeaths$Year==year))+length(which(lasts2[which(lasts==(year-1997))]<40))) ## how to fix because some bears die in last occ 
  # }   
  
  sumofMGMT[[paste0("y",year)]] <- sum(length(which(beardeaths$Year==year&beardeaths$`MORTALITY TYPE`=="MGMT"))+
                                         apply(mgmt2d,2,sum,na.rm=T)[counter] )
  counter=counter+1
}

plot(SWE[,-1], unlist(sumofMGMT))
# cor(as.vector(SWE[,-1]), as.numeric(sumofMGMT[1:nrow(firstfreeze1)]))
plot(SWE[,-1], unlist(sumofHBC))
?unlist
mod1 <- lm(unlist(sumofHBC)~SWE[,-1])
mod1log <- lm(log(unlist(sumofHBC))~SWE[,-1])
summary(mod1log)
summary(mod1)

## here is the lm for mgmt
# sumofMGMT1 <- 
mod2 <- lm(unlist(sumofMGMT)~SWE[,-1])
mod2log <- lm(log(unlist(sumofMGMT)+1)~SWE[,-1])    ### here we have a 0 in 1999 so we added a 1 to all values 
summary(mod2log)
summary(mod2)
#### now we are doing the sum of all remove to SWE

cor(as.vector(SWE[,-1]), as.numeric(sumofHBC[1:nrow(firstfreeze1)]))
plot(SWE[,-1], unlist(sumofMGMT) + unlist(sumofHBC) )

cor(SWE[,-1], unlist(sumofMGMT) + unlist(sumofHBC) )

mod3 = lm( unlist(sumofMGMT) + unlist(sumofHBC) ~ SWE[,-1])
summary(mod3)
mod3log = lm( log(unlist(sumofMGMT) + unlist(sumofHBC)) ~ SWE[,-1])
summary(mod3log)

## now we can look at the last freeze date due to morts 
##### here we will look at the plots first 
plot(firstfreeze1$jday, unlist(sumofMGMT))
cor(firstfreeze1$jday, unlist(sumofMGMT))
# cor(as.vector(firstfreeze1$jday), as.numeric(sumofMGMT[1:nrow(firstfreeze1)]))
plot(firstfreeze1$jday, unlist(sumofHBC))
cor(as.vector(firstfreeze1$jday), as.numeric(sumofHBC[1:nrow(firstfreeze1)]))
plot(firstfreeze1$jday, unlist(sumofMGMT) + unlist(sumofHBC))
cor(as.vector(firstfreeze1$jday), as.numeric(sumofallmort[1:nrow(firstfreeze1)]))
### bellow will be the LM 

### First it is HBC
mod4 = lm( unlist(sumofHBC) ~ firstfreeze1$jday)
summary(mod4)
mod4log = lm( log( unlist(sumofHBC)) ~ firstfreeze1$jday)
summary(mod4log)         ### this is starting to look interesting 
### Now MGMT
mod5 = lm( unlist(sumofMGMT) ~ firstfreeze1$jday)
summary(mod5)
mod5log = lm( log(unlist(sumofMGMT)+1) ~ firstfreeze1$jday)   ### here is where we add 1 value to all removes 
summary(mod5log)
## this is for MGMT ++ HBC  
mod6 = lm( unlist(sumofMGMT) + unlist(sumofHBC) ~ firstfreeze1$jday)
summary(mod6)
mod6log = lm( log(unlist(sumofMGMT) + unlist(sumofHBC)) ~ firstfreeze1$jday)
summary(mod6log)

# cor(as.vector(SWE[,-1]), as.numeric(sumofallmort[1:nrow(firstfreeze1)]))






### why is this happening maybe Josh messed it up?? 


# Save caphist and corresponding covariates -------------------------

filenm <- sprintf("caphist_etc_cjs_%s.RData",Sys.Date())
save(caphist3d,caphist2d,sexarray,agearray,bcsarray,huntarray,deadarray,hbcarray,mgmtarray,allids,deadarray1,file=filenm)
# load("C:/Users/jvasquez/Box/2022 Bear Project/caphist_etc_cjs_2023-01-29.RData")

nrow(caphist3d[,,1])    ## this is 699 now!


