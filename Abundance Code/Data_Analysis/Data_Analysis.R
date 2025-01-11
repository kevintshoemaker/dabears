rm(list=ls())
set.seed(2024)

################################################################################
################################################################################
### Prepare Bear Data For Analysis
################################################################################
################################################################################

# setwd("~/Dropbox/GitHub/Nevada-Black-Bear-Abundance")

###
### Load required packages
###

required.packages=c("MuMIn",
                    "openCR",
                    "readxl",
                    "splines",
                    "xtable")

#install.packages(required.packages)
lapply(required.packages,library,character.only=TRUE)

###
### Required Functions
###

pasty<-function(x){
  k<-ncol(x)
  n<-nrow(x)
  out<-array(dim=n)
  out<-apply(x, 1 , paste , collapse = "" )
  return(out)
}

#####################################################################
### Set analysis options
#####################################################################

###
### years to be analyzed
###

years=1997:2022
T=length(years)
seasons=1:4 
J=length(seasons)

###
### Remove cubs that were never seen again
###

remove.cubs=TRUE

###
### Analyze male data only, female data only, or both males and females
###

sex="Both"

###
### Seasonal or annual primary sessions
###

primary.sessions="annual"

###
### Which data sets to use ("all","marked bears")
###

data.set="marked bears"

###
### Number of basis functions (4 minimum)
###

basis.functions=5

################################################################################
### Load original data sent by Carl Lackey
################################################################################

###
### Load Data set 1 ('marked bears' sheet)
###

# data1=read_xlsx(paste0('https://github.com/perrywilliamsunr/',
#                        'Nevada-Black-Bear-Abundance/raw/main/',
#                        'Data/CaptureHistories.xlsx'),
#                 sheet="marked bears")
# 

###
### This won't work once repository is made private
### You can just download "CaptureHistories.xlsx" and 
### source it from the download on your machine (locally)
###

url <- paste0("https://raw.githubusercontent.com/perrywilliamsunr/",
              "Nevada-Black-Bear-Abundance/main/Data/",
              "CaptureHistories.xlsx")

# Temporary file path to save the downloaded file
temp_file <- tempfile(fileext = ".xlsx")

# Download the file from GitHub
response <- GET(url, write_disk(temp_file, overwrite = TRUE))

# Check if the download was successful
if (response$status_code != 200) {
  stop("Failed to download file. Status code: ", response$status_code)
}

# Verify that the downloaded file is not empty and is a valid Excel file
if (file.size(temp_file) == 0) {
  stop("Downloaded file is empty.")
}

# Read the Excel file
data1 <- read_xlsx(temp_file, sheet = "marked bears")


##  Set the first row of required data.
##  First two rows in Excel sheet were used for naming.
start.row=3  

## Create an object that describes the range of years data were collected
all.years.of.data=1997:((ncol(data1)-16)/4+1996)
all.years.ind=rep(all.years.of.data,each=length(seasons))

###
### Remove last season in last year because there is no data during that time
###

all.years.ind=all.years.ind[-length(all.years.of.data)*J]


################################################################################
### Create capture history matrix for data set 1
################################################################################

## total number of bears captured (all years)
n1.total=length(start.row:nrow(data1)) 

## 17 is column capture histories start,-1 removes column with #capture events
ch1.all.tmp=data1[start.row:n1.total,17:(ncol(data1)-1)] 

## Convert to Matrix                                        
ch1.all.m.tmp=as.matrix(ch1.all.tmp,nrow(ch1.all.tmp),ncol(ch1.all.tmp)) 

# add last season of recapture data
#ch1.all.m=cbind(ch1.all.m.tmp,0)                
ch1.all.m=cbind(ch1.all.m.tmp)                

################################################################################
### Include sex information
################################################################################

sex1.all=data1[start.row:n1.total,8] # gender is in the 8th column
sex1.all.f=factor(ifelse(sex1.all=="M","Male","Female"))

################################################################################
### Include age information
################################################################################

age1.all.tmp=data.frame(data1[start.row:n1.total,7]) # age is in the 7th column

## Fix data entry errors
age1.all.tmp[age1.all.tmp=="7 WK"]=7/52
age1.all.tmp[age1.all.tmp=="8 WEEKS"]=8/52
age1.all.tmp[age1.all.tmp=="4 MO"]=4/12
age1.all.tmp[age1.all.tmp=="4 MOS"]=4/12
age1.all.tmp[age1.all.tmp=="5 MO"]=5/12
age1.all.tmp[age1.all.tmp=="6 MO"]=6/12
age1.all.tmp[age1.all.tmp=="7 MO"]=7/12
age1.all.tmp[age1.all.tmp=="7 MOS"]=7/12
age1.all.tmp[age1.all.tmp=="8 MO"]=8/12
age1.all.tmp[age1.all.tmp=="8 MOS"]=8/12
age1.all.tmp[age1.all.tmp=="9 MO"]=9/12
age1.all.tmp[age1.all.tmp=="9 MOS"]=9/12
age1.all.tmp[age1.all.tmp=="10 MO"]=10/12
age1.all.tmp[age1.all.tmp=="10 MOS"]=10/12
age1.all.tmp[age1.all.tmp=="11 MO"]=11/12
age1.all.tmp[age1.all.tmp=="11 MOS"]=11/12
age1.all.tmp[age1.all.tmp=="15 MOS"]=15/12
names(age1.all.tmp)="age"
age1.all=as.numeric(age1.all.tmp$age)

###
### Indicator vectors to subset years of analysis
###

col.ind=all.years.ind%in%years
row.ind1=apply(ch1.all.m[,col.ind],1,sum)>0

###
### Remove bears < 16 months that are never captured again
###

if(remove.cubs){
  keep=(apply(ch1.all.m,1,sum)>1 | age1.all>(16/12))
}
if(!remove.cubs)keep=rep(TRUE,nrow(ch1.all.m))

ch1.m=ch1.all.m[row.ind1&keep,col.ind]
n=nrow(ch1.m)

###
### subset covariate data
###

sex1=sex1.all.f[row.ind1&keep]
age1=age1.all[row.ind1&keep]
length(age1)

###
### Session-specific covariates
###

season=factor(rep(1:4,T))
season=season[-(T*J)]
year=rep(1:T,each=J)
year=year[-(T*J)]
year.f=factor(year)
year2=year^2
bs=bs(year,basis.functions,intercept=TRUE)

sescov=data.frame(season=season,
                  year=year,
                  year.f=year.f,
                  year2=year2)

for(i in 1:basis.functions) {
  new=bs[,i]                     
  sescov[,ncol(sescov)+1]=new  
  colnames(sescov)[ncol(sescov)]=paste0("bs",i) 
}

if(data.set=="marked bears"){
  
  ###
  ### Format data for package analysis
  ###
  
  ch.m=ch1.m
  sex=sex1
  age=age1
  
  bear.df=data.frame(ch=pasty(ch1.m),sex=sex1)
  bearsCH=suppressWarnings(unRMarkInput(bear.df))
}

# nrow(bear1CH) #197 in manuscript for 2008 (192 here)
# sum(sex1=="Male") #123 
# sum(sex1=="Female") #74 
# sum(bear1CH) # 546 in manuscript  for 2008 

################################################################################
### Data set 2: dead on first event bears
################################################################################

if(data.set=="all"){
  data2 <- read_xlsx(temp_file, 
                     sheet = "dead on 1st event bears 97-2018")
  
  ###
  ### Data set 2: create capture history matrix
  ###
  
  n2.total=length(start.row:nrow(data2)) #total number of bears captured (all years)
  ch2.all.tmp=data2[start.row:n2.total,17:(ncol(data2)-1)] ## 17 is column capture histories start
  ## -1 removes column with #capture events
  ch2.all.m.tmp=as.matrix(ch2.all.tmp,nrow(ch2.all.tmp),ncol(ch2.all.tmp)) # convert to matrix
  ch2.all.m=matrix(as.numeric(ch2.all.m.tmp), nrow(ch2.all.m.tmp),ncol(ch2.all.m.tmp))
  
  ###
  ### Include sex information
  ###
  
  sex2.all=data2[start.row:n2.total,8] # gender is in the 8th column
  sex2.all.f=factor(ifelse(sex2.all=="M","Male","Female"))
  
  ###
  ### Include age information
  ###
  
  age2.all.tmp=data.frame(data2[start.row:n2.total,7]) # age is in the 7th column
  
  ## Fix data entry errors
  age2.all.tmp[age2.all.tmp=="UNK"]=NA
  age2.all.tmp[age2.all.tmp=="ADULT"]=NA
  age2.all.tmp[age2.all.tmp=="UK"]=NA
  names(age2.all.tmp)="age"
  age2.all=as.numeric(age2.all.tmp$age)
  
  ###
  ### Remove bears < 16 months that are never captured again
  ###
  
  row.ind2=apply(ch2.all.m[,col.ind],1,sum)>0
  if(remove.cubs){
    keep=(apply(ch2.all.m,1,sum)>1 | age2.all>(16/12))
  }
  if(!remove.cubs)keep=rep(TRUE,nrow(ch2.all.m))
  keep[is.na(keep)]=FALSE
  
  ch2.m=(ch2.all.m[row.ind2&keep,col.ind])
  n2=nrow(ch2.m)
  
  ###
  ### subset covariate data
  ###
  
  sex2=sex2.all.f[row.ind2&keep]
  age2=age2.all[row.ind2&keep]
  
  ch2.m=ch2.m
  ch.m=rbind(ch1.m,ch2.m)
  sex=c(sex1,sex2)
  
  bears.df=data.frame(ch=pasty(ch.m),sex=sex)
  bearsCH=unRMarkInput(bears.df)
}

################################################################################
### Dissolve to annual data
################################################################################

## Don't run this twice, it won't work the second time due to changing ch.m
if(primary.sessions=="annual"){
  
  ## Indicators to show the months that start and end of each year
  ind1=seq(1,ncol(ch.m),J) 
  ind2=seq(J,ncol(ch.m)+1,J)
  ind2[T]=ind2[T]-1
  
  ch.annual=matrix(NA,nrow(ch.m),T)
  for(t in 1:T){
    ch.annual[,t]=apply(ch.m[,ind1[t]:ind2[t]],1,sum)
  }
  ch.m=ifelse(ch.annual>0,1,0)
  bears.df=data.frame(ch=pasty(ch.m),sex=sex)
  bearsCH=suppressWarnings(unRMarkInput(bears.df))
  
  year=1:T
  year.f=factor(year)
  year2=year^2
  bs=bs(year,basis.functions,intercept=TRUE)
  sescov=data.frame(year=year,
                    year.f=year.f,
                    year2=year2)
  for(i in 1:basis.functions){
    new=bs[,i]                     
    sescov[,ncol(sescov)+1]=new  
    colnames(sescov)[ncol(sescov)]=paste0("bs",i) 
  }
}


################################################################################
### Model selection
################################################################################

if(basis.functions==6){
  phi=c(1,"sex","bs1 + bs2 + bs3 + bs4 + bs5 + bs6",
        "year","sex + bs1 + bs2 + bs3 + bs4 + bs5 + bs6","sex+year")
  p=c(1,"sex","bs1 + bs2 + bs3 + bs4 + bs5 + bs6",
      "year","sex + bs1 + bs2 + bs3 + bs4 + bs5 + bs6","sex+year")
  N=c(1,"sex","bs1 + bs2 + bs3 + bs4 + bs5 + bs6","year",
      "sex + bs1 + bs2 + bs3 + bs4 + bs5 + bs6","sex+year")
  models=expand.grid(phi,p,N)
}

if(basis.functions==5){
  phi=c(1,"sex","bs1 + bs2 + bs3 + bs4 + bs5",
        "year","sex + bs1 + bs2 + bs3 + bs4 + bs5","sex+year")
  p=c(1,"sex","bs1 + bs2 + bs3 + bs4 + bs5",
      "year","sex + bs1 + bs2 + bs3 + bs4 + bs5","sex+year")
  N=c(1,"sex","bs1 + bs2 + bs3 + bs4 + bs5","year",
      "sex + bs1 + bs2 + bs3 + bs4 + bs5","sex+year")
  models=expand.grid(phi,p,N)
}

if(basis.functions==4){
  phi=c(1,"sex","bs1 + bs2 + bs3 + bs4 ",
        "year","sex + bs1 + bs2 + bs3 + bs4 ","sex+year")
  p=c(1,"sex","bs1 + bs2 + bs3 + bs4 ",
      "year","sex + bs1 + bs2 + bs3 + bs4","sex+year")
  N=c(1,"sex","bs1 + bs2 + bs3 + bs4","year",
      "sex + bs1 + bs2 + bs3 + bs4","sex+year")
  models=expand.grid(phi,p,N)
}

args=list()
for(i in 1:nrow(models)){
  args[[i]]=list(bearsCH,type='JSSAN',
                 model=list(formula(noquote(paste0("phi~",models[i,1]))),
                            formula(noquote(paste0("p~",models[i,2]))),
                            formula(noquote(paste0("N~",models[i,3])))),
                 sessioncov=sescov,
                 method="Nelder-Mead",details=list(control=list(maxit=5000))
  )
}




###
### Fit model
###

suppressWarnings(fits <- par.openCR.fit(args, ncores = 16))
save(fits,file=paste0(getwd(),"/Results/fits_",basis.functions,"bs",head(years,1),"_",tail(years,1),".RData"))
load(paste0(getwd(),"/Results/fits_",basis.functions,"bs",head(years,1),"_",tail(years,1),".RData"))
remove=which(unname(sapply(fits, function(x) x$fit$convergence))==1)
if(length(remove)>0) fits2=fits[-remove]
if(length(remove)==0) fits2=fits
aic.table=AIC(fits2, criterion = 'AIC', use.rank = TRUE)

xtable((aic.table[1:10,]))
(top.model=aic.table$model[1])
(top.model.no=as.numeric(rownames(aic.table[1,])))
# top.model.no=77 # 5 bs
# top.model.no=150 # 4 bs
top.model=list(formula(noquote(paste0("phi~",models[top.model.no,1]))),
               formula(noquote(paste0("p~",models[top.model.no,2]))),
               formula(noquote(paste0("N~",models[top.model.no,3]))))


## replicate top model
m=fits[[top.model.no]]
# m=openCR.fit(bearsCH,type='JSSAN',
#              model=top.model,
#              sessioncov=sescov,
#              method='Nelder-Mead',
#              details=list(control=list(maxit=5000)
#              )
# )

pred=predict(m,all=TRUE)
## derived(m)

est=pred$N$estimate
estM=pred$N$estimate[pred$N$sex=="Male"]
lclM=pred$N$lcl[pred$N$sex=="Male"]
uclM=pred$N$ucl[pred$N$sex=="Male"]
estF=pred$N$estimate[pred$N$sex=="Female"]
lclF=pred$N$lcl[pred$N$sex=="Female"]
uclF=pred$N$ucl[pred$N$sex=="Female"]

mean=tail(unique(estF)+unique(estM),1)
lcl=tail(unique(lclF)+unique(lclM),1)
ucl=tail(unique(uclF)+unique(uclM),1)


###
### Phi
###

est.phi=pred$phi$estimate
estM.phi=pred$phi$estimate[pred$phi$sex=="Male"]
lclM.phi=pred$phi$lcl[pred$phi$sex=="Male"]
uclM.phi=pred$phi$ucl[pred$phi$sex=="Male"]
estF.phi=pred$phi$estimate[pred$phi$sex=="Female"]
lclF.phi=pred$phi$lcl[pred$phi$sex=="Female"]
uclF.phi=pred$phi$ucl[pred$phi$sex=="Female"]

mean.phi=tail(unique(estF.phi)+unique(estM.phi),1)
lcl.phi=tail(unique(lclF.phi)+unique(lclM.phi),1)
ucl.phi=tail(unique(uclF.phi)+unique(uclM.phi),1)

###
### lambda
###

lambda=derived(m)$estimates$lambda
save.image(file = paste0(getwd(),"/Results/Data_Analysis_Workspace.RData"))




