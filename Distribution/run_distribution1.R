
# This code provides several anomaly checks and outlier detection for eDNA primers. Of particular intrest are eONKI4 and eFISH1.

# Set parameters and load in libraries.
knitr::opts_chunk$set(fig.width=12, fig.height=10) 

suppressMessages(library(tidyverse))  #For general usage and piping

suppressMessages(library(fitdistrplus)) #This package allows us to fit models when there is no clear closed form

suppressMessages(library(EnvStats)) #This package provides qqPlot a slight alternative to qqplot

suppressMessages(library(actuar)) #We use this package to model Log-Logistic distributions

library(goftest) # This package contains the ad.test and cvm.test functions

suppressMessages(library(kableExtra)) #For building kable plots

suppressMessages(library(OutlierDetection)) #For outlier detection


# Read in Data and Wrangle.

Dat=read.csv('Helbing-LOD_LOQ_DataZero48-040319MJAML.csv',header=TRUE,stringsAsFactors = TRUE) #Read in the data set

Dat$SQf=as.factor(Dat$SQ) #Add a column that is the factor version of SQ

Dat=Dat %>%
  filter(SQ>=0.8) #We only wish to consider samples with an SQ of 0.8 or larger.

invisible(droplevels(Dat$SQf)) #Remove factor levels not used

Dat$SQf=factor(Dat$SQf) #Add a column(the SQ as a factor)




Dat2=Dat %>%
  filter(Target!='eFISH1') %>%
  filter(Target!='eONKI4') #We only wish to consider samples with an SQ of 0.8 or larger.

invisible(droplevels(Dat2$SQf)) #Remove factor levels not used


Dat2$SQf=as.factor(Dat2$SQ) #Add a column that is the factor version of SQ

invisible(droplevels(Dat2$SQf))
invisible(droplevels(Dat2$Target))


Dat2$SQf=factor(Dat2$SQf) #Add a column(the SQ as a factor)
Dat2$Target=factor(Dat2$Target)


# eFISH1

Dat.eFISH= Dat %>%
  filter(Target=='eFISH1') #Restrict attention to eFISH1


coef_var=rep(0,length(levels(Dat$SQf))) #Initalize a vector to store the regular cvs

coef_var_unbiased=rep(0,length(levels(Dat$SQf))) #Initalize a vector to store the unbiased cvs

p_val.cvm=rep(0,length(levels(Dat$SQf))) #Initalize a vector to store the p values from the cvm test
p_val.ad=rep(0,length(levels(Dat$SQf))) #Initalize a vector to store the p values from the ad test

entry=1 # this is used to loop over and populate the vectors


for(i in levels(Dat$SQf)){ # Loop over all SQ of intrest
  
  sub.quant= Dat.eFISH %>%
    filter(SQf==i)

  sub.quant=sub.quant[!(is.na(sub.quant$Cq)),]
  
  coef_var[entry]=round(sd(sub.quant$Cq)/mean(sub.quant$Cq),4) #Calculate the coefficient of variation
  
  coef_var_unbiased[entry]=round((1+1/(4*length(levels(Dat$SQf))))*coef_var[entry],4) #Unbiased CV for normal data
  
  
  #cvm.test and ad.test allow us to test our hypothesis of the fitted model of choice, we use the MLE for parameters
  # The null hypothesis is that the distribution we included is legitimate, if the p-value is small we would reject this
  #null hypothesis, so we go through all and see if it ever is rejected, if so we print stop. If stop is never printed that means we cannot reject the null that the distribution fits.
  
  cvm=cvm.test(sub.quant$Cq, "pnorm", mean=mean(sub.quant$Cq), sd=sd(sub.quant$Cq)) #the cvm test with supplied model being normal
  
  ad=ad.test(sub.quant$Cq, "pnorm", mean=mean(sub.quant$Cq), sd=sd(sub.quant$Cq)) #the ad test with supplied model being normal
  
  
  p_val.ad[entry]=round(ad$p.value,4) 
  p_val.cvm[entry]=round(cvm$p.value,4)
  
  entry=entry+1
  
 # if(ad$p.value<0.05 | cvm$p.value<0.05){ #If the p-value is below 0.05 for either test, we should reject the null hypothesis
 #   print("stop")
 # }
  cat("\n\n")
}




coef_var=rep(0,length(levels(Dat$SQf))) #Using MLE to calculate

coef_var_lognormal=rep(0,length(levels(Dat$SQf))) #Alternate form this is CVa

coef_test=rep(0)

coef_new_test=rep(0)

p_val.cvm=rep(0,length(levels(Dat$SQf)))
p_val.ad=rep(0,length(levels(Dat$SQf)))

entry=1


for(i in levels(Dat$SQf)){
  
  sub.quant= Dat.eFISH %>%
    filter(SQf==i)
  
  sub.quant=sub.quant[!(is.na(sub.quant$Cq)),] #Remove missing entries
  
 # qqnorm(log(sub.quant$Cq),main=paste("SQ",i))
 # qqline(log(sub.quant$Cq))
  
  
  n=length(sub.quant$Cq)
  
  mle.mu=sum(log(sub.quant$Cq))/(length(sub.quant$Cq)) #Calculate the MLE for the mu parameter of Log-Normal
  
  mle.sigmasq=sum((log(sub.quant$Cq)-mle.mu)^2)/(length(sub.quant$Cq)) #Calculate the MLE for the sigma squared parameter of the Log-Normal
  
  mle.sigma=sqrt(mle.sigmasq)
  

  
  mle.mean=exp(mle.mu+(mle.sigmasq/2)) #Estimate mean using MLE parameters
  
  
  mle.var=(exp(mle.sigmasq)-1)*exp(2*mle.mu+mle.sigmasq) #Estimate variance using MLE parameters
  
  mle.sd=sqrt(mle.var)
  
  
  log.mean=mean(log(sub.quant$Cq)) #Calculate the sample mean of the log cq value (ie y bar)
  
  log.sd=sd(log(sub.quant$Cq)) #this is s, the standard deviation of log data
  
  
  coef_test[entry]=round(sqrt(exp((mle.sigmasq))-1),4)                        
  
  sln=var(log(sub.quant$Cq)) #Calculate sln used for our alternatve cv calculation
  
  coef_var[entry]=round(log.sd/log.mean,4)  #Calculate CV using estimate sd and mean , 
  
  coef_var[entry]=round(sd(sub.quant$Cq)/mean(sub.quant$Cq),4)
  
  coef_var_lognormal[entry]=round(sqrt(exp((mle.sigmasq))-1),4)
  
  cvm=cvm.test(sub.quant$Cq, "plnorm", meanlog=log.mean, sdlog=log.sd)
  
  ad=ad.test(sub.quant$Cq, "plnorm", meanlog=log.mean, sdlog=log.sd)
  
  p_val.ad[entry]=round(ad$p.value,4)
  p_val.cvm[entry]=round(cvm$p.value,4)
  
  entry=entry+1
  
  
 # if(ad$p.value<0.05 | cvm$p.value<0.05){ #If the p-value is below 0.05 we should reject the null hypothesis.
  #  print(paste("Stop: Failure at SQ:",i))
#  }
  
  cat("\n\n")
  
}
  
  
  
  
  coef_var=rep(0,length(levels(Dat$SQf))) 
  coef_var.standard=rep(0,length(levels(Dat$Sqf)))
  p_val.cvm=rep(0,length(levels(Dat$SQf)))
  p_val.ad=rep(0,length(levels(Dat$SQf)))
  
  entry=1
  #pdf(file='eFISH1WeibullQQ.pdf')
  par(mfrow=c(3,2))
  
  for(i in levels(Dat$SQf)){
    
    sub.quant= Dat.eFISH %>%
      filter(SQf==i) %>%
      filter(!is.na(Cq))
    
    sub.quant=sub.quant[!(is.na(sub.quant$Cq)),]
    
    dat.weib=fitdist(sub.quant$Cq,"weibull") #Fit a Weibull distribution using fitdist to estimate the parameters
    
#    qqPlot(qweibull(ppoints(length(sub.quant$Cq)), shape = dat.weib$estimate[1], 
   #                 scale =  dat.weib$estimate[2]), sub.quant$Cq,add.line=TRUE,xlab="Theoretical Quantiles",ylab="Sample Quantiles",main=paste("SQ",i))
    
    
    cvm=cvm.test(sub.quant$Cq, "pweibull", shape=dat.weib$estimate[1], scale=dat.weib$estimate[2])
    
    ad=ad.test(sub.quant$Cq, "pweibull", shape=dat.weib$estimate[1], scale=dat.weib$estimate[2])
    
    
    mle.shape=dat.weib$estimate[1] #alpha
    
    mle.scale=dat.weib$estimate[2] #lambda
    
    weibull.mean=as.numeric(mle.scale*gamma(1+(1/mle.shape)))
    
    weibull.var=as.numeric(mle.scale^2*(gamma(1+(2/mle.shape))-gamma(1+(1/mle.shape)))^2)
    
    weibull.sd=sqrt(weibull.var)  
    
    

    coef_var[entry]=round(sd(sub.quant$Cq)/mean(sub.quant$Cq),4)
    
    
    coef_var.standard[entry]=round(sd(sub.quant$Cq)/mean(sub.quant$Cq),4)
    
    p_val.ad[entry]=round(ad$p.value,4)
    p_val.cvm[entry]=round(cvm$p.value,4)
    
    entry=entry+1
    
    
  #  if(ad$p.value<0.05 | cvm$p.value<0.05){
     # print("stop")
#}
    }
    
    

    coef_var_ll=rep(0,length(levels(Dat$SQf)))
    coef_var=rep(0,length(levels(Dat$SQf)))
    p_val.cvm=rep(0,length(levels(Dat$SQf)))
    p_val.ad=rep(0,length(levels(Dat$SQf)))
    entry=1
    

    
    for(i in levels(Dat$SQf)){
      
      sub.quant= Dat.eFISH %>%
        filter(SQf==i) %>%
        filter(!is.na(Cq))
      
      sub.quant=sub.quant[!(is.na(sub.quant$Cq)),]
      
      
      dat.loglogis=fitdist(sub.quant$Cq,"llogis")
      
     # qqPlot(qllogis(ppoints(length(sub.quant$Cq)), shape = dat.loglogis$estimate[1], 
     #                scale =  dat.loglogis$estimate[2]), log(sub.quant$Cq),add.line=TRUE,xlab="Theoretical Quantiles",ylab="Sample Quantiles",main=paste("SQ",i))
      
      
      
      cvm=cvm.test(sub.quant$Cq, "pllogis", shape = dat.loglogis$estimate[1], scale =  dat.loglogis$estimate[2])
      
      ad=ad.test(sub.quant$Cq, "pllogis", shape = dat.loglogis$estimate[1], scale =  dat.loglogis$estimate[2])
      
      
      mle.shape=as.numeric(dat.loglogis$estimate[1]) #Beta
      
      mle.scale=as.numeric(dat.loglogis$estimate[2]) #Alpha
      
      
      b=as.numeric(pi/mle.shape)
      
      loglogis.mean=(mle.scale*b)/sin(b)
      
      loglogis.var=(mle.scale^2)*((2*b/sin(2*b))-(b^2/sin(b)^2))
      
      loglogis.sd=sqrt(loglogis.var)
      
      
      coef_var_ll[entry]=round(loglogis.sd/loglogis.mean,4)
      coef_var[entry]=round(sd(sub.quant$Cq)/mean(sub.quant$Cq),4)
      
      p_val.ad[entry]=round(ad$p.value,4)
      p_val.cvm[entry]=round(cvm$p.value,4)
      
      entry=entry+1

      #if(ad$p.value<0.05 | cvm$p.value<0.05){
      #  print(paste("Log-Logistic does not fit at SQ:",i))
     # }
    #  cat("\n\n")
      
}

    
    
## eONKI4

coef_var=rep(0,length(levels(Dat$SQf)))

coef_var_unbiased=rep(0,length(levels(Dat$SQf)))


p_val.cvm=rep(0,length(levels(Dat$SQf)))
p_val.ad=rep(0,length(levels(Dat$SQf)))

entry=1

Dat.eONKI= Dat %>%
  filter(Target=='eONKI4')


for(i in levels(Dat$SQf)){
  
  sub.quant= Dat.eONKI %>%
    filter(SQf==i)
  
 # qqnorm(sub.quant$Cq,main=paste("SQ",i))
  #qqline(sub.quant$Cq)
  
  sub.quant=sub.quant[!(is.na(sub.quant$Cq)),]
  
  coef_var[entry]=round(sd(sub.quant$Cq)/mean(sub.quant$Cq),4)
  
  coef_var_unbiased[entry]=round((1+1/(4*length(levels(Dat$SQf))))*coef_var[entry],4)
  
  
  cvm=cvm.test(sub.quant$Cq, "pnorm", mean=mean(sub.quant$Cq), sd=sd(sub.quant$Cq))
  
  ad=ad.test(sub.quant$Cq, "pnorm", mean=mean(sub.quant$Cq), sd=sd(sub.quant$Cq))
  
  p_val.ad[entry]=round(ad$p.value,4)
  p_val.cvm[entry]=round(cvm$p.value,4)
  
  entry=entry+1
  
  #if(ad$p.value<0.05 | cvm$p.value<0.05){
 #   print("stop")
 # }
  
 # cat("\n\n")
}
#dev.off()








coef_var=rep(0,length(levels(Dat$SQf))) 
coef_var.standard=rep(0,length(levels(Dat$Sqf)))

p_val.cvm=rep(0,length(levels(Dat$SQf)))
p_val.ad=rep(0,length(levels(Dat$SQf)))

entry=1


for(i in levels(Dat$SQf)){
  
  sub.quant= Dat.eONKI %>%
    filter(SQf==i) %>%
    filter(!is.na(Cq)) #We don't include missing NA values.
  
  sub.quant=sub.quant[!(is.na(sub.quant$Cq)),]
  
  
  quant.mean=mean(sub.quant$Cq)
  quant.sd=sd(sub.quant$Cq)
  
  
  dat.weib=fitdist(sub.quant$Cq,"weibull")
  
  #qqPlot(qweibull(ppoints(length(sub.quant$Cq)), shape = dat.weib$estimate[1], 
   #               scale =  dat.weib$estimate[2]), sub.quant$Cq,add.line=TRUE,xlab="Theoretical Quantiles",ylab="Sample Quantiles",main=paste("SQ",i))
  
  
  mle.shape=dat.weib$estimate[1] #alpha
  
  mle.scale=dat.weib$estimate[2] #lambda
  
  
  weibull.mean=as.numeric(mle.scale*gamma(1+(1/mle.shape)))
  
  weibull.var=as.numeric(mle.scale^2*(gamma(1+(2/mle.shape))-gamma(1+(1/mle.shape)))^2)
  
  weibull.sd=sqrt(weibull.var)  
  
  
  
  #coef_var[entry]=round(weibull.sd/weibull.mean,4)
  coef_var[entry]=round(sd(sub.quant$Cq)/mean(sub.quant$Cq),4)
  
  #coef_var.standard[entry]=round(sd(sub.quant$Cq)/mean(sub.quant$Cq),4)
  
  
  cvm=cvm.test(sub.quant$Cq, "pweibull", shape=dat.weib$estimate[1], scale=dat.weib$estimate[2])
  
  ad=ad.test(sub.quant$Cq, "pweibull", shape=dat.weib$estimate[1], scale=dat.weib$estimate[2])
  
  p_val.ad[entry]=round(ad$p.value,4)
  p_val.cvm[entry]=round(cvm$p.value,4)
  
  entry=entry+1
  
  
  
 # if(ad$p.value<0.05 | cvm$p.value<0.05){
  #  print(paste("Weibull does not fit at SQ:",i))
 # }
 # cat("\n\n")
  
}
#dev.off()


coef_var_a=rep(0,length(levels(Dat$SQf)))
coef_var=rep(0,length(levels(Dat$SQf)))

p_val.cvm=rep(0,length(levels(Dat$SQf)))
p_val.ad=rep(0,length(levels(Dat$SQf)))

entry=1
par(mfrow=c(3,2))

for(i in levels(Dat$SQf)){
  
  sub.quant= Dat.eONKI %>%
    filter(SQf==i)
  
  sub.quant=sub.quant[!(is.na(sub.quant$Cq)),]
  
  
  dat.loglogis=fitdist(sub.quant$Cq,"llogis",lower=c(0,0)) #Fit a Log-logistic distribution using fitdist
  
  
  
  #qqPlot(qllogis(ppoints(length(sub.quant$Cq)), shape = dat.loglogis$estimate[1], 
    #             scale =  dat.loglogis$estimate[2]), log(sub.quant$Cq),add.line=TRUE,xlab="Theoretical Quantiles",ylab="Sample Quantiles",main=paste("SQ",i))
  
  
  mle.shape=as.numeric(dat.loglogis$estimate[1]) #Beta
  
  
  mle.scale=as.numeric(dat.loglogis$estimate[2]) #Alpha
  
  
  
  b=as.numeric(pi/mle.shape)
  
  loglogis.mean=(mle.scale*b)/sin(b)
  
  
  loglogis.var=(mle.scale^2)*((2*b/sin(2*b))-(b^2/sin(b)^2))
  
  
  loglogis.sd=sqrt(loglogis.var)
  
  
  
  coef_var_a[entry]=round(loglogis.sd/loglogis.mean,4)
  
  coef_var[entry]=round(sd(sub.quant$Cq)/mean(sub.quant$Cq),4)
  

  
  cvm=cvm.test(sub.quant$Cq, "pllogis", shape = dat.loglogis$estimate[1], scale = dat.loglogis$estimate[2])
  ad=ad.test(sub.quant$Cq, "pllogis", shape = dat.loglogis$estimate[1], scale = dat.loglogis$estimate[2])
  
  
  p_val.ad[entry]=round(ad$p.value,4)
  p_val.cvm[entry]=round(cvm$p.value,4)
  
  entry=entry+1
  
  
  
  #if(ad$p.value<0.05 | cvm$p.value<0.05){
 #   print(paste("Log-Logistic does not fit at SQ:",i))
 # }
 # cat("\n\n")
  
}


# All Targets

full_targets=c("eASMO9","eASTR4","eLIPI1", "eMIDO1","eMISA2","eRAAU1", "eRACA2","eRALU2", "eRAPR2")


Dat.new=Dat%>%
  filter(Target!='eONKI4')%>%
  filter(Target!='eFISH1')%>%
  filter(SQf!=1)%>%
  filter(SQf!=2)%>%
  filter(SQf!=3)


invisible(droplevels(Dat.new$SQf))
Dat.new$SQf=factor(Dat.new$SQf)

for(t in full_targets){
  
  Dat.e= Dat.new %>%
    filter(Target==t) #Restrict attention to each target
  
  
  coef_var=rep(0,length(levels(Dat$SQf))) #Initalize a vector to store the regular cvs
  
  coef_var_unbiased=rep(0,length(levels(Dat$SQf))) #Initalize a vector to store the unbiased cvs
  
  p_val.cvm=rep(0,length(levels(Dat$SQf))) #Initalize a vector to store the p values from the cvm test
  p_val.ad=rep(0,length(levels(Dat$SQf))) #Initalize a vector to store the p values from the ad test
  
  entry=1 # this is used to loop over and populate the vectors
  
  

  # SQ of  1, 2 and 3 were only for eONKI4 and eFISH1 so we can't loop over these.
  # Dat
  
  for(i in levels(Dat.new$SQf)){ # Loop over all SQ of intrest
    
    sub.quant=Dat.e %>% filter(SQf==i)

  }
  cat("\n\n")
}





