

# This script produces the various qqplot under several distributions for multiple eDNA targets. We focus on eFISH1 and eONKI4 but include
#normal qqplots for other targets as well.


# eFISH1

par(mfrow=c(3,2))

for(i in levels(Dat$SQf)){ # Loop over all SQ of intrest
  
  sub.quant= Dat.eFISH %>%
    filter(SQf==i)
  
   qqnorm(sub.quant$Cq,main=paste("Normal QQ-Plot:","SQ",i),sub='Target: eFISH1') # Make the Q-Q plot
   qqline(sub.quant$Cq) #Add the Q-Q line
}


par(mfrow=c(3,2))

for(i in levels(Dat$SQf)){
  
  sub.quant= Dat.eFISH %>%
    filter(SQf==i)
  
  sub.quant=sub.quant[!(is.na(sub.quant$Cq)),] #Remove missing entries
  
  qqnorm(log(sub.quant$Cq),main=paste("Log-Normal QQ-Plot:","SQ",i),sub='Target: eFISH1')
  qqline(log(sub.quant$Cq))
  
}



par(mfrow=c(3,2))

for(i in levels(Dat$SQf)){
  
  sub.quant= Dat.eFISH %>%
    filter(SQf==i) %>%
    filter(!is.na(Cq))
  
  sub.quant=sub.quant[!(is.na(sub.quant$Cq)),]
  
  dat.weib=fitdist(sub.quant$Cq,"weibull") #Fit a Weibull distribution using fitdist to estimate the parameters
  
   qqPlot(qweibull(ppoints(length(sub.quant$Cq)), shape = dat.weib$estimate[1], 
              scale =  dat.weib$estimate[2]), sub.quant$Cq,add.line=TRUE,xlab="Theoretical Quantiles",ylab="Sample Quantiles",main=paste("Weibull QQ-Plots:","SQ",i),sub='Target: eFISH1')
  
}




par(mfrow=c(3,2))

for(i in levels(Dat$SQf)){
  
  sub.quant= Dat.eFISH %>%
    filter(SQf==i) %>%
    filter(!is.na(Cq))
  
  sub.quant=sub.quant[!(is.na(sub.quant$Cq)),]
  
  
  dat.loglogis=fitdist(sub.quant$Cq,"llogis")
  
  qqPlot(qllogis(ppoints(length(sub.quant$Cq)), shape = dat.loglogis$estimate[1], 
      scale =  dat.loglogis$estimate[2]), log(sub.quant$Cq),add.line=TRUE,xlab="Theoretical Quantiles",ylab="Sample Quantiles",main=paste("Log-Logistic QQ-Plots:","SQ",i),sub='Target: eFISH1')
  
}



# eONKI4

Dat.eONKI= Dat %>%
  filter(Target=='eONKI4')

par(mfrow=c(3,2))

for(i in levels(Dat$SQf)){
  
  sub.quant= Dat.eONKI %>%
    filter(SQf==i)
  
  qqnorm(sub.quant$Cq,main=paste("Normal QQ-Plot:","SQ",i),sub="Target: eONKI4")
  qqline(sub.quant$Cq)
  
}


par(mfrow=c(3,2))

for(i in levels(Dat$SQf)){
  
  sub.quant= Dat.eONKI %>%
    filter(SQf==i)
  
  qqnorm(log(sub.quant$Cq),main=paste("Log-Normal QQ-Plot:","SQ",i),sub="Target: eONKI4")
  qqline(log(sub.quant$Cq))
  
}



par(mfrow=c(3,2))


for(i in levels(Dat$SQf)){
  
  sub.quant= Dat.eONKI %>%
    filter(SQf==i) %>%
    filter(!is.na(Cq)) #We don't include missing NA values.
  
  sub.quant=sub.quant[!(is.na(sub.quant$Cq)),]
  
  
  quant.mean=mean(sub.quant$Cq)
  quant.sd=sd(sub.quant$Cq)
  
  
  dat.weib=fitdist(sub.quant$Cq,"weibull")
  
  qqPlot(qweibull(ppoints(length(sub.quant$Cq)), shape = dat.weib$estimate[1], 
                  scale =  dat.weib$estimate[2]), sub.quant$Cq,add.line=TRUE,xlab="Theoretical Quantiles",ylab="Sample Quantiles",main=paste("Weibull QQ-Plot:","SQ",i),sub="Target eONKI4")
}


par(mfrow=c(3,2))

for(i in levels(Dat$SQf)){
  
  sub.quant= Dat.eONKI %>%
    filter(SQf==i)
  
  sub.quant=sub.quant[!(is.na(sub.quant$Cq)),]
  
  
  dat.loglogis=fitdist(sub.quant$Cq,"llogis",lower=c(0,0)) #Fit a Log-logistic distribution using fitdist
  
  
  
  qqPlot(qllogis(ppoints(length(sub.quant$Cq)), shape = dat.loglogis$estimate[1], 
                 scale =  dat.loglogis$estimate[2]), log(sub.quant$Cq),add.line=TRUE,xlab="Theoretical Quantiles",ylab="Sample Quantiles",main=paste('Log-Logistic QQ-Plot',"SQ",i),sub="Target: eONKI4")
  

}


# Other Targets

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
  
  
  par(mfrow=c(3,3)) #Set plotting parameters
  
  # SQ of  1, 2 and 3 were only for eONKI4 and eFISH1 so we can't loop over these.
  # Dat
  
  for(i in levels(Dat.new$SQf)){ # Loop over all SQ of intrest
    
    sub.quant=Dat.e %>% filter(SQf==i)
    #sub.quant=Dat.e%>%filter(!(is.na(Cq)))
    
    qqnorm(sub.quant$Cq,main=paste("SQ",i,"","Target",t),sub="Normal QQ Plot") # Make the Q-Q plot
    qqline(sub.quant$Cq) #Add the Q-Q line
  }
}


