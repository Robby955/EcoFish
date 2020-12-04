# This R script runs the density data analysis.


options(digits=2)

#Set your working directory to where you have "DensityDataUpdated.csv" and "densitybiomassoriginal.csv" saved.


# Load packages and libraries.
library(knitr) #For rendering output to pdf
library(dplyr) #For piping
library(kableExtra) #For making kable plots
library(Matrix) #For some matrix operations
library(fit.models) #For fitting some models
library(robust) #For robust analysis
library(tinytex) #For latex 
library(lme4) #For fitting linear mixed effect models
library(ggplot2)
library(ggthemes)
library(investr) #For calibration
library(latex2exp) #For latex in ggplots
library(broom) # plotting residuals in ggplot
library(patchwork) #Common Titles in ggplot


#Read in the initial dataset.
eco.dat<-read.csv("DensityDataUpdated.csv",fileEncoding="UTF-8-BOM")


#Check out the original dimension.
print(paste("The original dimension of the dataset is", nrow(eco.dat),"by",ncol(eco.dat)))

#Remove unimportant/blank columns that appeared in the raw data.
eco.dat <- eco.dat[,-(22:39)]
eco.dat$TCT<-eco.dat$Transformed.Ct..50.001.Ct.
eco.dat$CT<-eco.dat$Adjusted.Ct...N.A.50.

eco.dat=eco.dat %>%   #Remove sort code 128 since it did not pass the integritE DNA test in the lab.
  filter(Sort.Code != 128) 

# Use grep to extract number of fish and tank.
fish.ind0<-grep("^0 Fish,", eco.dat$Site.ID)
fish.ind1<-grep("^1 Fish,", eco.dat$Site.ID)
fish.ind2<-grep("^2 Fish,", eco.dat$Site.ID)
fish.ind4<-grep("^4 Fish,", eco.dat$Site.ID)
fish.ind8<-grep("^8 Fish,", eco.dat$Site.ID)
fish.ind16<-grep("^16 Fish,", eco.dat$Site.ID)
fish.ind32<-grep("^32 Fish,", eco.dat$Site.ID)
fish.ind65<-grep("^65 Fish,", eco.dat$Site.ID)

# Careful to make sure we have exact matches, we again use group to subset our four main tanks.

tank.ind1<-grep(" Tank 1$", eco.dat$Site.ID)
tank.ind2<-grep(" Tank 2", eco.dat$Site.ID)
tank.ind3<-grep(" Tank 3", eco.dat$Site.ID)
tank.ind4<-grep(" Tank 4",eco.dat$Site.ID)
tank.ind19<-grep(" Tank 19$", eco.dat$Site.ID)
tank.ind20<-grep(" Tank 20$", eco.dat$Site.ID)
tank.ind21<-grep(" Tank 21$", eco.dat$Site.ID)
tank.ind24<-grep(" Tank 24$", eco.dat$Site.ID)

# Add Tank and Fish id as new seperate columns to our original data frame.

eco.dat$Tank <- 1
eco.dat$Tank[tank.ind2] <- 2
eco.dat$Tank[tank.ind3] <- 3
eco.dat$Tank[tank.ind4] <- 4
eco.dat$Tank[tank.ind19] <- 19
eco.dat$Tank[tank.ind20] <- 20
eco.dat$Tank[tank.ind21] <- 21
eco.dat$Tank[tank.ind24] <- 24

eco.dat$Fish <- 0
eco.dat$Fish[fish.ind1] <- 1
eco.dat$Fish[fish.ind2] <- 2
eco.dat$Fish[fish.ind4] <- 4
eco.dat$Fish[fish.ind8] <- 8
eco.dat$Fish[fish.ind16] <- 16
eco.dat$Fish[fish.ind32] <- 32
eco.dat$Fish[fish.ind65] <- 65

eco.dat$TankF<-as.factor(eco.dat$Tank) #Add new column, tank as factor. 
eco.dat$FishF<-as.factor(eco.dat$Fish) #Add new column, Fish as factor.

# Manually assign dates according to experiment.

dp='08-12' # Pilot experiment date.
d1='08-19' # The official experiment began on August 19, 2015 and continued through that week on a daily basis.
d2='08-20'
d3='08-21'
d4='08-22'
d5='08-23'
d6='08-24'
d7='08-25'

eco.dat$Date<-"NA"

eco.dat$Date[fish.ind1] <- d1 # The number of fish was doubled each day, starting at 1 fish.
eco.dat$Date[fish.ind2] <- d2
eco.dat$Date[fish.ind4] <- d3
eco.dat$Date[fish.ind8] <- d4
eco.dat$Date[fish.ind16] <- d5
eco.dat$Date[fish.ind32] <- d6
eco.dat$Date[fish.ind65] <- d7

#Sort codes 141-161 were taken three a day starting on August 19.

for(i in 1:nrow(eco.dat)){
  if(eco.dat[i,]$Sort.Code %in% c(141:143)){
    eco.dat[i,]$Date=d1
  }
  else if(eco.dat[i,]$Sort.Code %in% c(144:146)){
    eco.dat[i,]$Date=d2
  }
  else if(eco.dat[i,]$Sort.Code %in% c(147:149)){
    eco.dat[i,]$Date=d3
  }
  else if(eco.dat[i,]$Sort.Code %in% c(150:152)){
    eco.dat[i,]$Date=d4
  }
  else if(eco.dat[i,]$Sort.Code %in% c(153:155)){
    eco.dat[i,]$Date=d5
  }
  else if(eco.dat[i,]$Sort.Code %in% c(156:158)){
    eco.dat[i,]$Date=d6
  }
  else if(eco.dat[i,]$Sort.Code %in% c(159:161)){
    eco.dat[i,]$Date=d7
  }
  else if(eco.dat[i,]$Sort.Code %in% c(162:173)){
    eco.dat[i,]$Date=dp
  }
}


mass.original <- read.csv('densitybiomassoriginal.csv', stringsAsFactors = FALSE, 
                          na.strings=c("NA", "NaN", "", 0),fileEncoding="UTF-8-BOM")

mass.original.sub <- mass.original[10:16, 8:14]
mass.original.sub <- mass.original.sub[,-c(2, 4, 6)] #Remove unneeded columns
iTank.original <- unique(mass.original$Tank) #Find the tank numbers
rep.tank <- rep(iTank.original, 7) #Repeat each of the tanks for each of the seven days

#There were seven days, 1:7 , for each day we wish to collapse over the four tanks 19, 20 , 21 and 24.
rep.fish <- rep(c(1, 2, 4, 8, 16, 32 ,65), each=4) #Each of the number of fish.
rep.day <- rep(1:7, each=4) # Each day of the week.
biomass <- as.numeric(t(as.matrix(mass.original.sub))) # We want biomass as a numeric.
eco.biomass <- cbind(Day=rep.day, Tank=rep.tank, Biomass=biomass, Fish=rep.fish) #Set up the dataframe with names.
eco.biomass<-data.frame(eco.biomass)


options(knitr.kable.NA = '') # We want NA to be blank.

# A function that computes unique elements of a vector.
numrep <- function(x){ length(unique(x))} 

#Apply the function to our dataframe by Sort Code
reps<- tapply(eco.dat$Sort.Code, list(eco.dat$Fish, eco.dat$Tank), numrep)

#Create the first kable plot.

k1=kable(eco.biomass,caption="Summary of Biomass (gram)",format="latex",booktabs=T,label="Table 1")%>%
  kable_styling(latex_options = 'striped')

k1=kable(eco.biomass,format="latex",booktabs=T,label="Table 1")%>%
  kable_styling(latex_options = 'striped')

reps<-as.data.frame(reps)

#Create the second table, sample replicate by number of fish and tank.

k2=kable(reps,caption="Number of Sample.replicates by number of fish and tank number.",booktabs=T,format="latex")%>%
  kable_styling(latex_options = "striped")%>%
  add_header_above(c("Fish"=1,"Tank"=8))

k2=kable(reps,booktabs=T,format="latex")%>%
  kable_styling(latex_options = "striped")%>%
  add_header_above(c("Fish"=1,"Tank"=8))



# We create several useful vectors that we will use to create our plots.
iFish <- c(1, 2, 4, 8, 16, 32, 65) # The number of fish
iTank <- c(19, 20, 21, 24) #The tanks used in the main experiment.
xvals<-c(1,2,3,4,5) #The sample replicates used.

# Initialize counter variables.
count1<-1
count2<-1
count3<-1

init.frame<-matrix(0,nrow=length(iFish)*length(iTank),ncol=length(iTank)+2) # An initlized dataframe that we will use to later store values.


for (j in iFish){
  for (k in iTank){
    trctall <- eco.dat$Transformed.Ct..50.001.Ct.[eco.dat$Fish==j & eco.dat$Tank==k]
    line.fill<-c(j,k,min(trctall),max(trctall),median(trctall),round(mean(trctall),2))
    init.frame[count1, ]<-line.fill
    
    count1<-count1+1
  }   
}


check.sort=c(141:173) #These correspond to the sort codes of 0 fish in tank 1 taken during the pilot experiment
check.sort2=c(174:193) #These sort codes refer to 0 fish taken during the second set of experiments on 0 fish.

new.frame=matrix(0,nrow=length(check.sort),ncol=8)
new.frame2=matrix(0,nrow=length(check.sort2),ncol=7)


#Populate the first frame. We wish to know the tank number, the max, minimum, median and mean TCT.

for(i in check.sort){
  
  trctall <- eco.dat$Transformed.Ct..50.001.Ct.[eco.dat$Sort.Code==i]
  tk=unique(eco.dat$Tank[eco.dat$Sort.Code==i])
  td=unique(eco.dat$Date[eco.dat$Sort.Code==i])
  line.fill<-c(td,0,tk,i,min(trctall),max(trctall),median(trctall),round(mean(trctall),2))
  new.frame[count2, ]<-line.fill
  count2<-count2+1
}

#Populate the second frame. We wish to know the tank number, the max, minimum, median and mean TCT.

for(i in check.sort2){
  trctall <- eco.dat$Transformed.Ct..50.001.Ct.[eco.dat$Sort.Code==i]
  tankN<-unique(eco.dat$Tank[eco.dat$Sort.Code==i])
  line.fill<-c(0,tankN,i,min(trctall),max(trctall),median(trctall),round(mean(trctall),2))
  new.frame2[count3, ]<-line.fill
  count3<-count3+1
}



#Table of results for fish,tank and corresponding Min TCT, Max TCT and Median/Mean TCT.
init.frame<-data.frame(init.frame)

init.frame <- setNames(init.frame, c("Fish","Tank","Min TCT","Max TCT","Median TCT","Mean TCT"))

k3=kable(init.frame, caption="Summary of TCTs by Tank and Fish over samples and techincal replicates.",booktabs=T,format='latex')%>%
  kable_styling(latex_options = 'striped')


k3=kable(init.frame,booktabs=T,format='latex')%>%
  kable_styling(latex_options = 'striped')


new.frame<-data.frame(new.frame)

new.frame<-setNames(new.frame,c("Date (2015)","Fish","Tank","Sort Code","Min TCT","Max TCT","Median TCT","Mean TCT"))

k4=kable(new.frame,caption="Zero Fish Summary-Pilot study.",format='latex',booktabs=T)%>%
  kable_styling(latex_options = 'striped',full_width = TRUE)

k4=kable(new.frame,format='latex',booktabs=T)%>%
  kable_styling(latex_options = 'striped',full_width = TRUE)


new.frame2<-data.frame(new.frame2)

new.frame2<-setNames(new.frame2,c("Fish","Tank","Sort Code","Min TCT","Max TCT","Median TCT","Mean TCT"))


k5=kable(new.frame2,caption="Zero Fish Summary",format='latex',booktabs=T)%>%
  kable_styling(latex_options = 'striped',full_width = TRUE)

k5=kable(new.frame2,format='latex',booktabs=T)%>%
  kable_styling(latex_options = 'striped',full_width = TRUE)

#Initilize a new column that we will fill with the tank value.

eco.dat$numtank=0

# Populate the new column.

for(i in 1:nrow(eco.dat)){
  
  if(eco.dat[i,'TankF']=='19'){
    
    eco.dat[i,'numtank']='Tank 19'
  }
  if(eco.dat[i,'TankF']=='20'){
    
    eco.dat[i,'numtank']='Tank 20'
  }
  if(eco.dat[i,'Tank']=='21'){
    
    eco.dat[i,'numtank']='Tank 21'
  }
  else if(eco.dat[i,'TankF']=='24'){
    
    eco.dat[i,'numtank']='Tank 24'
  }
  else if(eco.dat[i,'TankF']=='1'){
    
    eco.dat[i,'numtank']='Tank 1'
  }
  else if(eco.dat[i,'TankF']=='2'){
    eco.dat[i,'numtank']='Tank 2'
  }
  else if(eco.dat[i,'TankF']=='3'){
    eco.dat[i,'numtank']='Tank 3'
  }
  else if(eco.dat[i,'TankF']=='4'){
    eco.dat[i,'numtank']='Tank 4'
  }
}


#Init a new column for Fish.

eco.dat$numFish=0

#Populate the new column.

for(i in 1:nrow(eco.dat)){
  if(eco.dat[i,'Fish']=='1'){
    eco.dat[i,'numFish']='1 Fish'
  }
  if(eco.dat[i,'Fish']=='2'){
    eco.dat[i,'numFish']='2 Fish'
  }
  if(eco.dat[i,'Fish']=='4'){
    eco.dat[i,'numFish']='4 Fish'
  }
  else if(eco.dat[i,'Fish']=='8'){
    eco.dat[i,'numFish']='8 Fish'
  }
  else if(eco.dat[i,'Fish']=='16'){
    eco.dat[i,'numFish']='16 Fish'
  }
  else if(eco.dat[i,'Fish']=='32'){
    eco.dat[i,'numFish']='32 Fish'
  }
  else if(eco.dat[i,'Fish']=='65'){
    eco.dat[i,'numFish']='65 Fish'
  }
  else if(eco.dat[i,'Fish']=='0'){
    eco.dat[i,'numFish']='0 Fish'
  }
}

#Turn new columns into factors for plotting.

eco.dat$numFish <- factor(eco.dat$numFish,levels=c('0 Fish','1 Fish','2 Fish','4 Fish','8 Fish','16 Fish','32 Fish','65 Fish'))

eco.dat$numT<- factor(eco.dat$numtank,levels=c('Tank 1','Tank 2','Tank 3','Tank 4','Tank 19','Tank 20','Tank 21','Tank 24'))


jitter <- position_jitter(width = 0.15, height =0)

#Set up a labelling parameter.

new_labels <- c("1" = "Tank 19", "2" = "Tank 20", "3" = "Tank 21", "4" = "Tank 24")

gz_fish=ggplot(data=eco.dat%>%filter(Fish!=0))+
  geom_point(mapping=aes(x=Sample.replicate,y=TCT),alpha=0.5,position=jitter)+
  facet_grid(numFish~numtank)+
  ylim(0,25)+
  labs(y="TCT", x="Sample Replicate")+theme_bw()+ scale_y_continuous(breaks=c(0,12.5,25),
                                                                     labels=c("0", "12.5", "25"),limits=c(0,25))+ theme(strip.background =element_rect(fill="lightblue"))
#ggsave('gt.pdf',gz_fish,width=6,height=6,dpi=400)


rtt=range(141,161) # Select range of sort code for which we which to exclude.

`%notin%` <- Negate(`%in%`) #Define a negation operator.


gz_zerofish=ggplot(data=eco.dat%>%filter(Fish==0)%>%filter(Sort.Code %notin% 141:161))+
  geom_point(mapping=aes(x=Sample.replicate,y=TCT),position=jitter,alpha=0.5)+
  facet_wrap(numFish~numT,nrow=2)+
  ylim(0,20)+
  labs(y="TCT", x="Sample Replicate")+theme_bw()+ theme(strip.background =element_rect(fill="lightblue"))+theme( axis.line = element_line(colour = "black", 
                                                                                                                                                                                size = 0.5, linetype = "solid"))
#ggsave('gz2.pdf',gz_zerofish,width=6,height=6,dpi=400)


gz_zerofish2=ggplot(data=eco.dat%>%filter(Fish==0)%>%filter(Sort.Code %notin% 141:161))+
  geom_point(mapping=aes(x=Sample.replicate,y=TCT),position=jitter,alpha=0.5)+
  facet_wrap(.~numT,nrow=2)+
  ylim(0,20)+
  labs(y="TCT", x="Sample Replicate")+theme_bw()+ theme(strip.background =element_rect(fill="lightblue"))+theme( axis.line = element_line(colour = "black", 
                                                                                                                                                                                
                                                                                                                                                                                size = 0.5, linetype = "solid"))
#ggsave('gz3.pdf',gz_zerofish2,width=6,height=6,dpi=400)


eco.dat$NewCode=0

for(i in 1:nrow(eco.dat)){
  if(eco.dat[i,'Lab.Code'] %in% c('0.01.0.1','0.01.0.2', '0.01.0.3' )){
    eco.dat[i,'NewCode']='Set 1'
  }
  else if(eco.dat[i,'Lab.Code'] %in% c('0.02.0.1','0.02.0.2', '0.02.0.3' )){
    eco.dat[i,'NewCode']='Set 2'
  }
  else if(eco.dat[i,'Lab.Code'] %in% c('0.04.0.1','0.04.0.2', '0.04.0.3' )){
    eco.dat[i,'NewCode']='Set 3'
  }
  else    if(eco.dat[i,'Lab.Code'] %in% c('0.08.0.1','0.08.0.2', '0.08.0.3' )){
    eco.dat[i,'NewCode']='Set 4'
  }
  else  if(eco.dat[i,'Lab.Code'] %in% c('0.16.0.1','0.16.0.2', '0.16.0.3' )){
    eco.dat[i,'NewCode']='Set 5'
  }
  else  if(eco.dat[i,'Lab.Code'] %in% c('0.32.0.1','0.32.0.2', '0.32.0.3' )){
    eco.dat[i,'NewCode']='Set 6'
  }
  else  if(eco.dat[i,'Lab.Code'] %in% c('0.65.0.1','0.65.0.2', '0.65.0.3' )){
    eco.dat[i,'NewCode']='Set 7'
  }
}


gz_morezeros=ggplot(data=eco.dat%>%filter(Fish==0)%>%filter(Sort.Code %in% 141:161))+
  geom_point(mapping=aes(x=Sample.replicate,y=TCT),position=jitter,alpha=0.5)+
  facet_wrap(NewCode~.,nrow=2)+
  ylim(0,20)+ scale_x_continuous(breaks=c(1,2,3),
                                 labels=c('1', '2', '3'))+
  labs(y="TCT", x="Sample Replicate")+theme_bw()+
  theme(strip.background =element_rect(fill="lightblue"))+theme( axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"))+
  theme(panel.grid.major = element_line(colour="white"),panel.grid.minor = element_line(colour="white"))

#ggsave('gz3.pdf',gz_morezeros,width=6,height=6,dpi=400)

eco.dats <- eco.dat[order(eco.dat$Sort.Code),] ##Re-sorts our data in order of increasing sortcode 

#picks off first occurrence of a distinct Sort.Code, assume data sorted by Sort.Code
jind <- c(TRUE, eco.dats$Sort.Code[-length(eco.dats$Sort.Code)]!=eco.dats$Sort.Code[-1])

#Use tapply and factor with the function 'median', factor(eco.dats$Sort.Code) groups into the distinct sort codes.
eco.sum.dat <- data.frame(eco.dats[jind,])

#Add a column giving the median TCT.
eco.sum.dat$TCTmed <- tapply(eco.dats$Transformed.Ct..50.001.Ct., factor(eco.dats$Sort.Code), median) 

#Add a column giving the mean TCT.
eco.sum.dat$TCTmean <- tapply(eco.dats$Transformed.Ct..50.001.Ct., factor(eco.dats$Sort.Code), mean) 

#Add a column giving the standard deviation of TCT.
eco.sum.dat$TCTsd <- tapply(eco.dats$Transformed.Ct..50.001.Ct., factor(eco.dats$Sort.Code), sd)  

eco.sum.dat <- eco.sum.dat[eco.sum.dat$Fish!=0,] # Ignore zero fish for now.


for(i in iFish){
  for(j in iTank){
    eco.dat$TbioM[eco.dat$Fish==i & eco.dat$Tank==j] <- 
      eco.biomass$Biomass[eco.biomass$Fish==i & eco.biomass$Tank==j]
    eco.sum.dat$TbioM[eco.sum.dat$Fish==i & eco.sum.dat$Tank==j] <- 
      eco.biomass$Biomass[eco.biomass$Fish==i & eco.biomass$Tank==j]
  }
}


eco.sum.out.dat<-eco.sum.dat

eco.med<-eco.sum.out.dat$TCTmed # median transformed CT column corresponding to our data.
eco.mean<-eco.sum.out.dat$TCTmean
eco.tot<-eco.sum.out.dat$TbioM  # Total biomass column corresponding to our data.
eco.tanks<-eco.sum.out.dat$Tank #Non factor version of tank number.

eco.sum.out.dat$TankF=as.factor(eco.sum.out.dat$Tank) #Create a new column,TankF, which gives tank as a factor variable.
eco.sum.dat$TankF=as.factor(eco.sum.dat$Tank) #Create a new column,TankF, which gives tank as a factor variable.


eco.sum.out.dat$l2biom=log2(eco.sum.out.dat$TbioM) #Add as a new column



meanTCT=eco.sum.out.dat$TCTmean

sdTCT=eco.sum.out.dat$TCTsd

#Factor version of tank column. 
tankF=eco.sum.out.dat$TankF 

eco.sum.out.dat$tankF=tankF 


tankFm<-as.matrix(tankF) 


options(digits=2)


l.one.line<-lm(TCTmed~l2biom,data=eco.sum.out.dat) #Simple linear model based on just biomass.
lmparallel.tfac<-lm(TCTmed~l2biom+tankF,data=eco.sum.out.dat) #Biomass and Tank as a factor.
lfull.tfac<-lm(TCTmed~l2biom+tankF+l2biom*tankF,data=eco.sum.out.dat) #Include interaction and tank effect as a factor.





l.one.linem<-lm(TCTmean~l2biom,data=eco.sum.out.dat) #Simple linear model based on just biomass.
lmparallel.tfacm<-lm(TCTmean~l2biom+tankF,data=eco.sum.out.dat) #Biomass and Tank as a factor.
lfull.tfacm<-lm(TCTmean~l2biom+tankF+l2biom*tankF,data=eco.sum.out.dat) #Include interaction and tank effect as a factor.

az=anova(l.one.linem,lmparallel.tfacm,lfull.tfacm)

term_name=c("Intercept","Log2(Biomass")
kcoef=c(round(coef(l.one.linem)[[1]],2),round(coef(l.one.linem)[[2]],2))
kse=c(0.40,0.07)
tval=c(29.5,14.3)
pval=c('<2e-16','<2e-16')


dff2=data.frame(term_name,kcoef,kse,tval,pval)
names(dff2)=c("Term","Estimate",'Std Error','t value','Pr(>|t|)')

ktz2=kable(dff2,format="latex",booktabs=T,escape = FALSE)%>%
  kable_styling("striped")





term_name=c("Intercept","Log2(Biomass)","Tank 20","Tank 21","Tank 24")
kcoef=c(round(coef(lmparallel.tfacm)[[1]],2),round(coef(lmparallel.tfacm)[[2]],2),round(coef(lmparallel.tfacm)[[3]],2),round(coef(lmparallel.tfacm)[[4]],2),
        round(coef(lmparallel.tfacm)[[5]],2))
kse=c(0.39,0.06,0.32,0.32,0.32)
tval=c(31.1,17.5,1.8,-5.3,2.0)
pval=c('<2e-16','<2e-16',"0.07","4.6e-07","0.04")


dff2=data.frame(term_name,kcoef,kse,tval,pval)
names(dff2)=c("Term","Estimate",'Std Error','t value','Pr(>|t|)')

ktz2=kable(dff2,format="latex",booktabs=T,escape = FALSE)%>%
  kable_styling("striped")



term_name=c("Intercept","Log2(Biomass)","Tank 20","Tank 21","Tank 24","Log2(Biomass):Tank20","Log2(Biomass):Tank21","Log2(Biomass):Tank24")
kcoef=c(round(coef(lfull.tfacm)[[1]],2),round(coef(lfull.tfacm)[[2]],2),round(coef(lfull.tfacm)[[3]],2),round(coef(lfull.tfacm)[[4]],2),
        round(coef(lfull.tfacm)[[5]],2),round(coef(lfull.tfacm)[[6]],2),round(coef(lfull.tfacm)[[7]],2),round(coef(lfull.tfacm)[[8]],2))
kse=c(0.65,0.11,0.93,0.92,0.92,0.16,0.16,0.16)
tval=c(17.7,10.1,1.7,-2.2,2.7,-1.1,0.41,-2.1)
pval=c('<2e-16','<2e-16',"0.10","0.03","0.008","0.28","0.69","0.04")


dff2=data.frame(term_name,kcoef,kse,tval,pval)
names(dff2)=c("Term","Estimate",'Std Error','t value','Pr(>|t|)')

ktz2=kable(dff2,format="latex",booktabs=T,escape = FALSE)%>%
  kable_styling("striped")









coef_names=c("Intercept","Log2(Biomass)","Tank 20","Tank 21","Tank 24","Log2(Biomass):Tank20","Log2(Biomass):Tank21","Log2(Biomass):Tank24")
coef_terms=c(round(coef(lfull.tfac)[[1]],2),round(coef(lfull.tfac)[[2]],2),round(coef(lfull.tfac)[[3]],1),round(coef(lfull.tfac)[[4]],2),round(coef(lfull.tfac)[[5]],2),
             round(coef(lfull.tfac)[[6]],2),round(coef(lfull.tfac)[[7]],2),round(coef(lfull.tfac)[[8]],2))
d_terms=c(0.58,0.01,0.84,0.83,0.83,0.15,0.14,0.14)
tterms=c(20.12,10.60,1.40,-1.59,2.49,-0.80,-0.27,-1.90)
pterm=c("<2e-16","<2e-16","0.17","0.11","0.014","0.42","0.79","0.06")


dff=data.frame(coef_names,coef_terms,d_terms,tterms,pterm)
names(dff)=c("Term","Estimate",'Std Error','t value','Pr(>|t|)')

ktz=kable(dff,format="latex",booktabs=T,escape = FALSE)%>%
  kable_styling(latex_options = 'striped')

aa=anova(l.one.line,lmparallel.tfac,lfull.tfac)

coef_names=c("Res DF","RSS","DF","SQ","F",'Pr(>F)')
res=c(137,134,131)
RSS=c(294,192,186)
DF=c("","3","3")
SSQ=c("","102","6")
FF=c("","23.87","1.41")
PR=c("","2.2e-12","0.24")


RSS=c(366,240,227)
SSQ=c("","125.1","13.3")
FF=c("","24.1","2.6")
PR=c("","1.8e-12","0.058")


dff=data.frame(res,RSS,DF,SSQ,FF,PR)
names(dff)=c("Res DF","RSS","DF","SQ","F",'Pr(>F)')


ktz=kable(dff,format="latex",booktabs=T,escape = FALSE)%>%
  kable_styling(latex_options = 'striped')





coef_names=c("Intercept","Log2(Biomass)","Tank 20","Tank 21","Tank 24")
coef_terms=c(round(coef(lmparallel.tfac)[[1]],1),round(coef(lmparallel.tfac)[[2]],2),round(coef(lmparallel.tfac)[[3]],1),round(coef(lmparallel.tfac)[[4]],2),round(coef(lmparallel.tfac)[[5]],2))
d_terms=c(0.35,0.05,0.29,0.29,0.29)
tterms=c(35.70,18.65,1.84,-5.35,2.03)
pterm=c("<2e-16","<2e-16","0.068","3.7e-7","0.044")



dff=data.frame(coef_names,coef_terms,d_terms,tterms,pterm)
names(dff)=c("Term","Estimate",'Std Error','t value','Pr(>|t|)')

ktz=kable(dff,format="latex",booktabs=T,escape = FALSE)%>%
  kable_styling(latex_options = 'striped')




options(digits=3)
# Make kables

coef_names=c("Intercept","Log2(Biomass)")
lonelineinter=coef(l.one.line)[1]
lonelineslope=coef(l.one.line)[2]
coef_terms=c(lonelineinter,lonelineslope)

kcoef=c(round(coef(l.one.line)[[1]],1),round(coef(l.one.line)[[2]],2))
kse=c(0.36,0.06)
tval=c(33.90,15.22)
pval=c('<2e-16','<2e-16')


dff=data.frame(coef_names,kcoef,kse,tval,pval)
names(dff)=c("Term","Estimate",'Std Error','t value','Pr(>|t|)')

ktz=kable(dff,format="latex",booktabs=T,escape = FALSE)%>%
  kable_styling(latex_options = 'striped')

library(jtools)
#Extract linear coefficients.
lmcoef=coef(lfull.tfac)
t19coef=c(lmcoef[1],lmcoef[2])
t20coef=c(lmcoef[1]+lmcoef[3],lmcoef[2]+lmcoef[6])
t21coef=c(lmcoef[1]+lmcoef[4],lmcoef[2]+lmcoef[7])
t24coef=c(lmcoef[1]+lmcoef[5],lmcoef[2]+lmcoef[8])




l2biom<-as.matrix(eco.sum.out.dat$l2biom,nrow=nrow(eco.sum.out.dat),ncol=1)
eco.md<-as.matrix(eco.sum.out.dat$TCTmed,nrow=1,ncol=nrow(eco.sum.out.dat))

anova_med=anova(l.one.line,lmparallel.tfac,lfull.tfac,test='F') #anova to compare the three models. Additional sum of squares prinicpal and the F test indicates significant of our l4parallel.

lr<-robust::lmRob(eco.md~l2biom,data=eco.sum.out.dat) #Simple robust model based on just biomass.
lrparallel.tfac<-robust::lmRob(eco.md~l2biom+tankF,data=eco.sum.out.dat) #Robust model with biomass and tank as a factor.
lrfull.tfac<-robust::lmRob(eco.md~l2biom+tankF+l2biom*tankF,data=eco.sum.out.dat) #Include interaction and tank effect as a factor.
anova.lmRob(lr,lrparallel.tfac,lrfull.tfac)  #anova on robust models. #Model with tank included is highly significant. adding in the interaction does not appear to be of value in modelling median transformed CT.


#Extract coefficients from these three robust models.



coef_names=c("Terms","DF","RobustF","Pr(F)")
TERM=c("1","2","3")
DF=c("","3","3")
FF=c("","25.87","3.46")
PR=c("","2.2e-12","0.24")

dff=data.frame(TERM,DF,FF,PR)
names(dff)=c("Terms","Df","RobustF","Pr(F)")


ktz=kable(dff,format="latex",booktabs=T,escape = FALSE)%>%
  kable_styling(latex_options = 'striped')





lrcoef=coef(lrfull.tfac)

t19coef.robust=c(lrcoef[1],lrcoef[2])
t20coef.robust=c(lrcoef[1]+lrcoef[3],lrcoef[2]+lrcoef[6])
t21coef.robust=c(lrcoef[1]+lrcoef[4],lrcoef[2]+lrcoef[7])
t24coef.robust=c(lrcoef[1]+lrcoef[5],lrcoef[2]+lrcoef[8])


term_name=c("Intercept","Log2(Biomass)")
kcoef=c(round(coef(lr)[[1]],1),round(coef(lr)[[2]],2))
kse=c(0.38,0.06)
tval=c(32.90,14.30)
pval=c('<2e-16','<2e-16')


dff2=data.frame(term_name,kcoef,kse,tval,pval)
names(dff2)=c("Term","Estimate",'Std Error','t value','Pr(>|t|)')

ktz2=kable(dff2,format="latex",booktabs=T,escape = FALSE)%>%
  kable_styling("striped")


term_name=c("Intercept","Log2(Biomass)","Tank20","Tank21","Tank24")
kcoef=c(round(coef(lrparallel.tfac)[[1]],2),round(coef(lrparallel.tfac)[[2]],2),round(coef(lrparallel.tfac)[[3]],2),round(coef(lrparallel.tfac)[[4]],2),round(coef(lrparallel.tfac)[[5]],2))
kse=c(0.331,0.049,0.274,0.271,0.269)
tval=c(35.71,21.03,2.45,-6.12,2.22)
pval=c('<2e-16','<2e-16',"0.016","9.6e-09","0.028")


dff2=data.frame(term_name,kcoef,kse,tval,pval)
names(dff2)=c("Term","Estimate",'Std Error','t value','Pr(>|t|)')

ktz2=kable(dff2,format="latex",booktabs=T,escape = FALSE)%>%
  kable_styling("striped")






term_name=c("Intercept","Log2(Biomass)","Tank20","Tank21","Tank24","Log2(Biomass):Tank20",'Log2(Biomass):Tank21','Log2(Biomass):Tank24')
kcoef=c(round(coef(lrfull.tfac)[[1]],2),round(coef(lrfull.tfac)[[2]],2),round(coef(lrfull.tfac)[[3]],2),round(coef(lrfull.tfac)[[4]],2),
        round(coef(lrfull.tfac)[[5]],2),round(coef(lrfull.tfac)[[6]],2),round(coef(lrfull.tfac)[[7]],2),round(coef(lrfull.tfac)[[8]],2))
kse=c(0.85,0.14,1.19,1.21,1.14,0.21,0.20,0.19)
tval=c(13.72,7.43,0.62,-2.07,2.05,-0.06,0.73,-1.53)
pval=c('<2e-16','1.2e-11','0.54','0.04','0.04','0.94','0.47','0.13')


dff2=data.frame(term_name,kcoef,kse,tval,pval)
names(dff2)=c("Term","Estimate",'Std Error','t value','Pr(>|t|)')

ktz2=kable(dff2,format="latex",booktabs=T,escape = FALSE)%>%
  kable_styling("striped")



eco.dat$l2biom<-log2(eco.dat$TbioM)

eco.dat.zeroremoved=eco.dat %>% #eco.dat already had sort.code 128 removed so only need to remove zero fish for this portion.
  filter(Fish !=0)


lineartank19<-lm(TCTmed~l2biom, data=eco.sum.out.dat, subset=(TankF==19))
lineartank20<-lm(TCTmed~l2biom, data=eco.sum.out.dat, subset=(TankF==20))
lineartank21<-lm(TCTmed~l2biom, data=eco.sum.out.dat, subset=(TankF==21))
lineartank24<-lm(TCTmed~l2biom, data=eco.sum.out.dat, subset=(TankF==24))



kt.intercept=c(coef(lineartank19)[1],coef(lineartank20)[1],coef(lineartank21)[1],coef(lineartank24)[1])
kt.intercept_std=c(0.53,0.66,0.64,0.54)
  
kt.slope=c(coef(lineartank19)[2],coef(lineartank20)[2],coef(lineartank21)[2],coef(lineartank24)[2])
kt.slope_std=c(0.09,0.12,0.11,0.09)

kt.r2=round(c(summary(lineartank19)$r.squared,summary(lineartank20)$r.squared,summary(lineartank21)$r.squared,summary(lineartank24)$r.squared),3)
itank=c(19,20,21,24)

dkk=data.frame(itank,round(kt.intercept,3),kt.intercept_std,round(kt.slope,3),kt.slope_std,kt.r2)
names(dkk)=c("Tank","Intercept","Std Error","Slope","Std Error","$R^{2}$")

kzz=kable(dkk,format="latex",booktabs=T,escape = FALSE)%>%
  kable_styling(latex_options = 'striped')


col.vec=c("black","red","green","magenta")

gframe=data.frame(l2biom,eco.med,eco.sum.dat$Tank-18,eco.sum.dat$Tank-18)
names(gframe)=c("l2biom","eco.med","Tank","tank_shape")

gframe$Tank[gframe$Tank==1]="black"
gframe$Tank[gframe$Tank==2]="red"
gframe$Tank[gframe$Tank==3]="green"
gframe$Tank[gframe$Tank==6]="magenta"

gframe$tank_shape[gframe$Tank=='black']=1
gframe$tank_shape[gframe$Tank=='red']=2
gframe$tank_shape[gframe$Tank=='green']=3
gframe$tank_shape[gframe$Tank=='magenta']=6


gframe$tankf=as.factor(gframe$tank_shape)


gg_simple_median=ggplot(data=gframe)+
  geom_point(mapping=aes(x=l2biom,y=eco.med,colour=tankf,shape=tankf))+
  scale_shape_manual(name="Tank",values=c(1,2,3,6),labels=c("T19","T20","T21","T24"))+
  scale_color_manual(name="Tank",values=c("black", "red", "blue","magenta"),labels=c("T19","T20","T21","T24"))+ geom_smooth(mapping=aes(x=l2biom,y=eco.med),color='black',method="lm")+
  xlim(2,9)+
  ylim(10,25)+
  xlab(TeX("Fish Biomass ($log_{2}$ g)"))+
  ylab("Median Transformed CT")+ 
  theme(panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))+ theme(
          legend.position = c(1, 0.05),
          legend.justification = c("right", "bottom"))
# We can fit the plots where intercept is allowed to vary, but constant slope among tanks.
common_int=lmparallel.tfac$coefficients[1]
common_slope=lmparallel.tfac$coefficients[2]


#ggsave("ggplotnew.pdf",gg_simple_median,width=6,height=6,dpi=400)

gg_par_median=ggplot(data=gframe)+
  geom_point(mapping=aes(x=l2biom,y=eco.med,colour=tankf,shape=tankf))+
  scale_shape_manual(name="Tank",values=c(1,2,3,6),labels=c("T19","T20","T21","T24"))+
  scale_color_manual(name="Tank",values=c("black", "red", "blue","magenta"),labels=c("T19","T20","T21","T24"))+
  geom_abline(intercept=common_int,slope=common_slope,color="black")+
  geom_abline(intercept=common_int+lmparallel.tfac$coefficients[3],slope=common_slope,color="red")+
  geom_abline(intercept=common_int+lmparallel.tfac$coefficients[4],slope=common_slope,color="blue")+
  geom_abline(intercept=common_int+lmparallel.tfac$coefficients[5],slope=common_slope,color="magenta")+
  xlim(2,9)+
  ylim(10,25)+
  xlab(TeX("Fish Biomass ($log_{2}$ g)"))+
  ylab("Median Transformed CT")+ 
  theme(panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))+scale_y_continuous(breaks=10:25)


gg_par_median=gg_par_median+ theme(
  legend.position = c(1, 0.05),
  legend.justification = c("right", "bottom"))



#ggsave("parfits.pdf",gg_par_median,width=6,height=6,dpi=400)







gg_tank_linear=ggplot(data=gframe)+
  geom_point(mapping=aes(x=l2biom,y=eco.med,colour=tankf,shape=tankf))+
  scale_shape_manual(name="Tank",values=c(1,2,3,6),labels=c("T19","T20","T21","T24"))+
  scale_color_manual(name="Tank",values=c("black", "red", "blue","magenta"),labels=c("T19","T20","T21","T24"))+
  geom_abline(intercept=t19coef[1],slope=t19coef[2],color="black")+
  geom_abline(intercept=t20coef[1],slope=t20coef[2],color="red")+
  geom_abline(intercept=t21coef[1],slope=t21coef[2],color="blue")+
  geom_abline(intercept=t24coef[1],slope=t24coef[2],color="magenta")+
  xlim(2,9)+
  ylim(10,25)+
  xlab(TeX("Fish Biomass ($log_{2}$ g)"))+
  ylab("Median Transformed CT")+ 
  theme(panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))+scale_y_continuous(breaks=10:25)+ theme(legend.position = c(1, 0.05),legend.justification = c("right", "bottom"))

#ggsave("ggplotnew2.png",gg_tank_linear,width=6,height=6,dpi=400)

col_bluevec=eco.sum.dat$Tank-18
col_bluevec[col_bluevec==3]=4

l.one.line.mean<-lm(TCTmean~l2biom,data=eco.sum.out.dat)
lmparallel.mean<-lm(TCTmean~l2biom+tankF,data=eco.sum.out.dat) 

#Extract coefficients for plotting.
common_mean_int=lmparallel.mean$coefficients[1]
common_mean_slope=lmparallel.mean$coefficients[2]

gframe.mean=data.frame(l2biom,eco.mean,eco.sum.dat$Tank-18,eco.sum.dat$Tank-18)
names(gframe)=c("l2biom","eco.mean","Tank","tank_shape")


gframe.mean$Tank[gframe$Tank==1]="black"
gframe.mean$Tank[gframe$Tank==2]="red"
gframe.mean$Tank[gframe$Tank==3]="green"
gframe.mean$Tank[gframe$Tank==6]="magenta"

gframe.mean$tank_shape[gframe$Tank=='black']=1
gframe.mean$tank_shape[gframe$Tank=='red']=2
gframe.mean$tank_shape[gframe$Tank=='green']=3
gframe.mean$tank_shape[gframe$Tank=='magenta']=6

gframe.mean$tankf=as.factor(gframe.mean$tank_shape)

gg_par_mean=ggplot(data=gframe.mean)+
  geom_point(mapping=aes(x=l2biom,y=eco.mean,colour=tankf,shape=tankf))+
  scale_shape_manual(name="Tank",values=c(1,2,3,6),labels=c("T19","T20","T21","T24"))+
  scale_color_manual(name="Tank",values=c("black", "red", "blue","magenta"),labels=c("T19","T20","T21","T24"))+
  geom_abline(intercept=common_mean_int,slope=common_mean_slope,color="black")+
  geom_abline(intercept=common_mean_int+lmparallel.mean$coefficients[3],slope=common_mean_slope,color="red")+
  geom_abline(intercept=common_mean_int+lmparallel.mean$coefficients[4],slope=common_mean_slope,color="blue")+
  geom_abline(intercept=common_mean_int+lmparallel.mean$coefficients[5],slope=common_mean_slope,color="magenta")+
  xlim(2,9)+
  ylim(10,25)+
  xlab(TeX("Fish Biomass ($log_{2}$ g)"))+
  ylab(TeX("Mean Transformed CT"))+ 
  theme(panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))+scale_y_continuous(breaks=10:25)+
  theme(legend.position = c(1, 0.05),
legend.justification = c("right", "bottom"))

#ggsave('parmean.pdf',gg_par_mean,width=6,height=6,dpi=500)



gg_par_mean=gg_par_mean+labs(title="Mean Transformed CT by Log2 Fish Biomass (gram) \nwith differing intercepts")+ theme(
  legend.position = c(1, 0.05),
  legend.justification = c("right", "bottom"))



lfull.mean<-lm(TCTmean~l2biom+tankF+l2biom*tankF,data=eco.sum.out.dat) #Include interaction and tank effect as a factor.
gframe.mean=data.frame(l2biom,eco.mean,eco.sum.dat$Tank-18,eco.sum.dat$Tank-18)

names(gframe)=c("l2biom","eco.mean","Tank","tank_shape")

gframe.mean$Tank[gframe$Tank==1]="black"
gframe.mean$Tank[gframe$Tank==2]="red"
gframe.mean$Tank[gframe$Tank==3]="green"
gframe.mean$Tank[gframe$Tank==6]="magenta"

gframe.mean$tank_shape[gframe$Tank=='black']=1
gframe.mean$tank_shape[gframe$Tank=='red']=2
gframe.mean$tank_shape[gframe$Tank=='green']=3
gframe.mean$tank_shape[gframe$Tank=='magenta']=6


gframe.mean$tankf=as.factor(gframe.mean$tank_shape)


gg_tank_mean=ggplot(data=gframe.mean)+
  geom_point(mapping=aes(x=l2biom,y=eco.mean,colour=tankf,shape=tankf))+
  scale_shape_manual(name="Tank",values=c(1,2,3,6),labels=c("T19","T20","T21","T24"))+
  scale_color_manual(name="Tank",values=c("black", "red", "blue","magenta"),labels=c("T19","T20","T21","T24"))+ geom_smooth(mapping=aes(x=l2biom,y=eco.mean),color="black",method="lm")+
  xlim(2,9)+
  ylim(8,25)+
  xlab(TeX("Fish Biomass ($log_{2}$ g)"))+
  ylab("Mean Transformed CT")+ 
  theme(panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))+theme(
          legend.position = c(1, 0.05),
          legend.justification = c("right", "bottom"))

#ggsave('ggplotnew3.png',gg_tank_mean,height=6,width=6,dpi=500)

gg_tank_mean=gg_tank_mean+labs(title="Mean Transformed CT by Log2 Fish Biomass (gram)")+theme(
  legend.position = c(1, 0.05),
  legend.justification = c("right", "bottom"))





meanV=init.frame[,6] #Get the mean values from init.frame
medV=init.frame[,5]

l2biomTank=c(rep(l2biom[1],4),rep(l2biom[21],4),rep(l2biom[41],4),rep(l2biom[61],4),rep(l2biom[81],4),rep(l2biom[101],4),rep(l2biom[121],4))

new.tank=data.frame(l2biomTank,meanV) #Form a useful dataframe
new.tankmed=data.frame(l2biomTank,medV)


unique.bio<-unique(l2biomTank) # Extract the unique values of log2(biomass).

iframe<-data.frame(unique.bio,iFish) #Set up a frame of the biomass and fish.

labels.minor <- c("1","2","4","8","16","32","65") #String version of number of fish
label.pos<-unique.bio


names(iframe)<-c("x","Fish")


l.tankregression<-lm(meanV~l2biomTank,data=new.tank) # This is the model that ggplot fits. We have 28 observations. 4 for each of the seven number of fish. This models the mean TCT over each tank and specific number of fish.

l.tankregression.med<-lm(medV~l2biomTank,data=new.tankmed)


coef_names=c("Intercept","Log2(BiomPerTank)")
lonelineinter=coef(l.one.line)[1]
lonelineslope=coef(l.one.line)[2]
coef_terms=c(lonelineinter,lonelineslope)

kcoef=c(round(coef(l.tankregression)[[1]],2),round(coef(l.tankregression)[[2]],2))
kse=c(0.66,0.11)
tval=c(18.02,8.78)
pval=c('3.3e-16','2.9e-9')


dff=data.frame(coef_names,kcoef,kse,tval,pval)
names(dff)=c("Term","Estimate",'Std Error','t value','Pr(>|t|)')

ktz=kable(dff,format="latex",booktabs=T,escape = FALSE)%>%
  kable_styling(latex_options = 'striped')




options(digits=3)

nn=c('l.one.line','lmparallel.tfac','lfull.tfac','l.tankregression.med','l.one.line.mean','lmparallel.tfac.mean','lfull.tfac.mean','l.tankregression')
n2=c('Median TCT','Median TCT','Median TCT','Median TCT','Mean TCT','Mean TCT','Mean TCT','Mean TCT')
r2=c(0.626,0.756,0.763,0.748,0.692,0.735,0.750,0.748)
adjr2=c(0.624,0.748,0.751,0.738,0.688,0.728,0.737,0.739)
dz=data.frame(nn,n2,r2,adjr2)
names(dz)=c("Model","Response","$R^{2}$","adj$R^{2}$")

ktz=kable(dz,format="latex",booktabs=T,escape = FALSE)%>%
  kable_styling(latex_options = 'striped')


coef_names=c("Intercept","Log2(BiomPerTank)")


kcoef=c(round(coef(l.tankregression.med)[[1]],2),round(coef(l.tankregression.med)[[2]],2))
kse=c(0.64,0.11)
tval=c(19.13,8.79)
pval=c("<2e-16","2.9e-09")

dff=data.frame(coef_names,kcoef,kse,tval,pval)
names(dff)=c("Term","Estimate",'Std Error','t value','Pr(>|t|)')

ktz=kable(dff,format="latex",booktabs=T,escape = FALSE)%>%
  kable_styling(latex_options = 'striped')













#summary(l.tankregression)$r.squared #The R^2 value of 0.748 adjusted R squared of 0.7383

gn=ggplot(data=new.tank)+
  geom_point(mapping=aes(x=l2biomTank,y=meanV),shape=0)+
  geom_smooth(mapping=aes(x=l2biomTank,y=meanV),color="black",method="lm")+
  xlim(2,9)+
  ylim(10,25)+
  xlab("Number of Fish per 10,000L Tank")+
  ylab("Mean Transformed CT")+
  theme(panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  scale_x_continuous(position = "bottom") +
  scale_x_continuous(sec.axis=sec_axis(trans=~ . * 1, name=TeX("Coho biomass ($log_{2}$ g)")),labels=iFish,breaks=label.pos,position='top')

#ggsave('ggplotnew5.png',gn,width=6,height=6,dpi=500)



gn_mean=ggplot(data=new.tank)+
  geom_point(mapping=aes(x=l2biomTank,y=meanV),shape=0)+
  geom_smooth(mapping=aes(x=l2biomTank,y=meanV),color="black",method="lm")+
  xlim(2,9)+
  ylim(10,25)+
  xlab("Number of Fish per 10,000L Tank")+
  ylab("Mean Transformed CT")+
  theme(panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  scale_x_continuous(position = "bottom") +
  scale_x_continuous(sec.axis=sec_axis(trans=~ . * 1, name=TeX("Coho biomass ($log_{2}$ g)")),labels=iFish,breaks=label.pos,position='top')

#ggsave('ggplotnew4.png',gn_mean,width=6,height=6,dpi=500)


modular=augment(l.tankregression)

modular$TankNum=rep(c('19','20','21','24'),7)

gresid=ggplot(modular, aes(x = .fitted, y = .resid,col=factor(TankNum),shape=factor(TankNum))) + geom_point()+theme_bw()+ylab("Residuals")+geom_hline(yintercept=0)+scale_color_manual(name="Tank",values=c("black", "red", "blue","magenta"),labels=c("19","20","21","24"))+scale_shape_manual(name="Tank",values=c(1,2,3,6),labels=c("19","20","21","24"))+xlab("Fitted")+ggtitle("Residuals for l.tankregression")




gn2=ggplot(data=new.tankmed)+
  geom_point(mapping=aes(x=l2biomTank,y=medV),shape=0)+
  geom_smooth(mapping=aes(x=l2biomTank,y=medV),color="black",method="lm")+
  xlim(2,9)+
  ylim(10,25)+
  xlab("Number of Fish per 10,000L Tank")+
  ylab("Median Transformed CT")+
  theme(panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  scale_x_continuous(position = "bottom") +
  scale_x_continuous(sec.axis=sec_axis(trans=~ . * 1, name=TeX("Coho biomass ($log_{2}$ g)")),labels=iFish,breaks=label.pos,position='top')


#ggsave('ggplotnew5.png',gn2,width=6,height=6,dpi=500)

lineartank19.mean<-lm(TCTmean~l2biom, data=eco.sum.out.dat, subset=(TankF==19))
lineartank20.mean<-lm(TCTmean~l2biom, data=eco.sum.out.dat, subset=(TankF==20))
lineartank21.mean<-lm(TCTmean~l2biom, data=eco.sum.out.dat, subset=(TankF==21))
lineartank24.mean<-lm(TCTmean~l2biom, data=eco.sum.out.dat, subset=(TankF==24))


kt2.intercept=c(coef(lineartank19.mean)[1],coef(lineartank20.mean)[1],coef(lineartank21.mean)[1],coef(lineartank24.mean)[1])
kt2.slope=c(coef(lineartank19.mean)[2],coef(lineartank20.mean)[2],coef(lineartank21.mean)[2],coef(lineartank24.mean)[2])
kt2.r2=round(c(summary(lineartank19.mean)$r.squared,summary(lineartank20.mean)$r.squared,summary(lineartank21.mean)$r.squared,summary(lineartank24.mean)$r.squared),3)
itank=c(19,20,21,24)

intse=c(0.59,0.68,0.79,0.53)
slopese=c(0.10,0.12,0.14,0.09)


dkk2=data.frame(itank,round(kt2.intercept,3),intse,round(kt2.slope,3),slopese,kt2.r2)
names(dkk2)=c("Tank","Intercept","Std Error","Slope","Std Error","$R^{2}$")

kzz2=kable(dkk2,format="latex",booktabs=T,escape = FALSE)%>%
  kable_styling(latex_options = 'striped')

kzz2=kable(dkk2,format="latex",booktabs=T,escape = FALSE)%>%
  kable_styling(latex_options = 'striped')

# Extract the coefficeints.
t19coef.mean=coef(lineartank19.mean)
t20coef.mean=coef(lineartank20.mean)
t21coef.mean=coef(lineartank21.mean)
t24coef.mean=coef(lineartank24.mean)


gg_mean_full=ggplot(data=gframe.mean)+
  geom_point(mapping=aes(x=l2biom,y=eco.mean,colour=tankf,shape=tankf))+
  scale_shape_manual(name="Tank",values=c(1,2,3,6),labels=c("T19","T20","T21","T24"))+
  scale_color_manual(name="Tank",values=c("black", "red", "blue","magenta"),labels=c("T19","T20","T21","T24"))+
  geom_abline(intercept=t19coef.mean[1],slope=t19coef.mean[2],color="black")+
  geom_abline(intercept=t20coef.mean[1],slope=t20coef.mean[2],color="red")+
  geom_abline(intercept=t21coef.mean[1],slope=t21coef.mean[2],color="blue")+
  geom_abline(intercept=t24coef.mean[1],slope=t24coef.mean[2],color="magenta")+
  xlim(2,9)+
  ylim(8,25)+
  xlab(TeX("Fish Biomass (log_{2} g)"))+
  ylab("Mean Transformed CT")+ 
  theme(panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))+theme(
          legend.position = c(1, 0.05),
          legend.justification = c("right", "bottom"))

#ggsave('ggplotnew7.png',gg_mean_full,width=6,height=6,dpi=500)


gg_mean_full=gg_mean_full+labs(title="Mean Transformed CT by Log2 Fish Biomass (gram) \nwith Tank Specific Linear fits")+ theme(
  legend.position = c(1, 0.05),
  legend.justification = c("right", "bottom"))


## Show some of the plots.
## Note some of the above plots are commented out so won't show when you run the script, you can
## show them by typing them like as below or running them in the console.

#gz_fish
#gz_zerofish2
#gz_morezeros

# Plot Median Models

#gg_simple_median
#gg_par_median
#gg_tank_linear
#gn2

# Plot Mean Models

#gg_tank_mean
#gg_par_mean
#gg_mean_full
#gn


modular=augment(l.tankregression)

modular$TankNum=rep(c('19','20','21','24'),7)

gresid=ggplot(modular, aes(x = .fitted, y = .resid,col=factor(TankNum),shape=factor(TankNum))) + 
  geom_point()+theme_bw()+ylab("Residuals")+geom_hline(yintercept=0)+
  scale_color_manual(name="Tank",values=c("black", "red", "blue","magenta"),labels=c("19","20","21","24"))+
  scale_shape_manual(name="Tank",values=c(1,2,3,6),labels=c("19","20","21","24"))+xlab("Fitted")


# Color Blue Scheme
col_bluevec=eco.sum.dat$Tank-18
col_bluevec[col_bluevec==3]=4


#plot(residuals(lmparallel.tfac),col = eco.sum.dat$Tank-18,pch=eco.sum.dat$Tank-18,ylab="Residuals",las=1)
#abline(h=0)

#plot(residuals(lmparallel.tfac),col =col_bluevec,pch=eco.sum.dat$Tank-18,las=1,ylab="Residuals")
#abline(h=0)
#title("Residuals for lmparallel.tfac")





col_bluevec=eco.sum.dat$Tank-18
col_bluevec[col_bluevec==3]=4


#pdf(file='residuals2.pdf')
par(mfrow=c(2,2))

plot(residuals(l.one.line),col =col_bluevec,pch=eco.sum.dat$Tank-18,las=1,ylab="Residuals")
abline(h=0)
title("Residuals for linear fit l.one.line")

plot(residuals(lmparallel.tfac),col =col_bluevec,pch=eco.sum.dat$Tank-18,las=1,ylab="Residuals")
abline(h=0)
title("Residuals for lmparallel.tfac")


plot(residuals(lfull.tfac),col = col_bluevec,pch=eco.sum.dat$Tank-18,las=1,ylab="Residuals")
abline(h=0)
title("Residuals for linear fit lfull.tfac \nwith Tank and Interaction")
#dev.off()





#ggsave('ltankregressionresid.pdf',width=6,height=6,dpi=500)
