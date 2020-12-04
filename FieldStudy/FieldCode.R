# This is the R code for Analysis of Field Data set (Chapter 5).

#Set your working directory to where you have "EcoFieldUp.csv" saved.


# Load libraries and packages.
library(tidyverse) # For creating tidy dataframes and piping
library(ggplot2) #For plotting
library(gdata) #For unrolling matrices into vectors
library(kableExtra) # For making Kable plots
library(MASS) #Model selection
library(magrittr) # For formatting
library(lme4) # Mixed Models
library(MuMIn) #Model Averaging
library(caret) 
library(leaps) #Best Subset Method
library(factoextra) #Principal Component Analysis


#Read in the original data

fieldData=read.csv("EcoFieldUp.csv") #The original file name is Ecofish run of river salmonids - Combined raw ct scores and fish data. Changed up EcoFieldUp to read in.

# We remove observations from Stream DDD in 2018 as they failed integritE tests.
fieldData$StreamYear=paste(fieldData$Stream.Code,fieldData$Year) #Create a column that combines Stream.Code and Year.

fieldData=fieldData%>%                    
  filter(StreamYear!='DDD 2018')

fieldData$SSRS=paste(fieldData$Stream.Code,fieldData$Sample,fieldData$Reach,fieldData$Site.number) #Creates a unique identifier for each set of tech replicates Stream/Sample/Reach/Site
#Create new columns to assist in coding and plotting.

fieldData$ID=fieldData$ï..ID #Create a column ID to replace the odd symbol name.

fieldData$Site.numberF=as.factor(fieldData$Site.number) #Create a site number as a factor column.

fieldData$StreamReach=paste(fieldData$Stream.Code,fieldData$Reach) #Combine Stream and Reach to get a new column.

fieldData$StreamReachYear=paste(fieldData$Stream.Code,fieldData$Reach,fieldData$Year)

#Create Volume Columns for each species and fish general.

fieldData$Transect.Flow.ms=fieldData$Transect.Flow.cms*0.01 # Meter per second isntead of cm/s

fieldData$CO.Biomass.g.m3=fieldData$CO.Biomass.g.m2/fieldData$Transect.Depth.m #Create Volume by accounting for Transect Depth.

fieldData$CT.Biomass.g.m3=fieldData$CT.Biomass.g.m2/fieldData$Transect.Depth.m #Create Volume by accounting for Transect Depth.

fieldData$RB.Biomass.g.m3=fieldData$RB.Biomass.g.m2/fieldData$Transect.Depth.m #Create Volume by accounting for Transect Depth.

fieldData$Fish.Biomass.g.m3=fieldData$Fish.Biomass.g.m2/fieldData$Transect.Depth.m #Create Volume by accounting for Transect Depth.


#We collapse over each set of eight techinical replicates.
jind<-rep(c(TRUE,rep(FALSE,7)),nrow(fieldData)/8) 

field.collapse=fieldData[jind, ]

field.collapse$MeanTCTCt=tapply(fieldData$Transformed.cl,factor(fieldData$SSRS),mean,na.rm=T)

field.collapse$MeanTCTCo=tapply(fieldData$Transformed.ki,factor(fieldData$SSRS),mean,na.rm=T)

field.collapse$MeanTCTRb=tapply(fieldData$Transformed.my,factor(fieldData$SSRS),mean,na.rm=T)

field.collapse$MeanTCTEf=tapply(fieldData$Transformed.ef,factor(fieldData$SSRS),mean,na.rm=T)

field.collapse$CO.Total.Biomass.g[37:54]=0 # If we fix mistake, no coho were caught at stream DDD.
fieldData[fieldData$Stream.Code=='DDD','CO.Total.Biomass.g']=0
fieldData[fieldData$Stream.Code=='DDD','CO.Biomass.g.m2']=0
fieldData[fieldData$Stream.Code=='DDD','CO.Biomass.g.m3']=0
field.collapse$CO.Biomass.g.m2[37:54]=0 # If we fix mistake, no coho were caught at stream DDD.
field.collapse$CO.Biomass.g.m3[37:54]=0 # If we fix mistake, no coho were caught at stream DDD.


field.removeef=field.collapse%>%
  filter(SSRS != 'DDD C Diversion 2') #Remove the observation of the low mean TCT Ef as it is an outlier.


# Create dataframes we will use for plotting.

reach.vec_ef=c(rep("Diversion",12),rep("Upstream",12))

fish.vec_ef=rep(c("Cutthroat","Coho Salmon","Rainbow Trout","Fish"),3) #Will be used for plotting

site.vector_ef=c(rep(1,4),rep(2,4),rep(3,4))

new.site.vector_ef=rep(c(rep(1,4),rep(2,4),rep(3,4)),2)

new.fish.vec_ef=rep(fish.vec_ef,2) # We have twice as many samples now (18 instead of 9).


# Stream AAA

fieldData.AAA=fieldData%>%
  filter(Stream.Code=='AAA')



Stream_AAA_tidy=fieldData.AAA %>%
  group_by(Site.number) %>%
  summarize(Cutthroat=unique(CT.Abundance),Coho=unique(CO.Abundance),Rainbow=unique(RB.Abundance), Efish=unique(Fish.Abundance),
            Cutthroat.b=unique(CT.Total.Biomass.g),Coho.b=unique(CO.Total.Biomass.g),Rainbow.b=unique(RB.Total.Biomass.g),Efish.b=unique(Fish.Total.Biomass.g),
            Cutthroat.l=unique(CT.Avg.Length.mm),Coho.l=unique(CO.Avg.Length.mm),Rainbow.l=unique(RB.Avg.Length.mm),Efish.l=unique(Fish.Avg.Length.mm))


# Unroll our tidy data frame into a vector for plotting, making sure to turn into numeric vectors.

ar= unmatrix(Stream_AAA_tidy[ ,2:5],byrow=T) # This accounts for the Abundance.
ar=as.vector(as.numeric(ar))
ar2= unmatrix(Stream_AAA_tidy[ ,6:9],byrow=T) #This accounts for the Biomass.
ar2=as.vector(as.numeric(ar2))
ar3= unmatrix(Stream_AAA_tidy[ ,10:13],byrow=T) #This accounts for the Average Length.
ar3=as.vector(as.numeric(ar3))

AAA.frame=data.frame(fish.vec_ef,site.vector_ef,ar,ar2,ar3) # Combine into a dataframe
names(AAA.frame)=c("Species","Site","SiteAbundance","Biomass","Length") #Assign names to data frame

AAA.frame$Site=as.factor(AAA.frame$Site) #Make sure site is as a factor.

fieldData.AAA=fieldData%>%
  filter(Stream.Code=='AAA')

g1=ggplot(data=fieldData.AAA)+
  geom_point(mapping=aes(x=jitter(Site.number,0.2),y=Transformed.cl,color=Sample ),shape=19)+
  ylim(0,20)+
  xlab("Site Number")+
  ylab("Transformed CT")+
  ggtitle("Transformed CT for Cutthroat Trout at Stream AAA")+
  scale_x_discrete(limits = c("1","2","3"))+
  theme(legend.position = "none")+
  theme_minimal()

g2=ggplot(data=fieldData.AAA)+
  geom_point(mapping=aes(x=jitter(Site.number,0.2),y=Transformed.ki,color=Sample))+
  ylim(0,20)+
  xlab("Site Number")+
  ylab("Transformed CT")+
  scale_x_discrete(limits = c("1","2","3"))+
  theme(legend.position = "none")+
  theme_minimal()
#ggsave('AAA_co_tct.pdf',g2,width=6,height=3,dpi=400)

g3= ggplot(data=fieldData.AAA)+
  geom_point(mapping=aes(x=jitter(Site.number,0.2),y=Transformed.my,color=Sample))+
  ylim(0,20)+
  xlab("Site Number")+
  ylab("Transformed CT")+
  ggtitle("Transformed CT for Rainbow Trout at Stream AAA")+
  scale_x_discrete(limits = c("1","2","3"))+
  theme(legend.position = "none")+
  theme_minimal()



g4= ggplot(data=fieldData.AAA)+
  geom_point(mapping=aes(x=jitter(Site.number,0.2),y=Transformed.ef,color=Sample))+
  ylim(0,20)+
  xlab("Site Number")+
  ylab("Transformed CT")+
  scale_x_discrete(limits = c("1","2","3"))+
  theme(legend.position = "none")+
  theme_minimal()

#ggsave('AAA_ef_tct.pdf',g4,width=6,height=3)


gg_AAA_Biomass=ggplot(data=AAA.frame,mapping=aes(x=Species,y=Biomass,fill=Site,label=Biomass))+
  geom_bar(stat='identity')+
  geom_text(size = 3, position = position_stack(vjust = 0.5))+
  xlab("Species")+
  ylab("Total Biomass (g)")+
  theme_minimal()

#ggsave('AAA_Ef_new.pdf',gg_AAA_Biomass,width=6,height=3)

gg_AAA_Coho=ggplot(AAA.frame[AAA.frame$Species=="Coho Salmon", ],mapping=aes(x=Site,y=Biomass,label=Biomass))+
  geom_bar(stat='identity')+
  ylim(0,100)+
  geom_text(size = 3, position = position_stack(vjust = 0.5))+
  xlab("Site Number")+
  ylab("Total Biomass (g)")+
  theme_minimal()
#ggsave('AAA_Co_new.pdf',gg_AAA_Coho,width=6,height=3)


gg_AAA_Ct=ggplot(AAA.frame[AAA.frame$Species=="Cutthroat", ],mapping=aes(x=Site,y=Biomass,label=Biomass,fill=2))+
  geom_bar(stat='identity')+
  geom_text(size = 3, position = position_stack(vjust =1.05))+
  ylim(0,600)+
  xlab("Site Number")+
  ylab("Total Biomass (g)")+
  ggtitle("Total Biomass (gram) of Cutthroat for each site at Stream AAA")+
  theme_minimal()+
  theme(legend.position = "none")


gg_AAA_Co=ggplot(AAA.frame[AAA.frame$Species=="Coho Salmon", ],mapping=aes(x=Site,y=Biomass,label=Biomass,fill=2))+
  geom_bar(stat='identity')+
  geom_text(size = 3, position = position_stack(vjust =1.05))+
  xlab("Site Number")+
  ylim(0,600)+
  ylab("Total Biomass (g)")+
  ggtitle("Total Biomass (gram) of Coho Salmon for each site at Stream AAA")+
  theme_minimal()+
  theme(legend.position = "none")



gg_AAA_Rb=ggplot(AAA.frame[AAA.frame$Species=="Rainbow Trout", ],mapping=aes(x=Site,y=Biomass,label=Biomass,fill=2))+
  geom_bar(stat='identity')+
  ylim(0,600)+
  geom_text(size = 3, position = position_stack(vjust =1.05))+
  xlab("Site Number")+
  ylab("Total Biomass (g)")+
  ggtitle("Total Biomass (gram) of Rainbow Trout for each site at Stream AAA")+
  theme_minimal()+
  theme(legend.position = "none")


gg_AAA_Ef=ggplot(AAA.frame[AAA.frame$Species=="Fish", ],mapping=aes(x=Site,y=Biomass,label=Biomass,fill=2))+
  geom_bar(stat='identity')+
  geom_text(size = 3, position = position_stack(vjust =1.05))+
  ylim(0,600)+
  xlab("Site Number")+
  ylab("Total Biomass (g)")+
  theme_minimal()+
  theme(legend.position = "none")
ggsave('AAA_Ef_new.pdf',gg_AAA_Ef,width=6,height=3)

# Stream BBB

fieldData.BBB=fieldData%>%
  filter(Stream.Code=='BBB')


Stream_BBB_tidy=fieldData.BBB %>%
  group_by(Site.number) %>%
  summarize(Cutthroat=unique(CT.Abundance),Coho=unique(CO.Abundance),Rainbow=unique(RB.Abundance),Efish=unique(Fish.Abundance),
            Cutthroat.b=unique(CT.Total.Biomass.g),Coho.b=unique(CO.Total.Biomass.g),Rainbow.b=unique(RB.Total.Biomass.g),Efish.b=unique(Fish.Total.Biomass.g),
            Cutthroat.l=unique(CT.Avg.Length.mm),Coho.l=unique(CO.Avg.Length.mm),Rainbow.l=unique(RB.Avg.Length.mm),Efish.l=unique(Fish.Avg.Length.mm))


# Unroll our tidy data frame into a vector for plotting, making sure to turn into numeric vectors.
ar= unmatrix(Stream_BBB_tidy[ ,2:5],byrow=T) # This accounts for the Abundance.
ar=as.vector(as.numeric(ar))
ar2= unmatrix(Stream_BBB_tidy[ ,6:9],byrow=T) #This accounts for the Biomass.
ar2=as.vector(as.numeric(ar2))
ar3= unmatrix(Stream_BBB_tidy[ ,10:13],byrow=T) #This accounts for the Average Length.
ar3=as.vector(as.numeric(ar3))

BBB.frame=data.frame(fish.vec_ef,site.vector_ef,ar,ar2,ar3) # Combine into a dataframe
names(BBB.frame)=c("Species","Site","SiteAbundance","Biomass","Length") #Assign names to data frame

BBB.frame$Site=as.factor(BBB.frame$Site) #Make sure site is as a factor.




g1_BBB=ggplot(data=fieldData.BBB)+
  geom_point(mapping=aes(x=jitter(Site.number,0.2),y=Transformed.cl,color=Sample))+
  ylim(0,20)+
  xlab("Site Number")+
  ylab("Transformed CT")+
  ggtitle("Transformed CT for Cutthroat Trout at Stream BBB")+
  scale_x_discrete(limits = c("1","2","3"))+
  theme(legend.position = "none")+
  theme_minimal()


g2_BBB=ggplot(data=fieldData.BBB)+
  geom_point(mapping=aes(x=jitter(Site.number,0.2),y=Transformed.ki,color=Sample))+
  ylim(0,20)+
  xlab("Site Number")+
  ylab("Transformed CT")+
  scale_x_discrete(limits = c("1","2","3"))+
  theme(legend.position = "none")+
  theme_minimal()


#ggsave('BBB_co_tct.pdf',g2_BBB,width=6,height=3,dpi=400)

g3_BBB= ggplot(data=fieldData.BBB)+
  geom_point(mapping=aes(x=jitter(Site.number,0.2),y=Transformed.my,color=Sample))+
  ylim(0,20)+
  xlab("Site Number")+
  ylab("Transformed CT")+
  ggtitle("Transformed CT for Rainbow Trout at Stream BBB")+
  scale_x_discrete(limits = c("1","2","3"))+
  theme(legend.position = "none")+
  theme_minimal()


g4_BBB= ggplot(data=fieldData.BBB)+
  geom_point(mapping=aes(x=jitter(Site.number,0.2),y=Transformed.ef,color=Sample))+
  ylim(0,20)+
  xlab("Site Number")+
  ylab("Transformed CT")+
  scale_x_discrete(limits = c("1","2","3"))+
  theme(legend.position = "none")+
  theme_minimal()

#ggsave('BBB_ef_tct.pdf',g4_BBB,width=6,height=3,dpi=400)



gg_BBB_Biomass=ggplot(data=BBB.frame,mapping=aes(x=Species,y=Biomass,fill=Site,label=Biomass))+
  geom_bar(stat='identity')+
  geom_text(size = 2, position = position_stack(vjust = 0.5))+
  ylim(0,1000)+
  xlab("Species")+
  ylab("Total Biomass (g)")+
  theme_minimal()



gg_BBB_Ct=ggplot(BBB.frame[BBB.frame$Species=="Cutthroat", ],mapping=aes(x=Site,y=Biomass,label=Biomass,fill=2))+
  geom_bar(stat='identity')+
  geom_text(size = 3, position = position_stack(vjust =1.2))+
  ylim(0,400)+
  xlab("Site Number")+
  ylab("Total Biomass (g)")+
  ggtitle("Total Biomass (gram) of Cutthroat for each site at Stream BBB")+
  theme_minimal()+
  theme(legend.position = "none")


gg_BBB_Co=ggplot(BBB.frame[BBB.frame$Species=="Coho Salmon", ],mapping=aes(x=Site,y=Biomass,label=Biomass,fill=2))+
  geom_bar(stat='identity')+
  geom_text(size = 3, position = position_stack(vjust =1.2))+
  xlab("Site Number")+
  ylim(0,100)+
  ylab("Total Biomass (g)")+
  theme_minimal()+
  theme(legend.position = "none")

#ggsave('BBB_Co_new.pdf',gg_BBB_Co,width=6,height=3,dpi=400)


gg_BBB_Rb=ggplot(BBB.frame[BBB.frame$Species=="Rainbow Trout", ],mapping=aes(x=Site,y=Biomass,label=Biomass,fill=2))+
  geom_bar(stat='identity')+
  ylim(0,400)+
  geom_text(size = 3, position = position_stack(vjust =1.2))+
  xlab("Site Number")+
  ylab("Total Biomass (g)")+
  ggtitle("Total Biomass (gram) of Rainbow Trout for each site at Stream BBB")+
  theme_minimal()+
  theme(legend.position = "none")


gg_BBB_Ef=ggplot(BBB.frame[BBB.frame$Species=="Fish", ],mapping=aes(x=Site,y=Biomass,label=Biomass,fill=2))+
  geom_bar(stat='identity')+
  geom_text(size = 3, position = position_stack(vjust =1.05))+
  ylim(0,400)+
  xlab("Site Number")+
  ylab("Total Biomass (g)")+
  theme_minimal()+
  theme(legend.position = "none")

#ggsave('BBB_Ef_new.pdf',gg_BBB_Ef,width=6,height=3,dpi=400)

# Stream CCC

fieldData.CCC=fieldData%>%
  filter(Stream.Code=='CCC')

Stream_CCC_tidy=fieldData.CCC %>%
  group_by(Site.number,Reach) %>%
  summarize(Cutthroat=unique(CT.Abundance),Coho=unique(CO.Abundance),Rainbow=unique(RB.Abundance),Efish=unique(Fish.Abundance),
            Cutthroat.b=unique(CT.Total.Biomass.g),Coho.b=unique(CO.Total.Biomass.g),Rainbow.b=unique(RB.Total.Biomass.g),Efish.b=unique(Fish.Total.Biomass.g),
            Cutthroat.l=unique(CT.Avg.Length.mm),Coho.l=unique(CO.Avg.Length.mm),Rainbow.l=unique(RB.Avg.Length.mm),Efish.l=unique(Fish.Avg.Length.mm))


# Unroll our tidy data frame into a vector for plotting, making sure to turn into numeric vectors.
ar= unmatrix(Stream_CCC_tidy[ ,3:6],byrow=T) # This accounts for the Abundance.
ar=as.vector(as.numeric(ar))
ar2= unmatrix(Stream_CCC_tidy[ ,6:9],byrow=T) #This accounts for the Biomass.
ar2=as.vector(as.numeric(ar2))
ar3= unmatrix(Stream_CCC_tidy[ ,9:12],byrow=T) #This accounts for the Average Length.
ar3=as.vector(as.numeric(ar3))

CCC.frame=data.frame(new.fish.vec_ef,reach.vec_ef,new.site.vector_ef,ar,ar2,ar3) # Combine into a dataframe
names(CCC.frame)=c("Species","Reach","Site","SiteAbundance","Biomass","Length") #Assign names to data frame


CCC.frame$Site=as.factor(CCC.frame$Site) #Make sure site is as a factor.

g1_CCC=ggplot(data=fieldData.CCC)+
  geom_point(mapping=aes(x=jitter(Site.number,0.2),y=Transformed.cl,color=Sample))+
  facet_wrap(~Reach,ncol=2)+
  ylim(0,23)+
  xlab("Site Number")+
  ylab("Transformed CT")+
  ggtitle("Transformed CT for Cutthroat Trout at Stream CCC")+
  scale_x_discrete(limits = c("1","2","3"))+
  theme(legend.position = "none")+
  theme(axis.text.y = element_text( size = 8 ))+
  theme_minimal()


g2_CCC=ggplot(data=fieldData.CCC)+
  geom_point(mapping=aes(x=jitter(Site.number,0.2),y=Transformed.ki,color=Sample))+
  facet_wrap(~Reach,ncol=2)+
  ylim(0,23)+
  xlab("Site Number")+
  ylab("Transformed CT")+
  scale_x_discrete(limits = c("1","2","3"))+
  theme(legend.position = "none")+
  theme(axis.text.y = element_text( size = 8 ))+
  theme_minimal()

#ggsave('CCC_co_tct.pdf',g2_CCC,width=6,height=3,dpi=400)

g3_CCC=ggplot(data=fieldData.CCC)+
  geom_point(mapping=aes(x=jitter(Site.number,0.2),y=Transformed.my,color=Sample))+
  facet_wrap(~Reach,ncol=2)+
  ylim(0,23)+
  xlab("Site Number")+
  ylab("Transformed CT")+
  ggtitle("Transformed CT for Rainbow Trout at Stream CCC")+
  scale_x_discrete(limits = c("1","2","3"))+
  theme(legend.position = "none")+
  theme(axis.text.y = element_text( size = 8 ))+
  theme_minimal()


g4_CCC= ggplot(data=fieldData.CCC)+
  geom_point(mapping=aes(x=jitter(Site.number,0.2),y=Transformed.ef,color=Sample))+
  facet_wrap(~Reach,ncol=2)+
  xlab("Site Number")+
  ylab("Transformed CT")+
  ylim(0,23)+
  scale_x_discrete(limits = c("1","2","3"))+
  theme(legend.position = "none")+
  theme(axis.text.y = element_text( size = 8 ))+
  theme_minimal()

#ggsave('CCC_ef_tct.pdf',g4_CCC,width=6,height=3,dpi=400)

gg_CCC_Biomass=ggplot(data=CCC.frame,mapping=aes(x=Species,y=Biomass,fill=Site,label=Biomass))+
  geom_bar(stat='identity')+
  facet_wrap(~Reach)+
  geom_text(size = 2, position = position_stack(vjust = 0.6))+
  ylim(0,1000)+
  xlab("Species")+
  ylab("Total Biomass (g)")+
  ggtitle("Total Biomass (gram) of Fish for each site at Stream CCC")+
  theme(axis.text.x = element_text(size=1))+
  theme_minimal()+theme(axis.text.x = element_text(angle = 50, hjust = 1))


gg_CCC_Ct=ggplot(CCC.frame[CCC.frame$Species=="Cutthroat", ],mapping=aes(x=Site,y=Biomass,label=Biomass,fill=2))+
  geom_bar(stat='identity')+
  facet_wrap(~Reach,ncol=2)+
  geom_text(size = 3, position = position_stack(vjust =1.5))+
  ylim(0,400)+
  xlab("Site Number")+
  ylab("Total Biomass (g)")+
  ggtitle("Total Biomass (gram) of Cutthroat for each site at Stream CCC")+
  theme_minimal()+
  theme(legend.position = "none")


gg_CCC_Co=ggplot(CCC.frame[CCC.frame$Species=="Coho Salmon", ],mapping=aes(x=Site,y=Biomass,label=Biomass,fill=2))+
  facet_wrap(~Reach,ncol=2)+
  geom_bar(stat='identity')+
  geom_text(size = 3, position = position_stack(vjust =1.2))+
  xlab("Site Number")+
  ylim(0,400)+
  ylab("Total Biomass (g)")+
  theme_minimal()+
  theme(legend.position = "none")


#ggsave('CCC_Co_new.pdf',gg_CCC_Co,width=6,height=3,dpi=400)

gg_CCC_Rb=ggplot(CCC.frame[CCC.frame$Species=="Rainbow Trout", ],mapping=aes(x=Site,y=Biomass,label=Biomass,fill=2))+
  geom_bar(stat='identity')+
  ylim(0,400)+
  facet_wrap(~Reach,ncol=2)+
  geom_text(size = 3, position = position_stack(vjust =1.2))+
  xlab("Site Number")+
  ylab("Total Biomass (g)")+
  ggtitle("Total Biomass (gram) of Rainbow Trout for each site at Stream CCC")+
  theme_minimal()+
  theme(legend.position = "none")


gg_CCC_Ef=ggplot(CCC.frame[CCC.frame$Species=="Fish", ],mapping=aes(x=Site,y=Biomass,label=Biomass,fill=2))+
  geom_bar(stat='identity')+
  geom_text(size = 3, position = position_stack(vjust =1.1))+
  ylim(0,600)+
  xlab("Site Number")+
  ylab("Total Biomass (g)")+
  facet_wrap(~Reach,ncol=2)+
  theme_minimal()+
  theme(legend.position = "none")

#ggsave('CCC_Ef_new.pdf',gg_CCC_Ef,width=6,height=3,dpi=400)



# Stream DDD

fieldData.DDD=fieldData%>%
  filter(Stream.Code=='DDD')


#Fix Excel mistake. Zero Coho were found at Site DDD!
fieldData.DDD$CO.Abundance=0
fieldData.DDD$CO.Total.Biomass.g=0
fieldData.DDD$CO.Avg.Length.mm=0


reach.vec_ef_DDD=c(rep("Diversion",12),rep("Downstream",12))

Stream.DDD_tidy<-fieldData.DDD%>%
  group_by(Site.number,Reach)%>%
  summarize(Cutthroat=unique(CT.Abundance),Coho=unique(CO.Abundance),Rainbow=unique(RB.Abundance),EFish=unique(Fish.Abundance),
            Cutthroat.b=unique(CT.Total.Biomass.g),Coho.b=unique(CO.Total.Biomass.g),Rainbow.b=unique(RB.Total.Biomass.g),EFish.b=unique(Fish.Total.Biomass.g),
            Cutthroat.l=unique(CT.Avg.Length.mm),Coho.l=unique(CO.Avg.Length.mm),Rainbow.l=unique(RB.Avg.Length.mm),EFish.l=unique(Fish.Avg.Length.mm))%>%
  arrange(Reach)


ar= unmatrix(Stream.DDD_tidy[ ,3:6],byrow=T)
ar=as.vector(as.numeric(ar))
ar2= unmatrix(Stream.DDD_tidy[ ,6:9],byrow=T)
ar2=as.vector(as.numeric(ar2))
ar3= unmatrix(Stream.DDD_tidy[ ,9:12],byrow=T)
ar3=as.vector(as.numeric(ar3))

DDD_df=data.frame(new.fish.vec_ef,reach.vec_ef_DDD,new.site.vector_ef,ar,ar2,ar3)
names(DDD_df)=c("Species","Reach","Site","Abundance","Biomass","Length")


#DDD_df$Site[DDD_df$Reach=='Downstream']=DDD_df$Site+3 # The three sites Downstream are NOT the same as the three sites at the Diversion. So we scale so we can differentiate. 

DDD_df$Site=as.factor(DDD_df$Site) #Make sure site is as a factor.


DDD.frame=DDD_df

g1_DDD=ggplot(data=fieldData.DDD)+
  geom_point(mapping=aes(x=jitter(Site.number,0.2),y=Transformed.cl,color=Sample))+
  facet_wrap(~Reach,ncol=2)+
  ylim(0,20)+
  xlab("Site Number")+
  ylab("Transformed CT")+
  ggtitle("Transformed CT for Cutthroat Trout at Stream DDD")+
  scale_x_discrete(limits = c("1","2","3"))+
  theme(legend.position = "none")+
  theme(axis.text.y = element_text( size = 8 ))+
  theme_minimal()


g2_DDD=ggplot(data=fieldData.DDD)+
  geom_point(mapping=aes(x=jitter(Site.number,0.2),y=Transformed.ki,color=Sample))+
  facet_wrap(~Reach,ncol=2)+
  ylim(0,20)+
  xlab("Site Number")+
  ylab("Transformed CT")+
  scale_x_discrete(limits = c("1","2","3"))+
  theme(legend.position = "none")+
  theme(axis.text.y = element_text( size = 8 ))+
  theme_minimal()


#ggsave('DDD_co_tct.pdf',g2_DDD,width=6,height=3,dpi=400)


g3_DDD=ggplot(data=fieldData.DDD)+
  geom_point(mapping=aes(x=jitter(Site.number,0.2),y=Transformed.my,color=Sample))+
  facet_wrap(~Reach,ncol=2)+
  ylim(0,20)+
  xlab("Site Number")+
  ylab("Transformed CT")+
  ggtitle("Transformed CT for Rainbow Trout at Stream DDD")+
  scale_x_discrete(limits = c("1","2","3"))+
  theme(legend.position = "none")+
  theme(axis.text.y = element_text( size = 8 ))+
  theme_minimal()


g4_DDD= ggplot(data=fieldData.DDD)+
  geom_point(mapping=aes(x=jitter(Site.number,0.2),y=Transformed.ef,color=Sample))+
  facet_wrap(~Reach,ncol=2)+
  xlab("Site Number")+
  ylab("Transformed CT")+
  ylim(0,20)+
  scale_x_discrete(limits = c("1","2","3"))+
  theme(legend.position = "none")+
  theme(axis.text.y = element_text( size = 8 ))+
  theme_minimal()



#ggsave('DDD_ef_tct.pdf',g4_DDD,width=6,height=3,dpi=400)



gdd_Biomass=ggplot(data=DDD_df,mapping=aes(x=Species,y=Biomass,fill=Site,label=Biomass))+
  geom_bar(stat='identity')+
  facet_grid(~Reach)+
  xlab("Species")+
  ylab("Total Biomass (g)")+
  ylim(0,1000)+
  geom_text(size = 2, position = position_stack(vjust = 0.5))+
  ggtitle("Total Biomass (gram) of Fish for Stream DDD")+
  theme(axis.text.x = element_text(size=6))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 50, hjust = 1))



gg_DDD_Ct=ggplot(DDD.frame[DDD.frame$Species=="Cutthroat", ],mapping=aes(x=Site,y=Biomass,label=Biomass,fill=2))+
  geom_bar(stat='identity')+
  facet_wrap(~Reach,ncol=2)+
  geom_text(size = 3, position = position_stack(vjust =1.5))+
  ylim(0,400)+
  xlab("Site Number")+
  ylab("Total Biomass (g)")+
  ggtitle("Total Biomass (gram) of Cutthroat for each site at Stream DDD")+
  theme_minimal()+
  theme(legend.position = "none")


gg_DDD_Co=ggplot(DDD.frame[DDD.frame$Species=="Coho Salmon", ],mapping=aes(x=Site,y=Biomass,label=Biomass,fill=2))+
  facet_wrap(~Reach,ncol=2)+
  geom_bar(stat='identity')+
  geom_text(size = 3, position = position_stack(vjust =1.2))+
  xlab("Site Number")+
  ylim(0,400)+
  ylab("Total Biomass (g)")+
  theme_minimal()+
  theme(legend.position = "none")

#ggsave('DDD_Co_new.pdf',gg_DDD_Co,width=6,height=3,dpi=500)

gg_DDD_Rb=ggplot(DDD.frame[DDD.frame$Species=="Rainbow Trout", ],mapping=aes(x=Site,y=Biomass,label=Biomass,fill=2))+
  geom_bar(stat='identity')+
  ylim(0,400)+
  facet_wrap(~Reach,ncol=2)+
  geom_text(size = 3, position = position_stack(vjust =1.2))+
  xlab("Site Number")+
  ylab("Total Biomass (g)")+
  ggtitle("Total Biomass (gram) of Rainbow Trout for each site at Stream DDD")+
  theme_minimal()+
  theme(legend.position = "none")


gg_DDD_Ef=ggplot(DDD.frame[DDD.frame$Species=="Fish", ],mapping=aes(x=Site,y=Biomass,label=Biomass,fill=2))+
  geom_bar(stat='identity')+
  geom_text(size = 3, position = position_stack(vjust =1.1))+
  ylim(0,600)+
  xlab("Site Number")+
  ylab("Total Biomass (g)")+
  facet_wrap(~Reach,ncol=2)+
  theme_minimal()+
  theme(legend.position = "none")


#ggsave('DDD_Ef_new.pdf',gg_DDD_Ef,height=3,width=6,dpi=400)



#We collapse over each set of eight techinical replicates.
jind<-rep(c(TRUE,rep(FALSE,7)),nrow(fieldData)/8) 

# Only keep the first from each set.
field.collapse=fieldData[jind, ]

field.collapse$MeanTCTCt=tapply(fieldData$Transformed.cl,factor(fieldData$SSRS),mean,na.rm=T)

field.collapse$MeanTCTCo=tapply(fieldData$Transformed.ki,factor(fieldData$SSRS),mean,na.rm=T)

field.collapse$MeanTCTRb=tapply(fieldData$Transformed.my,factor(fieldData$SSRS),mean,na.rm=T)

field.collapse$MeanTCTEf=tapply(fieldData$Transformed.ef,factor(fieldData$SSRS),mean,na.rm=T)

field.collapse$CO.Total.Biomass.g[37:54]=0 # If we fix mistake, no coho were caught at stream DDD.
field.collapse$CO.Biomass.g.m2[37:54]=0 # If we fix mistake, no coho were caught at stream DDD.
field.collapse$CO.Biomass.g.m3[37:54]=0 # If we fix mistake, no coho were caught at stream DDD.


Stream=field.collapse$StreamReach

field.collapse$MedTCTCt=tapply(fieldData$Transformed.cl,factor(fieldData$SSRS),median,na.rm=T)

field.collapse$MedTCTCo=tapply(fieldData$Transformed.ki,factor(fieldData$SSRS),median,na.rm=T)

field.collapse$MedTCTRb=tapply(fieldData$Transformed.my,factor(fieldData$SSRS),median,na.rm=T)

field.collapse$MedTCTEf=tapply(fieldData$Transformed.ef,factor(fieldData$SSRS),median,na.rm=T)


#Simple Linear Models Biomass vs Mean TCT

model_co=lm(MeanTCTCo~CO.Total.Biomass.g, data=field.collapse) 
#R^2 of 0.600 and adjusted R^2 of 0.592


model_EF3=lm(MeanTCTEf~Fish.Total.Biomass.g, data=field.collapse) 


model_ct=lm(MeanTCTCt~CT.Total.Biomass.g, data=field.collapse) 
#R^2 of 0.716 and adjusted R^2 of 0.711

model_rb=lm(MeanTCTRb~RB.Total.Biomass.g, data=field.collapse) 
#R^2 of 0.108 and adjusted R^2 of 0.091. highly depedent on stream.


model_co.m2=lm(MeanTCTCo~CO.Biomass.g.m2, data=field.collapse) 
#R^2 of 0.600 and adjusted R^2 of 0.592

model_ct.m2=lm(MeanTCTCt~CT.Biomass.g.m2, data=field.collapse) 
#R^2 of 0.716 and adjusted R^2 of 0.711

model_rb.m2=lm(MeanTCTRb~RB.Biomass.g.m2, data=field.collapse) 
#R^2 of 0.108 and adjusted R^2 of 0.091. highly depedent on stream.




model_co.m3=lm(MeanTCTCo~CO.Biomass.g.m3, data=field.collapse) 
#R^2 of 0.600 and adjusted R^2 of 0.592

model_ct.m3=lm(MeanTCTCt~CT.Biomass.g.m3, data=field.collapse) 
#R^2 of 0.716 and adjusted R^2 of 0.711

model_rb.m3=lm(MeanTCTRb~RB.Biomass.g.m3, data=field.collapse) 
#R^2 of 0.108 and adjusted R^2 of 0.091. highly depedent on stream.


model_co_med=lm(MedTCTCo~CO.Total.Biomass.g, data=field.collapse) 
#R^2 of 0.7 and adjusted R^2 of 0.694

model_ct_med=lm(MedTCTCt~CT.Total.Biomass.g, data=field.collapse) 
#R^2 of 0.701 and adjusted R^2 of 0.695

model_rb_med=lm(MedTCTRb~RB.Total.Biomass.g, data=field.collapse) 
#R^2 of 0.0986 and adjusted R^2 of 0.081. highly depedent on stream.




# GGplots

gg_CO2=ggplot(data=field.collapse)+
  geom_point(mapping=aes(x=CO.Total.Biomass.g,y=MeanTCTCo,color=Stream,shape=Stream))+
  geom_smooth(mapping=aes(x=CO.Total.Biomass.g,y=MeanTCTCo),method="lm",color='black')+
  xlab("Total Biomass (grams)")+
  ylab("Mean TCT")+
  ylim(0,17)+
  scale_shape_manual(values = 0:7)+
  theme_classic()



#ggsave('model_CO.pdf',gg_CO2,height=6,width=6,dpi=500)



gg_CT2=ggplot(data=field.collapse)+
  geom_point(mapping=aes(x=CT.Total.Biomass.g,y=MeanTCTCt,color=Stream,shape=Stream))+
  geom_smooth(mapping=aes(x=CT.Total.Biomass.g,y=MeanTCTCt),method="lm",color='black')+
  xlab("Total Biomass (grams)")+
  ylab("Mean TCT")+
  ylim(0,17)+
  ggtitle("Mean TCT versus Total Biomass for Cutthroat")+
  scale_shape_manual(values = 0:7)+
  theme_classic()



gg_RB2=ggplot(data=field.collapse)+
  geom_point(mapping=aes(x=RB.Total.Biomass.g,y=MeanTCTRb,color=Stream,shape=Stream))+
  geom_smooth(mapping=aes(x=RB.Total.Biomass.g,y=MeanTCTRb),method="lm",color='black')+
  xlab("Total Biomass (grams)")+
  ylab("Mean TCT")+
  ylim(0,17)+
  ggtitle("Mean TCT versus Total Biomass for Rainbow Trout")+
  scale_shape_manual(values = 0:7)+
  theme_classic()






gg_EF_Outlier=ggplot(data=field.collapse)+
  geom_point(mapping=aes(x=Fish.Total.Biomass.g,y=MeanTCTEf,color=Stream,shape=Stream))+
  geom_smooth(mapping=aes(x=Fish.Total.Biomass.g,y=MeanTCTEf),method="lm",color='black')+
  xlab("Total Biomass (grams)")+
  ylab("Mean TCT")+
  scale_shape_manual(values = 0:7)+
  theme_classic()

#ggsave('gg_ef.pdf',gg_EF_Outlier,width=6,height=6,dpi=500)


field.removeef=field.collapse%>%
  filter(SSRS != 'DDD C Diversion 2') #Remove the observation of the low mean TCT Ef as it is an outlier.


model_coho=lm(MeanTCTCo~CO.Total.Biomass.g+Transect.Flow.cms+CO.Total.Biomass.g*Transect.Flow.cms+Water.Temperature.C+pH,data=field.collapse)



model_coho_full=lm(MeanTCTCo~CO.Total.Biomass.g+Transect.Flow.cms+CO.Total.Biomass.g*Transect.Flow.cms+Water.Temperature.C+pH+Wetted.Width.m+Site.Area.m2+eDNA.Distance.from.Shore.m+eDNA.Total.Water.Depth.m,data=field.collapse)


model_coho_fin=lm(MeanTCTCo~CO.Biomass.g.m3+Transect.Flow.cms+CO.Biomass.g.m3*Transect.Flow.cms+
                    Water.Temperature.C+pH+eDNA.Distance.from.Shore.m+eDNA.Total.Water.Depth.m,
                  data=field.collapse)




redmodel=lm(MeanTCTCo~pH+CO.Biomass.g.m3+Transect.Flow.cms+CO.Biomass.g.m3*Transect.Flow.cms,data=field.collapse)




# Set options needed for MuMIn package.
options(na.action = "na.fail")
options(digits=3)
ma_coho <- MuMIn::dredge(model_coho_fin)
mma=model.avg(ma_coho, subset = delta < 4)
confset.95p <- get.models(ma_coho, cumsum(weight) <= .95)
avgmod.95p <- model.avg(confset.95p)

# Create subset models
models_sub <- regsubsets(MeanTCTCo~CO.Biomass.g.m3+Transect.Flow.cms+CO.Biomass.g.m3*Transect.Flow.cms+Water.Temperature.C+pH+eDNA.Distance.from.Shore.m+eDNA.Total.Water.Depth.m, data = field.collapse, nvmax = 7)


res.sum <- summary(models_sub)
data.frame(
  Adj.R2 = which.max(res.sum$adjr2),
  CP = which.min(res.sum$cp),
  BIC = which.min(res.sum$bic)
)

compare_models=data.frame(1:7,round(res.sum$adjr2,2),round(res.sum$cp,2),round(res.sum$bic,2))

names(compare_models)=c("p","AdjR^{2}","Cp","BIC")




k_subsets=kable(compare_models,booktabs=T,format="latex",escape = FALSE)%>%
  kable_styling(latex_options = 'striped')%>%
  add_header_above(c("p"=1,"Adj$R^{2}$"=1,"Cp"=1,"BIC"=1),escape=FALSE)%>%
  row_spec(row=3,col='red')




pred_frame=c("Biomass","Flow","Bio*Flow","Temp","pH","Shore","Depth")

best1=c("*","","",""," ",""," ")
best2=c("*","","*","","","","")
best3=c("*","","*","","*","","")
best4=c("*","","*","*","*","","")
best5=c("*","","*","*","","*","")
best6=c("*","*","*","*","*","*","")

best=rbind(best1,best2,best3,best4,best5,best6)
names(best)=1:7



k_best=kable(best,caption="Best Subset Method",booktabs=T,format="latex")%>%
  kable_styling(latex_options = 'striped')%>%
  add_header_above(c("p"=1,"Biomass"=1,"Flow"=1,"Biomass*Flow"=1,"Temp"=1,"pH"=1,"Shore"=1,"Depth"=1))%>%
  row_spec(row=3,col='red')


# Make a dataframe summarizing key covariates.

temp_fieldall2=fieldData %>% # This allows us to view the mean temperature at each stream for a given year.
  group_by(Stream.Code,Year,Site.number,Reach) %>%
  summarize(temp=round(mean(Water.Temperature.C),3),ph=round(mean(pH),3),flow=round(mean(Transect.Flow.cms),3),depth=round(mean(eDNA.Total.Water.Depth.m),3),shore=round(mean(eDNA.Distance.from.Shore.m),3))


temp_fieldall2=data.frame(temp_fieldall2)

names(temp_fieldall2)=c("Stream","Year","Site","Reach","Temperature","pH","Flow (cm/s)","Depth","Shore")





# Remove Identifying variables such as Stream, Site , Reach and Year

temp_fieldall_pca=temp_fieldall2 %>% dplyr::select(-Stream,-Site,-Reach,-Year)

p=paste(temp_fieldall2$Stream,temp_fieldall2$Reach) # Get stream and reach

ll=as.numeric(as.factor(paste(temp_fieldall2$Stream,temp_fieldall2$Reach)))

streams=c('AA1','AA2','AA3','BB1','BB2','BB3','CC1','CC2','CC3','CC4','CC5','CC6','DD1','DD2','DD3','DD4','DD5','DD6')

row.names(temp_fieldall_pca)=streams

# Create PCA
res.pca_new<-prcomp(temp_fieldall_pca,scale=TRUE)



# Get coefficients
res.pca_new


#Scree Plot
f1=fviz_eig(res.pca_new,addlabels = TRUE,barfill='blueviolet')


#Contributions
f2=fviz_pca_var(res.pca_new,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("blue", "purple", "red"),
             
             repel = TRUE     # Avoid text overlapping
)


# Plot by first two components by stream/reach.

f3=fviz_pca_ind(res.pca_new,
             col.ind = p, # color by groups
             palette = c("red","blue","green","purple","black","orange"),
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Stream",
             repel = TRUE,geom.ind = "point")



# Create custom table
code1= c('CO.Total.Biomass.g','CO.Biomass.g.m2','CO.Biomass.g.m3','eDNA.Distance.from.Shore.m','eDNA.Total.Water.Depth.m','pH','Transect.Flow.cms',
        'Water.Temperature.C','CO.Biomass.g.m3:Transect.Flow.cms')

exp=c("Total Coho biomass (gram) caught in the transect","Coho biomass per metre squared of transect (g/m2)","Coho biomass per metre cubed of transect (g/m3)",
      "Distance from shore in meters for which the sample was taken","Water depth in meters for which the sample was taken",
      "the pH value","The rate of flow of the water in cm/s","The water temperature in Celcius",'The interaction between Flow and Biomass per meter cubed')

dz=data.frame(code1,exp)
names(dz)=c("Term","Explanation")

ktz2=kable(dz,format="latex",booktabs=T,escape = FALSE)%>%
  kable_styling("striped")



code= c('CO.Biomass.g.m3: 1','eDNA.Distance.from.Shore.m: 2','eDNA.Total.Water.Depth.m: 3','pH: 4','Transect.Flow.cms: 5',
'Water.Temperature.C: 6','CO.Biomass.g.m3:Transect.Flow.cms: 7')

df=data.frame(code)
names(df)='Terms'

ktz2=kable(df,format="latex",booktabs=T,escape = FALSE)%>%
  kable_styling("striped")



temp_fieldall2=fieldData %>% # This allows us to view the mean temperature at each stream for a given year.
  group_by(Stream.Code,Year,Site.number,Reach) %>%
  summarize(temp=round(mean(Water.Temperature.C),3),ph=round(mean(pH),3),flow=round(mean(Transect.Flow.cms),3),depth=round(mean(eDNA.Total.Water.Depth.m),3),shore=round(mean(eDNA.Distance.from.Shore.m),3))              


temp_fieldall2=data.frame(temp_fieldall2)
names(temp_fieldall2)=c("Stream","Year","Site","Reach","Temperature","pH","Flow (cm/s)","Depth","Shore")

k_new2=kable(temp_fieldall2,booktabs=T,format="latex")%>%
  kable_styling(latex_options = 'striped')%>%
  add_header_above(c("Stream"=1,"Year"=1,"Site"=1,"Reach"=1,"Temperature "=1,"pH"=1,"Flow"=1,"Depth","Shore Distance"))

#pdf(file="updatedkable2.pdf",width=6,height=6)
#temp_fieldall2
#dev.off()



field.check=field.collapse%>%
  filter(!(Stream.Code=='DDD' & Reach=='Diversion' & Site.number==3 & Sample.replicate=='C'))

fieldData.check=fieldData%>%
  filter(!(Stream.Code=='DDD' & Reach=='Diversion' & Site.number==3 & Sample=='C'))


field.collapse$MeanTCTEf=tapply(fieldData$Transformed.ef,factor(fieldData$SSRS),mean,na.rm=T)


field.removeef=field.collapse%>%
  filter(SSRS != 'DDD C Diversion 2') #Remove the observation of the low mean TCT Ef as it is an outlier.

field.pairsubset_co=field.collapse%>%
  dplyr::select(MeanTCTCo,CO.Total.Biomass.g,CO.Biomass.g.m2,CO.Biomass.g.m3,pH,Water.Temperature.C,Transect.Flow.cms,Stream.Code)


field.pairsubset_ef=field.removeef%>%
  dplyr::select(MeanTCTEf,Fish.Total.Biomass.g,Fish.Biomass.g.m2,Fish.Biomass.g.m3,pH,Water.Temperature.C,Transect.Flow.cms,Stream.Code)


field.pairsubset_ct=field.collapse%>%
  dplyr::select(MeanTCTCt,CT.Total.Biomass.g,CT.Biomass.g.m2,CT.Biomass.g.m3,pH,Water.Temperature.C,Transect.Flow.cms,Stream.Code)


field.pairsubset_rb=field.collapse%>%
  dplyr::select(MeanTCTRb,RB.Total.Biomass.g,RB.Biomass.g.m2,RB.Biomass.g.m3,pH,Water.Temperature.C,Transect.Flow.cms,Stream.Code)
#pdf(file="Cohopairs.pdf",width=7,height=7)

pairs(field.collapse[, c('MeanTCTCo','CO.Total.Biomass.g','CO.Biomass.g.m2','CO.Biomass.g.m3','Water.Temperature.C','pH','Transect.Flow.cms')], col=factor(field.collapse$Stream.Code),cex=0.8,cex.labels=0.5)

#dev.off()

#pdf(file="EFishpairs_outlierincluded.pdf",width=7,height=7)


pairs(field.collapse[, c('MeanTCTEf','Fish.Total.Biomass.g','Fish.Biomass.g.m2','Fish.Biomass.g.m3','Water.Temperature.C','pH','Transect.Flow.cms')],col=factor(field.collapse$Stream.Code),cex=0.8,cex.labels=0.6)

#dev.off()

#pdf(file="EFishpairs.pdf",width=7,height=7)

pairs(field.removeef[, c('MeanTCTEf','Fish.Total.Biomass.g','Fish.Biomass.g.m2','Fish.Biomass.g.m3','Water.Temperature.C','pH','Transect.Flow.cms')], col = factor(field.removeef$Stream.Code),cex=0.8,cex.labels=0.6)


#dev.off()
