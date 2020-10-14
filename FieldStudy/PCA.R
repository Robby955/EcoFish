
# This code is for applying Principal Component Analysis on the Stream Field Data.

# Load Libraries

library(tidyverse) # For creating tidy dataframes and piping
library(ggplot2) #For plotting
library(magrittr)
library(dplyr)
library(factoextra)

#The original file name is Ecofish run of river salmonids - Combined raw ct scores and fish data. Changed up EcoFieldUp to read in.
fieldData=read.csv("EcoFieldUp.csv") 


# We remove observations from Stream DDD in 2018 as they failed integritE tests.
fieldData$StreamYear=paste(fieldData$Stream.Code,fieldData$Year) #Create a column that combines Stream.Code and Year.

fieldData=fieldData%>%                    
  filter(StreamYear!='DDD 2018')
fieldData$SSRS=paste(fieldData$Stream.Code,fieldData$Sample,fieldData$Reach,fieldData$Site.number) #Creates a unique identifier for each set of tech replicates Stream/Sample/Reach/Site



#Create new columns to assist in coding and plotting.

fieldData$ID=fieldData$ï..ID  #Create a column ID to replace the odd symbol name.

fieldData$Site.numberF=as.factor(fieldData$Site.number) #Create a site number as a factor column.

fieldData$StreamReach=paste(fieldData$Stream.Code,fieldData$Reach) #Combine Stream and Reach to get a new column.

fieldData$StreamReachYear=paste(fieldData$Stream.Code,fieldData$Reach,fieldData$Year)

#Create Volume Columns for each species and fish general.

fieldData$Transect.Flow.ms=fieldData$Transect.Flow.cms*0.01 # Meter per second isntead of cm/s

fieldData$CO.Biomass.g.m3=fieldData$CO.Biomass.g.m2/fieldData$Transect.Depth.m #Create Volume by accounting for Transect Depth.

fieldData$CT.Biomass.g.m3=fieldData$CT.Biomass.g.m2/fieldData$Transect.Depth.m #Create Volume by accounting for Transect Depth.

fieldData$RB.Biomass.g.m3=fieldData$RB.Biomass.g.m2/fieldData$Transect.Depth.m #Create Volume by accounting for Transect Depth.

fieldData$Fish.Biomass.g.m3=fieldData$Fish.Biomass.g.m2/fieldData$Transect.Depth.m #Create Volume by accounting for Transect Depth.


jind<-rep(c(TRUE,rep(FALSE,7)),nrow(fieldData)/8) #We collapse over each set of eight techinical replicates.

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
fviz_eig(res.pca_new,addlabels = TRUE,barfill='blueviolet')


#Contributions
fviz_pca_var(res.pca_new,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("blue", "purple", "red"),
             
             repel = TRUE     # Avoid text overlapping
)


# Plot by first two components by stream/reach.

fviz_pca_ind(res.pca_new,
             col.ind = p, # color by groups
             palette = c("red","blue","green","purple","black","orange"),
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Stream",
             repel = TRUE,geom.ind = "point")