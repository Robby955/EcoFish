# This document contains the R script used for the Dilution Analysis (Chapter 4).

# Set the working directory to location where "Flow_experiment_complete-070419MJA.csv" is saved.

# Note: The Bent Cable fitting searches a large grid of possible values, so don't be alarmed if it takes 30+ seconds to fit.
# Saving plots is commented out.

#Load libraries.
library(knitr) # For rendering
library(dplyr)  #For piping and general use
library(kableExtra) #For making kable plots
library(ggplot2) #ggplot for plotting
library(ggthemes) # themes for ggplots
library(nlstools) #For nonlinear regression
library(SiZer) #For bent.cable function

# For sending plots to file
sink.indicator<-TRUE

# Read in the dilution data
flow.dat<-read.csv("Flow-4July2019.csv") 

#The original name of the excel/csv file was "Flow_experiment_complete-070419MJA.csv"

#Remove unimportant/blank columns
flow.dat <- flow.dat[1:864,-c(2,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22:ncol(flow.dat))]

#Rename the columns for easy acess.
names(flow.dat)<-c("Sort.Code","Site.ID","Sample.replicate","Lab.Code","Transformed.Ct")  

#We need to remove samples which did not pass the integrity score. These correspond to sort.code=84 and sort.code=31
sort.index<-c(31,84)

#Define negation function to use
'%ni%' <- Negate('%in%')

#Use negation to remove bad sort codes
flow.dat<-flow.dat[flow.dat$Sort.Code %ni% sort.index, ] #Remove those from our data set.

#Use grep, these indictors will tell us what rows contain that kL of flow and tank number.
# ^ in grep means 'starts with' and $ in grep means ends with.

tank.ind19<-grep(" Tank 19$", flow.dat$Site.ID)
tank.ind20<-grep(" Tank 20$", flow.dat$Site.ID)
tank.ind21<-grep(" Tank 21$", flow.dat$Site.ID)
tank.ind24<-grep(" Tank 24$", flow.dat$Site.ID)

flow.ind10<-grep("^10kL Flow,",flow.dat$Site.ID)
flow.ind20<-grep("^20kL Flow,",flow.dat$Site.ID)
flow.ind40<-grep("^40kL Flow,",flow.dat$Site.ID)
flow.ind80<-grep("^80kL Flow,",flow.dat$Site.ID)
flow.ind160<-grep("^160kL Flow,",flow.dat$Site.ID)
flow.ind1000<-grep("^1000kL Flow,",flow.dat$Site.ID)
flow.indSink<-grep("^sink",flow.dat$Site.ID)
flow.indPond<-grep("^pond",flow.dat$Site.ID)

flow.dat$Tank <- 1
flow.dat$Tank[tank.ind19] <- 19
flow.dat$Tank[tank.ind20] <- 20
flow.dat$Tank[tank.ind21] <- 21
flow.dat$Tank[tank.ind24] <- 24

flow.dat$Flow<-0
flow.dat$Flow[flow.ind10]<-10
flow.dat$Flow[flow.ind20] <-20
flow.dat$Flow[flow.ind40] <- 40
flow.dat$Flow[flow.ind80] <- 80
flow.dat$Flow[flow.ind160] <- 160
flow.dat$Flow[flow.ind1000] <- 1000
flow.dat$Flow[flow.indSink] <- "sink"
flow.dat$Flow[flow.indPond] <- "pond"

flow.dat$Lab.Code.char<-as.character(flow.dat$Lab.Code)
g1=grep("10.1a" ,flow.dat$Lab.Code.char)
g2=grep("10.2a" ,flow.dat$Lab.Code.char)
g3=grep("10.3a" ,flow.dat$Lab.Code.char)
flow.dat$Lab.Code.bool<-TRUE
flow.dat$Lab.Code.bool[c(g1, g2, g3)] <- FALSE

flow.cut<- flow.dat %>%           #New data frame which only contains those we wish to keep.
  filter(Lab.Code.bool==TRUE)  

# Set knitr options for kables , NA show as blank
options(knitr.kable.NA = '')

# Define a counting function that counts number of unique elements.
numrep <- function(x){ length(unique(x))}

#Define some vectors that we will use for kables.
tf=c("",0,10,20,40,80,160,1000,'pond','sink')
tf2=c("","",round(log2(10),2),round(log2(20),2),round(log2(40),2),round(log2(80),2),round(log2(160),2),round(log2(1000),2),"","")

tf3=c(1,"","","","","","","",3,2)
tf4=c(19,3,4,2,3,3,3,4,"","")
tf5=c(20,3,4,3,3,3,3,3,"","")
tf6=c(21,3,3,3,3,3,3,4,"","")
tf7=c(24,3,3,3,3,3,3,4,"","")

reps3=cbind(tf,tf2,tf3,tf4,tf5,tf6,tf7)

names(reps3)=c("",1,19,20,21,24)


#Apply the function over each sort code and level of flow/tank.
reps2 <- tapply(flow.cut$Sort.Code, list(flow.cut$Flow, flow.cut$Tank), numrep)

# Create first flow kable, need to use latex if want to re-render.

k_sample_rep1=kable(reps2, caption="Number of Sample.replicates by Flow/Tank",format='latex',booktabs=T)%>%
  kable_styling("striped") %>%
  add_header_above(c("Flow (kL)" = 1, "Tank" = 5))


k_sample_rep2=kable(reps3,format='latex',booktabs=T)%>%
  kable_styling("striped") %>%
  add_header_above(c("Flow (kL)" = 1,"Log2(Flow)"=1, "Tank" = 5))


#initialize blank columns.
flow.dat$tanknum=0
flow.dat$flownum=0

#Populate new columns.

for(i in 1:nrow(flow.dat)){
  
  if(flow.dat[i,'Tank']==19){
    flow.dat[i,'tanknum']='Tank 19'
  }
  if(flow.dat[i,'Tank']==20){
    flow.dat[i,'tanknum']='Tank 20'
  }
  if(flow.dat[i,'Tank']==21){
    flow.dat[i,'tanknum']='Tank 21'
  }
  if(flow.dat[i,'Tank']==24){
    flow.dat[i,'tanknum']='Tank 24'
  }
  if(flow.dat[i,'Tank']==1 & flow.dat[i,'Flow']=='pond'){
    flow.dat[i,'tanknum']='Tank 1: Pond'
  }
  if(flow.dat[i,'Tank']==1 & flow.dat[i,'Flow']=='sink'){
    flow.dat[i,'tanknum']='Tank 1: Sink'
  }
}


for(i in 1:nrow(flow.dat)){
  
  if(flow.dat[i,'Flow']=='10'){
    flow.dat[i,'flownum']='10 kL'
  }
  if(flow.dat[i,'Flow']=='20'){
    flow.dat[i,'flownum']='20 kL'
  }
  if(flow.dat[i,'Flow']=='40'){
    flow.dat[i,'flownum']='40 kL'
  }
  if(flow.dat[i,'Flow']=='80'){
    flow.dat[i,'flownum']='80 kL'
  }
  if(flow.dat[i,'Flow']=='160'){
    flow.dat[i,'flownum']='160 kL'
  }
  if(flow.dat[i,'Flow']=='1000'){
    flow.dat[i,'flownum']='1000 kL'
  }
}


#Turn new columns into factors and assign an order.

flow.dat$flownum=factor(flow.dat$flownum)

# Possible flow values
flow.dat$flownum <- factor(flow.dat$flownum,levels=c('10 kL',"20 kL","40 kL","80 kL","160 kL","1000 kL"))

# Jitter x axis for plotting
jitter <- position_jitter(width = 0.15, height =0)


# Plot of zero fish
gg_zerofish=ggplot(data=flow.dat%>%filter(Flow=='0'))+
  geom_point(mapping=aes(x=Sample.replicate,y=Transformed.Ct),alpha=0.5,position=jitter)+
  facet_wrap(.~tanknum)+theme_bw()+ylab("TCT")+
  xlab("Sample Replicate")+
  theme(strip.text.x = element_text(margin = margin( b = 1.5, t = 3),size=8))+
  scale_x_continuous(breaks=c(1,2,3),
                     labels=c("1", "2", "3"))+
  theme(strip.background =element_rect(fill="tomato1"))

  #pdf('..\\EcoFish-master\\Thesis\\Chapter4Images\\zerofishtct.pdf')
  #gg_zerofish
  #dev.off()
  


#ggsave('zerofishtct.pdf',gg_zerofish,height=6,width=6,dpi=400)

# Sink and Pond Plot

gg_pond_sink=ggplot(data=flow.dat%>%filter(Flow==c('pond','sink')))+
  geom_point(mapping=aes(x=Sample.replicate,y=Transformed.Ct),alpha=0.5,position=jitter)+
  facet_wrap(.~tanknum)+theme_bw()+ylab("TCT")+
  xlab("Sample Replicate")+theme(strip.text.x = element_text(margin = margin( b = 1.5, t = 3),size=8))+
  theme(strip.background =element_rect(fill="tomato1"))+
  theme( axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"))


  #pdf('..\\EcoFish-master\\Thesis\\Chapter4Images\\pondandsink.pdf')
  #gg_pond_sink
  #dev.off()

#ggsave('pondandsink.pdf',gg_pond_sink,height=6,width=6,dpi=400)

# TCT by Tank and Flow

gz_final_flow=ggplot(data=flow.dat%>%filter(Flow!='0')%>%filter(Flow!='sink')%>%filter(Flow!='pond')%>%filter(Lab.Code.bool==TRUE))+
  geom_point(mapping=aes(x=Sample.replicate,y=Transformed.Ct),alpha=0.5,position=jitter)+ylim(0,22.5)+
  facet_grid(flownum~tanknum)+theme_bw()+ylab("TCT")+
  xlab("Sample Replicate")+theme(strip.text.x = element_text(margin = margin( b = 1.5, t = 3),size=8))+
  scale_x_continuous(breaks=c(1,2,3,4),labels=c("1", "2", "3","4"))+
  theme(strip.background =element_rect(fill="tomato1"))+
  theme( axis.line = element_line(colour = "black", 
                                  size = 0.5, linetype = "solid"))


  #pdf('..\\EcoFish-master\\Thesis\\Chapter4Images\\pondandsink.pdf')
  #gg_pond_sink
  #dev.off()

#ggsave('finalflowplots.pdf',gz_final_flow,height=6,width=6,dpi=500)

# Numeric versions of Flow, Tank and Log2 Flow

iFlowN <- c(10, 20, 40, 80, 160, 1000)
iTank <- c(19, 20, 21, 24)
flow.init=c(rep(3.322,4),rep(4.322,4),rep(5.322,4),rep(6.322,4),rep(7.322,4),rep(9.965,4))


#Init a matrix that we will use to store values for more summary statistics.
init.frame<-matrix(0,nrow=length(iFlowN)*length(iTank),ncol=length(iTank)+2)
count1=1

# Populate
for (j in iFlowN){
  for (k in iTank){
    flow.sub<-flow.cut %>%
      filter(Flow==j,Tank==k)
    line.fill<-c(j,k,min(flow.sub$Transformed.Ct),max(flow.sub$Transformed.Ct),median(flow.sub$Transformed.Ct),round(mean(flow.sub$Transformed.Ct),2))
    
    init.frame[count1, ]<-line.fill
    count1<-count1+1
  }
}

#Flow values of interest , except zero flow.
iFlow<-c(10,20,40,80,160,1000)

#Round to two decimal places.
ynum<-round(log2(iFlow),2)

# Turn to characters instead of integer.
ynum.char=as.character(ynum) 


#initialize a matrix.
num.matrix<-matrix(0,ncol=1,nrow=length(iFlow)-1) 

reps3=reps2[1:7,]
reps3=reps3[c('10','20','40','80','160','1000'), ] #use rep 2 from Table 1 to give total number of replicates for each flow.

for(i in 1:(length(iFlow))){ #Populate the matrix
  num.matrix[i]= sum(reps3[i, ],na.rm=TRUE) #na.rm makes sure there is not an error when we sum over the data.
}

num.reps<-as.vector(num.matrix)  

log.dataframe<-data.frame(iFlow,ynum,num.reps)

names(log.dataframe)<-c("Flow","log2(Flow)","Total number of Sample.replicates")

# Second flow kable
kf=kable(log.dataframe,caption="Values of log2(Flow) for associated Flow value in kL",format="latex",booktabs=T) 



flow.cut.order<-flow.cut[order(flow.cut$Sort.Code), ] #Form a new data frame which is our original flow.cut , but ordered with respect to increasing Sort.Code.


#flow.new removes any 0 flow, sink and pond from our ordered data frame.
flow.new<- flow.cut.order %>% 
  filter(Flow!=0 & Flow !='sink' & Flow!='pond' )


#flow.new.zero is a subset of data only containing a Flow of 0.
flow.new.zero<-flow.cut.order %>%
  filter(Flow==0)

# Start forming summary statistics over sort codes

jmed_flow.test<-tapply(flow.new$Transformed.Ct, factor(flow.new$Sort.Code), median) #Gives the median for each sort code.

jmean_flow<-tapply(flow.new$Transformed.Ct,factor(flow.new$Sort.Code),mean) #Gives the mean for each sort code.

jmed_flow.test.zero<-tapply(flow.new.zero$Transformed.Ct, factor(flow.new.zero$Sort.Code), median) #Gives the median for each sort code.

jmean_flow.zero<-tapply(flow.new.zero$Transformed.Ct,factor(flow.new.zero$Sort.Code),mean) #Gives the mean for each sort code.


#picks off first occurrence of a distinct Sort.Code, assume data sorted by Sort.Code

jind <- c(TRUE, flow.new$Sort.Code[-length(flow.new$Sort.Code)]!=flow.new$Sort.Code[-1]) #Only the 10kL of interest is included.


jind.zero<-c(TRUE, flow.new.zero$Sort.Code[-length(flow.new.zero$Sort.Code)]!=flow.new.zero$Sort.Code[-1])

flow.new.sum.dat<-data.frame(flow.new[jind, ]) #Only the 10kL of intrerst is included.

flow.new.sum.zero<-data.frame(flow.new.zero[jind.zero, ]) #This contains the 12 sort codes for 0 fish.


flow.new.sum.dat$TCTmed<-jmed_flow.test 
flow.new.sum.dat$TCTmean<-jmean_flow  
flow.new.sum.zero$TCTmean<-jmean_flow.zero  


flow.new.sum.dat$FlowN<-as.numeric(flow.new.sum.dat$Flow)
flow.new.sum.dat$l2Flow<-log2(flow.new.sum.dat$FlowN) #make a column called l2Flow for log2(flow)
flow.new.sum.dat$TankF<-as.factor(flow.new.sum.dat$Tank)


medTCT.red<-as.matrix(flow.new.sum.dat$TCTmed,nrow=1,ncol=sum(num.reps)) #Refers to cut data (only the 10kl we want included)
logFlow.red<-as.matrix(log2(flow.new.sum.dat$FlowN),nrow=sum(num.reps),ncol=1) #Refers to full data (only the 10KL we want included)


msg.UCV<-"" #This definition allows the code to run without an error

l.one.line=lm(TCTmed~l2Flow,data=flow.new.sum.dat)#Only included the 10kl we want. A simple linear fit.

options(digits=3)




term_name=c("Intercept","Log2(Flow)")
kcoef=c(round(coef(l.one.line)[[1]],2),round(coef(l.one.line)[[2]],2))
kse=c(1.46,0.22)
tval=c(16.75,-12.62)
pval=c('<2e-16','<2e-16')


dff2=data.frame(term_name,kcoef,kse,tval,pval)
names(dff2)=c("Term","Estimate",'Std Error','t value','Pr(>|t|)')

ktz2=kable(dff2,format="latex",booktabs=T,escape = FALSE)%>%
  kable_styling("striped")

term_name=c("Intercept","Log2(Flow)")
kcoef=c(round(coef(l.one.line.mean)[[1]],2),round(coef(l.one.line.mean)[[2]],2))
kse=c(1.36,0.21)
tval=c(17.55,-12.99)
pval=c('<2e-16','<2e-16')


dff2=data.frame(term_name,kcoef,kse,tval,pval)
names(dff2)=c("Term","Estimate",'Std Error','t value','Pr(>|t|)')

ktz2=kable(dff2,format="latex",booktabs=T,escape = FALSE)%>%
  kable_styling("striped")


flow.l.one.line=l.one.line
l.one.line.mean=lm(formula = TCTmean ~ l2Flow, data = flow.new.sum.dat) #Same thing as l.one.lone but for mean.


col.vec=as.factor(flow.new.sum.dat$Tank-18)
col_bluevec=flow.new.sum.dat$Tank-18
col_bluevec[col_bluevec==3]=4


gframe_flow.med=data.frame(log2(flow.new.sum.dat$FlowN),flow.new.sum.dat$TCTmed,col.vec,col.vec)
gframe_flow.mean=data.frame(log2(flow.new.sum.dat$FlowN),flow.new.sum.dat$TCTmean,col.vec,col.vec)

names(gframe_flow.med)=c("l2Flow","TCTmed","Tank","tank_shape")
names(gframe_flow.mean)=c("l2Flow","TCTmean","Tank","tank_shape")


gframe_flow.med$tankf=as.factor(gframe_flow.med$tank_shape)

med.coef=coef(l.one.line)


p2flow=ggplot(data=gframe_flow.med)+
  geom_point(mapping=aes(x=l2Flow,y=TCTmed,color=tankf,shape=tankf))+
  scale_shape_manual(name="Tank",values=c(1,2,3,6),labels=c("T19","T20","T21","T24"))+
  scale_color_manual(name="Tank",values=c("black", "red", "blue","magenta"),labels=c("T19","T20","T21","T24"))+ 
  #geom_abline(intercept=med.coef[1],slope=med.coef[2])+  
  geom_smooth(mapping=aes(x=l2Flow,y=TCTmed),color='black',method="lm")+
  xlab("Log2 Flow (KL)")+
  ylab("Median Transformed Ct")+ 
  xlim(3,10)+
  ylim(-6,20)+
  geom_hline(aes(yintercept=0),lty='dashed')+
  theme(panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))+

  #pdf('..\\EcoFish-master\\Thesis\\Chapter4Images\\flowmedtctgg.pdf')
  #p2flow
  #dev.off()


#ggsave('flowmedtctgg.pdf',p2flow,width=6,height=6,dpi=500)




# Linear Regression on median
l.one.line <- lm(TCTmed~l2Flow,data=flow.new.sum.dat)

flow.new.sum.dat$TankF<-as.factor(flow.new.sum.dat$Tank)
l.tank<-lm(TCTmed~l2Flow+TankF,data=flow.new.sum.dat) 
l.four.line<-lm(TCTmed~l2Flow+TankF+l2Flow*TankF,data=flow.new.sum.dat)


dff2=data.frame(term_name,kcoef,kse,tval,pval)
names(dff2)=c("Term","Estimate",'Std Error','t value','Pr(>|t|)')

ktz2=kable(dff2,format="latex",booktabs=T,escape = FALSE)%>%
  kable_styling("striped")


dff2=data.frame(term_name,kcoef,kse,tval,pval)
names(dff2)=c("Term","Estimate",'Std Error','t value','Pr(>|t|)')

ktz2=kable(dff2,format="latex",booktabs=T,escape = FALSE)%>%
  kable_styling("striped")



a1=anova(l.one.line,l.tank,l.four.line)

# Linear Regressions on Mean TCT
l.one.line.mean=lm(formula = TCTmean ~ l2Flow, data = flow.new.sum.dat)
l.tank.mean=lm(formula = TCTmean ~ l2Flow + TankF, data = flow.new.sum.dat)
l.four.line.mean=lm(formula = TCTmean ~ l2Flow + TankF + l2Flow * TankF, data = flow.new.sum.dat)


a2=anova(l.one.line.mean,l.tank.mean,l.four.line.mean)

options(digits=3)

n2=c("Broken-Stick","Bent-Cable","Hyperbolic Tangent","Lowess")
n3=c(614.76,614.77,595.24,589.70)
n4=c(-189.24,-189.24,-187.99,-185.34)
n5=c("Mean TCT","Mean TCT","Mean TCT","Mean TCT")

dz=data.frame(n2,n5,n3,n4)
names(dz)=c("Model","Response","RSS","Log-Likelihood")

kt=kable(dz,format="latex",booktabs=T,escape = FALSE)%>%
  kable_styling("striped")


# Broken Stick
bstick.lm.mean <- nls(TCTmean ~ cbind("intercept"=1, "l2Flow"=l2Flow, 
                                      "l2FlowBr"=ifelse(l2Flow > Br, l2Flow - Br, 0)),
                      start=list(Br=6), algorithm="plinear", data=flow.new.sum.dat)

pp=predict(bstick.lm.mean,flow.new.sum.dat$l2Flow)
rss1=(flow.new.sum.dat$TCTmean-pp)

# Hyperbolic Tangent Models
x1 <- seq(min(flow.new.sum.dat$TCTmean), max(flow.new.sum.dat$TCTmean), 0.1)
xbs <- with(flow.new.sum.dat, seq(min(l2Flow), max(l2Flow), length=100))

TCTmean.htan<-nls(TCTmean~cbind("intercept"=1,"B1"=l2Flow-br,"B2"=(l2Flow-br)*tanh((l2Flow-br)/gamma)), start=list(br=6.5,gamma=0.5),algorithm='plinear',nls.control(maxiter=10000),data=flow.new.sum.dat)


# Lowess
TCTmean.lowess <- lowess(flow.new.sum.dat$l2Flow, flow.new.sum.dat$TCTmean,f=3/4,iter=3)

meanlowess.lm=lm(TCTmean.lowess$y~TCTmean.lowess$x)


# Bent Cable

model.bc.mean=bent.cable(flow.new.sum.dat$l2Flow,flow.new.sum.dat$TCTmean,grid.size =200)

x.grid <- seq(min(flow.new.sum.dat$TCTmean), max(flow.new.sum.dat$TCTmean), length=200)


new.seq4=seq(min(flow.new.sum.dat$l2Flow),max(flow.new.sum.dat$l2Flow),length=200)
new.seq4=data.frame(new.seq4)
names(new.seq4)='l2Flow'
new.seq5=new.seq4
names(new.seq5)='x'


a.hat200=model.bc.mean$alpha #MLE of alpha
g.hat200=model.bc.mean$gamma #MLE of gamma

#build our q matrix

q.hat200=rep(0) #initialize

for(i in 1:length(new.seq5)){
  I.middle <- 1 * (abs(new.seq5[i] - a.hat200) < g.hat200)
  I.left <- 1 * (new.seq5[i] >= a.hat200 + g.hat200)
  q.hat200[i] <- (new.seq5[i] - a.hat200 + g.hat200)^2/(4 * g.hat200) * I.middle + 
    (new.seq5 - a.hat200) * I.left
}

nf200=data.frame(new.seq5,q.hat200)
names(nf200)=c('x','q')

pz200=predict(model.bc.mean$model,newdata=nf200,interval='confidence')


new.tank=data.frame(flow.init,init.frame[,6]) #extract from earlier
names(new.tank)=c("Flow","Mean")

#The S.E of break point estimate is 0.5417, times 1.96 is 1.062
se=0.5417
ci=se*1.96

meanflows=data.frame(new.tank) #Our points, mean over tanks

flowfit200=data.frame(new.seq5$x,pz200[,1],pz200[,2],pz200[,3]) # Build a data frame that contains the x values and fits+C.I
names(flowfit200)=c('x','fit','lower','upper')


breaks=c(3,4,5,6,model.bc.mean$alpha,7.5,8,9,10)

gz4=ggplot(data=flowfit200)+
  geom_point(data=new.tank,aes(x=Flow,y=Mean),shape=0)+
  ylim(-2,20)+
  geom_line(aes(x=x,y=fit),size=1)+
  geom_ribbon(aes(x=x,ymin=lower, ymax=upper), alpha=0.3)+
  geom_line(aes(x=x,y=upper),linetype='dashed')+
  geom_line(aes(x=x,y=lower),linetype='dashed')+
  geom_segment(aes(x=model.bc.mean$alpha-ci, y=-2, xend=model.bc.mean$alpha-ci, yend=6.6),linetype='dashed')+
  geom_segment(aes(x=model.bc.mean$alpha+ci, y=-2, xend=model.bc.mean$alpha+ci, yend=1.3),linetype='dashed')+
  geom_segment(aes(x=model.bc.mean$alpha, y=-2, xend=model.bc.mean$alpha, yend=2.2),linetype='solid')+
  xlab("Log2 Water Volume (KL) diluted")+
  ylab("Mean Transformed CT")+
  theme(panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))+ 
  scale_x_continuous(name="Log2 Water Volume (KL) diluted", 
                     breaks=breaks,labels=c('3','4','5','6',round(model.bc.mean$alpha,2),'7.5','8','9','10'))+ scale_y_continuous(expand = c(0,0))

gz4

  #pdf('..\\EcoFish-master\\Thesis\\Chapter4Images\\flowggplot1.pdf')
  #gz4
  #dev.off()

#ggsave('flowggplot1.pdf',gz4,width=6,height=6,dpi=500)

# For plotting the niche models

x <- seq(min(flow.new.sum.dat$TCTmed), max(flow.new.sum.dat$TCTmed), 0.1)

TCTmean.htan<-nls(TCTmean~cbind("intercept"=1,"B1"=l2Flow-br,"B2"=(l2Flow-br)*tanh((l2Flow-br)/gamma)), start=list(br=6.5,gamma=0.5),algorithm='plinear',nls.control(maxiter=10000),data=flow.new.sum.dat)

p2<-coef(TCTmean.htan)


#gz_final_flow
#gg_pond_sink
#gg_zerofish

p.tanh<-coef(TCTmean.htan)
p <- coef(bstick.lm.mean)


 #pdf('..\\EcoFish-master\\Thesis\\Chapter4Images\\meanTCTmodelcomparison.pdf')
plot(flow.new.sum.dat$l2Flow,flow.new.sum.dat$TCTmean,xlab="Log2 Flow (KL)",ylab="Mean TCT",las=1)
lines(TCTmean.lowess,col='green')
lines(x1,p.tanh[3]+p.tanh[4]*(x1-p.tanh[1])+p.tanh[5]*(x1-p.tanh[1])*tanh((x1-p.tanh[1])/p.tanh[2]),lwd=1,lty='solid',col='purple')
lines(x.grid, predict(model.bc.mean, x.grid),lty='solid',col='blue')
lines(xbs, p[2] + p[3]*xbs + ifelse(xbs>p[1], p[4]*(xbs-p[1]), 0), lwd=2, lty='solid', col='red')
legend("topright",c("Lowess","Tanh","Bent Cable","Broken Stick"),col=c("green","purple","blue","red"),lty=rep('solid',4))
#dev.off()

gz_final_flow
