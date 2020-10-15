gz_final_flow
gg_pond_sink
gg_zerofish



p.tanh<-coef(TCTmean.htan)
p <- coef(bstick.lm.mean)
p2<-coef(TCTmean.htan)



p.tanh<-coef(TCTmean.htan)
p <- coef(bstick.lm.mean)

plot(flow.new.sum.dat$l2Flow,flow.new.sum.dat$TCTmean,xlab="Log2 Flow (KL)",ylab="Mean TCT",main="Models for Mean TCTs",las=1)
lines(TCTmean.lowess,col='green')
lines(x1,p.tanh[3]+p.tanh[4]*(x1-p.tanh[1])+p.tanh[5]*(x1-p.tanh[1])*tanh((x1-p.tanh[1])/p.tanh[2]),lwd=1,lty='solid',col='purple')
lines(x.grid, predict(model.bc.mean, x.grid),lty='solid',col='blue')
lines(xbs, p[2] + p[3]*xbs + ifelse(xbs>p[1], p[4]*(xbs-p[1]), 0), lwd=2, lty='solid', col='red')
legend("topright",c("Lowess","Tanh","Bent Cable","Broken Stick"),col=c("green","purple","blue","red"),lty=rep('solid',4))
##LogLikilihoods
paste("Log-Likilihoods")
paste("Lowess: ",round(logLik(meanlowess.lm),3))
paste("Tanh: ",round(logLik(TCTmean.htan),3))
paste("Broken Stick: ",round(logLik(bstick.lm.mean),3))
paste("Bent Cable :", round(logLik(model.bc.mean),3))



gz4
