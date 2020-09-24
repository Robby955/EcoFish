gz_final_flow
gg_pond_sink
gg_zerofish



p.tanh<-coef(TCTmean.htan)
p <- coef(bstick.lm.mean)
plot(flow.new.sum.dat$l2Flow,flow.new.sum.dat$TCTmean,xlab="Log2 Flow (KL)",ylab="Mean TCT",main="Models for Mean TCTs",las=1)
lines(TCTmean.lowess,col='green')
lines(x,p.tanh[3]+p.tanh[4]*(x-p.tanh[1])+p.tanh[5]*(x-p.tanh[1])*tanh((x-p.tanh[1])/p.tanh[2]),lwd=1,lty=2,col='blue')
lines(x.grid, predict(model.bc.mean, x.grid),lty=10,col='firebrick4')
lines(xbs, p[2] + p[3]*xbs + ifelse(xbs>p[1], p[4]*(xbs-p[1]), 0), lwd=2, lty=8, col='darkorchid4')
legend("topright",c("Lowess","Tanh","Bent Cable","Broken Stick"),col=c("green","blue","firebrick4","darkorchid4"),lty=c(1,2,10,8))



gz3
