plot.regr <-
function(x,...){
	
dev.new(width=5,height=2.5*length(x$constituents));
par(mfrow=c(length(x$constituents),2));
for(i in 1:length(x$constituents)){
mi<-min(na.omit(c(x$cal.mea.pre[,x$constituents[i]],x$cal.mea.pre[,paste("pred.",x$constituents[i],sep="")],x$val.mea.pre[,x$constituents[i]],x$val.mea.pre[,paste("pred.",x$constituents[i],sep="")])));
if(mi<0){mi<-0};
ma<-max(na.omit(c(x$cal.mea.pre[,x$constituents[i]],x$cal.mea.pre[,paste("pred.",x$constituents[i],sep="")],x$val.mea.pre[,x$constituents[i]],x$val.mea.pre[,paste("pred.",x$constituents[i],sep="")])));
plot(x$cal.mea.pre[,x$constituents[i]]~x$cal.mea.pre[,paste("pred.",x$constituents[i],sep="")],xlab=paste("Predicted ",x$constituents[i],sep=""),ylab=paste("Measured ",x$constituents[i],sep=""),main="Calibration",xlim=c(mi,ma),ylim=c(mi,ma));
plot(x$val.mea.pre[,x$constituents[i]]~x$val.mea.pre[,paste("pred.",x$constituents[i],sep="")],xlab=paste("Predicted ",x$constituents[i],sep=""),ylab=paste("Measured ",x$constituents[i],sep=""),main="Validation",xlim=c(mi,ma),ylim=c(mi,ma));
}

	
	}
