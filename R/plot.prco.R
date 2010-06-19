plot.prco <-
function(x,...){
	if(x$i.pc==0){pairs(x$prcomp$x[,1:5])};
	if(x$i.pc==1){pairs(x$prcomp$x[,1:5])};
	if(x$i.pc>1){pairs(x$prco[,1:x$i.pc])}
	
	x11();
plot(x$prcomp$rotation[,1]~as.numeric(rownames(x$prcomp$rotation)),ylab="Loading values",xlab="Wavebands (cm^-1)",type="l",ylim=c(min(x$prcomp$rotation[,1:2],na.rm=T),max(x$prcomp$rotation[,1:2],na.rm=T)));
lines(x$prcomp$rotation[,2]~as.numeric(rownames(x$prcomp$rotation)),col="red")
legend(x="bottomright",legend=c("First PC","Second PC"),col=c("black","red"),lty=1)

	}

