view.spectra <-
function(spec){
	if(nrow(spec)>0)
		
		k<-ncol(spec)

spec<-as.matrix(spec[,2:k])
wave<-as.numeric(substr(colnames(spec),2,11))
#Check which columns are blank
spn.1<-as.numeric(spec[1,])
s<-which(spn.1!="NA")

plot(wave[s],spec[1,s],type="l",col=1,xlim=c(max(wave),min(wave)),pch="n",ylim=c(min(na.omit(spec[,s])),max(na.omit(spec[,s]))),ylab="Infrared measurements",xlab=expression("Wavenumbers cm"^-1))
for (p in 1:nrow(spec))
{
lines(wave[s],spec[p,s],col=p)}
			#quartz();
}
