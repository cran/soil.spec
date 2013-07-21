average<-function(spec, type=""){
ssn<-as.character(spec[,1])
uss<-c()
for (i in 1:length(ssn)){
zl<-nchar(ssn[i])
z<-ifelse(zl>9,9,zl)
uss[i]<-substr(ssn[i],1,z)}

spec[,1]<-uss
ussn<-unique(uss)
frr<-c()
avg<-c()
avd<-c()
	for (k in 1:length(ussn)){
	frr<-subset(spec[,2:ncol(spec)],spec[,1]==ussn[k])
 	avg<-t(as.matrix(apply(frr,2,mean)))
	avd<-rbind(avd,avg)}
	
rownames(avd)<-ussn

#Check number of spectra averaged
dim(avd)
k<-ncol(avd)

if(type=="spectra"){
spec.avd<-as.matrix(avd[,1:k])
wave<-as.numeric(substr(colnames(spec.avd),2,15))
	
	#Check which columns are blank
spn.1<-as.numeric(spec.avd[1,])
s<-which(spn.1!="NA")

plot(wave[s],spec.avd[1,s],type="l",col=1,xlim=c(max(wave),min(wave)),pch="n",ylim=c(min(na.omit(spec.avd[,s])),max(na.omit(spec.avd[,s]))),ylab="Infrared measurements",xlab=expression("Wavenumbers cm"^-1),main="Averaged spectra")
for (p in 1:nrow(spec.avd))
{
lines(wave[s],spec.avd[p,s],col=p)}

test <- menu(c("Spectra looks ok - please continue", "Spectrum(s) not ok - please stop"), 
                graphics = T)
            graphics.off()
            if (test == 2) 
                stop("Spectra averaged does not look ok; investigate!")
}
#################################################
avds<-cbind(rownames(avd),avd)
colnames(avds)<-c("SSN",colnames(avd))
avds<-as.data.frame(avds)


output<-list(rep.spec=spec,averaged=avds)
return(output)}
