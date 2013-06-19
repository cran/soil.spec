average <-
function(spec){
	spec<-raw
ssn<-as.character(spec[,1])
uss<-c()
for (i in 1:length(ssn)){
zl<-nchar(ssn[i])
z<-zl-0
uss[i]<-substr(ssn[i],1,z)}

spec[,1]<-uss
ussn<-unique(uss)
frr<-c()
av<-c()
avd<-c()
#for (j in 1:length(ssn)){
	for (k in 1:length(ussn)){
		#for (k in 1:243){
	frr<-subset(spec[,2:ncol(spec)],spec[,1]==ussn[k])
 	av<-t(as.matrix(apply(frr,2,mean)))
	avd<-rbind(avd,av)}
	
	rownames(avd)<-ussn

#Check number of spectra averaged
dim(avd)
#avd[1:3,1:4]
k<-ncol(avd)

spec<-as.matrix(avd[,1:k])
wave<-as.numeric(substr(colnames(spec),2,11))
#Check which columns are blank
spn.1<-as.numeric(spec[1,])
s<-which(spn.1!="NA")

plot(wave[s],spec[1,s],type="l",col=1,xlim=c(max(wave),min(wave)),pch="n",ylim=c(min(na.omit(spec[,s])),max(na.omit(spec[,s]))),ylab="Infrared measurements",xlab=expression("Wavenumbers cm"^-1))
for (p in 1:nrow(spec))
{
lines(wave[s],spec[p,s],col=p)}


#################################################
avds<-cbind(rownames(avd),avd)
colnames(avds)<-c("SSN",colnames(avd))

output<-list(rep.spec=raw,averaged=as.data.frame(avds))
class(output)<-"average"
return(output)
}
