explore.data <-
function(datafile){
size<-c()
minim<-c()
maxim<-c()
mn<-c()
dev<-c()
for (j in 1:ncol(datafile)){
	
siz<-ifelse(is.numeric(datafile[,j])=="FALSE",length(datafile[,j]),length(na.omit(datafile[,j])))
mi<-ifelse(is.numeric(datafile[,j])=="FALSE","NA",min(na.omit(datafile[,j])))
mx<-ifelse(is.numeric(datafile[,j])=="FALSE","NA",max(na.omit(datafile[,j])))
av<-ifelse(is.numeric(datafile[,j])=="FALSE","NA",round(mean(na.omit(datafile[,j])),2))
sdev<-ifelse(is.numeric(datafile[,j])=="FALSE","NA",round(sd(na.omit(datafile[,j])),2))

size<-c(size,siz)
minim<-c(minim,mi)
maxim<-c(maxim,mx)
mn<-c(mn,av)
dev<-c(dev,sdev)
summa<-cbind(as.vector(colnames(datafile)),size,mn,minim,maxim,dev)
colnames(summa)<-c("Variable(s)","(n)","Mean","Min","Max","Sd")
}
as.data.frame(summa)
}
