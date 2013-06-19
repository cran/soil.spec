copy <-
function(source,destination,ddmonyyyy){
start.period<-des.di<-NULL # Setting the date and source to NULL
finf<-file.info(dir(source))
fol<-as.vector(rownames(finf))

lo<-c()
all.filedetails<-c()
for(u in 1:length(fol)){
folu<-fol[u]
#########################
#Specify the start date
#######################



z <- as.Date(start.period, "%d%b%Y")
finp<-paste(source,folu,sep="/")
setwd(finp)
tmpf<-list.files()

std<-tmpf

stds<-c()
for(l in 1:length(std)) {
	stds<-c(stds,length(std[[l]]))}

l.stds<-which(stds>1)

hts.stds<-c()
for (f in 1:length(l.stds)){
hts.stds<-c(hts.stds,tmpf[l.stds[f]])}

#
file.copy(hts.stds,des.di,recursive=TRUE,overwrite=TRUE)


all<-list.files(des.di)
z<-length(all)
setwd(des.di)

new<-c()

	for (a in 1:length(hts.stds)){
		new[a]<-which(all==hts.stds[a])
		}
	
#new<-which(new>10)


a<-length(new)



ifelse(a>0,for(i in 1:length(new)){file.rename(all[new[i]],paste(folu,new[i]," Mua",".0",sep=""))},"none")
}
}
