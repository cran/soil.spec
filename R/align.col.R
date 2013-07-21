rm(list=ls())
align.col<-function(x,y){
#Set file 1to be the master table
y.file1.all<-c()
	for (i in 1:ncol(y)){
y.file1<-which(colnames(y)[i]==colnames(x))
y.file1.all<-c(y.file1.all,y.file1)}

if(length(y.file1.all)<1){stop("The files are incompatible: at least one column MUST be found in both files")};

#Map the columns to match those in the master
blank.file2_master<-matrix("NA",nrow(y),ncol(x))
colnames(blank.file2_master)<-colnames(x)

for (i in 1:nrow(blank.file2_master)){
	for (j in 1:length(y.file1.all)){
		blank.file2_master[i,y.file1.all[j]]<-as.vector(y[i,j])
	}
}

x.y<-rbind(x,as.data.frame(blank.file2_master))
return(x.y)
}
