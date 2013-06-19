align <-
function(file1,file2){
#Set file 1to be the master table
file2.file1.all<-c()
	for (i in 1:ncol(file2)){
file2.file1<-which(colnames(file2)[i]==colnames(file1))
file2.file1.all<-c(file2.file1.all,file2.file1)}

if(length(file2.file1.all)<1){stop("The files are incompatible: at least one column should be found in both files")};

#Map the columns to match those in the master
blank.file2_master<-matrix("NA",nrow(file2),ncol(file1))
colnames(blank.file2_master)<-colnames(file1)

for (i in 1:nrow(blank.file2_master)){
	for (j in 1:length(file2.file1.all)){
		blank.file2_master[i,file2.file1.all[j]]<-as.vector(file2[i,j])
	}
}

file1.file2<-rbind(file1,as.data.frame(blank.file2_master))
combined<-file1.file2
return(combined)
}
