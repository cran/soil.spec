make.comp <-
function(spec,new,repl.out="F",sav="F",path="",out="Csm",save.as="ws"){
	
# Check function body:

if(class(spec)!="data.frame" & class(spec)!="matrix"){stop("Invalid argument: class of 'spec' has to be either 'data.frame' or 'matrix'.")};
if(class(spec)=="data.frame"){spec<-as.matrix(spec)};
if(class(new)!="numeric"){stop("Invalid argument: class of 'new' has to be 'numeric'.")};
if(is.na(match(repl.out,c("F","T")))){stop("Invalid argument: 'repl.out' has to be either 'T' or 'F'.")};
if(is.na(match(sav,c("F","T")))){stop("Invalid argument: 'sav' has to either 'F' or 'T'.")}
if(class(out)!="character"){stop("Invalid argument: 'out' has to be of class 'character'.")}
if(is.na(match(save.as,c("ws","csv.file")))){stop("invalid argument: 'save.as' has to be either 'ws' or 'csv.file'.")}

# Make compatible:

comp<-matrix(nrow=nrow(spec),ncol=ncol(spec),dimnames=list(rownames(spec),colnames(spec)));
old.wa<-as.numeric(colnames(spec));
for(i in 1:nrow(spec)){
	comp[i,]<-round(spline(old.wa,spec[i,],method="natural",xout=new)[[2]],6)
if(repl.out=="T"){
		for(l in 1:length(old.wa)){
	if((new[l]>old.wa[1])==T)break
	if(new[l]<old.wa[1])
	{
	comp[i,l]<-NA	
	}}
		for(m in (length(new)):(length(new)-		length(old.wa))){
	if((new[m]<old.wa[length(old.wa)])==T)break;
	if(new[m]>old.wa[length(old.wa)])
	{
	comp[i,m]<-NA	
	}}
}
}
colnames(comp)<-new;

# Prepare output:

output<-list(original.spectra=spec,compatible.spectra=comp);
class(output)<-"make.comp";

if(sav=="T"){
if(path!=""){
		setwd(path);}
		if(save.as=="ws"){save(output,file=out)}
		if(save.as=="csv.file"){write.table(comp,file=paste(out,".csv",sep=""),dec=".",sep=",")}
		}
return(output);
	

	}
