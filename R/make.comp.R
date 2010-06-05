make.comp <-
function(spec,new.waveb,repl.out.range="FALSE",sav="FALSE",save.path="NULL",output.name="Compatible spectral matrix",save.as="workspace"){
	
# Check function body:

if(class(spec)!="data.frame" & class(spec)!="matrix"){stop("Invalid argument: class of 'spec' has to be either 'data.frame' or 'matrix'.")};
if(class(spec)=="data.frame"){spec<-as.matrix(spec)};
if(class(new.waveb)!="numeric"){stop("Invalid argument: class of 'new.waveb' has to be 'numeric'.")};
if(is.na(match(repl.out.range,c("FALSE","TRUE")))){stop("Invalid argument: 'repl.out.range' has to be either 'TRUE' or 'FALSE'.")};
if(is.na(match(sav,c("FALSE","TRUE")))){stop("Invalid argument: 'sav' has to either 'FALSE' or 'TRUE'.")}
if(class(output.name)!="character"){stop("Invalid argument: 'output.name' has to be of class 'character'.")}
if(is.na(match(save.as,c("workspace","csv.file")))){stop("invalid argument: 'save.as' has to be either 'workspace' or 'csv.file'.")}

# Make compatible:

comp<-matrix(nrow=nrow(spec),ncol=ncol(spec),dimnames=list(rownames(spec),colnames(spec)));
old.wa<-as.numeric(colnames(spec));
for(i in 1:nrow(spec)){
	comp[i,]<-round(spline(old.wa,spec[i,],method="natural",xout=new.waveb)[[2]],6)
if(repl.out.range=="TRUE"){
		for(l in 1:length(old.wa)){
	if((new.waveb[l]>old.wa[1])==T)break
	if(new.waveb[l]<old.wa[1])
	{
	comp[i,l]<-NA	
	}}
		for(m in (length(new.waveb)):(length(new.waveb)-		length(old.wa))){
	if((new.waveb[m]<old.wa[length(old.wa)])==T)break;
	if(new.waveb[m]>old.wa[length(old.wa)])
	{
	comp[i,m]<-NA	
	}}
}
}
colnames(comp)<-new.waveb;

# Prepare output:

output<-list(original.spectra=spec,compatible.spectra=comp);
class(output)<-"make.comp";

if(sav=="TRUE"){
if(save.path!="NULL"){
		setwd(save.path);}
		if(save.as=="workspace"){save(output,file=output.name)}
		if(save.as=="csv.file"){write.table(comp,file=paste(output.name,".csv",sep=""),dec=".",sep=",")}
		}
return(output);
	

	}

