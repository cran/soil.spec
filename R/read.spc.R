read.spc <-
function(loa="NULL",save.path="NULL",sav="FALSE",output.name="Spectral matrix",save.as="workspace",wavenumber="first.sample"){
	
## Check function body arguments:

tmp<-c() ;
if(loa!="NULL"){
	setwd(loa);
	}
	tmp<-list.files();	
	a<-nchar(tmp);
	a<-substr(tmp,a-2,a);
	a<-which(a!="spc" & a!="SPC");
	if(length(a)>0){tmp<-tmp[-a]};
	rm(a);
	if(length(tmp)==0)stop("Invalid argument: the path is either incorrect or the folder empty.");
if(is.na(match(sav,c("FALSE","TRUE")))){stop("Invalid argument: 'sav' has to either 'FALSE' or 'TRUE'.")}
if(class(output.name)!="character"){stop("Invalid argument: 'output.name' has to be of class 'character'.")}
if(is.na(match(save.as,c("workspace","csv.file")))){stop("invalid argument: 'save.as' has to be either 'workspace' or 'csv.file'.")}

if(is.na(match(wavenumber,c("first.sample","ICRAF")))){stop("Invalid argument: 'wavenumber' has to be either 'first.sample' or 'ICRAF'.")};
if(wavenumber=="ICRAF"){
	n.mpa<-2307
Flast.mpa<-12493.20;
Ffirst.mpa<-3598.69;
bla<-c(rep(NA,n.mpa-2));
for(k in 1:(n.mpa-2)){
		bla[k]<-k*(Flast.mpa-Ffirst.mpa)/(n.mpa-1);
		}
	w<-Ffirst.mpa+bla;
	waveb.mpa<-round(c(Ffirst.mpa,w,Flast.mpa),6);
n.hts<-3578
Flast.hts<-7497.964;
Ffirst.hts<-599.760;
bla<-c(rep(NA,n.hts-2));
for(k in 1:(n.hts-2)){
		bla[k]<-k*(Flast.hts-Ffirst.hts)/(n.hts-1);
		}
	w<-Ffirst.hts+bla;
	waveb.hts<-round(c(Ffirst.hts,w,Flast.hts),6);
	}


# Check which range the spectrum covers (Mir or Nir):

rang<-c();
for(i in 1:length(tmp)){
require(hexView,quietly=T);# 9836 14920 
bla<-readRaw(tmp[i],offset=4,nbytes=4,human="int",size=4,endian="little")[[5]]
test<-readRaw(tmp[i],offset=(bla*4+544+64),nbytes= 431,human="char",endian="little");
test<-blockString(test);
test<-strsplit(test,"..",fixed=T);
test<-substr(test[[1]][17],nchar(test[[1]][17])-2,nchar(test[[1]][17]));
ifelse(test==c("NIR") & nchar(test)==3,rang[i]<-test,rang[i]<-c("MIR"))	}
rm(test,i);

mir<-which(rang=="MIR");
nir<-which(rang=="NIR");

# Check if the spectral are scanned using the same method:

if(length(mir)>0){

# Get the Batch_Labid from each sample name:

test<-c();
for(i in 1:length(tmp[mir])){
	test[i]<-nchar(tmp[mir][i],type="chars")
	}

test.1<-c(rep(NA,length(tmp[mir])));
for(i in 1:length(tmp[mir])){
	test.1[i]<-substr(tmp[mir][i],1,(test[i]-6))
	}
	
test.2<-c(rep(NA,length(tmp[mir])));
for(i in 1:length(tmp[mir])){
	test.2[i]<-substr(tmp[mir][i],1,(test[i]-4))
	}

# Get the positions of the measurement replications for each sample:

a<-unique(test.1);
b<-list();
for(i in 1:length(a)){
	b[[i]]<-which(test.1==a[i]);
	}

# Create aditional information data frame for the replicates:

mir.ad.tmp<-as.data.frame(matrix(nrow=length(mir),ncol=7,dimnames=list(c(1:length(tmp[mir])),c("Batch_Labid","Wavebands","Material","Spectral_range","Scan_date","Resolution","Zero_filling"))));
mir.ad.tmp[,"Batch_Labid"]<-test.1;
mir.ad.tmp[,"Spectral_range"]<-c("MIR");

for(i in 1:length(tmp[mir])){
	n<-readRaw(tmp[mir][i],offset=4,nbytes=4,human="int",size=4,endian="little")[[5]];
	mir.ad.tmp[i,"Wavebands"]<-n;
Ffirst<-readRaw(tmp[mir][1],offset=8,nbytes=8,human="real",size=8,endian="little")[[5]];
Flast<-readRaw(tmp[mir][1],offset=16,nbytes=8,human="real",size=8,endian="little")[[5]];
if(Ffirst>Flast){c<-Ffirst;Ffirst<-Flast;Flast<-c}
bla<-c(rep(NA,n-2));
	for(k in 1:(n-2)){
		bla[k]<-k*(Flast-Ffirst)/(n-1);
		}
	w<-Ffirst+bla;
	waveb.hts<-round(c(Ffirst,w,Flast),6);	
	test<-as.character(readRaw(tmp[mir][i],offset=32,nbytes=4,human="int",size=1.5,machine="binary",signed=F,endian="little"));
test.a<-paste(substr(test,35,42),substr(test,26,33),substr(test,17,21),sep="")
year<-substr(test.a,1,12);
c<-nchar(year)
d<-c(rep(NA,c));
k<-0;
for(j in c:1){
	d[j]<-as.numeric(substr(year,j,j))*2^k;
	k<-k+1
	}
year <-sum(d)
month<-substr(test.a,13,16);
c<-nchar(month)
d<-c(rep(NA,c));
k<-0;
for(j in c:1){
	d[j]<-as.numeric(substr(month,j,j))*2^k;
	k<-k+1
	}
month<-sum(d)
day<-substr(test.a,17,21);
c<-nchar(day)
d<-c(rep(NA,c));
k<-0;
for(j in c:1){
	d[j]<-as.numeric(substr(day,j,j))*2^k;
	k<-k+1
	}
day <-sum(d)
mir.ad.tmp[i,"Scan_date"]<-paste(month,"/",day,"/",year,sep="")

test<-readRaw(tmp[mir][i],offset=(length(waveb.hts)*4+544+64),nbytes= 431,human="char",endian="little");
test<-blockString(test);
test<-strsplit(test,"..",fixed=T);
if(grepl("Plant",test[[1]][8])){mir.ad.tmp[i,"Material"]<-"Plant"};
if(grepl("Soil",test[[1]][8])){mir.ad.tmp[i,"Material"]<-"Soil"};
mir.ad.tmp[i,"Resolution"]<-as.integer(substr(test[[1]][7],nchar(test[[1]])[7],nchar(test[[1]])[7]));
mir.ad.tmp[i,"Zero_filling"]<-as.integer(substr(test[[1]][20],nchar(test[[1]])[20],nchar(test[[1]])[20]));
	}

# Make checks if parameters are the same for all samples:

if(length(summary(as.factor(mir.ad.tmp[,"Material"])))!=1 |  length(summary(as.factor(mir.ad.tmp[,"Resolution"])))!=1 | length(summary(as.factor(mir.ad.tmp[,"Zero_filling"])))!=1 | length(summary(as.factor(mir.ad.tmp[,"Wavebands"])))!=1){test<-menu("Scan settings for mir-spectra not uniform. Find a txt-file with detailed scaning information named 'mir.ad.tmp.txt' in 'loa'.",graphics=T);if(test==1){write.table(mir.ad.tmp,file="mir.ad.tmp.txt",sep=",",dec=".",row.names=FALSE)}}

if(length(summary(as.factor(mir.ad.tmp[,"Material"])))!=1){ma<-names(summary(as.factor(mir.ad.tmp[,"Material"])));mat<-c();mat[1]<-paste("Different materials for mir-spectra. Shall only ",ma[1]," be saved?",sep="");for(i in 2:length(ma)){mat[i]<-paste("Different materials for mir-spectra. Shall only ",ma[i]," be saved?",sep="")};test<-menu(mat,graphics=T);for(i in 1:length(ma)){if(test==i){bla<-which(mir.ad.tmp[,"Material"]==ma[i]);mir.ad.tmp<-mir.ad.tmp[bla,]}}}


if(length(summary(as.factor(mir.ad.tmp[,"Resolution"])))!=1){re<-as.integer(names(summary(as.factor(mir.ad.tmp[,"Resolution"]))));res<-c();res[1]<-paste("Different resolutions for 'mir' spectra. Shall only ",re[1]," be saved?",sep="");for(i in 2:length(re)){res[i]<-paste("Different resolutions for 'mir' spectra. Shall only ",re[i]," be saved?",sep="")};test<-menu(res,graphics=T);for(i in 1:length(re)){if(test==i){bla<-which(mir.ad.tmp[,"Resolution"]==re[i]);mir.ad.tmp<-mir.ad.tmp[bla,]}}}

if(length(summary(as.factor(mir.ad.tmp[,"Zero_filling"])))!=1){ze<-as.integer(names(summary(as.factor(mir.ad.tmp[,"Zero_filling"]))));zer<-c();zer[1]<-paste("Different Zero_fillings for 'mir' spectra. Shall only ",ze[1]," be saved?",sep="");for(i in 2:length(ze)){zer[i]<-paste("Different Zero_fillings. Shall only ",ze[i]," be saved?",sep="")};test<-menu(zer,graphics=T);for(i in 1:length(ze)){if(test==i){bla<-which(mir.ad.tmp[,"Zero_filling"]==ze[i]);mir.ad.tmp<-mir.ad.tmp[bla,]}}}

if(length(summary(as.factor(mir.ad.tmp[,"Wavebands"])))!=1){test<-menu(c("Differing number of wavebands  for 'mir' spectra - continue according to function argument 'wavenumber'.","Differing number of wavebands  for 'mir' spectra - stop."),graphics=T);if(test==2){stop("You stopped reading spc-files.")}}

# Create additional data.frame:

a<-unique(mir.ad.tmp[,"Batch_Labid"]);
b<-list();
for(i in 1:length(a)){
	b[[i]]<-which(mir.ad.tmp[,"Batch_Labid"]==a[i]);
	}

mir.ad<-as.data.frame(matrix(nrow=length(a),ncol=7,dimnames=list(a,c("Batch_Labid","Wavebands","Material","Spectral_range","Scan_date","Resolution","Zero_filling"))));
mir.ad[,"Batch_Labid"]<-unique(mir.ad.tmp[,"Batch_Labid"]);
mir.ad[,"Wavebands"]<-mir.ad.tmp[1,"Wavebands"];
mir.ad[,"Material"]<-mir.ad.tmp[1,"Material"];
mir.ad[,"Spectral_range"]<-mir.ad.tmp[1,"Spectral_range"];
mir.ad[,"Resolution"]<-mir.ad.tmp[1,"Resolution"];
mir.ad[,"Zero_filling"]<-mir.ad.tmp[1,"Zero_filling"];
for(i in 1:length(a)){
if(length(which(mir.ad.tmp[b[[i]],"Scan_date"]==mir.ad.tmp[b[[i]][1],"Scan_date"]))==length(b[[i]])){mir.ad[i,"Scan_date"]<-mir.ad.tmp[b[[i]][1],"Scan_date"]}
}

mir<-mir[as.integer(rownames(mir.ad.tmp))];

test.2<-c(rep(NA,length(tmp[mir])));
for(i in 1:length(tmp[mir])){
	test.2[i]<-substr(tmp[mir][i],1,(test[i]-4))
	}


}

# Check if the spectral are scanned using the same method:

if(length(nir)>0){

# Get the Batch_Labid from each sample name:

test<-c();
for(i in 1:length(tmp[nir])){
	test[i]<-nchar(tmp[nir][i],type="chars")
	}

test.1<-c(rep(NA,length(tmp[nir])));
for(i in 1:length(tmp[nir])){
	test.1[i]<-substr(tmp[nir][i],1,(test[i]-6))
	}

# Create aditional information data frame:

nir.ad<-as.data.frame(matrix(nrow=length(test.1),ncol=7,dimnames=list(test.1,c("Batch_Labid","Wavebands","Material","Spectral_range","Scan_date","Resolution","Zero_filling"))));
nir.ad[,"Batch_Labid"]<-test.1;
nir.ad[,"Spectral_range"]<-c("NIR");

for(i in 1:length(tmp[nir])){
	n<-readRaw(tmp[nir][i],offset=4,nbytes=4,human="int",size=4,endian="little")[[5]];
	nir.ad[i,"Wavebands"]<-n;
Ffirst<-readRaw(tmp[nir][1],offset=8,nbytes=8,human="real",size=8,endian="little")[[5]];
Flast<-readRaw(tmp[nir][1],offset=16,nbytes=8,human="real",size=8,endian="little")[[5]];
if(Ffirst>Flast){c<-Ffirst;Ffirst<-Flast;Flast<-c}
bla<-c(rep(NA,n-2));
	for(k in 1:(n-2)){
		bla[k]<-k*(Flast-Ffirst)/(n-1);
		}
	w<-Ffirst+bla;
	waveb.mpa<-round(c(Ffirst,w,Flast),6);	
	test<-as.character(readRaw(tmp[nir][i],offset=32,nbytes=4,human="int",size=1.5,machine="binary",signed=F,endian="little"));
test.a<-paste(substr(test,35,42),substr(test,26,33),substr(test,17,21),sep="")
year<-substr(test.a,1,12);
c<-nchar(year)
d<-c(rep(NA,c));
k<-0;
for(j in c:1){
	d[j]<-as.numeric(substr(year,j,j))*2^k;
	k<-k+1
	}
year <-sum(d)
month<-substr(test.a,13,16);
c<-nchar(month)
d<-c(rep(NA,c));
k<-0;
for(j in c:1){
	d[j]<-as.numeric(substr(month,j,j))*2^k;
	k<-k+1
	}
month<-sum(d)
day<-substr(test.a,17,21);
c<-nchar(day)
d<-c(rep(NA,c));
k<-0;
for(j in c:1){
	d[j]<-as.numeric(substr(day,j,j))*2^k;
	k<-k+1
	}
day <-sum(d)
nir.ad[i,"Scan_date"]<-paste(month,"/",day,"/",year,sep="")

test<-readRaw(tmp[nir][i],offset=(length(waveb.mpa)*4+544+64),nbytes= 431,human="char",endian="little");
test<-blockString(test);
test<-strsplit(test,"..",fixed=T);
if(grepl("Plant",test[[1]][8])){nir.ad[i,"Material"]<-"Plant"};
if(grepl("Soil",test[[1]][8])){nir.ad[i,"Material"]<-"Soil"};
nir.ad[i,"Resolution"]<-as.integer(substr(test[[1]][7],nchar(test[[1]])[7],nchar(test[[1]])[7]));
nir.ad[i,"Zero_filling"]<-as.integer(substr(test[[1]][20],nchar(test[[1]])[20],nchar(test[[1]])[20]));
	}

# Make checks if parameters are the same for all samples:

if(length(summary(as.factor(nir.ad[,"Material"])))!=1 |  length(summary(as.factor(nir.ad[,"Resolution"])))!=1 | length(summary(as.factor(nir.ad[,"Zero_filling"])))!=1 | length(summary(as.factor(nir.ad[,"Wavebands"])))!=1){test<-menu("Scan settings for nir-spectra not uniform. Find a txt-file with detailed scaning information named 'nir.ad.txt' in 'loa'.",graphics=T);if(test==1){write.table(nir.ad,file="nir.ad.txt",sep=",",dec=".",row.names=FALSE)}}

if(length(summary(as.factor(nir.ad[,"Material"])))!=1){ma<-names(summary(as.factor(nir.ad[,"Material"])));mat<-c();mat[1]<-paste("Different materials for nir-spectra. Shall only ",ma[1]," be saved?",sep="");for(i in 2:length(ma)){mat[i]<-paste("Different materials for nir-spectra. Shall only ",ma[i]," be saved?",sep="")};test<-menu(mat,graphics=T);for(i in 1:length(ma)){if(test==i){bla<-which(nir.ad[,"Material"]==ma[i]);nir.ad<-nir.ad[bla,]}}}


if(length(summary(as.factor(nir.ad[,"Resolution"])))!=1){re<-as.integer(names(summary(as.factor(nir.ad[,"Resolution"]))));res<-c();res[1]<-paste("Different resolutions for nir-spectra. Shall only ",re[1]," be saved?",sep="");for(i in 2:length(re)){res[i]<-paste("Different resolutions. Shall only ",re[i]," be saved?",sep="")};test<-menu(res,graphics=T);for(i in 1:length(re)){if(test==i){bla<-which(nir.ad[,"Resolution"]==re[i]);nir.ad<-nir.ad[bla,]}}}

if(length(summary(as.factor(nir.ad[,"Zero_filling"])))!=1){ze<-as.integer(names(summary(as.factor(nir.ad[,"Zero_filling"]))));zer<-c();zer[1]<-paste("Different Zero_fillings for nir-spectra. Shall only ",ze[1]," be saved?",sep="");for(i in 2:length(ze)){zer[i]<-paste("Different Zero_fillings. Shall only ",ze[i]," be saved?",sep="")};test<-menu(zer,graphics=T);for(i in 1:length(ze)){if(test==i){bla<-which(nir.ad[,"Zero_filling"]==ze[i]);nir.ad<-nir.ad[bla,]}}}

if(length(summary(as.factor(nir.ad[,"Wavebands"])))!=1){test<-menu(c("Differing number of wavebands  for nir-spectra - continue according to function argument 'wavenumber'.","Differing number of wavebands  for nir-spectra - stop."),graphics=T);if(test==2){stop("You stopped reading spc-files.")}}

nir<-nir[match(rownames(nir.ad),test.1)];

}

# Further checks:

if(wavenumber=="first.sample"){
if(length(mir)>0){n<-readRaw(tmp[mir][1],offset=4,nbytes=4,human="int",size=4,endian="little")[[5]];
Ffirst<-readRaw(tmp[mir][1],offset=8,nbytes=8,human="real",size=8,endian="little")[[5]];
Flast<-readRaw(tmp[mir][1],offset=16,nbytes=8,human="real",size=8,endian="little")[[5]];
if(Ffirst>Flast){c<-Ffirst;Ffirst<-Flast;Flast<-c}
bla<-c(rep(NA,n-2));
	for(k in 1:(n-2)){
		bla[k]<-k*(Flast-Ffirst)/(n-1);
		}
	w<-Ffirst+bla;
	waveb.hts<-round(c(Ffirst,w,Flast),6);
	}
if(length(nir)>0){n<-readRaw(tmp[nir][1],offset=4,nbytes=4,human="int",size=4,endian="little")[[5]];
Ffirst<-readRaw(tmp[nir][1],offset=8,nbytes=8,human="real",size=8,endian="little")[[5]];
Flast<-readRaw(tmp[nir][1],offset=16,nbytes=8,human="real",size=8,endian="little")[[5]];
if(Ffirst>Flast){c<-Ffirst;Ffirst<-Flast;Flast<-c}
bla<-c(rep(NA,n-2));
	for(k in 1:(n-2)){
		bla[k]<-k*(Flast-Ffirst)/(n-1);
		}
	w<-Ffirst+bla;
	waveb.mpa<-round(c(Ffirst,w,Flast),6);
	}
	}
	
# Make for MIR:

if(length(mir)>0){

# Create temporary and final mir averaged spectra matrix:

mir.spec<-matrix(nrow=nrow(mir.ad),ncol=length(waveb.hts),dimnames=list(rownames(mir.ad),waveb.hts));
mir.spec.tmp<-matrix(nrow=length(mir),ncol=length(waveb.hts),dimnames=list(test.2,waveb.hts));

# Read spc-file into R:

for(i in 1:length(tmp[mir])){

n<-readRaw(tmp[mir][i],offset=4,nbytes=4,human="int",size=4,endian="little")[[5]];
Ffirst<-readRaw(tmp[mir][i],offset=8,nbytes=8,human="real",size=8,endian="little")[[5]];
Flast<-readRaw(tmp[mir][i],offset=16,nbytes=8,human="real",size=8,endian="little")[[5]];
spec<-readRaw(tmp[mir][i],offset=544,nbytes=(4*n),human="real",size=4,endian="little")[[5]];
if(Ffirst>Flast){c<-Ffirst;Ffirst<-Flast;Flast<-c;spec<-rev(spec)};
# Check if the wavebands have the same first waveband and number of wavebands compared to filemaker:
if(round(Ffirst,6)==waveb.hts[1] & n==length(waveb.hts))
# If yes, read the spectrum:
{	mir.spec.tmp[i,]<-spec;
}
# If the wavebands are different
else
{	# Read and calculate the differeing wavebands:
	
	bla<-c(rep(NA,n-2));
	for(k in 1:(n-2)){
		bla[k]<-k*(Flast-Ffirst)/(n-1);
		}
	w<-Ffirst+bla;
	wa<-round(c(Ffirst,w,Flast),6);
	# Make the spectrum compatible:
	mir.spec.tmp[i,]<-round(spline(wa,spec,method="natural",xout=waveb.hts)[[2]],6);
	# Replace interpolations outside the measured range 		by missing values:
	for(l in 1:length(wa)){
	if((waveb.hts[l]>wa[1])==T)break
	if(waveb.hts[l]<wa[1])
	{
	mir.spec.tmp[i,l]<-NA	
	}}
		for(m in (length(waveb.hts)):(length(waveb.hts)-		length(wa))){
	if((waveb.hts[m]<wa[length(wa)])==T)break;
	if(waveb.hts[m]>wa[length(wa)])
	{
	mir.spec.tmp[i,m]<-NA	
	}}
	
	}
	}
	
# Average spectra:

for(i in 1:length(a)){
	if(length(mir.spec.tmp[b[[i]]])>1){
	mir.spec[i,]<-colMeans(mir.spec.tmp[b[[i]],]);}
	if(length(mir.spec.tmp[b[[i]]])==1){
	mir.spec[i,]<-mir.spec.tmp[b[[i]],]};	
	}
	
# Check if spectra look ok:

if(dim(mir.spec)[1]<1000){f<-1}else{f<-round(dim(mir.spec)[1]/1000,0)+1}
g<-1;
for(i in 1:f){
plot(mir.spec[g,]~waveb.hts,type="l",ylim=c(min(mir.spec,na.rm=T),max(mir.spec,na.rm=T)),xlab="Wavebands (cm^-1)",ylab="Absorption");
	
	if(dim(mir.spec)[1]>1){
	for(j in (g+1):(if((dim(mir.spec)[1]-g)<1000){g+dim(mir.spec)[1]-g}else{(g+999)})){
		lines(mir.spec[j,]~waveb.hts);
		}}
	# Ask for user command y/n:
	test<-menu(c("Spectra ok - please continue","Spectra not ok - please stop"),graphics=T);
	graphics.off();
	if(test==2)stop("You stopped importing those spectra");
	g<-g+1000;
	}
	}

# Make for NIR:

if(length(nir)>0){

# Create nir spectra matrix:

nir.spec<-matrix(nrow=nrow(nir.ad),ncol=length(waveb.mpa),dimnames=list(rownames(nir.ad),waveb.mpa));

# Read spc-file into R:

for(i in 1:length(tmp[nir])){

n<-readRaw(tmp[nir][i],offset=4,nbytes=4,human="int",size=4,endian="little")[[5]];
Ffirst<-readRaw(tmp[nir][i],offset=8,nbytes=8,human="real",size=8,endian="little")[[5]];
Flast<-readRaw(tmp[nir][i],offset=16,nbytes=8,human="real",size=8,endian="little")[[5]];
spec<-readRaw(tmp[nir][i],offset=544,nbytes=(4*n),human="real",size=4,endian="little")[[5]];
if(Ffirst>Flast){c<-Ffirst;Ffirst<-Flast;Flast<-c;spec<-rev(spec)};
# Check if the wavebands have the same first waveband and number of wavebands compared to filemaker:
if(round(Ffirst,6)==waveb.mpa[1] & n==length(waveb.mpa))
# If yes, read the spectrum:
{	nir.spec[i,]<-spec;
}
# If the wavebands are different
else
{	# Read and calculate the differeing wavebands:
	
	bla<-c(rep(NA,n-2));
	for(k in 1:(n-2)){
		bla[k]<-k*(Flast-Ffirst)/(n-1);
		}
	w<-Ffirst+bla;
	wa<-round(c(Ffirst,w,Flast),6);
	# Make the spectrum compatible:
	nir.spec[i,]<-round(spline(wa,spec,method="natural",xout=waveb.mpa)[[2]],6);
	# Replace interpolations outside the measured range 		by missing values:
	for(l in 1:length(wa)){
	if((waveb.mpa[l]>wa[1])==T)break
	if(waveb.mpa[l]<wa[1])
	{
	nir.spec[i,l]<-NA	
	}}
		for(m in (length(waveb.mpa)):(length(waveb.mpa)-		length(wa))){
	if((waveb.mpa[m]<wa[length(wa)])==T)break;
	if(waveb.mpa[m]>wa[length(wa)])
	{
	nir.spec[i,m]<-NA	
	}}
	
	}
	}
		
# Check if spectra look ok:

if(dim(nir.spec)[1]<1000){f<-1}else{f<-round(dim(nir.spec)[1]/1000,0)+1}
g<-1;
for(i in 1:f){
plot(nir.spec[g,]~waveb.mpa,type="l",ylim=c(min(nir.spec,na.rm=T),max(nir.spec,na.rm=T)),xlab="Wavebands (cm^-1)",ylab="Absorption");
	
	if(dim(nir.spec)[1]>1){
	for(j in (g+1):(if((dim(nir.spec)[1]-g)<1000){g+dim(nir.spec)[1]-g}else{(g+999)})){
		lines(nir.spec[j,]~waveb.mpa);
		}}
	# Ask for user command y/n:
	test<-menu(c("Spectra ok - please continue","Spectra not ok - please stop"),graphics=T);
	graphics.off();
	if(test==2)stop("You stopped importing those spectra");
	g<-g+1000;
	}
}

# Prepare function output:

# When mir and nir files were read:

if(length(mir)>0 & length(nir)>0){
	round(nir.spec,6);
	round(mir.spec,6);
	output<-list(mir.spectra=mir.spec,mir.additional.information=mir.ad,nir.spectra=nir.spec,nir.additional.information=nir.ad);
	class(output)<-"read.spc";
	if(sav=="TRUE"){
if(save.path!="NULL"){
		setwd(save.path);}
		if(save.as=="workspace"){save(output,file=output.name)}
		if(save.as=="csv.file"){write.table(mir.spec,file=paste("Mir.",output.name,".csv",sep=""),dec=".",sep=",");write.table(nir.spec,file=paste("Nir.",output.name,".csv",sep=""),dec=".",sep=",")}
		}
return(output);

	}

# When mir files were read:

	
if(length(mir)>0 & length(nir)==0){
	round(mir.spec,6);
	output<-list(spectra=mir.spec,additional.information=mir.ad);
class(output)<-"read.spc";

if(sav=="TRUE"){
if(save.path!="NULL"){
		setwd(save.path);}
		if(save.as=="workspace"){save(output,file=output.name)}
		if(save.as=="csv.file"){write.table(mir.spec,file=paste(output.name,".csv",sep=""),dec=".",sep=",")}
		}
return(output);
	}
		
# When nir files were read:

if(length(nir)>0 & length(mir)==0){
	round(nir.spec,);
	output<-list(spectra=nir.spec,additional.information=nir.ad);
class(output)<-"read.spc";

if(sav=="TRUE"){
if(save.path!="NULL"){
		setwd(save.path);}
		if(save.as=="workspace"){save(output,file=output.name)}
		if(save.as=="csv.file"){write.table(nir.spec,file=paste(output.name,".csv",sep=""),dec=".",sep=",")}
		}
return(output);

	}
		
}

