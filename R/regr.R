#This requires revision it has more parameters that extends beyond the limit of 95 characters line width
regr <-
function(x,y,sw="NL",sp="wt",r="svm",p="T",pn=0.3,n,md="test",drv=1,bandwidth=21,val="CV",ft="h",lv=3,dis="g",nt=1000,s=0.01,k="ra"){
#Set the two global variables: distribution AND validation to NL
distribution<-validation<-NULL

#Set other global variables: 
box.cox.powers<-bc<-MSEP<-gbm.fit<-gbm.perf<-NULL
# Check on function body:

a<-c(rownames(x),rownames(y));
if(nrow(x)==nrow(y) & length(unique(a))!=nrow(y)){
	b<-a[which(duplicated(a)=="T")];
	x<-x[b,];
	y<-y[b,];
	rm(b);
	}
test<-list();
for(i in 1:dim(y)[2]){
	test[[i]]<-which(y[,i]!="NA")
	}
c<-unique(test[[1]]);
if(dim(y)[2]>1){
	for(j in 2:dim(y)[2]){
	c<-unique(c(c,test[[j]]));
		}
	}
y<-y[c,];
x<-x[rownames(y),];
rm(a,test,c,i,j);
if(class(x)!="data.frame" & class(x)!="matrix"){stop("Invalid argument: 'x' must be of class 'data.frame' or 'matrix'.")};
if(class(x)=="data.frame"){x<-as.matrix(x)};
if(class(as.numeric(colnames(x)))!="numeric"){stop("Invalid argument: the colnames of 'x', which should be the waveband positions, are not coercible to class 'numeric'.")}
if(as.numeric(colnames(x)[1])>as.numeric(colnames(x)[2])){
	test<-x;
	for(i in 1:nrow(x)){
		test[i,]<-rev(test[i,]);
		}
		colnames(test)<-rev(colnames(test));
		x<-test;
		rm(test);
	}
if(class(y)!="data.frame" & class(y)!="matrix"){stop("Invalid argument: 'y' must be of class 'data.frame' or 'matrix'.")};
if(class(y)=="matrix"){y<-as.data.frame(y)};
if(sw=="NL"){sw<-getwd()};
if(sp!="raw" & sp!="derivative" & sp!="continuum removed" & sp!="wt"){stop("Invalid argument: 'sp' has to be either 'raw', 'derivative', 'continuum removed' or 'wt'.")};
if(is.na(match(r,c("pls","brt","svm")))){stop("Invalid argument: 'r' has to be either 'pls', 'brt' or 'svm'.")};
if(r=="brt" & nrow(x)<70){stop("Invalid argument: when 'r' is equl to 'brt' the number of samples has to be more than 70.")};
if(class(p)!="character"){stop("Invalid argument: 'p' has to be of class 'character'.")}
if(is.na(match(p,c("T","FALSE")))){stop("Invalid argument: 'p' has to be either 'T' or 'FALSE'.")}
if(class(pn)!="numeric"){stop("Invalid argument: 'p' has to be of class 'numeric'.")};
if(pn<0 | pn>1){stop("Invalid argument: 'p' has to be between 0 and 1.")}
if(p=="T"){
	n<-c();
	for(i in 1:ncol(y)){
	n[i]<-round(pn*length(na.omit(y[,i])),0);
		}
}
if(p=="FALSE"){
	if(exists("n")=="FALSE"){stop("Invalid argument: object 'n' not found.")}
	if(class(as.integer(n))!="integer"){stop("Invalid argument: 'n' has to be of class 'integer'.")};
	if(n<1 | n>=nrow(x)){stop("Invalid argument: 'n' has to be between 1 and the number of samples minus one.")};
	n<-c();
	for(i in 1:ncol(y)){
	n[i]<-n;	
		}
}
if(class(md)!="character"){stop("Invalid argument: 'md' has to be of class 'character'.")}

# Check additional objects for changing other functions (locpoly, pls, brt, svm):

# locpoly function:
if(class(drv)!="integer" & class(drv)!="numeric"){stop("Invalid argument: 'drv' has to be of class 'integer' or 'numeric'.")}
if(class(drv)=="numeric"){drv<-as.integer(drv)};
if(drv<0 | drv > 3){stop("Invalid argument: 'drv' has to be between 0 and 3.")}
if(class(bandwidth)!="integer" & class(bandwidth)!="numeric"){stop("Invalid argument: 'bandwidth' has to be of class 'integer' or 'numeric'.")}
if(class(bandwidth)=="numeric"){bandwidth <-as.integer(bandwidth)};
if(bandwidth<1 | bandwidth > 30){stop("Invalid argument: 'bandwidth' has to be between 1 and 30.")}
# mvr function:
if(class(validation)!="character"){stop("Invalid argument: 'validation' has to be of class 'character'.")}
if(validation!="CV" & validation!="none" & validation!="LOO"){stop("Invalid argument: 'validation' has to be either 'none', 'CV' or 'LOO'.")}
# dwt function:
if(class(ft)!="character"){stop("Invalid argument: 'ft' has to be of class 'character'.")}
if(ft!="h" &ft!="d4" &ft!="d6" &ft!="d8" &ft!="d10" &ft!="d12" &ft!="d14" &ft!="d16" &ft!="d18" &ft!="d20" &ft!="la8" &ft!="la10" &ft!="la12" &ft!="la14" &ft!="la16" &ft!="la18" &ft!="la20" &ft!="bl14" &ft!="bl18" &ft!="bl20" &ft!="c6" &ft!="c12" &ft!="c18" &ft!="c24" &ft!="c30"){stop("Invalid argument: ft' has to be either 'h', 'd4', 'd6', 'd8', 'd10', 'd12', 'd14', 'd16', 'd18', 'd20', 'la8', 'la10', 'la12', 'la14', 'la16', 'la18', 'la20', 'bl14', 'bl18', 'bl20', 'c6', 'c12', 'c18', 'c24', or 'c30'.")}
if(class(lv)!="numeric"){stop("Invalid argument: 'lv' has to be of class 'numeric'.")}
if(is.na(match(lv,c(1:9)))){stop("Invalid argument: 'lv' has to be between 1 and 9.")}
slo<-10-lv;
# gmb.fit:
if(class(distribution)!="character"){stop("Invalid argument: 'distribution' has to be of class 'character'.")}
if(is.na(match(distribution,c("g","laplace","bernoulli","adaboost","poisson","coxph")))){stop("Invalid argument: 'distribution' has to be either 'g', 'laplace', 'bernoulli', 'adaboost', 'poisson' or 'coxph'.")}
if(class(nt)!="numeric" & class(nt)!="integer"){stop("Invalid argument: 'nt' has to be of class 'numeric' or 'integer'.")}
if(class(nt)=="numeric"){nt <-as.integer(nt)};
if(class(s)!="numeric"){stop("Invalid argument: 's' has to be of class 'numeric'.")}
# svm:
if(class(k)!="character"){stop("Invalid argument: 'k' has to be of class 'character'.")}
if(is.na(match(k,c("polynomial","linear","sigmoid","ra")))){stop("Invalid argument: 'k' has to be either 'polynomial', 'linear','sigmoid' or 'ra'.")}

# Spectral transformation:

if(sp=="raw"){
	x.tr<-x;
	}
	
if(sp=="derivative"){
x.tr<-matrix(nrow=nrow(x),ncol=ncol(x),dimnames=list(rownames(x),colnames(x)));
waveb<-as.numeric(colnames(x));
#require(KernSmooth,quietly=T);
for(i in 1:nrow(x)){
	x.tr[i,]<-locpoly(waveb,x[i,],drv=drv,bandwidth=bandwidth,gridsize=ncol(x))[[2]]
	}
#detach(package:KernSmooth);
	}
	
if(sp=="continuum removed"){
x.tr<-matrix(nrow=nrow(x),ncol=ncol(x),dimnames=list(rownames(x),colnames(x)));
waveb<-as.numeric(colnames(x));
#require(KernSmooth,quietly=T);	
test<-x;
for(i in 1:nrow(x)){
	test.1<-cbind(waveb,test[i,]);
	test.1<-sortedXyData(test.1[,1],test.1[,2]);
	ch<-chull(test.1);
	ch.1<-ch;
	ch<-ch[1:(which(ch==1))]
	ch<-sort(ch);
	ch<-c(ch,ncol(x));
	appr.ch<-approx(test.1[ch,],xout=test.1[,1],method="linear",ties="mean");
	cr<-test.1[[2]]-appr.ch[[2]];
	x.tr[i,]<-cr;
	}
x.tr<-x.tr[,2:(ncol(x)-2)];
#detach(package:KernSmooth);
	}
	
if(sp=="wt"){
waveb<-as.numeric(colnames(x));
waveb.1024.up<-round(max(waveb));
waveb.1024.down<-round(min(waveb));
waveb.1024.n<-1023;
waveb.1024.step<-(waveb.1024.up-waveb.1024.down)/waveb.1024.n;
waveb.1024<-c();
waveb.1024[1]<-waveb.1024.down;
for(i in 2:1024){
	waveb.1024[i]<-round(waveb.1024.down+(i-1)*waveb.1024.step,5)
	}
x.comp<-matrix(nrow=nrow(x),ncol=length(waveb.1024),dimnames=list(rownames(x),waveb.1024));
for(i in 1:nrow(x)){
	x.comp[i,]<-round(spline(waveb,x[i,],method="natural",xout=waveb.1024)[[2]],6);
}
#require(wavelets,quietly=T);
x.tr<-matrix(nrow=nrow(x.comp),ncol=2^slo,dimnames=list(rownames(x.comp),paste("WC_",c(1:2^slo),sep="")));
for(i in 1:nrow(x.tr)){
	blub<-dwt(x.comp[i,],filter=ft);
	x.tr[i,]<-slot(blub,"W")[[lv]]
	}
#detach(package:wavelets);
	}

# Reference data transformation:

# Transform yerence values:

y.n.n<-colnames(y);
y.bc.n<-paste("bc.",colnames(y),sep="");
y.sqrt.n<-paste("sqrt.",colnames(y),sep="");
y.log.n<-paste("log.",colnames(y),sep="");

y.bc<-matrix(nrow=nrow(y),ncol=length(y.bc.n),dimnames=list(rownames(y),y.bc.n));
y.sqrt<-matrix(nrow=nrow(y),ncol=length(y.sqrt.n),dimnames=list(rownames(y),y.sqrt.n));
y.log<-matrix(nrow=nrow(y),ncol=length(y.log.n),dimnames=list(rownames(y),y.log.n));

#require(car,quietly=T);
lambda<-c();
y.1<-y+1;
for(i in 1:ncol(y.1)){
lambda[i]<-box.cox.powers(na.omit(y.1[,i]))[[1]];
y.bc[,i]<-bc(y.1[,i],lambda[i]);
	};
y.bc<-as.data.frame(y.bc);
#detach(package:car);
y.sqrt<-sqrt(y);
y.log<-log(y+1);

tr<-c("Untransformed","Box-Cox transformed","Square root transformed","Log transformed");
tr.n<-c();
for(i in 1:ncol(y)){
	dev.new(width=8,height=8);
	par(mfrow=c(2,2));
	plot(density(na.omit(y[,i])),main=paste("Untransformed",y.n.n[i]));
	plot(density(na.omit(y.bc[,i])),main=paste("Box-Cox transformed",y.n.n[i]));
	plot(density(na.omit(y.sqrt[,i])),main=paste("Square root transformed",y.n.n[i]));
	plot(density(na.omit(y.log[,i])),main=paste("Log transformed",y.n.n[i]));
	# Ask for user command 1 to 4:
	test<-menu(c("Untransformed","Box-Cox transformed","Square root transformed","Log transformed"),graphics=T,title="Which graph represent the closest normal data distribution?");
	graphics.off();
	if(test==1)tr.n[i]<-tr[1];
	if(test==2)tr.n[i]<-tr[2];
	if(test==3)tr.n[i]<-tr[3];
	if(test==4)tr.n[i]<-tr[4];
	}

test<-list();
test.n<-c();
for(i in 1:length(tr.n)){
	if(tr.n[i]==tr[1]){test[[i]]<-y[,i];test.n[i]<-y.n.n[i]};
	if(tr.n[i]==tr[2]){test[[i]]<-y.bc[,i];test.n[i]<-y.bc.n[i]};
	if(tr.n[i]==tr[3]){test[[i]]<-y.sqrt[,i];test.n[i]<-y.sqrt.n[i]};
	if(tr.n[i]==tr[4]){test[[i]]<-y.log[,i];test.n[i]<-y.log.n[i]};
	}
y.tr<-matrix(nrow=nrow(y),ncol=ncol(y),dimnames=list(rownames(y),test.n));
for(i in 1:length(test)){
	y.tr[,i]<-test[[i]];
	}
y.tr<-as.data.frame(y.tr);

# Name lambda according to its use:

names(lambda)<-y.n.n;
a<-which(substr(colnames(y.tr),1,2)=="bc");
lambda<-lambda[a];

# Sample selection:

# Get validation samples for all soil properties:

cal<-list();
cal.n<-list();
val<-list();
val.n<-list();
for(i in 1:dim(y)[2]){

# Make a Principal component analysis and get the important principal components:

bl<-rownames(y.tr[which(y.tr[,colnames(y.tr)[i]]!="NA"),]);
pca<-prcomp(x.tr[bl,],scale=T);
prco<-pca$x[,1:20];
cpv<-summary(pca)[[6]][3,1:20];
zzz<-matrix(nrow=1,ncol=length(cpv)-4);
for (j in 1:16){
	e<-(cpv[j]+0.04)<cpv[j+3];
	zzz[j]<-e
	}
pc<-(which(zzz==FALSE)-1)[1];
if(pc==1){pc<-2};
prco<-prco[,1:pc];	


# Get the number of validation samples:

va.n<-n[i];

# Get the most extreme points (min and max) for each important principal component:

ca.min<-c(rep(1,ncol(prco)));
ca.max<-c(rep(1,ncol(prco)));
for (j in 1:ncol(prco)){
	blub<-which(prco[,j]==min(prco[,j]));
	ca.min[j]<-blub[1];
	bla<-which(prco[,j]==max(prco[,j]));
	ca.max[j]<-bla[1];
	}
ca.min<-rownames(prco)[ca.min];
ca.max<-rownames(prco)[ca.max];
ca.start<-unique(c(ca.min,ca.max));
ca.start.n<-match(ca.start,rownames(prco));

# Get the second most extreme points (min and max) for each PC and assign them for the val set:

va.min<-c(rep(1,ncol(prco)));
va.max<-c(rep(1,ncol(prco)));
for (j in 1:ncol(prco)){
	blub<-which(prco[-ca.start.n,j]==min(prco[-ca.start.n,j]));
	va.min[j]<-blub[sample(length(blub),1)];
	bla<-which(prco[-ca.start.n,j]==max(prco[-ca.start.n,j]));
	va.max[j]<-bla[sample(length(bla),1)];
	}
va.min<-rownames(prco[-ca.start.n,])[va.min];
va.max<-rownames(prco[-ca.start.n,])[va.max];
va.start<-unique(c(va.min,va.max));
va.start.n<-match(va.start,rownames(prco));
ca.va<-c(ca.start,va.start);
ca.va.start<-match(c(ca.start,va.start),rownames(prco));

# Calculate the Euclidean distance matrix:

euc<-as.matrix(dist(prco));

# Get the remaining validation samples up to the desired number and retrieve recursive the calibration samples:

start<-rownames(prco)[-ca.va.start];
start.b<-start
va<-va.start;
for(k in 1:(va.n-length(va.start))){
test<-apply(euc[start.b,va],1,min);
bla<-names(which(test==max(test)));
va<-c(va,bla);
start.b<-start.b[-(which(match(start.b,bla)!="NA"))];
}
va.n<-match(va,rownames(prco));
ca.n<-c(1:nrow(prco))[-va.n];
ca<-rownames(prco)[ca.n];
cal[[i]]<-ca;
cal.n[[i]]<-ca.n;
val[[i]]<-va;
val.n[[i]]<-va.n;
}

# Model calculation:

# Prepare needed objects:

im<-c();
for(i in 1:length(tr.n)){
	im[i]<-which(tr.n[i]==tr);
	}
a<-unique(c(y.n.n,colnames(y.tr)));
b<-paste("pred.",a,sep="");
variables.all<-c(c(a,b));
variables<-colnames(y);

# Create statistics ojbect for cal-set:

cal.stat.variables<-c("n","r2","a","bias","RMSEC","RPD","n LV","n bc out","n trees");
cal.stat<-matrix(nrow=length(variables),ncol=length(cal.stat.variables),dimnames=list(variables,cal.stat.variables));
cal.n.bl<-cal[[1]];
for(i in 2:length(variables)){
	cal.n.bl<-unique(c(cal.n.bl,cal[[i]]));
	}
cal.mp<-matrix(nrow=length(cal.n.bl),ncol=(length(variables.all)),dimnames=list(cal.n.bl,variables.all));
for(i in 1:length(variables)){
	cal.mp[cal[[i]],colnames(y)[i]]<-y[cal[[i]],colnames(y)[i]];
	cal.mp[cal[[i]],colnames(y.tr)[i]]<-y.tr[cal[[i]],colnames(y.tr)[i]];
	}
cal.mp<-as.data.frame(cal.mp);

# Create statistics ojbect for the val set:

val.stat.variables<-c("n","r2","a","bias","RMSEP","RPD","n bc out");
val.stat<-matrix(nrow=length(variables),ncol=length(val.stat.variables),dimnames=list(variables,val.stat.variables));
val.n.bl<-val[[1]];
for(i in 2:length(variables)){
	val.n.bl<-unique(c(val.n.bl,val[[i]]));
	}
val.mp<-matrix(nrow=length(val.n.bl),ncol=(length(variables.all)),dimnames=list(val.n.bl,variables.all));
for(i in 1:length(variables)){
	val.mp[val[[i]],colnames(y)[i]]<-y[val[[i]],colnames(y)[i]];
	val.mp[val[[i]],colnames(y.tr)[i]]<-y.tr[val[[i]],colnames(y.tr)[i]];
	}
val.mp<-as.data.frame(val.mp);

if(r=="pls"){

# PLS:

#require(pls,quietly=T);

model<-list();
bc.cal<-list();
bc.val<-list();
bla<-c();
for(j in 1:ncol(y.tr)){
	bla[j]<-length(na.omit(y.tr[cal[[j]],j]))
	}
	
for(i in 1:length(variables)){

# Set F-distribution alpha value and calculate according F-limit:

if(bla[i]<=200){alpha<-0.95};
if(bla[i]>200){alpha<-0.975};
f.value<-qf(alpha,bla[i],bla[i]);

cal.pls<-mvr(cal.mp[cal[[i]],colnames(y.tr)[i]] ~x.tr[cal[[i]],],method="kernelpls",validation=validation,ncomp=20);
press.cal.pls<-MSEP(cal.pls)[[1]][1,1,2:21] *length(cal[[i]]);
min.cal.pls<-which(press.cal.pls==min(press.cal.pls));
fh.cal.pls<-c();
for(j in 1:min.cal.pls){
	fh.cal.pls[j]<-press.cal.pls[j]/press.cal.pls[min.cal.pls];
	}
cal.pls.lv<-which(fh.cal.pls<f.value)[1];
cal.pls$ncomp<-cal.pls.lv
# [length(cal[[i]])]
# Calculate cal-subset statistics:
if(im[i]==1){
	cal.mp[cal[[i]],paste("pred.",variables[i],sep="")]<-predict(cal.pls,newdata=x.tr[cal[[i]],],ncomp=cal.pls.lv);
	}
if(im[i]==2){
	cal.mp[cal[[i]],paste("pred.",colnames(y.tr)[i],sep="")]<-predict(cal.pls,newdata=x.tr[cal[[i]],],ncomp=cal.pls.lv);
	cal.mp[cal[[i]],paste("pred.",variables[i],sep="")]<-(cal.mp[cal[[i]],paste("pred.",colnames(y.tr)[i],sep="")]*lambda[variables[i]]+1)^(1/lambda[variables[i]])-1;
	}
if(im[i]==3){
	cal.mp[cal[[i]],paste("pred.",colnames(y.tr)[i],sep="")]<-predict(cal.pls,newdata=x.tr[cal[[i]],],ncomp=cal.pls.lv);
	cal.mp[cal[[i]],paste("pred.",variables[i],sep="")]<-cal.mp[cal[[i]],paste("pred.",colnames(y.tr)[i],sep="")]^2	}
if(im[i]==4){
	cal.mp[cal[[i]],paste("pred.",colnames(y.tr)[i],sep="")]<-predict(cal.pls,newdata=x.tr[cal[[i]],],ncomp=cal.pls.lv);
	cal.mp[cal[[i]],paste("pred.",variables[i],sep="")]<-exp(cal.mp[cal[[i]],paste("pred.",colnames(y.tr)[i],sep="")])-1	}
bc.cal[[i]]<-which(cal.mp[cal[[i]],paste("pred.",variables[i],sep="")]=="NaN");

model[[i]]<-cal.pls
model[[i]]$constituent<-variables[i];
# n LV:
cal.stat[variables[i],"n LV"]<-cal.pls.lv
# Calculate val-subset statistics:
if(im[i]==1){
	val.mp[val[[i]],paste("pred.",variables[i],sep="")]<-predict(cal.pls,newdata=x.tr[val[[i]],],ncomp=cal.pls.lv);
	}
if(im[i]==2){
	val.mp[val[[i]],paste("pred.",colnames(y.tr)[i],sep="")]<-predict(cal.pls,newdata=x.tr[val[[i]],],ncomp=cal.pls.lv);
	val.mp[val[[i]],paste("pred.",variables[i],sep="")]<-(val.mp[val[[i]],paste("pred.",colnames(y.tr)[i],sep="")]*lambda[variables[i]]+1)^(1/lambda[variables[i]])-1;
	}
if(im[i]==3){
	val.mp[val[[i]],paste("pred.",colnames(y.tr)[i],sep="")]<-predict(cal.pls,newdata=x.tr[val[[i]],],ncomp=cal.pls.lv);
	val.mp[val[[i]],paste("pred.",variables[i],sep="")]<-val.mp[val[[i]],paste("pred.",colnames(y.tr)[i],sep="")]^2	}
if(im[i]==4){
	val.mp[val[[i]],paste("pred.",colnames(y.tr)[i],sep="")]<-predict(cal.pls,newdata=x.tr[val[[i]],],ncomp=cal.pls.lv);
	val.mp[val[[i]],paste("pred.",variables[i],sep="")]<-exp(val.mp[val[[i]],paste("pred.",colnames(y.tr)[i],sep="")])-1	}
bc.val[[i]]<-which(val.mp[val[[i]],paste("pred.",variables[i],sep="")]=="NaN");
}
#detach(package:pls);
	}

if(r=="brt"){

# GBM:

#require(splines,quietly=T);
#require(survival,quietly=T);
#require(lattice,quietly=T);
#require(gbm,quietly=T);

model<-list();
bc.cal<-list();
bc.val<-list();

for(i in 1:length(variables)){

cal.gbm<-gbm.fit(x.tr[cal[[i]],],cal.mp[cal[[i]],colnames(y.tr)[i]],distribution = distribution,s=s,nt=nt,verbose=FALSE);
cal.stat[variables[i],"n trees"]<-gbm.perf(cal.gbm,method="OOB",plot.it=FALSE)
# Calculate cal-subset statistics:
if(im[i]==1){
	cal.mp[cal[[i]],paste("pred.",variables[i],sep="")]<-predict(cal.gbm,newdata=x.tr[cal[[i]],],nt=cal.stat[variables[i],"n trees"],type="response");
	}
if(im[i]==2){
	cal.mp[cal[[i]],paste("pred.",colnames(cal.mp)[i],sep="")]<-predict(cal.gbm,newdata=x.tr[cal[[i]],],nt=cal.stat[variables[i],"n trees"],type="response");
	cal.mp[cal[[i]],paste("pred.",variables[i],sep="")]<-(cal.mp[cal[[i]],paste("pred.",colnames(cal.mp)[i],sep="")]*lambda[variables[i]]+1)^(1/lambda[variables[i]])-1;
	}
if(im[i]==3){
	cal.mp[cal[[i]],paste("pred.",colnames(cal.mp)[i],sep="")]<-predict(cal.gbm,newdata=x.tr[cal[[i]],],nt=cal.stat[variables[i],"n trees"],type="response");
	cal.mp[cal[[i]],paste("pred.",variables[i],sep="")]<-cal.mp[cal[[i]],paste("pred.",colnames(cal.mp)[i],sep="")]^2	}
if(im[i]==4){
	cal.mp[cal[[i]],paste("pred.",colnames(cal.mp)[i],sep="")]<-predict(cal.gbm,newdata=x.tr[cal[[i]],],nt=cal.stat[variables[i],"n trees"],type="response");
	cal.mp[cal[[i]],paste("pred.",variables[i],sep="")]<-exp(cal.mp[cal[[i]],paste("pred.",colnames(cal.mp)[i],sep="")])-1	}
bc.cal[[i]]<-which(cal.mp[cal[[i]],paste("pred.",variables[i],sep="")]=="NaN");

model[[i]]<-cal.gbm
model[[i]]$constituent<-variables[i];
# Calculate val-subset statistics:
if(im[i]==1){
	val.mp[val[[i]],paste("pred.",variables[i],sep="")]<-predict(cal.gbm,newdata=x.tr[val[[i]],],nt=cal.stat[variables[i],"n trees"],type="response");
	}
if(im[i]==2){
	val.mp[val[[i]],paste("pred.",colnames(cal.mp)[i],sep="")]<-predict(cal.gbm,newdata=x.tr[val[[i]],],nt=cal.stat[variables[i],"n trees"],type="response");
	val.mp[val[[i]],paste("pred.",variables[i],sep="")]<-(val.mp[val[[i]],paste("pred.",colnames(cal.mp)[i],sep="")]*lambda[variables[i]]+1)^(1/lambda[variables[i]])-1;
	}
if(im[i]==3){
	val.mp[val[[i]],paste("pred.",colnames(cal.mp)[i],sep="")]<-predict(cal.gbm,newdata=x.tr[val[[i]],],nt=cal.stat[variables[i],"n trees"],type="response");
	val.mp[val[[i]],paste("pred.",variables[i],sep="")]<-val.mp[val[[i]],paste("pred.",colnames(cal.mp)[i],sep="")]^2	}
if(im[i]==4){
	val.mp[val[[i]],paste("pred.",colnames(cal.mp)[i],sep="")]<-predict(cal.gbm,newdata=x.tr[val[[i]],],nt=cal.stat[variables[i],"n trees"],type="response");
	val.mp[val[[i]],paste("pred.",variables[i],sep="")]<-exp(val.mp[val[[i]],paste("pred.",colnames(cal.mp)[i],sep="")])-1	}
bc.val[[i]]<-which(val.mp[val[[i]],paste("pred.",variables[i],sep="")]=="NaN");
}
#detach(package:gbm);
#detach(package:splines);
#detach(package:survival);
#detach(package:lattice);
	}
	
if(r=="svm"){
#require(class,quietly=T);
#require(e1071,quietly=T);

model<-list();
bc.cal<-list();
bc.val<-list();

for(i in 1:length(variables)){

test<-tune.svm(type="eps",kernel=k,x.tr[cal[[i]],],cal.mp[cal[[i]],colnames(y.tr)[i]],cost=c(2,3,4,5,6),epsilon=c(0.05,0.1,0.2),gamma=c(0.01,0.001,0.0001,0.00001),tunecontrol=tune.control(sampling="fix",fix=0.7));
test<-tune.svm(type="eps",kernel=k,x.tr[cal[[i]],],cal.mp[cal[[i]],colnames(y.tr)[i]],cost=test[[1]][2][1,1],epsilon=test[[1]][3][1,1],gamma=c(test[[1]][1][1,1]-test[[1]][1][1,1]*0.6,test[[1]][1][1,1]-test[[1]][1][1,1]*0.4,test[[1]][1][1,1]-test[[1]][1][1,1]*0.2,test[[1]][1][1,1],test[[1]][1][1,1]+test[[1]][1][1,1]*2,test[[1]][1][1,1]+test[[1]][1][1,1]*4,test[[1]][1][1,1]+test[[1]][1][1,1]*6),tunecontrol=tune.control(sampling="fix",fix=0.7));
c<-test[[1]][2][1,1];
e<-test[[1]][3][1,1];
g<-test[[1]][1][1,1];
cal.svm<-svm(x.tr[cal[[i]],],cal.mp[cal[[i]],colnames(y.tr)[i]],type="eps",kernel=k,cost=c,epsilon=e,gamma=g);
# Calculate cal-subset statistics:
if(im[i]==1){
	cal.mp[cal[[i]],paste("pred.",variables[i],sep="")]<-predict(cal.svm,newdata=x.tr[cal[[i]],]);
	}
if(im[i]==2){
	cal.mp[cal[[i]],paste("pred.",colnames(y.tr)[i],sep="")]<-predict(cal.svm,newdata=x.tr[cal[[i]],]);
	cal.mp[cal[[i]],paste("pred.",variables[i],sep="")]<-(cal.mp[cal[[i]],paste("pred.",colnames(y.tr)[i],sep="")]*lambda[variables[i]]+1)^(1/lambda[variables[i]])-1;
	}
if(im[i]==3){
	cal.mp[cal[[i]],paste("pred.",colnames(y.tr)[i],sep="")]<-predict(cal.svm,newdata=x.tr[cal[[i]],]);
	cal.mp[cal[[i]],paste("pred.",variables[i],sep="")]<-cal.mp[cal[[i]],paste("pred.",colnames(y.tr)[i],sep="")]^2	}
if(im[i]==4){
	cal.mp[cal[[i]],paste("pred.",colnames(y.tr)[i],sep="")]<-predict(cal.svm,newdata=x.tr[cal[[i]],]);
	cal.mp[cal[[i]],paste("pred.",variables[i],sep="")]<-exp(cal.mp[cal[[i]],paste("pred.",colnames(y.tr)[i],sep="")])-1	}
bc.cal[[i]]<-which(cal.mp[cal[[i]],paste("pred.",variables[i],sep="")]=="NaN");
model[[i]]<-cal.svm
model[[i]]$constituent<-variables[i];
# Calculate val-subset statistics:
if(im[i]==1){
	val.mp[val[[i]],paste("pred.",variables[i],sep="")]<-predict(cal.svm,newdata=x.tr[val[[i]],]);
	}
if(im[i]==2){
	val.mp[val[[i]],paste("pred.",colnames(y.tr)[i],sep="")]<-predict(cal.svm,newdata=x.tr[val[[i]],]);
	val.mp[val[[i]],paste("pred.",variables[i],sep="")]<-(val.mp[val[[i]],paste("pred.",colnames(y.tr)[i],sep="")]*lambda[variables[i]]+1)^(1/lambda[variables[i]])-1;
	}
if(im[i]==3){
	val.mp[val[[i]],paste("pred.",colnames(y.tr)[i],sep="")]<-predict(cal.svm,newdata=x.tr[val[[i]],]);
	val.mp[val[[i]],paste("pred.",variables[i],sep="")]<-val.mp[val[[i]],paste("pred.",colnames(y.tr)[i],sep="")]^2	}
if(im[i]==4){
	val.mp[val[[i]],paste("pred.",colnames(y.tr)[i],sep="")]<-predict(cal.svm,newdata=x.tr[val[[i]],]);
	val.mp[val[[i]],paste("pred.",variables[i],sep="")]<-exp(val.mp[val[[i]],paste("pred.",colnames(y.tr)[i],sep="")])-1	}
bc.val[[i]]<-which(val.mp[val[[i]],paste("pred.",variables[i],sep="")]=="NaN");
}
#detach(package:e1071);
#detach(package:class);
	}

# Compute statistics independend of regression method:

for(i in 1:length(variables)){

# n bc out:
cal.stat[variables[i],"n bc out"]<-length(bc.cal[[i]])
# n:
cal.stat[variables[i],"n"]<-length(na.omit(cal.mp[,paste("pred.",variables[i],sep="")]));
if(length(bc.cal[[i]])!=0){
# R2:
cal.stat[variables[i],"r2"]<-
round(cor(cal.mp[cal[[i]][-bc.cal[[i]]],variables[i]],cal.mp[cal[[i]][-bc.cal[[i]]],paste("pred.",variables[i],sep="")])^2,2);
# a:
cal.stat[variables[i],"a"]<-round(lm(cal.mp[cal[[i]],variables[i]]~cal.mp[cal[[i]],paste("pred.",variables[i],sep="")])$coefficients[2],2);
# Bias:
residual<-na.omit(cal.mp[,paste("pred.",variables[i],sep="")]-cal.mp[,variables[i]]);
cal.stat[variables[i],"bias"]<-round(sum(residual)/length(cal.mp[,variables[i]]),2)
# RMSEC:
cal.stat[variables[i],"RMSEC"]<-
round(sqrt(sum((cal.mp[cal[[i]][-bc.cal[[i]]],paste("pred.",variables[i],sep="")]-cal.mp[cal[[i]][-bc.cal[[i]]],variables[i]])^2)/(length(cal.mp[cal[[i]][-bc.cal[[i]]],variables[i]]))),2);
# RPD:
cal.stat[variables[i],"RPD"]<-round(sd(cal.mp[cal[[i]],variables[i]])/cal.stat[variables[i],"RMSEC"],2);
}
if(length(bc.cal[[i]])==0){
# R2:
cal.stat[variables[i],"r2"]<-
round(cor(cal.mp[cal[[i]],variables[i]],cal.mp[cal[[i]],paste("pred.",variables[i],sep="")])^2,2);
# a:
cal.stat[variables[i],"a"]<-round(lm(cal.mp[cal[[i]],variables[i]]~cal.mp[cal[[i]],paste("pred.",variables[i],sep="")])$coefficients[2],2);
# Bias:
residual<-na.omit(cal.mp[,paste("pred.",variables[i],sep="")]-cal.mp[,variables[i]]);
cal.stat[variables[i],"bias"]<-round(sum(residual)/length(cal.mp[,variables[i]]),2)
# RMSEC:
cal.stat[variables[i],"RMSEC"]<-
round(sqrt(sum((cal.mp[cal[[i]],paste("pred.",variables[i],sep="")]-cal.mp[cal[[i]],variables[i]])^2)/(length(cal.mp[cal[[i]],variables[i]]))),2);
# RPD:
cal.stat[variables[i],"RPD"]<-round(sd(cal.mp[cal[[i]],variables[i]])/cal.stat[variables[i],"RMSEC"],2);
}

# n bc out:
val.stat[variables[i],"n bc out"]<-length(bc.val[[i]])
# n:
val.stat[variables[i],"n"]<-length(na.omit(val.mp[,paste("pred.",variables[i],sep="")]));
if
(length(bc.val[[i]])!=0)
{
# R2:
val.stat[variables[i],"r2"]<-round(cor(val.mp[val[[i]][-bc.val[[i]]],variables[i]],val.mp[val[[i]][-bc.val[[i]]],paste("pred.",variables[i],sep="")])^2,2);
# a:
val.stat[variables[i],"a"]<-round(lm(val.mp[val[[i]],variables[i]]~val.mp[val[[i]],paste("pred.",variables[i],sep="")])$coefficients[2],2);
# Bias:
residual<-na.omit(val.mp[,paste("pred.",variables[i],sep="")]-val.mp[,variables[i]]);
val.stat[variables[i],"bias"]<-round(sum(residual)/length(val.mp[,variables[i]]),2)
# RMSEP:
val.stat[variables[i],"RMSEP"]<-round(sqrt(sum((val.mp[val[[i]][-bc.val[[i]]],paste("pred.",variables[i],sep="")]-val.mp[val[[i]][-bc.val[[i]]],variables[i]])^2)/(length(val.mp[val[[i]][-bc.val[[i]]],variables[i]]))),2);
# RPD:
val.stat[variables[i],"RPD"]<-round(sd(val.mp[val[[i]],variables[i]])/val.stat[variables[i],"RMSEP"],2);
}
if(length(bc.val[[i]])==0){
# R2:
val.stat[variables[i],"r2"]<-round(cor(val.mp[val[[i]],variables[i]],val.mp[val[[i]],paste("pred.",variables[i],sep="")])^2,2);
# a:
val.stat[variables[i],"a"]<-round(lm(val.mp[val[[i]],variables[i]]~val.mp[val[[i]],paste("pred.",variables[i],sep="")])$coefficients[2],2);
# Bias:
residual<-na.omit(val.mp[,paste("pred.",variables[i],sep="")]-val.mp[,variables[i]]);
val.stat[variables[i],"bias"]<-round(sum(residual)/length(val.mp[,variables[i]]),2)
# RMSEP:
val.stat[variables[i],"RMSEP"]<-round(sqrt(sum((val.mp[val[[i]],paste("pred.",variables[i],sep="")]-val.mp[val[[i]],variables[i]])^2)/(length(val.mp[val[[i]],variables[i]]))),2);
# RPD:
val.stat[variables[i],"RPD"]<-round(sd(val.mp[val[[i]],variables[i]])/val.stat[variables[i],"RMSEP"],2);
}
}	

# Create objects for prediction:

# Calculate the pca-space and mahalanobis-space of the calibration sets:

pca.l<-list();
for(i in 1:length(variables)){
	pca.l[[i]]<-prcomp(x.tr[cal[[i]],],scale.=T);
		}

# Calculate the mahalanobis-space of the calibration sets:

mah.l<-list();
if(dim(x.tr[cal[[i]],])[1]>=dim(x.tr[cal[[i]],])[2])
for(i in 1:length(variables)){
	if(dim(x.tr[cal[[i]],])[1]>=dim(x.tr[cal[[i]],])[2])
{mah.l[[i]]<-sqrt(mahalanobis(pca.l[[i]]$x,colMeans(pca.l[[i]]$x),cov(pca.l[[i]]$x)))}
	if(dim(x.tr[cal[[i]],])[1]<dim(x.tr[cal[[i]],])[2]){
		mah.l[[i]]<-NA
		}
		}

# Get the range of the calibration set:

ran<-list();
for(i in 1:length(variables)){
	ran[[i]]<-range(y[cal[[i]],variables[i]])
	}

# Get the rmsep and the sorted linear regression models  of the validation set:

rmsep.l<-list();
lm.l<-list();
for(i in 1:length(variables)){
	n<-length(na.omit(val.mp[,paste("pred.",variables[i],sep="")]));
	if(n>=250){d<-50};
	if(n<250){d<-round(n*0.2,-1)};
	rmsep.1<-c(rep(NA,n-d));
	if(length(bc.val[[i]])!=0){	p<-which(val.mp[,paste("pred.",variables[i],sep="")]!="NA")[-bc.val[[i]]]}
	if(length(bc.val[[i]])==0){p<-which(val.mp[,paste("pred.",variables[i],sep="")]!="NA")}
	b<-val.mp[p,paste("pred.",variables[i],sep="")];
	a<-sort(b,decreasing=F,index.return=T)$ix;
	for(j in 1:(n-d)){
b.1<-sort(b,decreasing=F,index.return=T)$ix[j:(j+d)];
rmsep.1[j]<-round(sqrt(median((val.mp[p[b.1],paste("pred.",variables[i],sep="")]-val.mp[p[b.1],variables[i]])^2)/0.4549364),2)
};
	rmsep<-c(rep(NA,n));
	rmsep[1:round(d/2,0)]<-rmsep.1[1];
	rmsep[(round(d/2,0)+1):(length(rmsep.1)+(d-round(d/2,0)))]<-rmsep.1;
	rmsep[(n-(d-round(d/2,0)-1)):n]<-rmsep.1[length(rmsep.1)]
	rmsep.l[[i]]<-rmsep;
	lm.l[[i]]<-lm(val.mp[,variables[i]]~val.mp[,paste("pred.",variables[i],sep="")])$fitted.values[a];
	}

# Plot linear regression of measured against predicted of calibration and validation sets:

dev.new(width=5,height=2.5*length(variables));
par(mfrow=c(length(variables),2));
for(i in 1:length(variables)){
mi<-min(na.omit(c(cal.mp[,variables[i]],cal.mp[,paste("pred.",variables[i],sep="")],val.mp[,variables[i]],val.mp[,paste("pred.",variables[i],sep="")])));
if(mi<0){mi<-0};
ma<-max(na.omit(c(cal.mp[,variables[i]],cal.mp[,paste("pred.",variables[i],sep="")],val.mp[,variables[i]],val.mp[,paste("pred.",variables[i],sep="")])));
plot(cal.mp[,variables[i]]~cal.mp[,paste("pred.",variables[i],sep="")],xlab=paste("Predicted ",variables[i],sep=""),ylab=paste("Measured ",variables[i],sep=""),main="Calibration",xlim=c(mi,ma),ylim=c(mi,ma));
plot(val.mp[,variables[i]]~val.mp[,paste("pred.",variables[i],sep="")],xlab=paste("Predicted ",variables[i],sep=""),ylab=paste("Measured ",variables[i],sep=""),main="Validation",xlim=c(mi,ma),ylim=c(mi,ma));
}

cat("Calibration set statistics:\n")
print(cal.stat);
cat("Validation set statistics:\n")
print(val.stat);

# Prepare output:

output<-list("model"=model,"model"=model,"x.tr"=x.tr,"spectral.transformation"=sp,"constituents"=variables,"constituents.transformation"=tr.n,"lambda"=lambda,"method"=r,"cal.samples"=cal,"val.samples"=val,"cal.statistics"=cal.stat,"cal.mea.pre"=cal.mp,"val.statistics"=val.stat,"val.mea.pre"=val.mp,"cal.pca"=pca.l,"mahalanobis"=mah.l,"cal.range"=ran,"rmsep"=rmsep.l,"lm"=lm.l,"wavebands"=waveb,"drv"=drv,"bandwidth"=bandwidth,"filter"=ft,"lv"=lv);
class(output)<-"regr";
if(sw!="NL"){
		setwd(sw);
		}
save(output,file=model)
return(output);

	}

