isric<-function(spec,chem,plo){

# Set percentage of validation samples:

per<-0.3;

# Define the variables used:

variables<-colnames(chem);
variables.bc<-paste("bc.",variables,sep="")[-1]
waveb<-as.numeric(colnames(spec));

# Get information about soil layers:

plo.sep<-list();
for(i in 1: length(unique(plo[,1]))){
	plo.sep[[i]]<-subset(plo,subset=plo[,1]==unique(plo[,1])[i]);
	}
# Number of profiles:
prof.n<-length(unique(plo[,1]));
# Number of layer per profile:
prof.lay<-c();
for(i in 1:prof.n){
	prof.lay[i]<-length(which(plo[,1]==unique(plo[,1])[i]))
	};
prof.lay.sum<-summary(prof.lay);

# Compute the Box-Cox transformation (except pH2):

library(car);
lambda<-c();
chem.bc <-matrix(nrow=nrow(chem),ncol=(ncol(chem)-1),dimnames=list(rownames(chem),variables.bc));
for(i in 1:(ncol(chem)-1)){
lambda[i]<-box.cox.powers(chem[,i+1])[[1]];
chem.bc[,i]<-bc(chem[,i+1],lambda[i]);
	};
chem.bc<-as.data.frame(chem.bc);
detach(package:car);
chem.ph.bc<-cbind(chem[,"pH2"],chem.bc);
colnames(chem.ph.bc)<-c("pH2",variables.bc);

# 1st derivative:

library(KernSmooth);
cols <- ncol(spec);
rows <- nrow(spec);
m <- 1;	
Deri.Smo<-matrix(nrow=rows,ncol=cols);

for (i in 1:rows) {
smoderi<-locpoly(waveb, spec[i,], drv = 1,bandwidth=21, gridsize = cols);
a<-smoderi[[2]];
Deri.Smo[m,]<-a;
m=m+1
}
	
rownames(Deri.Smo)<-rownames(spec); 
colnames(Deri.Smo)<-colnames(spec);
detach(package:KernSmooth);

# Make a PCA of the derivative spectra with reference values:

pca<-prcomp(Deri.Smo,scale=T);
prco<-pca$x[,1:20];
cpv<-summary(pca)[[6]][3,1:20];
zzz<-matrix(nrow=1,ncol=(length(cpv)-4));
for (i in 1:16){
	e<-(cpv[i]+0.04)<cpv[i+3];
	zzz[i]<-e
	}
pc<-(which(zzz==FALSE)-1)[1];	

# Get the cumulative proportion of explained variance (Table 2):

cum<-summary(pca)[[6]][,1:20]

# Get val samples:

val.n<-round(per*nrow(prco),0);

# Get the most extreme points (min and max) for each PC and assign them for the cal set:

cal.min<-c(rep(1,ncol(prco[,1:pc])));
cal.max<-c(rep(1,ncol(prco[,1:pc])));
for (i in 1:ncol(prco[,1:pc])){
	blub<-which(prco[,i]==min(prco[,i]));
	cal.min[i]<-blub[1];
	bla<-which(prco[,i]==max(prco[,i]));
	cal.max[i]<-bla[1];
	}
cal.min<-rownames(prco)[cal.min];
cal.max<-rownames(prco)[cal.max];
cal.start<-unique(c(cal.min,cal.max));
cal.start.n<-match(cal.start,rownames(chem));

# Retrieve all samples which belong to the  profiles of the starting calibration samples:

cal.prof<-unique(plo[cal.start,"Plot_code"]);
cal.prof.out<-list();
for(i in 1:length(cal.prof)){
	cal.prof.out[[i]]<-which(plo[,1]==cal.prof[i]);
	}
cal.prof.out.all<-c();
blub<-cal.prof.out[[1]];
for(i in 2:length(cal.prof)){
	blub<-c(blub,cal.prof.out[[i]]);
	cal.prof.out.all<-blub;
	}

# Exclude all soil profiles where cal samples were assigned to and get from the remaining samples the most extreme points (min and max) for each PC and assign them for the val set:

val.min<-c(rep(1,ncol(prco[-cal.prof.out.all,1:pc])));
val.max<-c(rep(1,ncol(prco[-cal.prof.out.all,1:pc])));
for (i in 1:ncol(prco[-cal.prof.out.all,1:pc])){
	blub<-which(prco[-cal.prof.out.all,i]==min(prco[-cal.prof.out.all,i]));
	val.min[i]<-blub[1];
	bla<-which(prco[-cal.prof.out.all,i]==max(prco[-cal.prof.out.all,i]));
	val.max[i]<-bla[1];
	}
val.min<-rownames(prco[-cal.prof.out.all,])[val.min];
val.max<-rownames(prco[-cal.prof.out.all,])[val.max];
val.start<-unique(c(val.min,val.max));
val.start.n<-match(val.start,rownames(chem));

# Retrieve all samples which belong to the  profiles of the starting validation samples:

val.prof<-unique(plo[val.start,"Plot_code"]);
val.prof.out<-list();
for(i in 1:length(val.prof)){
	val.prof.out[[i]]<-which(plo[,1]==val.prof[i]);
	}
val.prof.out.all<-c();
blub<-val.prof.out[[1]];
for(i in 2:length(val.prof)){
	blub<-c(blub,val.prof.out[[i]]);
	val.prof.out.all<-blub;
	}

# Retrieve all samples, which are already assigned to either cal or val set:

cal.val<-c(cal.prof.out.all,val.prof.out.all);
cal.val.bl<-rownames(chem)[cal.val];
cal.val.start<-match(cal.val.bl,rownames(chem));

# Calculate the inter point distance (euclidean distance):

euc<-as.matrix(dist(prco[,1:pc]));

# Calculate for each object not in the cal.val set the sum of distance to the points in the val.start set:

prco.start<-rownames(prco[-cal.val.start,]);
prco.start.b<-prco.start
val<-rownames(chem)[val.prof.out.all];
for(k in 1:(val.n-length(val.prof.out.all))){
	
	test<-apply(euc[prco.start.b,val],1,min);
	bla<-names(which(test==max(test)));

bla.pr<-plo[bla,"Plot_code"];
bla.pr.bl<-rownames(plo)[which(plo[,1]==bla.pr)];
val<-c(val,bla.pr.bl);

prco.start.b<-prco.start.b[-(which(match(prco.start.b,bla.pr.bl)!="NA"))];

if(length(val)>=val.n)break;
}
val.n<-match(val,rownames(chem));
cal.n<-c(1:nrow(chem))[-val.n];
cal<-rownames(chem)[cal.n];

# Get to know how many plots are only in the cal set or only in the val set:

plo.sep.cal<-matrix(nrow=length(unique(plo[,1])),ncol=3,dimnames=list(as.character(unique(plo[,1])),c("n per plot","n per plot in cal","ratio")));
plo.sep.cal[,1]<-prof.lay;
for(i in 1:length(unique(plo[,1]))){	
	plo.sep.cal[i,2]<-length(which(match(rownames(plo.sep[[i]]),rownames(chem[cal,]))!="NA"));
	}
plo.sep.cal[,3]<-plo.sep.cal[,1]/plo.sep.cal[,2];
plo.cal<-length(which(plo.sep.cal[,3]==1));
plo.val<-length(which(plo.sep.cal[,3]=="Inf"));
plo.n<-length(unique(plo[,1]));

# Get the statistics of calibration and validation samples (Table 1):

chem.stat.all<-summary(chem);
chem.stat.val<-summary(chem[val,]);
chem.stat.cal<-summary(chem[cal,]);
chem.stat.sd.all<-sd(chem);

# Get correlation between transformed soil properties including two new variables for the calibration set (Table 4):

bla<-chem.ph.bc[cal,];
chem.cor.cal<-list();
cor<-c();
for(i in 1:(ncol(bla)-1)){
	for(j in i:ncol(bla)){
		cor[j]<-cor(bla[,i],bla[,j]);
		}
	names(cor)<-colnames(bla)
	chem.cor.cal[[i]]<-cor;
	names(chem.cor.cal)[i]<-colnames(bla)[i];
	cor<-c();
	}

# Get correlation between transformed soil properties including two new variables for all samples:

bla<-c(1:ncol(chem.ph.bc));
chem.cor.all<-list();
cor.all<-c();
for(i in 1:(ncol(chem.ph.bc)-1)){
	for(j in i:ncol(chem.ph.bc)){
		cor.all[j]<-cor(chem.ph.bc[,i],chem.ph.bc[,j]);
		}
	names(cor.all)<-colnames(chem.ph.bc)
	chem.cor.all[[i]]<-cor.all;
	names(chem.cor.all)[i]<-colnames(chem.ph.bc)[i];
	cor.all<-c();
	}

# Make a PCA of the derivative spectra with reference values from the cal set:

cal.pca<-prcomp(Deri.Smo[cal,],scale=T);
cal.prco<-cal.pca$x[,1:20];
cpv<-summary(cal.pca)[[6]][3,1:20];
zzz<-matrix(nrow=1,ncol=(length(cpv)-4));
for (i in 1:16){
	e<-(cpv[i]+0.04)<cpv[i+3];
	zzz[i]<-e
	}
cal.pc<-(which(zzz==FALSE)-1)[1];		

# Plot 10 raw spectra covering the spectral space (Fig. 2):

x11(width=12,height=7);
plot(spec[val[1],]~ waveb,xlim=c(max(waveb),min(waveb)),ylim=c(min(spec[val[1:10],]),max(spec[val[1:10],])),ylab="Absorbance",xlab=expression(paste("Wavenumber", ,sep=" ")  (cm^-1)),type="l",lwd=0.5);
for(i in 2:10){
	lines(spec[val[i],]~waveb,lwd=0.5);
	}

# Set the F-distribution alpha-value:

alpha<-0.975

# Define the variables used:

variables.full.statistics<-c(variables,variables.bc,"pred.pH2","pred.bc.Ca","pred.bc.Mg","pred.bc.Corg","pred.bc.Clay","pred.bc.Sand","pred.bc.CEC.soil","pred.Ca","pred.Mg","pred.Corg","pred.Clay","pred.Sand","pred.CEC.soil")

# Make PLSR for the cal set (R2, bias, SEP, RPD):

# Create statistics ojbect for cal-set (Table 3):

library(pls);
cal.stat.variables<-c("r2","a","bias","RMSEC","RPD","n LV","n bc out");
cal.stat<-matrix(nrow=length(variables),ncol=length(cal.stat.variables),dimnames=list(variables,cal.stat.variables));
cal.mea.pre<-matrix(nrow=length(cal),ncol=(length(variables)*4-2),dimnames=list(cal,variables.full.statistics));
for(i in 1:length(variables)){
	cal.mea.pre[,variables[i]]<-chem[cal,variables[i]];
	}
for(i in 1:length(variables.bc)){
	cal.mea.pre[,variables.bc[i]]<-chem.bc[cal,variables.bc[i]];
	}
cal.mea.pre<-as.data.frame(cal.mea.pre);

# Create statistics ojbect for the val set (Table 3):

val.stat.variables<-c("r2","a","bias","RMSEP","RPD","n bc out");
val.stat<-matrix(nrow=length(variables),ncol=length(val.stat.variables),dimnames=list(variables,val.stat.variables));
val.mea.pre<-matrix(nrow=length(val),ncol=(length(variables)*4-2),dimnames=list(val,variables.full.statistics));
for(i in 1:length(variables)){
	val.mea.pre[,variables[i]]<-chem[val,variables[i]];
	}
for(i in 1:length(variables.bc)){
	val.mea.pre[,variables.bc[i]]<-chem.bc[val,variables.bc[i]];
	}
val.mea.pre<-as.data.frame(val.mea.pre);

# Calculate the F-distribution limits:

f.value<-c();
for(i in 1:nrow(spec)){
	f.value[i]<-qf(alpha,i,i);
	}

# PLS for pH:

cal.ph<-mvr(cal.mea.pre[,"pH2"]~Deri.Smo[cal,],method="kernelpls",validation="CV",ncomp=ifelse(dim(cal.pca$x)[2]<20,dim(cal.pca$x)[2],20));
press.cal.ph<-MSEP(cal.ph)[[1]][1,1,2:ifelse(dim(cal.pca$x)[2]<20,dim(cal.pca$x)[2],21)]*length(cal);
min.cal.ph<-which(press.cal.ph==min(press.cal.ph));
fh.cal.ph<-c();
for(i in 1:min.cal.ph){
	fh.cal.ph[i]<-press.cal.ph[i]/press.cal.ph[min.cal.ph];
	}
cal.ph.lv<-which(fh.cal.ph<f.value[length(cal)])[1];
#x11();
#plot(cal.ph,main="Derivatives cal subset pH",ncomp=cal.ph.lv);
# Calculate cal-subset statistics:
cal.mea.pre[,"pred.pH2"]<-predict.mvr(cal.ph,newdata=Deri.Smo[cal,],ncomp=cal.ph.lv);
#plot(cal.mea.pre[,"pred.pH2"]~cal.mea.pre[,"pH2"]);
# R2:
cal.stat["pH2","r2"]<-
round(cor(cal.mea.pre[,"pH2"],cal.mea.pre[,"pred.pH2"])^2,2);
# a:
cal.stat["pH2","a"]<-round(lm(pred.pH2~pH2,data=cal.mea.pre)$coefficients[2],2);
# Bias:
residual<-cal.mea.pre[,"pred.pH2"]-cal.mea.pre[,"pH2"];
cal.stat["pH2","bias"]<-round(sum(residual)/length(cal.mea.pre[,"pH2"]),2)
# RMSEC:
cal.stat["pH2","RMSEC"]<-
round(sqrt(sum((cal.mea.pre[,"pred.pH2"]-cal.mea.pre[,"pH2"])^2)/(length(cal.mea.pre[,"pH2"])-cal.ph.lv-1)),2);
# RPD:
cal.stat["pH2","RPD"]<-round(sd(cal.mea.pre[,"pH2"])/cal.stat["pH2","RMSEC"],2);
# n LV:
cal.stat["pH2","n LV"]<-cal.ph.lv;

# Calculate val-subset statistics:
val.mea.pre[,"pred.pH2"]<-predict.mvr(cal.ph,newdata=Deri.Smo[val,],ncomp=cal.ph.lv);
#plot(val.mea.pre[,"pred.pH2"]~val.mea.pre[,"pH2"]);
# R2:
val.stat["pH2","r2"]<-round(cor(val.mea.pre[,"pH2"],val.mea.pre[,"pred.pH2"])^2,2);
# a:
val.stat["pH2","a"]<-round(lm(pred.pH2~pH2,data=val.mea.pre)$coefficients[2],2);
# Bias:
residual<-val.mea.pre[,"pred.pH2"]-val.mea.pre[,"pH2"];
val.stat["pH2","bias"]<-round(sum(residual)/length(val.mea.pre[,"pH2"]),2)
# RMSEP:
val.stat["pH2","RMSEP"]<-round(sqrt(sum((val.mea.pre[,"pred.pH2"]-val.mea.pre[,"pH2"])^2)/(length(val.mea.pre[,"pH2"]))),2);
# RPD:
val.stat["pH2","RPD"]<-round(sd(val.mea.pre[,"pH2"])/val.stat["pH2","RMSEP"],2);

# Check on poor predictions in val set for pH>9.9:

ph.10<-plo[rownames(val.mea.pre)[which(val.mea.pre[,"pH2"]>10)],];

# PLS for Ca:

cal.ca<-mvr(cal.mea.pre[,"bc.Ca"] ~Deri.Smo[cal,],method="kernelpls",validation="CV",ncomp=ifelse(dim(cal.pca$x)[2]<20,dim(cal.pca$x)[2],20));
press.cal.ca<-MSEP(cal.ca)[[1]][1,1,2:ifelse(dim(cal.pca$x)[2]<20,dim(cal.pca$x)[2],21)] *length(cal);
min.cal.ca<-which(press.cal.ca==min(press.cal.ca));
fh.cal.ca<-c();
for(i in 1:min.cal.ca){
	fh.cal.ca[i]<-press.cal.ca[i]/press.cal.ca[min.cal.ca];
	}
cal.ca.lv<-which(fh.cal.ca<f.value[length(cal)])[1];
#plot(cal.ca,main="Derivatives cal subset Ca");
# Calculate cal-subset statistics:
cal.mea.pre[,"pred.bc.Ca"]<-predict.mvr(cal.ca,newdata=Deri.Smo[cal,],ncomp=cal.ca.lv);
cal.mea.pre[,"pred.Ca"]<-(cal.mea.pre[,"pred.bc.Ca"]*lambda[1]+1)^(1/lambda[1]);
cal.ca.bc.pred.out<-which(cal.mea.pre[,"pred.Ca"]=="NaN");
#plot(cal.mea.pre[,"pred.Ca"]~cal.mea.pre[,"Ca"]);
if(length(cal.ca.bc.pred.out)!=0)
{# R2:
cal.stat["Ca","r2"]<-round(cor(cal.mea.pre[-cal.ca.bc.pred.out,"Ca"],cal.mea.pre[-cal.ca.bc.pred.out,"pred.Ca"])^2,2);
# a:
cal.stat["Ca","a"]<-round(lm(pred.Ca~Ca,data=cal.mea.pre)$coefficients[2],2);
# Bias:
residual<-cal.mea.pre[-cal.ca.bc.pred.out,"pred.Ca"]-cal.mea.pre[-cal.ca.bc.pred.out,"Ca"];
cal.stat["Ca","bias"]<-round(sum(residual)/length(cal.mea.pre[-(length(cal.ca.bc.pred.out)),"Ca"]),2)
# RMSEC:
cal.stat["Ca","RMSEC"]<-
round(sqrt(sum((cal.mea.pre[-cal.ca.bc.pred.out,"pred.Ca"]-cal.mea.pre[-cal.ca.bc.pred.out,"Ca"])^2)/(length(cal.mea.pre[-cal.ca.bc.pred.out,"Ca"])-cal.ph.lv-1)),2);
# RPD:
cal.stat["Ca","RPD"]<-round(sd(cal.mea.pre[-cal.ca.bc.pred.out,"Ca"])/cal.stat["Ca","RMSEC"],2);
# n LV:
cal.stat["Ca","n LV"]<-cal.ca.lv;
# n bc out:
cal.stat["Ca","n bc out"]<-length(cal.ca.bc.pred.out);
}
if(length(cal.ca.bc.pred.out)==0)
{# R2:
cal.stat["Ca","r2"]<-round(cor(cal.mea.pre[,"Ca"],cal.mea.pre[,"pred.Ca"])^2,2);
# a:
cal.stat["Ca","a"]<-round(lm(pred.Ca~Ca,data=cal.mea.pre)$coefficients[2],2);
# Bias:
residual<-cal.mea.pre[,"pred.Ca"]-cal.mea.pre[,"Ca"];
cal.stat["Ca","bias"]<-round(sum(residual)/length(cal.mea.pre[,"Ca"]),2)
# RMSEC:
cal.stat["Ca","RMSEC"]<-
round(sqrt(sum((cal.mea.pre[,"pred.Ca"]-cal.mea.pre[,"Ca"])^2)/(length(cal.mea.pre[,"Ca"])-cal.ph.lv-1)),2);
# RPD:
cal.stat["Ca","RPD"]<-round(sd(cal.mea.pre[,"Ca"])/cal.stat["Ca","RMSEC"],2);
# n LV:
cal.stat["Ca","n LV"]<-cal.ca.lv;
# n bc out:
cal.stat["Ca","n bc out"]<-length(cal.ca.bc.pred.out);
}

# Calculate val-subset statistics:
val.mea.pre[,"pred.bc.Ca"]<-predict.mvr(cal.ca,newdata=Deri.Smo[val,],ncomp=cal.ca.lv);
val.mea.pre[,"pred.Ca"]<-(val.mea.pre[,"pred.bc.Ca"]*lambda[1]+1)^(1/lambda[1]);
val.ca.bc.pred.out<-which(val.mea.pre[,"pred.Ca"]=="NaN");
#plot(val.mea.pre[,"pred.Ca"]~val.mea.pre[,"Ca"]);
if(length(val.ca.bc.pred.out)!=0)
{# R2:
val.stat["Ca","r2"]<-round(cor(val.mea.pre[-val.ca.bc.pred.out,"Ca"],val.mea.pre[-val.ca.bc.pred.out,"pred.Ca"])^2,2);
# a:
val.stat["Ca","a"]<-round(lm(pred.Ca~Ca,data=val.mea.pre)$coefficients[2],2);
# Bias:
residual<-val.mea.pre[-val.ca.bc.pred.out,"pred.Ca"]-val.mea.pre[-val.ca.bc.pred.out,"Ca"];
val.stat["Ca","bias"]<-round(sum(residual)/length(val.mea.pre[-(length(val.ca.bc.pred.out)),"Ca"]),2)
# RMSEP:
val.stat["Ca","RMSEP"]<-
round(sqrt(sum((val.mea.pre[-val.ca.bc.pred.out,"pred.Ca"]-val.mea.pre[-val.ca.bc.pred.out,"Ca"])^2)/(length(val.mea.pre[-val.ca.bc.pred.out,"Ca"]))),2)
# RPD:
val.stat["Ca","RPD"]<-round(sd(val.mea.pre[-val.ca.bc.pred.out,"Ca"])/val.stat["Ca","RMSEP"],2);
# n bc out:
val.stat["Ca","n bc out"]<-length(val.ca.bc.pred.out);
}
if(length(val.ca.bc.pred.out)==0)
{# R2:
val.stat["Ca","r2"]<-round(cor(val.mea.pre[,"Ca"],val.mea.pre[,"pred.Ca"])^2,2);
# a:
val.stat["Ca","a"]<-round(lm(pred.Ca~Ca,data=val.mea.pre)$coefficients[2],2);
# Bias:
residual<-val.mea.pre[,"pred.Ca"]-val.mea.pre[,"Ca"];
val.stat["Ca","bias"]<-round(sum(residual)/length(val.mea.pre[,"Ca"]),2)
# RMSEP:
val.stat["Ca","RMSEP"]<-
round(sqrt(sum((val.mea.pre[,"pred.Ca"]-val.mea.pre[,"Ca"])^2)/(length(val.mea.pre[,"Ca"]))),2)
# RPD:
val.stat["Ca","RPD"]<-round(sd(val.mea.pre[,"Ca"])/val.stat["Ca","RMSEP"],2);
# n bc out:
val.stat["Ca","n bc out"]<-length(val.ca.bc.pred.out);
}

# PLS for Mg:

cal.mg<-mvr(cal.mea.pre[,"bc.Mg"] ~Deri.Smo[cal,],method="kernelpls",validation="CV",ncomp=ifelse(dim(cal.pca$x)[2]<20,dim(cal.pca$x)[2],20));
press.cal.mg<-MSEP(cal.mg)[[1]][1,1,2:ifelse(dim(cal.pca$x)[2]<20,dim(cal.pca$x)[2],21)] *length(cal);
min.cal.mg<-which(press.cal.mg==min(press.cal.mg));
fh.cal.mg<-c();
for(i in 1:min.cal.mg){
	fh.cal.mg[i]<-press.cal.mg[i]/press.cal.mg[min.cal.mg];
	}
cal.mg.lv<-which(fh.cal.mg<f.value[length(cal)])[1];
#plot(cal.mg,main="Derivatives cal subset Mg");
# Calculate cal-subset statistics:
cal.mea.pre[,"pred.bc.Mg"]<-predict.mvr(cal.mg,newdata=Deri.Smo[cal,],ncomp=cal.mg.lv);
cal.mea.pre[,"pred.Mg"]<-(cal.mea.pre[,"pred.bc.Mg"]*lambda[2]+1)^(1/lambda[2]);
cal.mg.bc.pred.out<-as.numeric(which(cal.mea.pre[,"pred.Mg"]=="NaN"));
#plot(cal.mea.pre[,"pred.Mg"]~cal.mea.pre[,"Mg"]);
if(length(cal.mg.bc.pred.out)!=0)
{# R2:
cal.stat["Mg","r2"]<-round(cor(cal.mea.pre[-cal.mg.bc.pred.out,"Mg"],cal.mea.pre[-cal.mg.bc.pred.out,"pred.Mg"])^2,2);
# a:
cal.stat["Mg","a"]<-round(lm(pred.Mg~Mg,data=cal.mea.pre)$coefficients[2],2);
# Bias:
residual<-cal.mea.pre[-cal.mg.bc.pred.out,"pred.Mg"]-cal.mea.pre[-cal.mg.bc.pred.out,"Mg"];
cal.stat["Mg","bias"]<-round(sum(residual)/length(cal.mea.pre[-(length(cal.mg.bc.pred.out)),"Mg"]),2)
# RMSEC:
cal.stat["Mg","RMSEC"]<-
round(sqrt(sum((cal.mea.pre[-cal.mg.bc.pred.out,"pred.Mg"]-cal.mea.pre[-cal.mg.bc.pred.out,"Mg"])^2)/(length(cal.mea.pre[-cal.mg.bc.pred.out,"Mg"])-cal.ph.lv-1)),2);
# RPD:
cal.stat["Mg","RPD"]<-round(sd(cal.mea.pre[-cal.mg.bc.pred.out,"Mg"])/cal.stat["Mg","RMSEC"],2);
# n LV:
cal.stat["Mg","n LV"]<-cal.mg.lv;
# n bc out:
cal.stat["Mg","n bc out"]<-length(cal.mg.bc.pred.out);
}
if(length(cal.mg.bc.pred.out)==0)
{# R2:
cal.stat["Mg","r2"]<-round(cor(cal.mea.pre[,"Mg"],cal.mea.pre[,"pred.Mg"])^2,2);
# a:
cal.stat["Mg","a"]<-round(lm(pred.Mg~Mg,data=cal.mea.pre)$coefficients[2],2);
# Bias:
residual<-cal.mea.pre[,"pred.Mg"]-cal.mea.pre[,"Mg"];
cal.stat["Mg","bias"]<-round(sum(residual)/length(cal.mea.pre[,"Mg"]),2)
# RMSEC:
cal.stat["Mg","RMSEC"]<-
round(sqrt(sum((cal.mea.pre[,"pred.Mg"]-cal.mea.pre[,"Mg"])^2)/(length(cal.mea.pre[,"Mg"])-cal.ph.lv-1)),2);
# RPD:
cal.stat["Mg","RPD"]<-round(sd(cal.mea.pre[,"Mg"])/cal.stat["Mg","RMSEC"],2);
# n LV:
cal.stat["Mg","n LV"]<-cal.mg.lv;
# n bc out:
cal.stat["Mg","n bc out"]<-length(cal.mg.bc.pred.out);
}

# Calculate val-subset statistics:
val.mea.pre[,"pred.bc.Mg"]<-predict.mvr(cal.mg,newdata=Deri.Smo[val,],ncomp=cal.mg.lv);
val.mea.pre[,"pred.Mg"]<-(val.mea.pre[,"pred.bc.Mg"]*lambda[2]+1)^(1/lambda[2]);
val.mg.bc.pred.out<-as.numeric(which(val.mea.pre[,"pred.Mg"]=="NaN"));
#plot(val.mea.pre[,"pred.Mg"]~val.mea.pre[,"Mg"]);
if(length(val.mg.bc.pred.out)!=0)
{
# R2:
val.stat["Mg","r2"]<-round(cor(val.mea.pre[-val.mg.bc.pred.out,"Mg"],val.mea.pre[-val.mg.bc.pred.out,"pred.Mg"])^2,2);
# a:
val.stat["Mg","a"]<-round(lm(pred.Mg~Mg,data=val.mea.pre)$coefficients[2],2);
# Bias:
residual<-val.mea.pre[-val.mg.bc.pred.out,"pred.Mg"]-val.mea.pre[-val.mg.bc.pred.out,"Mg"];
val.stat["Mg","bias"]<-round(sum(residual)/length(val.mea.pre[-(length(val.mg.bc.pred.out)),"Mg"]),2)
# RMSEP:
val.stat["Mg","RMSEP"]<-
round(sqrt(sum((val.mea.pre[-val.mg.bc.pred.out,"pred.Mg"]-val.mea.pre[-val.mg.bc.pred.out,"Mg"])^2)/(length(val.mea.pre[-val.mg.bc.pred.out,"Mg"]))),2);
# RPD:
val.stat["Mg","RPD"]<-round(sd(val.mea.pre[-val.mg.bc.pred.out,"Mg"])/val.stat["Mg","RMSEP"],2);
# n bc out:
val.stat["Mg","n bc out"]<-length(val.mg.bc.pred.out);
}
if(length(val.mg.bc.pred.out)==0)
{
# R2:
val.stat["Mg","r2"]<-round(cor(val.mea.pre[,"Mg"],val.mea.pre[,"pred.Mg"])^2,2);
# a:
val.stat["Mg","a"]<-round(lm(pred.Mg~Mg,data=val.mea.pre)$coefficients[2],2);
# Bias:
residual<-val.mea.pre[,"pred.Mg"]-val.mea.pre[,"Mg"];
val.stat["Mg","bias"]<-round(sum(residual)/length(val.mea.pre[,"Mg"]),2)
# RMSEP:
val.stat["Mg","RMSEP"]<-
round(sqrt(sum((val.mea.pre[,"pred.Mg"]-val.mea.pre[,"Mg"])^2)/(length(val.mea.pre[,"Mg"]))),2);# RPD:
val.stat["Mg","RPD"]<-round(sd(val.mea.pre[,"Mg"])/val.stat["Mg","RMSEP"],2);
# n bc out:
val.stat["Mg","n bc out"]<-length(val.mg.bc.pred.out);
}

# PLS for Corg:

cal.orgc<-mvr(cal.mea.pre[,"bc.Corg"] ~Deri.Smo[cal,],method="kernelpls",validation="CV",ncomp=ifelse(dim(cal.pca$x)[2]<20,dim(cal.pca$x)[2],20));
press.cal.orgc<-MSEP(cal.orgc)[[1]][1,1,2:ifelse(dim(cal.pca$x)[2]<20,dim(cal.pca$x)[2],21)] *length(cal);
min.cal.orgc<-which(press.cal.orgc==min(press.cal.orgc));
fh.cal.orgc<-c();
for(i in 1:min.cal.orgc){
	fh.cal.orgc[i]<-press.cal.orgc[i]/press.cal.orgc[min.cal.orgc];
	}
cal.orgc.lv<-which(fh.cal.orgc<f.value[length(cal)])[1];
#plot(cal.orgc,main="Derivatives cal subset Corg");
# Calculate cal-subset statistics:
cal.mea.pre[,"pred.bc.Corg"]<-predict.mvr(cal.orgc,newdata=Deri.Smo[cal,],ncomp=cal.orgc.lv);
cal.mea.pre[,"pred.Corg"]<-(cal.mea.pre[,"pred.bc.Corg"]*lambda[3]+1)^(1/lambda[3]);
cal.orgc.bc.pred.out<-as.numeric(which(cal.mea.pre[,"pred.Corg"]=="NaN"));
#plot(cal.mea.pre[,"pred.Corg"]~cal.mea.pre[,"Corg"]);
if(length(cal.orgc.bc.pred.out)!=0)
{
# R2:
cal.stat["Corg","r2"]<-round(cor(cal.mea.pre[-cal.orgc.bc.pred.out
,"Corg"],cal.mea.pre[-cal.orgc.bc.pred.out
,"pred.Corg"])^2,2);
# a:
cal.stat["Corg","a"]<-round(lm(pred.Corg~Corg,data=cal.mea.pre)$coefficients[2],2);
# Bias:
residual<-cal.mea.pre[-cal.orgc.bc.pred.out
,"pred.Corg"]-cal.mea.pre[-cal.orgc.bc.pred.out
,"Corg"];
cal.stat["Corg","bias"]<-round(sum(residual)/length(cal.mea.pre[-(length(cal.orgc.bc.pred.out)),"Corg"]),2)
# RMSEC:
cal.stat["Corg","RMSEC"]<-
round(sqrt(sum((cal.mea.pre[-cal.orgc.bc.pred.out,"pred.Corg"]-cal.mea.pre[-cal.orgc.bc.pred.out,"Corg"])^2)/(length(cal.mea.pre[-cal.orgc.bc.pred.out,"Corg"])-cal.ph.lv-1)),2);
# RPD:
cal.stat["Corg","RPD"]<-round(sd(cal.mea.pre[-cal.orgc.bc.pred.out
,"Corg"])/cal.stat["Corg","RMSEC"],2);
# n LV:
cal.stat["Corg","n LV"]<-cal.orgc.lv;
# n bc out:
cal.stat["Corg","n bc out"]<-length(cal.orgc.bc.pred.out);
}
if(length(cal.orgc.bc.pred.out)==0)
{
# R2:
cal.stat["Corg","r2"]<-round(cor(cal.mea.pre[,"Corg"],cal.mea.pre[,"pred.Corg"])^2,2);
# a:
cal.stat["Corg","a"]<-round(lm(pred.Corg~Corg,data=cal.mea.pre)$coefficients[2],2);
# Bias:
residual<-cal.mea.pre[,"pred.Corg"]-cal.mea.pre[,"Corg"];
cal.stat["Corg","bias"]<-round(sum(residual)/length(cal.mea.pre[,"Corg"]),2)
# RMSEC:
cal.stat["Corg","RMSEC"]<-
round(sqrt(sum((cal.mea.pre[,"pred.Corg"]-cal.mea.pre[,"Corg"])^2)/(length(cal.mea.pre[,"Corg"])-cal.ph.lv-1)),2);
# RPD:
cal.stat["Corg","RPD"]<-round(sd(cal.mea.pre[,"Corg"])/cal.stat["Corg","RMSEC"],2);
# n LV:
cal.stat["Corg","n LV"]<-cal.orgc.lv;
# n bc out:
cal.stat["Corg","n bc out"]<-length(cal.orgc.bc.pred.out);
}

# Calculate val-subset statistics:
val.mea.pre[,"pred.bc.Corg"]<-predict.mvr(cal.orgc,newdata=Deri.Smo[val,],ncomp=cal.orgc.lv);
val.mea.pre[,"pred.Corg"]<-(val.mea.pre[,"pred.bc.Corg"]*lambda[3]+1)^(1/lambda[3]);
val.orgc.bc.pred.out<-as.numeric(which(val.mea.pre[,"pred.Corg"]=="NaN"));
#plot(val.mea.pre[,"pred.Corg"]~val.mea.pre[,"Corg"]);
if(length(val.orgc.bc.pred.out)!=0)
{
# R2:
val.stat["Corg","r2"]<-round(cor(val.mea.pre[-val.orgc.bc.pred.out,"Corg"],val.mea.pre[-val.orgc.bc.pred.out,"pred.Corg"])^2,2);
val.stat["Corg","a"]<-round(lm(pred.Corg~Corg,data=val.mea.pre)$coefficients[2],2);
# Bias:
residual<-val.mea.pre[-val.orgc.bc.pred.out,"pred.Corg"]-val.mea.pre[-val.orgc.bc.pred.out,"Corg"];
val.stat["Corg","bias"]<-round(sum(residual)/length(val.mea.pre[-(length(val.orgc.bc.pred.out)),"Corg"]),2)
# RMSEP:
val.stat["Corg","RMSEP"]<-
round(sqrt(sum((val.mea.pre[-val.orgc.bc.pred.out,"pred.Corg"]-val.mea.pre[-val.orgc.bc.pred.out,"Corg"])^2)/(length(val.mea.pre[-val.orgc.bc.pred.out,"Corg"]))),2);
# RPD:
val.stat["Corg","RPD"]<-round(sd(val.mea.pre[-val.orgc.bc.pred.out,"Corg"])/val.stat["Corg","RMSEP"],2);
# n bc out:
val.stat["Corg","n bc out"]<-length(val.orgc.bc.pred.out);
}
if(length(val.orgc.bc.pred.out)==0)
{
# R2:
val.stat["Corg","r2"]<-round(cor(val.mea.pre[,"Corg"],val.mea.pre[,"pred.Corg"])^2,2);
val.stat["Corg","a"]<-round(lm(pred.Corg~Corg,data=val.mea.pre)$coefficients[2],2);
# Bias:
residual<-val.mea.pre[,"pred.Corg"]-val.mea.pre[,"Corg"];
val.stat["Corg","bias"]<-round(sum(residual)/length(val.mea.pre[,"Corg"]),2)
# RMSEP:
val.stat["Corg","RMSEP"]<-
round(sqrt(sum((val.mea.pre[,"pred.Corg"]-val.mea.pre[,"Corg"])^2)/(length(val.mea.pre[,"Corg"]))),2);
# RPD:
val.stat["Corg","RPD"]<-round(sd(val.mea.pre[,"Corg"])/val.stat["Corg","RMSEP"],2);
# n bc out:
val.stat["Corg","n bc out"]<-length(val.orgc.bc.pred.out);
}

# PLS for Clay:

cal.clay<-mvr(cal.mea.pre[,"bc.Clay"] ~Deri.Smo[cal,],method="kernelpls",validation="CV",ncomp=ifelse(dim(cal.pca$x)[2]<20,dim(cal.pca$x)[2],20));
press.cal.clay<-MSEP(cal.clay)[[1]][1,1,2:ifelse(dim(cal.pca$x)[2]<20,dim(cal.pca$x)[2],21)] *length(cal);
min.cal.clay<-which(press.cal.clay==min(press.cal.clay));
fh.cal.clay<-c();
for(i in 1:min.cal.clay){
	fh.cal.clay[i]<-press.cal.clay[i]/press.cal.clay[min.cal.clay];
	}
cal.clay.lv<-which(fh.cal.clay<f.value[length(cal)])[1];
#plot(cal.clay,main="Derivatives cal subset Clay");
# Calculate cal-subset statistics:
cal.mea.pre[,"pred.bc.Clay"]<-predict.mvr(cal.clay,newdata=Deri.Smo[cal,],ncomp=cal.clay.lv);
cal.mea.pre[,"pred.Clay"]<-(cal.mea.pre[,"pred.bc.Clay"]*lambda[4]+1)^(1/lambda[4]);
cal.clay.bc.pred.out<-as.numeric(which(cal.mea.pre[,"pred.Clay"]=="NaN"));
clay<-which(cal.mea.pre[,"pred.Clay"]>=100);
if(length(clay)>0)
{cal.mea.pre[clay,"pred.Clay"]<-100};
#plot(cal.mea.pre[,"pred.Clay"]~cal.mea.pre[,"Clay"]);
if(length(cal.clay.bc.pred.out)!=0)
{
# R2:
cal.stat["Clay","r2"]<-round(cor(cal.mea.pre[-cal.clay.bc.pred.out,"Clay"],cal.mea.pre[-cal.clay.bc.pred.out,"pred.Clay"])^2,2);
# a:
cal.stat["Clay","a"]<-round(lm(pred.Clay~Clay,data=cal.mea.pre)$coefficients[2],2);
# Bias:
residual<-cal.mea.pre[-cal.clay.bc.pred.out,"pred.Clay"]-cal.mea.pre[-cal.clay.bc.pred.out,"Clay"];
cal.stat["Clay","bias"]<-round(sum(residual)/length(cal.mea.pre[-(length(cal.clay.bc.pred.out)),"Clay"]),2)
# RMSEC:
cal.stat["Clay","RMSEC"]<-
round(sqrt(sum((cal.mea.pre[-cal.clay.bc.pred.out,"pred.Clay"]-cal.mea.pre[-cal.clay.bc.pred.out,"Clay"])^2)/(length(cal.mea.pre[-cal.clay.bc.pred.out,"Clay"])-cal.ph.lv-1)),2);
# RPD:
cal.stat["Clay","RPD"]<-round(sd(cal.mea.pre[-cal.clay.bc.pred.out,"Clay"])/cal.stat["Clay","RMSEC"],2);
# n LV:
cal.stat["Clay","n LV"]<-cal.clay.lv;
# n bc out:
cal.stat["Clay","n bc out"]<-length(cal.clay.bc.pred.out);
}
if(length(cal.clay.bc.pred.out)==0)
{
# R2:
cal.stat["Clay","r2"]<-round(cor(cal.mea.pre[,"Clay"],cal.mea.pre[,"pred.Clay"])^2,2);
# a:
cal.stat["Clay","a"]<-round(lm(pred.Clay~Clay,data=cal.mea.pre)$coefficients[2],2);
# Bias:
residual<-cal.mea.pre[,"pred.Clay"]-cal.mea.pre[,"Clay"];
cal.stat["Clay","bias"]<-round(sum(residual)/length(cal.mea.pre[,"Clay"]),2)
# RMSEC:
cal.stat["Clay","RMSEC"]<-
round(sqrt(sum((cal.mea.pre[,"pred.Clay"]-cal.mea.pre[,"Clay"])^2)/(length(cal.mea.pre[,"Clay"])-cal.ph.lv-1)),2);
# RPD:
cal.stat["Clay","RPD"]<-round(sd(cal.mea.pre[,"Clay"])/cal.stat["Clay","RMSEC"],2);
# n LV:
cal.stat["Clay","n LV"]<-cal.clay.lv;
# n bc out:
cal.stat["Clay","n bc out"]<-length(cal.clay.bc.pred.out);
}

# Calculate val-subset statistics:
val.mea.pre[,"pred.bc.Clay"]<-predict.mvr(cal.clay,newdata=Deri.Smo[val,],ncomp=cal.clay.lv);
val.mea.pre[,"pred.Clay"]<-(val.mea.pre[,"pred.bc.Clay"]*lambda[4]+1)^(1/lambda[4]);
val.clay.bc.pred.out<-as.numeric(which(val.mea.pre[,"pred.Clay"]=="NaN"));
clay<-which(val.mea.pre[,"pred.Clay"]>=100);
if(length(clay)>0)
{val.mea.pre[clay,"pred.Clay"]<-100};
#plot(val.mea.pre[,"pred.Clay"]~val.mea.pre[,"Clay"]);
if(length(val.clay.bc.pred.out)!=0)
{
# R2:
val.stat["Clay","r2"]<-round(cor(val.mea.pre[-val.clay.bc.pred.out,"Clay"],val.mea.pre[-val.clay.bc.pred.out,"pred.Clay"])^2,2);
# a:
val.stat["Clay","a"]<-round(lm(pred.Clay~Clay,data=val.mea.pre)$coefficients[2],2);
# Bias:
residual<-val.mea.pre[-val.clay.bc.pred.out,"pred.Clay"]-val.mea.pre[-val.clay.bc.pred.out,"Clay"];
val.stat["Clay","bias"]<-round(sum(residual)/length(val.mea.pre[-(length(val.clay.bc.pred.out)),"Clay"]),2)
# RMSEP:
val.stat["Clay","RMSEP"]<-
round(sqrt(sum((val.mea.pre[-val.clay.bc.pred.out,"pred.Clay"]-val.mea.pre[-val.clay.bc.pred.out,"Clay"])^2)/(length(val.mea.pre[-val.clay.bc.pred.out,"Clay"]))),2);
# RPD:
val.stat["Clay","RPD"]<-round(sd(val.mea.pre[-val.clay.bc.pred.out,"Clay"])/val.stat["Clay","RMSEP"],2);
# n bc out:
val.stat["Clay","n bc out"]<-length(val.clay.bc.pred.out);
}
if(length(val.clay.bc.pred.out)==0)
{
# R2:
val.stat["Clay","r2"]<-round(cor(val.mea.pre[,"Clay"],val.mea.pre[,"pred.Clay"])^2,2);
# a:
val.stat["Clay","a"]<-round(lm(pred.Clay~Clay,data=val.mea.pre)$coefficients[2],2);
# Bias:
residual<-val.mea.pre[,"pred.Clay"]-val.mea.pre[,"Clay"];
val.stat["Clay","bias"]<-round(sum(residual)/length(val.mea.pre[,"Clay"]),2)
# RMSEP:
val.stat["Clay","RMSEP"]<-
round(sqrt(sum((val.mea.pre[,"pred.Clay"]-val.mea.pre[,"Clay"])^2)/(length(val.mea.pre[,"Clay"]))),2);
# RPD:
val.stat["Clay","RPD"]<-round(sd(val.mea.pre[,"Clay"])/val.stat["Clay","RMSEP"],2);
# n bc out:
val.stat["Clay","n bc out"]<-length(val.clay.bc.pred.out);
}

# PLS for Sand:

cal.tots<-mvr(cal.mea.pre[,"bc.Sand"] ~Deri.Smo[cal,],method="kernelpls",validation="CV",ncomp=ifelse(dim(cal.pca$x)[2]<20,dim(cal.pca$x)[2],20));
press.cal.tots<-MSEP(cal.tots)[[1]][1,1,2:ifelse(dim(cal.pca$x)[2]<20,dim(cal.pca$x)[2],21)] *length(cal);
min.cal.tots<-which(press.cal.tots==min(press.cal.tots));
fh.cal.tots<-c();
for(i in 1:min.cal.tots){
	fh.cal.tots[i]<-press.cal.tots[i]/press.cal.tots[min.cal.tots];
	}
cal.tots.lv<-which(fh.cal.tots<f.value[length(cal)])[1];
#plot(cal.tots,main="Derivatives cal subset Sand");
# Calculate cal-subset statistics:
cal.mea.pre[,"pred.bc.Sand"]<-predict.mvr(cal.tots,newdata=Deri.Smo[cal,],ncomp=cal.tots.lv);
cal.mea.pre[,"pred.Sand"]<-(cal.mea.pre[,"pred.bc.Sand"]*lambda[5]+1)^(1/lambda[5]);
cal.tots.bc.pred.out<-as.numeric(which(cal.mea.pre[,"pred.Sand"]=="NaN"));
sand<-which(cal.mea.pre[,"pred.Sand"]>=100);
sand.cal<-which(cal.mea.pre[,"pred.Sand"]>=100);
if(length(sand)>0)
{cal.mea.pre[sand,"pred.Sand"]<-100};
#plot(cal.mea.pre[,"pred.Sand"]~cal.mea.pre[,"Sand"]);
if(length(cal.tots.bc.pred.out)!=0)
{
# R2:
cal.stat["Sand","r2"]<-round(cor(cal.mea.pre[-cal.tots.bc.pred.out,"Sand"],cal.mea.pre[-cal.tots.bc.pred.out,"pred.Sand"])^2,2);
# a:
cal.stat["Sand","a"]<-round(lm(pred.Sand~Sand,data=cal.mea.pre)$coefficients[2],2);
# Bias:
residual<-cal.mea.pre[-cal.tots.bc.pred.out,"pred.Sand"]-cal.mea.pre[-cal.tots.bc.pred.out,"Sand"];
cal.stat["Sand","bias"]<-round(sum(residual)/length(cal.mea.pre[-(length(cal.tots.bc.pred.out)),"Sand"]),2)
# RMSEC:
cal.stat["Sand","RMSEC"]<-
round(sqrt(sum((cal.mea.pre[-cal.sand.bc.pred.out,"pred.Sand"]-cal.mea.pre[-cal.sand.bc.pred.out,"Sand"])^2)/(length(cal.mea.pre[-cal.sand.bc.pred.out,"Sand"])-cal.ph.lv-1)),2);
# RPD:
cal.stat["Sand","RPD"]<-round(sd(cal.mea.pre[-cal.tots.bc.pred.out,"Sand"])/cal.stat["Sand","RMSEC"],2);
# n LV:
cal.stat["Sand","n LV"]<-cal.tots.lv;
# n bc out:
cal.stat["Sand","n bc out"]<-length(cal.tots.bc.pred.out);
}
if(length(cal.tots.bc.pred.out)==0)
{
# R2:
cal.stat["Sand","r2"]<-round(cor(cal.mea.pre[,"Sand"],cal.mea.pre[,"pred.Sand"])^2,2);
# a:
cal.stat["Sand","a"]<-round(lm(pred.Sand~Sand,data=cal.mea.pre)$coefficients[2],2);
# Bias:
residual<-cal.mea.pre[,"pred.Sand"]-cal.mea.pre[,"Sand"];
cal.stat["Sand","bias"]<-round(sum(residual)/length(cal.mea.pre[,"Sand"]),2)
# RMSEC:
cal.stat["Sand","RMSEC"]<-
round(sqrt(sum((cal.mea.pre[,"pred.Sand"]-cal.mea.pre[,"Sand"])^2)/(length(cal.mea.pre[,"Sand"])-cal.ph.lv-1)),2);
# RPD:
cal.stat["Sand","RPD"]<-round(sd(cal.mea.pre[,"Sand"])/cal.stat["Sand","RMSEC"],2);
# n LV:
cal.stat["Sand","n LV"]<-cal.tots.lv;
# n bc out:
cal.stat["Sand","n bc out"]<-length(cal.tots.bc.pred.out);
}

# Calculate val-subset statistics:
val.mea.pre[,"pred.bc.Sand"]<-predict.mvr(cal.tots,newdata=Deri.Smo[val,],ncomp=cal.tots.lv);
val.mea.pre[,"pred.Sand"]<-(val.mea.pre[,"pred.bc.Sand"]*lambda[5]+1)^(1/lambda[5]);
val.tots.bc.pred.out<-as.numeric(which(val.mea.pre[,"pred.Sand"]=="NaN"));
sand<-which(val.mea.pre[,"pred.Sand"]>=100);
sand.val<-which(val.mea.pre[,"pred.Sand"]>=100);
if(length(sand)>0)
{val.mea.pre[sand,"pred.Sand"]<-100};
#plot(val.mea.pre[,"pred.Sand"]~val.mea.pre[,"Sand"]);
if(length(val.tots.bc.pred.out)!=0)
{
# R2:
val.stat["Sand","r2"]<-round(cor(val.mea.pre[-val.tots.bc.pred.out,"Sand"],val.mea.pre[-val.tots.bc.pred.out,"pred.Sand"])^2,2);
# a:
val.stat["Sand","a"]<-round(lm(pred.Sand~Sand,data=val.mea.pre)$coefficients[2],2);
# Bias:
residual<-val.mea.pre[-val.tots.bc.pred.out,"pred.Sand"]-val.mea.pre[-val.tots.bc.pred.out,"Sand"];
val.stat["Sand","bias"]<-round(sum(residual)/length(val.mea.pre[-(length(val.tots.bc.pred.out)),"Sand"]),2)
# RMSEP:
val.stat["Sand","RMSEP"]<-
round(sqrt(sum((val.mea.pre[-val.tots.bc.pred.out,"pred.Sand"]-val.mea.pre[-val.tots.bc.pred.out,"Sand"])^2)/(length(val.mea.pre[-val.tots.bc.pred.out,"Sand"]))),2);
# RPD:
val.stat["Sand","RPD"]<-round(sd(val.mea.pre[-val.tots.bc.pred.out,"Sand"])/val.stat["Sand","RMSEP"],2);
# n bc out:
val.stat["Sand","n bc out"]<-length(val.tots.bc.pred.out);
}
if(length(val.tots.bc.pred.out)==0)
{
# R2:
val.stat["Sand","r2"]<-round(cor(val.mea.pre[,"Sand"],val.mea.pre[,"pred.Sand"])^2,2);
# a:
val.stat["Sand","a"]<-round(lm(pred.Sand~Sand,data=val.mea.pre)$coefficients[2],2);
# Bias:
residual<-val.mea.pre[,"pred.Sand"]-val.mea.pre[,"Sand"];
val.stat["Sand","bias"]<-round(sum(residual)/length(val.mea.pre[,"Sand"]),2)
# RMSEP:
val.stat["Sand","RMSEP"]<-
round(sqrt(sum((val.mea.pre[,"pred.Sand"]-val.mea.pre[,"Sand"])^2)/(length(val.mea.pre[,"Sand"]))),2);
# RPD:
val.stat["Sand","RPD"]<-round(sd(val.mea.pre[,"Sand"])/val.stat["Sand","RMSEP"],2);
# n bc out:
val.stat["Sand","n bc out"]<-length(val.tots.bc.pred.out);
}

# PLS for CEC:

cal.cecsoil<-mvr(cal.mea.pre[,"bc.CEC.soil"] ~Deri.Smo[cal,],method="kernelpls",validation="CV",ncomp=ifelse(dim(cal.pca$x)[2]<20,dim(cal.pca$x)[2],20));
press.cal.cecsoil<-MSEP(cal.cecsoil)[[1]][1,1,2:ifelse(dim(cal.pca$x)[2]<20,dim(cal.pca$x)[2],21)] *length(cal);
min.cal.cecsoil<-which(press.cal.cecsoil==min(press.cal.cecsoil));
fh.cal.cecsoil<-c();
for(i in 1:min.cal.cecsoil){
	fh.cal.cecsoil[i]<-press.cal.cecsoil[i]/press.cal.cecsoil[min.cal.cecsoil];
	}
cal.cecsoil.lv<-which(fh.cal.cecsoil<f.value[length(cal)])[1];
#plot(cal.cecsoil,main="Derivatives cal subset CEC");
# Calculate cal-subset statistics:
cal.mea.pre[,"pred.bc.CEC.soil"]<-predict.mvr(cal.cecsoil,newdata=Deri.Smo[cal,],ncomp=cal.cecsoil.lv);
cal.mea.pre[,"pred.CEC.soil"]<-(cal.mea.pre[,"pred.bc.CEC.soil"]*lambda[6]+1)^(1/lambda[6]);
cal.cecsoil.bc.pred.out<-as.numeric(which(cal.mea.pre[,"pred.CEC.soil"]=="NaN"));
#plot(cal.mea.pre[,"pred.CEC.soil"]~cal.mea.pre[,"CEC.soil"]);
if(length(cal.cecsoil.bc.pred.out)!=0)
{
# R2:
cal.stat["CEC.soil","r2"]<-round(cor(cal.mea.pre[-cal.cecsoil.bc.pred.out,"CEC.soil"],cal.mea.pre[-cal.cecsoil.bc.pred.out,"pred.CEC.soil"])^2,2);
# a:
cal.stat["CEC.soil","a"]<-round(lm(pred.CEC.soil~CEC.soil,data=cal.mea.pre)$coefficients[2],2);
# Bias:
residual<-cal.mea.pre[-cal.cecsoil.bc.pred.out,"pred.CEC.soil"]-cal.mea.pre[-cal.cecsoil.bc.pred.out,"CEC.soil"];
cal.stat["CEC.soil","bias"]<-round(sum(residual)/length(cal.mea.pre[-(length(cal.cecsoil.bc.pred.out)),"CEC.soil"]),2)
# RMSEC:
cal.stat["CEC.soil","RMSEC"]<-
round(sqrt(sum((cal.mea.pre[-cal.cecsoil.bc.pred.out,"pred.CEC.soil"]-cal.mea.pre[-cal.cecsoil.bc.pred.out,"CEC.soil"])^2)/(length(cal.mea.pre[-cal.cecsoil.bc.pred.out,"CEC.soil"])-cal.ph.lv-1)),2);
# RPD:
cal.stat["CEC.soil","RPD"]<-round(sd(cal.mea.pre[-cal.cecsoil.bc.pred.out,"CEC.soil"])/cal.stat["CEC.soil","RMSEC"],2);
# n LV:
cal.stat["CEC.soil","n LV"]<-cal.cecsoil.lv;
# n bc out:
cal.stat["CEC.soil","n bc out"]<-length(cal.cecsoil.bc.pred.out);
}
if(length(cal.cecsoil.bc.pred.out)==0)
{
# R2:
cal.stat["CEC.soil","r2"]<-round(cor(cal.mea.pre[,"CEC.soil"],cal.mea.pre[,"pred.CEC.soil"])^2,2);
# a:
cal.stat["CEC.soil","a"]<-round(lm(pred.CEC.soil~CEC.soil,data=cal.mea.pre)$coefficients[2],2);
# Bias:
residual<-cal.mea.pre[,"pred.CEC.soil"]-cal.mea.pre[,"CEC.soil"];
cal.stat["CEC.soil","bias"]<-round(sum(residual)/length(cal.mea.pre[,"CEC.soil"]),2)
# RMSEC:
cal.stat["CEC.soil","RMSEC"]<-
round(sqrt(sum((cal.mea.pre[,"pred.CEC.soil"]-cal.mea.pre[,"CEC.soil"])^2)/(length(cal.mea.pre[,"CEC.soil"])-cal.ph.lv-1)),2);
# RPD:
cal.stat["CEC.soil","RPD"]<-round(sd(cal.mea.pre[,"CEC.soil"])/cal.stat["CEC.soil","RMSEC"],2);
# n LV:
cal.stat["CEC.soil","n LV"]<-cal.cecsoil.lv;
# n bc out:
cal.stat["CEC.soil","n bc out"]<-length(cal.cecsoil.bc.pred.out);
}

# Calculate val-subset statistics:
val.mea.pre[,"pred.bc.CEC.soil"]<-predict.mvr(cal.cecsoil,newdata=Deri.Smo[val,],ncomp=cal.cecsoil.lv);
val.mea.pre[,"pred.CEC.soil"]<-(val.mea.pre[,"pred.bc.CEC.soil"]*lambda[6]+1)^(1/lambda[6]);
val.cecsoil.bc.pred.out<-as.numeric(which(val.mea.pre[,"pred.CEC.soil"]=="NaN"));
#plot(val.mea.pre[,"pred.CEC.soil"]~val.mea.pre[,"CEC.soil"]);
if(length(val.cecsoil.bc.pred.out)!=0)
{
# R2:
val.stat["CEC.soil","r2"]<-round(cor(val.mea.pre[-val.cecsoil.bc.pred.out,"CEC.soil"],val.mea.pre[-val.cecsoil.bc.pred.out,"pred.CEC.soil"])^2,2);
# a:
val.stat["CEC.soil","a"]<-round(lm(pred.CEC.soil~CEC.soil,data=val.mea.pre)$coefficients[2],2);
# Bias:
residual<-val.mea.pre[-val.cecsoil.bc.pred.out,"pred.CEC.soil"]-val.mea.pre[-val.cecsoil.bc.pred.out,"CEC.soil"];
val.stat["CEC.soil","bias"]<-round(sum(residual)/length(val.mea.pre[-(length(val.cecsoil.bc.pred.out)),"CEC.soil"]),2)
# RMSEP:
val.stat["CEC.soil","RMSEP"]<-
round(sqrt(sum((val.mea.pre[-val.cecsoil.bc.pred.out,"pred.CEC.soil"]-val.mea.pre[-val.cecsoil.bc.pred.out,"CEC.soil"])^2)/(length(val.mea.pre[-val.cecsoil.bc.pred.out,"CEC.soil"]))),2);
# RPD:
val.stat["CEC.soil","RPD"]<-round(sd(val.mea.pre[-val.cecsoil.bc.pred.out,"CEC.soil"])/val.stat["CEC.soil","RMSEP"],2);
# n bc out:
val.stat["CEC.soil","n bc out"]<-length(val.cecsoil.bc.pred.out);
}
if(length(val.cecsoil.bc.pred.out)==0)
{
# R2:
val.stat["CEC.soil","r2"]<-round(cor(val.mea.pre[,"CEC.soil"],val.mea.pre[,"pred.CEC.soil"])^2,2);
# a:
val.stat["CEC.soil","a"]<-round(lm(pred.CEC.soil~CEC.soil,data=val.mea.pre)$coefficients[2],2);
# Bias:
residual<-val.mea.pre[,"pred.CEC.soil"]-val.mea.pre[,"CEC.soil"];
val.stat["CEC.soil","bias"]<-round(sum(residual)/length(val.mea.pre[,"CEC.soil"]),2)
# RMSEP:
val.stat["CEC.soil","RMSEP"]<-
round(sqrt(sum((val.mea.pre[,"pred.CEC.soil"]-val.mea.pre[,"CEC.soil"])^2)/(length(val.mea.pre[,"CEC.soil"]))),2);
# RPD:
val.stat["CEC.soil","RPD"]<-round(sd(val.mea.pre[,"CEC.soil"])/val.stat["CEC.soil","RMSEP"],2);
# n bc out:
val.stat["CEC.soil","n bc out"]<-length(val.cecsoil.bc.pred.out);
}

# Make a table containing the prediction bc outlier:

cal.all.bc.pred.out<-list();
cal.all.bc.pred.out[[2]]<-cal.ca.bc.pred.out;
cal.all.bc.pred.out[[3]]<-cal.mg.bc.pred.out;
cal.all.bc.pred.out[[4]]<-cal.orgc.bc.pred.out;
cal.all.bc.pred.out[[5]]<-cal.clay.bc.pred.out;
cal.all.bc.pred.out[[6]]<-cal.tots.bc.pred.out;
cal.all.bc.pred.out[[7]]<-cal.cecsoil.bc.pred.out;
cal.all.bc.pred.out[[9]]<-0;

val.all.bc.pred.out<-list();
val.all.bc.pred.out[[2]]<-val.ca.bc.pred.out;
val.all.bc.pred.out[[3]]<-val.mg.bc.pred.out;
val.all.bc.pred.out[[4]]<-val.orgc.bc.pred.out;
val.all.bc.pred.out[[5]]<-val.clay.bc.pred.out;
val.all.bc.pred.out[[6]]<-val.tots.bc.pred.out;
val.all.bc.pred.out[[7]]<-val.cecsoil.bc.pred.out;
val.all.bc.pred.out[[9]]<-0;

# Change % values to g kg-1 (for cal and val set for Corg, Sand and Clay):

val.mea.pre.1<-val.mea.pre;
val.mea.pre.1[,"Corg"]<-val.mea.pre.1[,"Corg"]*10;
val.mea.pre.1[,"pred.Corg"]<-val.mea.pre.1[,"pred.Corg"]*10;

val.mea.pre.1[,"Clay"]<-val.mea.pre.1[,"Clay"]*10;
val.mea.pre.1[,"pred.Clay"]<-val.mea.pre.1[,"pred.Clay"]*10;

val.mea.pre.1[,"Sand"]<-val.mea.pre.1[,"Sand"]*10;
val.mea.pre.1[,"pred.Sand"]<-val.mea.pre.1[,"pred.Sand"]*10;

cal.mea.pre.1<-cal.mea.pre;
cal.mea.pre.1[,"Corg"]<-cal.mea.pre.1[,"Corg"]*10;
cal.mea.pre.1[,"pred.Corg"]<-cal.mea.pre.1[,"pred.Corg"]*10;

cal.mea.pre.1[,"Clay"]<-cal.mea.pre.1[,"Clay"]*10;
cal.mea.pre.1[,"pred.Clay"]<-cal.mea.pre.1[,"pred.Clay"]*10;

cal.mea.pre.1[,"Sand"]<-cal.mea.pre.1[,"Sand"]*10;
cal.mea.pre.1[,"pred.Sand"]<-cal.mea.pre.1[,"pred.Sand"]*10;

# Calculate the statistics (percentiles and sd for all soil variables):

ref<-rbind(cal.mea.pre.1[,c("pH2","Ca","Mg","Corg","Clay","Sand","CEC.soil")],val.mea.pre.1[,c("pH2","Ca","Mg","Corg","Clay","Sand","CEC.soil")]);

# Make Regression of measured against predicted for val samples (Fig. 3):

pred.variables<-c("pred.pH2","pred.Ca","pred.Mg","pred.Corg","pred.Clay","pred.Sand","pred.CEC.soil");
x.min<-c();
for(i in 1:length(variables)){
	x.min[i]<-min(val.mea.pre.1[,variables[i]])
	}
x.max<-c();
for(i in 1:length(variables)){
	x.max[i]<-max(val.mea.pre.1[,variables[i]])
	}
y.min<-c();
for(i in 1:length(variables)){
	if(length(val.all.bc.pred.out[[i]])!=0)
	{y.min[i]<-min(val.mea.pre.1[-val.all.bc.pred.out[[i]],pred.variables[i]])}
	else
	{y.min[i]<-min(val.mea.pre.1 	[,pred.variables[i]])}
	}
y.max<-c();
for(i in 1:length(variables)){
	if(length(val.all.bc.pred.out[[i]])!=0)
	{y.max[i]<-max(val.mea.pre.1[-val.all.bc.pred.out[[i]],pred.variables[i]])}
	else
	{y.max[i]<-max(val.mea.pre.1 	[,pred.variables[i]])}
	}
max<-ifelse(x.max>=y.max,x.max,y.max);
min<-ifelse(x.min<=y.min,x.min,y.min);

var.names<-c("pH value",expression(paste("Exch. Ca content (",cmol[c]," ",kg^-1,")")),expression(paste("Exch. Mg content (",cmol[c]," ",kg^-1,")")),expression(paste("Organic C content (g ",kg^-1,")")),expression(paste("Clay content (g ",kg^-1,")")),expression(paste("Sand content (g ",kg^-1,")")),expression(paste("CEC (",cmol[c]," ",kg^-1,")")));
stat.val.r2<-as.character(val.stat[,"r2"]);
stat.val.rpd<-val.stat[,"RPD"];
x11(width=6,height=12);
par(mfrow=c(4,2),mar=c(4,4.2,1,1));
for(i in 1:length(variables)){
	plot(val.mea.pre.1[,pred.variables[i]]~val.mea.pre.1[,variables[i]],xlab="Measured values",ylab="Predicted values",ylim=c(min[i],max[i]),xlim=c(min[i],max[i]),cex=0.3,cex.axis=1.4,cex.lab=1.4);
abline(lm(val.mea.pre.1[,pred.variables[i]]~val.mea.pre.1[,variables[i]]));
text(min[i]-((max[i]-min[i])*0.07),max[i]-((max[i]-min[i])*0.03),pos=4,labels=var.names[i],cex=1.4);
text(max[i]+((max[i]-min[i])*0.07),min[i]+((max[i]-min[i])*0.1),pos=2,labels=paste("R2 = ",stat.val.r2[i]),cex=1.4);
text(max[i]+((max[i]-min[i])*0.07),min[i]+((max[i]-min[i])*0.01),pos=2,labels=paste("RPD = ",stat.val.rpd[i]),cex=1.4);
	}

# Determine variable importance:

# Get the regression coefficients:

rc.ph<-coef(cal.ph,comps=cal.ph.lv);
rc.ca<-coef(cal.ca,comps=cal.ca.lv);
rc.orgc<-coef(cal.orgc,comps=cal.orgc.lv);
rc.clay<-coef(cal.clay,comps=cal.clay.lv);
rc.tots<-coef(cal.tots,comps=cal.tots.lv);
rc.cecsoil<-coef(cal.cecsoil,comps=cal.cecsoil.lv);
rc<-list();
rc[[1]]<-rc.ph;
rc[[2]]<-rc.orgc;
rc[[3]]<-rc.ca;
rc[[4]]<-rc.cecsoil;
rc[[5]]<-rc.clay;
rc[[6]]<-rc.tots;

k<-ncol(Deri.Smo)

## Get the VIP:

# pH value:

lw.ph<-loading.weights(cal.ph)[,cal.ph.lv];
ssa.ph<-sum((cal.ph$fitted.values[,1,cal.ph.lv]-cal.ph$Ymeans)^2);
sst.ph<-sum((cal.ph$fitted.values[,1,dim(cal.ph$fitted.values[,,])[2]]-cal.ph$Ymeans)^2)
vip.ph<-c();
for(i in 1:length(lw.ph)){
	vip.ph[i]<-k*sum(lw.ph[i]^2*(ssa.ph/sst.ph));
	}
	
i.w.ph<-which(abs(rc.ph)>=sd(rc.ph) & vip.ph>1);

# Exch. Ca content:

lw.ca<-loading.weights(cal.ca)[,cal.ca.lv];
ssa.ca<-sum((cal.ca$fitted.values[,1,cal.ca.lv]-cal.ca$Ymeans)^2);
sst.ca<-sum((cal.ca$fitted.values[,1,dim(cal.ca$fitted.values[,,])[2]]-cal.ca$Ymeans)^2)
vip.ca<-c();
for(i in 1:length(lw.ca)){
	vip.ca[i]<-k*sum(lw.ca[i]^2*(ssa.ca/sst.ca));
	}
	
i.w.ca<-which(abs(rc.ca)>=sd(rc.ca) & vip.ca>1);

# OrgC content:

lw.orgc<-loading.weights(cal.orgc)[,cal.orgc.lv];
ssa.orgc<-sum((cal.orgc$fitted.values[,1,cal.orgc.lv]-cal.orgc$Ymeans)^2);
sst.orgc<-sum((cal.orgc$fitted.values[,1,dim(cal.orgc$fitted.values[,,])[2]]-cal.orgc$Ymeans)^2)
vip.orgc<-c();
for(i in 1:length(lw.orgc)){
	vip.orgc[i]<-k*sum(lw.orgc[i]^2*(ssa.orgc/sst.orgc));
	}
	
i.w.orgc<-which(abs(rc.orgc)>=sd(rc.orgc) & vip.orgc>1);

# Clay content:

lw.clay<-loading.weights(cal.clay)[,cal.clay.lv];
ssa.clay<-sum((cal.clay$fitted.values[,1,cal.clay.lv]-cal.clay$Ymeans)^2);
sst.clay<-sum((cal.clay$fitted.values[,1,dim(cal.clay$fitted.values[,,])[2]]-cal.clay$Ymeans)^2)
vip.clay<-c();
for(i in 1:length(lw.clay)){
	vip.clay[i]<-k*sum(lw.clay[i]^2*(ssa.clay/sst.clay));
	}
	
i.w.clay<-which(abs(rc.clay)>=sd(rc.clay) & vip.clay>1);

# Sand content:

lw.sand<-loading.weights(cal.tots)[,cal.tots.lv];
ssa.sand<-sum((cal.tots$fitted.values[,1,cal.tots.lv]-cal.tots$Ymeans)^2);
sst.sand <-sum((cal.tots$fitted.values[,1,dim(cal.tots$fitted.values[,,])[2]]-cal.tots$Ymeans)^2)
vip.sand<-c();
for(i in 1:length(lw.sand)){
	vip.sand[i]<-k*sum(lw.sand[i]^2*(ssa.sand/sst.sand));
	}
	
i.w.sand<-which(abs(rc.tots)>=sd(rc.tots) & vip.sand>1);

# CEC content:

lw.cecsoil<-loading.weights(cal.cecsoil)[,cal.cecsoil.lv];
ssa.cecsoil<-sum((cal.cecsoil$fitted.values[,1,cal.cecsoil.lv]-cal.cecsoil$Ymeans)^2);
sst.cecsoil <-sum((cal.cecsoil$fitted.values[,1,dim(cal.cecsoil$fitted.values[,,])[2]]-cal.cecsoil$Ymeans)^2)
vip.cecsoil<-c();
for(i in 1:length(lw.cecsoil)){
	vip.cecsoil[i]<-k*sum(lw.cecsoil[i]^2*(ssa.cecsoil/sst.cecsoil));
	}
	
i.w.cecsoil<-which(abs(rc.cecsoil)>=sd(rc.cecsoil) & vip.cecsoil>1);

i.w<-list();
i.w[[1]]<-i.w.ph;
i.w[[2]]<-i.w.orgc;
i.w[[3]]<-i.w.ca
i.w[[4]]<-i.w.cecsoil;
i.w[[5]]<-i.w.clay;
i.w[[6]]<-i.w.sand

# Plot variable selection (Fig. 4):

sp<-spec[100,];
sp<-sp-min(sp);

x11(width=6,height=8);
par(mfrow=c(3,2),mar=c(4,4.5,1,1));

plot(rev(rc[[1]])~rev(waveb),type="l",xlim=c(max(waveb),min(waveb)),ylim=c(min(rc[[1]])*1.5,max(rc[[1]])*1.3),ylab=expression(paste("PLSR ",b[coeff]," pH")),cex.lab=1.3,xlab="");
ra.ph<-abs(min(rc[[1]]*1.5)-min(rc[[1]]))/max(sp);
sp.ph<-sp*ra.ph-abs(min(rc[[1]]*1.5));
lines((rev(sp.ph)~rev(waveb)));
points(rev(sp.ph[i.w[[1]]])~rev(waveb[i.w[[1]]]),pch=19,cex=0.5)
text(rev(waveb)[10],max(rc[[1]])*1.25,pos=4,labels="(a)",cex=1.3);
text(rev(waveb)[length(waveb)],max(rc[[1]])*1.25,pos=2,labels=paste("VIP > 1; ",round(sd(rc[[1]]),1)," < bcoeff < -",round(sd(rc[[1]]),1),sep=""),cex=1.2)

plot(rev(rc[[2]])~rev(waveb),type="l",xlim=c(max(waveb),min(waveb)),ylim=c(min(rc[[2]])*1.5,max(rc[[2]])*1.3),ylab=expression(paste("PLSR ",b[coeff]," Organic C")),cex.lab=1.3,xlab="");
ra.orgc<-abs(min(rc[[2]]*1.5)-min(rc[[2]]))/max(sp);
sp.orgc<-sp*ra.orgc-abs(min(rc[[2]]*1.5));
lines((rev(sp.orgc)~rev(waveb)));
points(rev(sp.orgc[i.w[[2]]])~rev(waveb[i.w[[2]]]),pch=19,cex=0.5)
text(rev(waveb)[10],max(rc[[2]])*1.25,pos=4,labels="(b)",cex=1.3);
text(rev(waveb)[length(waveb)],max(rc[[2]])*1.25,pos=2,labels=paste("VIP > 1; ",round(sd(rc[[2]]),1)," < bcoeff < -",round(sd(rc[[2]]),1),sep=""),cex=1.3)

plot(rev(rc[[3]])~rev(waveb),type="l",xlim=c(max(waveb),min(waveb)),ylim=c(min(rc[[3]])*1.5,max(rc[[3]])*1.3),ylab=expression(paste("PLSR ",b[coeff]," Exch. Ca")),cex.lab=1.3,xlab="");
ra.ca<-abs(min(rc[[3]]*1.5)-min(rc[[3]]))/max(sp);
sp.ca<-sp*ra.ca-abs(min(rc[[3]]*1.5));
lines((rev(sp.ca)~rev(waveb)));
points(rev(sp.ca[i.w[[3]]])~rev(waveb[i.w[[3]]]),pch=19,cex=0.5)
text(rev(waveb)[10],max(rc[[3]])*1.25,pos=4,labels="(c)",cex=1.3);
text(rev(waveb)[length(waveb)],max(rc[[3]])*1.25,pos=2,labels=paste("VIP > 1; ",round(sd(rc[[3]]),1)," < bcoeff < -",round(sd(rc[[3]]),1),sep=""),cex=1.3)

plot(rev(rc[[4]])~rev(waveb),type="l",xlim=c(max(waveb),min(waveb)),ylim=c(min(rc[[4]])*1.5,max(rc[[4]])*1.3),ylab=expression(paste("PLSR ",b[coeff]," CEC")),cex.lab=1.3,xlab="");
ra.cec<-abs(min(rc[[4]]*1.5)-min(rc[[4]]))/max(sp);
sp.cec<-sp*ra.cec-abs(min(rc[[4]]*1.5));
lines((rev(sp.cec)~rev(waveb)));
points(rev(sp.cec[i.w[[4]]])~rev(waveb[i.w[[4]]]),pch=19,cex=0.5)
text(rev(waveb)[10],max(rc[[4]])*1.25,pos=4,labels="(d)",cex=1.3);
text(rev(waveb)[length(waveb)],max(rc[[4]])*1.25,pos=2,labels=paste("VIP > 1; ",round(sd(rc[[4]]),1)," < bcoeff < -",round(sd(rc[[4]]),1),sep=""),cex=1.3)

par(mar=c(4.5,4.5,1,1))

plot(rev(rc[[5]])~rev(waveb),type="l",xlim=c(max(waveb),min(waveb)),ylim=c(min(rc[[5]])*1.5,max(rc[[5]])*1.3),ylab=expression(paste("PLSR ",b[coeff]," Clay")),cex.lab=1.3,xlab=expression(paste("Wavenumber (",cm^-1,")")));
ra.clay<-abs(min(rc[[5]]*1.5)-min(rc[[5]]))/max(sp);
sp.clay<-sp*ra.clay-abs(min(rc[[5]]*1.5));
lines((rev(sp.clay)~rev(waveb)));
points(rev(sp.clay[i.w[[5]]])~rev(waveb[i.w[[5]]]),pch=19,cex=0.5)
text(rev(waveb)[10],max(rc[[5]])*1.25,pos=4,labels="(e)",cex=1.3);
text(rev(waveb)[length(waveb)],max(rc[[5]])*1.25,pos=2,labels=paste("VIP > 1; ",round(sd(rc[[5]]),1)," < bcoeff < -",round(sd(rc[[5]]),1),sep=""),cex=1.3)

plot(rev(rc[[6]])~rev(waveb),type="l",xlim=c(max(waveb),min(waveb)),ylim=c(min(rc[[6]])*1.5,max(rc[[6]])*1.3),ylab=expression(paste("PLSR ",b[coeff]," Sand")),cex.lab=1.3,xlab=expression(paste("Wavenumber (",cm^-1,")")));
ra.sand<-abs(min(rc[[6]]*1.5)-min(rc[[6]]))/max(sp);
sp.sand<-sp*ra.sand-abs(min(rc[[6]]*1.5));
lines((rev(sp.sand)~rev(waveb)));
points(rev(sp.sand[i.w[[6]]])~rev(waveb[i.w[[6]]]),pch=19,cex=0.5)
text(rev(waveb)[10],max(rc[[6]])*1.25,pos=4,labels="(f)",cex=1.3);
text(rev(waveb)[length(waveb)],max(rc[[6]])*1.25,pos=2,labels=paste("VIP > 1; ",round(sd(rc[[6]]),1)," < bcoeff < -",round(sd(rc[[6]]),1),sep=""),cex=1.3)

# Prepare output:

output<-list(Table_1=cum,Table_2_summary=chem.stat.all,Table_2_sd=chem.stat.sd.all,Table_3_calibration=cal.stat,Table_3_validation=val.stat,Table_4_all=chem.cor.all,Table_4_calibration=chem.cor.cal,Calibration_predicted=cal.mea.pre.1,Validation_predicted=val.mea.pre.1);
class(output)<-"isric";
return(output);

}