pred.regr <-
function(new,model,output.name="test",s="NL"){
	
# Check function body:

if(is.na(match(class(new),c("data.frame","matrix")))){stop("Invalid argument: 'new' has to be of class 'data frame' or 'matrix'.")}
if(class(new)=="data.frame"){new<-as.matrix(new)};
if(class(as.numeric(colnames(new)))!="numeric"){stop("Invalid argument: colnames of 'new' can't be coerced to class 'numeric'.")}
if(as.numeric(colnames(new)[1])>as.numeric(colnames(new)[2])){
	test<-new;
	for(i in 1:nrow(new)){
		test[i,]<-rev(test[i,]);
		}
		colnames(test)<-rev(colnames(test));
		new<-test;
		rm(test);
	}
if(as.numeric(colnames(new)[1])==model$wavebands[1] & length(colnames(new))==length(model$wavebands))
{	new<-new;}
if(as.numeric(colnames(new)[1])!=model$wavebands[1]){
	test<-matrix(nrow=nrow(new),ncol=ncol(new),dimnames=list(rownames(new),colnames(new)));
	for(i in 1:nrow(new)){
	test[i,]<-round(spline(as.numeric(colnames(new)),new[i,],method="natural",xout=model$wavebands)[[2]],6);
	}
	new<-test;
	rm(test);
}

if(is.na(match(c("x.tr"),names(model)[3]))){stop("Invalid argument: 'model' is not of class 'regr'.")};
class(model)<-"regr";
if(class(output.name)!="character"){stop("Invalid argument: 'output.name' has to be of class 'character'.")}

# Make spectral transformation:

if(model$spectral.transformation=="raw"){
	new.tr<-new;
	}

if(model$spectral.transformation=="derivative"){
	new.tr<-matrix(nrow=nrow(new),ncol=ncol(new),dimnames=list(rownames(new),colnames(new)));
waveb<-as.numeric(colnames(new));
#require(KernSmooth,quietly=T);
for(i in 1:nrow(new)){
	new.tr[i,]<-locpoly(waveb,new[i,],drv=model$drv,bandwidth=model$bandwidth,gridsize=ncol(new))[[2]]
	}
    #detach(package:KernSmooth);
	}

if(model$spectral.transformation=="continuum removed"){
	new.tr<-matrix(nrow=nrow(new),ncol=ncol(new),dimnames=list(rownames(new),colnames(new)));
waveb<-as.numeric(colnames(new));
#require(KernSmooth,quietly=T);	
test<-new;
for(i in 1:nrow(new)){
	test.1<-cbind(waveb,test[i,]);
	test.1<-sortedXyData(test.1[,1],test.1[,2]);
	ch<-chull(test.1);
	ch.1<-ch;
	ch<-ch[1:(which(ch==1))]
	ch<-sort(ch);
	ch<-c(ch,ncol(new));
	appr.ch<-approx(test.1[ch,],xout=test.1[,1],method="linear",ties="mean");
	cr<-test.1[[2]]-appr.ch[[2]];
	new.tr[i,]<-cr;
	}
new.tr<-new.tr[,2:(ncol(new)-2)];
#detach(package:KernSmooth);
	}

if(model$spectral.transformation=="wavelet transformed"){
	waveb<-as.numeric(colnames(new));
waveb.1024.up<-round(max(waveb));
waveb.1024.down<-round(min(waveb));
waveb.1024.n<-1023;
waveb.1024.step<-(waveb.1024.up-waveb.1024.down)/waveb.1024.n;
waveb.1024<-c();
waveb.1024[1]<-waveb.1024.down;
for(i in 2:1024){
	waveb.1024[i]<-round(waveb.1024.down+(i-1)*waveb.1024.step,5)
	}
new.comp<-matrix(nrow=nrow(new),ncol=length(waveb.1024),dimnames=list(rownames(new),waveb.1024));
for(i in 1:nrow(new)){
	new.comp[i,]<-round(spline(waveb,new[i,],method="natural",xout=waveb.1024)[[2]],6);
}
#require(wavelets,quietly=T);
new.tr<-matrix(nrow=nrow(new.comp),ncol=2^(10-model$level),dimnames=list(rownames(new.comp),paste("WC_",c(1:2^(10-model$level)),sep="")));
for(i in 1:nrow(new.tr)){
	blub<-dwt(new.comp[i,],filter=model$filte);
	new.tr[i,]<-slot(blub,"W")[[model$level]]
	}
#detach(package:wavelets);
	}
	
# Make preparations for predictions:

test<-c();
mah.test<-c();
new.tr.l<-list();
pred.l<-list();
tr<-c();
pred.variables<-c(c(paste("pred.",model$constituents[1],sep=""),paste(model$constituents[1],".lower.confidence",sep=""),paste(model$constituents[1],".upper.confidence",sep="")));
if(length(model$constituents)>1){
	for(i in 2:length(model$constituents)){
	bla<-c(paste("pred.",model$constituents[i],sep=""),paste(model$constituents[i],".lower.confidence",sep=""),paste(model$constituents[i],".upper.confidence",sep=""));
	pred.variables<-c(pred.variables,bla);
	}
	}
pred<-matrix(nrow=nrow(new.tr),ncol=length(pred.variables),dimnames=list(rownames(new.tr),pred.variables));

# Test if the spectra fall within the calibration library mahalanobis range:

for(i in 1:length(model$constituents)){
test<-predict(model$cal.pca[[i]],new.tr);
if(length(model$mahalanobis)>0){
mah.test<-sqrt(mahalanobis(test,colMeans(model$cal.pca[[i]]$x),cov(model$cal.pca[[i]]$x)));
out<-which(mah.test>max(model$mahalanobis[[i]]));
if(length(out)!=0){new.tr.l[[i]]<-new.tr;new.tr.l[[i]][out,]<-NA}
if(length(out)==0){new.tr.l[[i]]<-new.tr}
};
if(length(model$mahalanobis)==0){
	new.tr.l[[i]]<-new.tr;
	}


# Test if the samples have local neighbours (remove when there is no local neighbour within a distance of 5):

out<-c();
if(length(mah.test)>0){
for(j in 1:length(mah.test)){
	out[j]<-min(abs(mah.test[j]-model$mahalanobis[[i]]));
	}
out<-which(out>5);
if(length(out)!=0){new.tr.l[[i]][out,]<-NA};
if(dim(na.omit(new.tr.l[[i]]))[1]==0){stop("No sample is represented by the calibration set.")}
}

if(model$constituents.transformation[i]=="Untransformed"){tr[i]<-""};
if(model$constituents.transformation[i]=="Box-Cox transformed"){tr[i]<-"bc."};
if(model$constituents.transformation[i]=="Square root transformed"){tr[i]<-"sqrt."};
if(model$constituents.transformation[i]=="Log transformed"){tr[i]<-"log."};

}

# Make predictions:
	
if(model$method=="pls"){

#require(pls,quietly=T);

for(i in 1:length(model$constituents)){
	if(tr[i]==""){pred.l[[i]]<-matrix(nrow=nrow(new.tr.l[[i]]),ncol=1,dimnames=list(rownames(new.tr.l[[i]]),c(paste("pred.",model$constituents[i],sep=""))));}
if(tr[i]!=""){pred.l[[i]]<-matrix(nrow=nrow(new.tr.l[[i]]),ncol=2,dimnames=list(rownames(new.tr.l[[i]]),c(paste("pred.",tr[i],model$constituents[i],sep=""),paste("pred.",model$constituents[i],sep=""))));}
if(tr[i]==""){pred.l[[i]][rownames(na.omit(new.tr.l[[i]])),paste("pred.",tr[i],model$constituents[i],sep="")]<-predict(model$model[[i]],na.omit(new.tr.l[[i]]),ncomp=model$model[[i]]$ncomp)}
if(tr[i]=="bc."){pred.l[[i]][rownames(na.omit(new.tr.l[[i]])),paste("pred.",tr[i],model$constituents[i],sep="")]<-predict(model$model[[i]],na.omit(new.tr.l[[i]]),ncomp=model$model[[i]]$ncomp);
pred.l[[i]][,paste("pred.",model$constituents[i],sep="")]<-(pred.l[[i]][,paste("pred.",tr[i],model$constituents[i],sep="")]*model$lambda[model$constituents[i]][model$constituents[i]]+1)^(1/model$lambda[model$constituents[i]][model$constituents[i]])-1;}
if(tr[i]=="sqrt."){pred.l[[i]][rownames(na.omit(new.tr.l[[i]])),paste("pred.",tr[i],model$constituents[i],sep="")]<-predict(model$model[[i]],na.omit(new.tr.l[[i]]),ncomp=model$model[[i]]$ncomp);
pred.l[[i]][,paste("pred.",model$constituents[i],sep="")]<-pred.l[[i]][rownames(new.tr.l[[i]]),paste("pred.",tr[i],model$constituents[i],sep="")]^2;
}
if(tr[i]=="log."){pred.l[[i]][rownames(na.omit(new.tr.l[[i]])),paste("pred.",tr[i],model$constituents[i],sep="")]<-predict(model$model[[i]],na.omit(new.tr.l[[i]]),ncomp=model$model[[i]]$ncomp);
pred.l[[i]][,paste("pred.",model$constituents[i],sep="")]<-exp(pred.l[[i]][rownames(new.tr.l[[i]]),paste("pred.",tr[i],model$constituents[i],sep="")])-1}
	}

#detach(package:pls);	
	}
	
if(model$method=="brt"){
	
#require(splines,quietly=T);
#require(survival,quietly=T);
#require(lattice,quietly=T);
#require(gbm,quietly=T);



for(i in 1:length(model$constituents)){
if(tr[i]==""){pred.l[[i]]<-matrix(nrow=nrow(new.tr.l[[i]]),ncol=1,dimnames=list(rownames(new.tr.l[[i]]),c(paste("pred.",model$constituents[i],sep=""))));}
if(tr[i]!=""){pred.l[[i]]<-matrix(nrow=nrow(new.tr.l[[i]]),ncol=2,dimnames=list(rownames(new.tr.l[[i]]),c(paste("pred.",tr[i],model$constituents[i],sep=""),paste("pred.",model$constituents[i],sep=""))));}
if(tr[i]==""){pred.l[[i]][rownames(na.omit(new.tr.l[[i]])),paste("pred.",tr[i],model$constituents[i],sep="")]<-na.omit(predict(model$model[[i]],na.omit(new.tr.l[[i]]),n.trees=model$cal.stat[model$constituents[i],"n trees"],type.convert="response"))}
if(tr[i]=="bc."){pred.l[[i]][rownames(na.omit(new.tr.l[[i]])),paste("pred.",tr[i],model$constituents[i],sep="")]<-na.omit(predict(model$model[[i]],na.omit(new.tr.l[[i]]),n.trees=model$cal.stat[model$constituents[i],"n trees"],type="response"));
pred.l[[i]][,paste("pred.",model$constituents[i],sep="")]<-(pred.l[[i]][,paste("pred.",tr[i],model$constituents[i],sep="")]*model$lambda[model$constituents[i]][model$constituents[i]]+1)^(1/model$lambda[model$constituents[i]][model$constituents[i]])-1;}
if(tr[i]=="sqrt."){pred.l[[i]][rownames(na.omit(new.tr.l[[i]])),paste("pred.",tr[i],model$constituents[i],sep="")]<-na.omit(predict(model$model[[i]],na.omit(new.tr.l[[i]]),n.trees=model$cal.stat[model$constituents[i],"n trees"],type="response"));
pred.l[[i]][,paste("pred.",model$constituents[i],sep="")]<-pred.l[[i]][rownames(new.tr.l[[i]]),paste("pred.",tr[i],model$constituents[i],sep="")]^2;
}
if(tr[i]=="log."){pred.l[[i]][rownames(na.omit(new.tr.l[[i]])),paste("pred.",tr[i],model$constituents[i],sep="")]<-na.omit(predict(model$model[[i]],na.omit(new.tr.l[[i]]),n.trees=model$cal.stat[model$constituents[i],"n trees"],type="response"));
pred.l[[i]][,paste("pred.",model$constituents[i],sep="")]<-exp(pred.l[[i]][rownames(new.tr.l[[i]]),paste("pred.",tr[i],model$constituents[i],sep="")])-1}
	}

#detach(package:gbm);
#detach(package:splines);
#detach(package:survival);
#detach(package:lattice);
	
	}

if(model$method=="svm"){

#require(class,quietly=T);
#require(e1071,quietly=T);
for(i in 1:length(model$constituents)){
	if(tr[i]==""){pred.l[[i]]<-matrix(nrow=nrow(na.omit(new.tr.l[[i]])),ncol=1,dimnames=list(rownames(na.omit(new.tr.l[[i]])),c(paste("pred.",model$constituents[i],sep=""))));}
if(tr[i]!=""){pred.l[[i]]<-matrix(nrow=nrow(na.omit(new.tr.l[[i]])),ncol=2,dimnames=list(rownames(na.omit(new.tr.l[[i]])),c(paste("pred.",tr[i],model$constituents[i],sep=""),paste("pred.",model$constituents[i],sep=""))));}
if(tr[i]==""){pred.l[[i]][rownames(na.omit(new.tr.l[[i]])),paste("pred.",tr[i],model$constituents[i],sep="")]<-predict(model$model[[i]],na.omit(new.tr.l[[i]]))}
if(tr[i]=="bc."){pred.l[[i]][rownames(na.omit(new.tr.l[[i]])),paste("pred.",tr[i],model$constituents[i],sep="")]<-predict(model$model[[i]],na.omit(new.tr.l[[i]]));
pred.l[[i]][,paste("pred.",model$constituents[i],sep="")]<-(pred.l[[i]][,paste("pred.",tr[i],model$constituents[i],sep="")]*model$lambda[model$constituents[i]][model$constituents[i]]+1)^(1/model$lambda[model$constituents[i]][model$constituents[i]])-1;}
if(tr[i]=="sqrt."){pred.l[[i]][rownames(na.omit(new.tr.l[[i]])),paste("pred.",tr[i],model$constituents[i],sep="")]<-predict(model$model[[i]],na.omit(new.tr.l[[i]]));
pred.l[[i]][,paste("pred.",model$constituents[i],sep="")]<-pred.l[[i]][rownames(new.tr.l[[i]]),paste("pred.",tr[i],model$constituents[i],sep="")]^2;
}
if(tr[i]=="log."){pred.l[[i]][rownames(na.omit(new.tr.l[[i]])),paste("pred.",tr[i],model$constituents[i],sep="")]<-predict(model$model[[i]],na.omit(new.tr.l[[i]]));
pred.l[[i]][,paste("pred.",model$constituents[i],sep="")]<-exp(pred.l[[i]][rownames(new.tr.l[[i]]),paste("pred.",tr[i],model$constituents[i],sep="")])-1}
}

#detach(package:class);
#detach(package:e1071);
	
	}
	
# Test if the predictions are within the calibration set range:
## Test if the right Batch_Labid is used (due to na.omit):

for(i in 1:length(model$constituents)){
mi<-which(pred.l[[i]][,paste("pred.",model$constituents[i],sep="")]<model$cal.range[[i]][1])
ma<-which(pred.l[[i]][,paste("pred.",model$constituents[i],sep="")]>model$cal.range[[i]][2]);
out<-c(mi,ma);
if(length(out)!=0){new.tr.l[[i]][out,]<-NA};
if(dim(na.omit(new.tr.l[[i]]))[1]==0){stop("No sample is represented by the calibration set.")}

pred[rownames(pred.l[[i]]),colnames(pred.l[[i]])[max(length(colnames(pred.l[[i]])))]]<-pred.l[[i]][,max(length(colnames(pred.l[[i]])))]

# Uncertainty prediction:

a<-c();
mean.rmsep<-c();
for(j in 1:nrow(pred.l[[i]])){
a[j]<-which(model$val.mea.pre[as.numeric(names(model$lm[[i]])),model$constituents[i]]>pred.l[[i]][j,paste("pred.",model$constituents[i],sep="")])[1]
if(length(na.omit(a[j]))==0){mean.rmsep[j]<-NA}
if(length(na.omit(a[j]))!=0){
	if(a[j]==1){mean.rmsep[j]<-model$rmsep[[i]][a[j]]}
	if(a[j]!=1){mean.rmsep[j]<-mean(model$rmsep[[i]][a[j]-1],model$rmsep[[i]][a[j]],na.rm=T)}};
if(length(na.omit(a[j]))==0){pred[rownames(pred)[j],paste(model$constituents[i],".lower.confidence",sep="")]<-NA}
if(length(na.omit(a[j]))!=0){pred[rownames(pred)[j],paste(model$constituents[i],".lower.confidence",sep="")]<-pred[rownames(pred)[j],paste("pred.",model$constituents[i],sep="")]-mean.rmsep[j]};
if(length(na.omit(a[j]))==0){pred[rownames(pred)[j],paste(model$constituents[i],".upper.confidence",sep="")]<-NA}
if(length(na.omit(a[j]))!=0){pred[rownames(pred)[j],paste(model$constituents[i],".upper.confidence",sep="")]<-pred[rownames(pred)[j],paste("pred.",model$constituents[i],sep="")]+mean.rmsep[j];}
	}
}
pred<-round(pred,2);

# Write output to

if(s!="NULL"){
		setwd(s);
		}
test<-cbind(rownames(pred),pred);
colnames(test)[1]<-"Sample_ID";
write.table(test,file=paste("Predictions ",output.name,".csv",sep=""),sep=",",dec=".",row.names=F)

# Create output:

output<-list("output.name"=output.name,"predicted.values"=pred,"method"=model$method,"spectral.transformations"=model$spectral.transformation);
class(output)<-"pred.regr";
return(output);
		
		}
