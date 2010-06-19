prco <-
function(x,sav="FALSE",save.path="NULL",output.name="prco"){

if(is.na(match(sav,c("FALSE","TRUE")))){stop("Invalid argument: 'sav' has to either 'FALSE' or 'TRUE'.")};
if(class(output.name)!="character"){stop("Invalid argument: 'output.name' has to be of class 'character'.")};

x<-na.omit(x);
pca<-prcomp(x,scale=T);
cpv<-summary(pca)[[6]][3,1:20];
zzz<-matrix(nrow=1,ncol=length(cpv)-4);
for (i in 1:16){
	e<-(cpv[i]+0.04)<cpv[i+3];
	zzz[i]<-e
	}
pc<-(which(zzz==FALSE)-1)[1];
if(pc==0 | pc==1){prco<-pca$x};
if(pc>0){prco<-pca$x[,1:pc]};	

# Prepare output:

output<-list(prcomp=pca,prco=prco,i.pc=pc)
class(output)<-"prco"

# Plot PC spaces against each other:

x11();
if(output$i.pc==0){pairs(output$prcomp$x[,1:5])};
if(output$i.pc==1){pairs(output$prcomp$x[,1:5])};
if(output$i.pc>1){pairs(output$prco[,1:output$i.pc])}
x11();
plot(output$prcomp$rotation[,1]~as.numeric(rownames(output$prcomp$rotation)),ylab="Loading values",xlab="Wavebands (cm^-1)",type="l",ylim=c(min(output$prcomp$rotation[,1:2],na.rm=T),max(output$prcomp$rotation[,1:2],na.rm=T)));
lines(output$prcomp$rotation[,2]~as.numeric(rownames(output$prcomp$rotation)),col="red")
legend(x="bottomright",legend=c("First PC","Second PC"),col=c("black","red"),lty=1)


if(sav=="TRUE"){
	if(save.path!="NULL"){
		setwd(save.path)
		}
		save(output,file=output.name)
	}
return(output);
}

