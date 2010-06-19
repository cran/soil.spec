pc.space <-
function(base,new){
	
# Make PCA on base:

base<-na.omit(base);
pca<-prcomp(base,scale=T);
prco<-pca$x;

# Predict new in PC space of base:

ne<-predict(object=pca,newdata=new);

# Plot new in base:

x11(height=,width=9);
par(mfrow=c(5,5),mar=c(1,1,1,1));
plot(prco[,1]~prco[,1],type="n",xaxt="n",yaxt="n");
text(x=max((prco[,2]))-(max((prco[,2]))-min((prco[,2])))/2,y=max((prco[,1]))-(max((prco[,1]))-min((prco[,1])))/2,labels="PC1",cex=2,pos=4);
plot(prco[,1]~prco[,2],ylim=c(min(c(prco[,1],ne[,1])),max(c(prco[,1],ne[,1]))),xlim=c(min(c(prco[,2],ne[,2])),max(c(prco[,2],ne[,2]))));
points(ne[,1]~ne[,2],col="red")
plot(prco[,1]~prco[,3],ylim=c(min(c(prco[,1],ne[,1])),max(c(prco[,1],ne[,1]))),xlim=c(min(c(prco[,3],ne[,3])),max(c(prco[,3],ne[,3]))));
points(ne[,1]~ne[,3],col="red")
plot(prco[,1]~prco[,4],ylim=c(min(c(prco[,1],ne[,1])),max(c(prco[,1],ne[,1]))),xlim=c(min(c(prco[,4],ne[,4])),max(c(prco[,4],ne[,4]))));
points(ne[,1]~ne[,4],col="red")
plot(prco[,1]~prco[,5],ylim=c(min(c(prco[,1],ne[,1])),max(c(prco[,1],ne[,1]))),xlim=c(min(c(prco[,5],ne[,5])),max(c(prco[,5],ne[,5]))));
points(ne[,1]~ne[,5],col="red")
plot(prco[,2]~prco[,1],ylim=c(min(c(prco[,2],ne[,2])),max(c(prco[,2],ne[,2]))),xlim=c(min(c(prco[,1],ne[,1])),max(c(prco[,1],ne[,1]))));
points(ne[,2]~ne[,1],col="red")
plot(prco[,2]~prco[,2],type="n",xaxt="n",yaxt="n");
text(x=max((prco[,2]))-(max((prco[,2]))-min((prco[,2])))/2,y=max((prco[,2]))-(max((prco[,2]))-min((prco[,2])))/2,labels="PC2",cex=2,pos=4);

plot(prco[,2]~prco[,3],ylim=c(min(c(prco[,2],ne[,2])),max(c(prco[,2],ne[,2]))),xlim=c(min(c(prco[,3],ne[,3])),max(c(prco[,3],ne[,3]))));
points(ne[,2]~ne[,3],col="red")
plot(prco[,2]~prco[,4],ylim=c(min(c(prco[,2],ne[,2])),max(c(prco[,2],ne[,2]))),xlim=c(min(c(prco[,4],ne[,4])),max(c(prco[,4],ne[,4]))));
points(ne[,2]~ne[,4],col="red")
plot(prco[,2]~prco[,5],ylim=c(min(c(prco[,2],ne[,2])),max(c(prco[,2],ne[,2]))),xlim=c(min(c(prco[,5],ne[,5])),max(c(prco[,5],ne[,5]))));
points(ne[,2]~ne[,5],col="red")
plot(prco[,3]~prco[,1],ylim=c(min(c(prco[,3],ne[,3])),max(c(prco[,3],ne[,3]))),xlim=c(min(c(prco[,1],ne[,1])),max(c(prco[,1],ne[,1]))));
points(ne[,3]~ne[,1],col="red")
plot(prco[,3]~prco[,2],ylim=c(min(c(prco[,3],ne[,3])),max(c(prco[,3],ne[,3]))),xlim=c(min(c(prco[,2],ne[,2])),max(c(prco[,2],ne[,2]))));
points(ne[,3]~ne[,2],col="red")
plot(prco[,3]~prco[,3],type="n",xaxt="n",yaxt="n");
text(x=max((prco[,2]))-(max((prco[,2]))-min((prco[,2])))/2,y=max((prco[,3]))-(max((prco[,3]))-min((prco[,3])))/2,labels="PC3",cex=2,pos=4);
plot(prco[,3]~prco[,4],ylim=c(min(c(prco[,3],ne[,3])),max(c(prco[,3],ne[,3]))),xlim=c(min(c(prco[,4],ne[,4])),max(c(prco[,4],ne[,4]))));
points(ne[,3]~ne[,4],col="red")
plot(prco[,3]~prco[,5],ylim=c(min(c(prco[,3],ne[,3])),max(c(prco[,3],ne[,3]))),xlim=c(min(c(prco[,5],ne[,5])),max(c(prco[,5],ne[,5]))));
points(ne[,3]~ne[,5],col="red")
plot(prco[,4]~prco[,1],ylim=c(min(c(prco[,4],ne[,4])),max(c(prco[,4],ne[,4]))),xlim=c(min(c(prco[,1],ne[,1])),max(c(prco[,1],ne[,1]))));
points(ne[,4]~ne[,1],col="red")
plot(prco[,4]~prco[,2],ylim=c(min(c(prco[,4],ne[,4])),max(c(prco[,4],ne[,4]))),xlim=c(min(c(prco[,2],ne[,2])),max(c(prco[,2],ne[,2]))));
points(ne[,4]~ne[,2],col="red")
plot(prco[,4]~prco[,3],ylim=c(min(c(prco[,4],ne[,4])),max(c(prco[,4],ne[,4]))),xlim=c(min(c(prco[,3],ne[,3])),max(c(prco[,3],ne[,3]))));
points(ne[,4]~ne[,3],col="red")
plot(prco[,4]~prco[,4],type="n",xaxt="n",yaxt="n");
text(x=max((prco[,2]))-(max((prco[,2]))-min((prco[,2])))/2,y=max((prco[,4]))-(max((prco[,4]))-min((prco[,4])))/2,labels="PC4",cex=2,pos=4);
plot(prco[,4]~prco[,5],ylim=c(min(c(prco[,4],ne[,4])),max(c(prco[,4],ne[,4]))),xlim=c(min(c(prco[,5],ne[,5])),max(c(prco[,5],ne[,5]))));
points(ne[,4]~ne[,5],col="red")
plot(prco[,5]~prco[,1],ylim=c(min(c(prco[,5],ne[,5])),max(c(prco[,5],ne[,5]))),xlim=c(min(c(prco[,1],ne[,1])),max(c(prco[,1],ne[,1]))));
points(ne[,5]~ne[,1],col="red")
plot(prco[,5]~prco[,2],ylim=c(min(c(prco[,5],ne[,5])),max(c(prco[,5],ne[,5]))),xlim=c(min(c(prco[,2],ne[,2])),max(c(prco[,2],ne[,2]))));
points(ne[,5]~ne[,2],col="red")
plot(prco[,5]~prco[,3],ylim=c(min(c(prco[,5],ne[,5])),max(c(prco[,5],ne[,5]))),xlim=c(min(c(prco[,3],ne[,3])),max(c(prco[,3],ne[,3]))));
points(ne[,5]~ne[,3],col="red")
plot(prco[,5]~prco[,4],ylim=c(min(c(prco[,5],ne[,5])),max(c(prco[,5],ne[,5]))),xlim=c(min(c(prco[,4],ne[,4])),max(c(prco[,4],ne[,4]))));
points(ne[,5]~ne[,4],col="red")
plot(prco[,5]~prco[,5],type="n",xaxt="n",yaxt="n");
text(x=max((prco[,2]))-(max((prco[,2]))-min((prco[,2])))/2,y=max((prco[,5]))-(max((prco[,5]))-min((prco[,5])))/2,labels="PC5",cex=2,pos=4);

output<-list(prco.base=prco,prco.new=ne);
class(output)<-"pc.space";
return(output);
	}

