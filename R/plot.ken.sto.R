plot.ken.sto <-
function(x,...){
	if(x$"Calibration and validation set"=="FALSE"){
		x11(width=13,height=8);
if(x$"Number important PC"<=2){par(mfrow=c(2,2),mar=c(1,1,1,1))};
for(i in 3:5){
if(x$"Number important PC"==i){par(mfrow=c(i,i),mar=c(1,1,1,1))};
	}
if(x$"Number important PC">5){par(mfrow=c(5,5),mar=c(1,1,1,1))};
for(i in 1:if(x$"Number important PC"<=5){x$"Number important PC"}else{5}){
	for(j in 1:if(x$"Number important PC"<=5){x$"Number important PC"}else{5}){
		plot(x$"PC space important PC"[,i]~x$"PC space important PC"[,j],cex=0.3);
		points(x$"PC space important PC"[x$"Chosen row number",i]~x$"PC space important PC"[x$"Chosen row number",j],col="green",cex=0.3);
		}	
	}
		}
	
	if(x$"Calibration and validation set"=="TRUE"){
		x11(width=13,height=8);
if(x$"Number important PC"<=2){par(mfrow=c(2,2),mar=c(1,1,1,1))};
for(i in 3:5){
if(x$"Number important PC"==i){par(mfrow=c(i,i),mar=c(1,1,1,1))};
	}
if(x$"Number important PC">5){par(mfrow=c(5,5),mar=c(1,1,1,1))};
for(i in 1:if(x$"Number important PC"<=5){x$"Number important PC"}else{5}){
	for(j in 1:if(x$"Number important PC"<=5){x$"Number important PC"}else{5}){
		plot(x$"PC space important PC"[,i]~x$"PC space important PC"[,j],cex=0.3);
		points(x$"PC space important PC"[x$"Chosen validation row number",i]~x$"PC space important PC"[x$"Chosen validation row number",j],col="green",cex=0.3);
		}	
	}
		}	
	
	}

