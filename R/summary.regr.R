summary.regr <-
function(object,...){

output<-list("Model name"=object$model.name,"Spectral transformation"=object$spectral.transformation,"Constituents"=object$constituents,"Constituent transformation"=object$constituents.transformation,"Regression method"=object$method,"Calibration set statistics"=object$cal.statistics,"Validation set statistics"=object$val.statistics);
print(output);
	
	
	}
