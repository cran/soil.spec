summary.prco <-
function(object,...){
	cat("Summary prcomp function:\n")
	print(summary(object$prcomp)[[6]][2:3,1:10])
	cat("\n")
	cat("Important principal components:\n")
	print(object$i.pc)
	}

