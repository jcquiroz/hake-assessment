#' Read a ADMB fit file
#' 
#' \code{read.fit} returns the sum of all the values present in its arguments.
#' 
#' This is a generic function: methods can be defined for it directly or via the
#' \code{\link{read.psv}} group generic. For this to work properly, the
#' arguments \code{...} should be unnamed, and dispatch is on the first
#' argument.
#' 
#' @param ifile list of ADMB files after convergence
#' @return list of ADMB fits
#' 
#' @family aggregate functions
#' @seealso \code{\link{prod}} for products, \code{\link{cumsum}} for cumulative
#'   sums, and \code{\link{colSums}}/\code{\link{rowSums}} marginal sums over 
#'   high-dimensional arrays.


read.fit <-
function(ifile)
{
	# __Example:             
	#	file <-("~/admb/simple")
	#	A <- reptoRlist(file)
	#	Note there is no extension on the file name.
	
	## The following is a contribution from:
	## Anders Nielsen that reads the par & cor files.
	ret<-list() 
	parfile<-as.numeric(scan(paste(ifile,'.par', sep=''),   
	 what='', n=16, quiet=TRUE)[c(6,11,16)]) 
	ret$nopar<-as.integer(parfile[1]) 
	ret$nlogl<-parfile[2] 
	ret$maxgrad<-parfile[3] 
	file<-paste(ifile,'.cor', sep='')
	if(file.exists(file))
	{
		lin<-readLines(file) 
		ret$npar<-length(lin)-2 
		ret$logDetHess<-as.numeric(strsplit(lin[1], '=')[[1]][2]) 
		sublin<-lapply(strsplit(lin[1:ret$npar+2], ' '),function(x)x[x!='']) 
		ret$names<-unlist(lapply(sublin,function(x)x[2])) 
		ret$est<-as.numeric(unlist(lapply(sublin,function(x)x[3]))) 
		ret$std<-as.numeric(unlist(lapply(sublin,function(x)x[4]))) 
		ret$cor<-matrix(NA, ret$npar, ret$npar) 
		corvec<-unlist(sapply(1:length(sublin), function(i)sublin[[i]][5:(4+i)])) 
		ret$cor[upper.tri(ret$cor, diag=TRUE)]<-as.numeric(corvec) 
		ret$cor[lower.tri(ret$cor)] <- t(ret$cor)[lower.tri(ret$cor)] 
		ret$cov<-ret$cor*(ret$std%o%ret$std)		
	}
	return(ret)
}
