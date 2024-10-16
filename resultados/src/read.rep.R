#' Read a ADMB report file
#' 
#' \code{read.rep} returns the sum of all the values present in its arguments.
#' 
#' This is a generic function: methods can be defined for it directly or via the
#' \code{\link{read.fit}} group generic. For this to work properly, the
#' arguments \code{...} should be unnamed, and dispatch is on the first
#' argument.
#' 
#' @param fn list of ADMB files after convergence
#' @return list of ADMB fits
#' 
#' @family aggregate functions
#' @seealso \code{\link{prod}} for products, \code{\link{cumsum}} for cumulative
#'   sums, and \code{\link{colSums}}/\code{\link{rowSums}} marginal sums over 
#'   high-dimensional arrays.
#'   
#' @export


read.rep <-
function(fn)
{
	# The following reads a report file
	# Then the 'A' object contains a list structure
	# with all the elemements in the report file.
	# In the REPORT_SECTION of the AMDB template use 
	# the following format to output objects:
	#  	report<<"object \n"<<object<<endl;
	#
	# The part in quotations becomes the list name.
	# Created By Steven Martell

	# A slight problem with read.table for ragged objects. read.table only
	# scans the first 5 rows to determine the number of colums.  FIX developed Sep 5 2013

	options(warn=-1)  #Suppress the NA message in the coercion to double
	
	ifile=scan(fn,what="character",flush=TRUE,blank.lines.skip=FALSE,quiet=TRUE)
	idx=sapply(as.double(ifile),is.na)
	vnam=ifile[idx] #list names
	# cat(vnam)
	nv=length(vnam) #number of objects
	A=list()
	ir=0
	for(i in 1:nv)
	{
		ir=match(vnam[i],ifile)
		if(i!=nv) irr=match(vnam[i+1],ifile) else irr=length(ifile)+1 #next row
		dum=NA
		if(irr-ir==2) dum=as.double(scan(fn,skip=ir,nlines=1,quiet=TRUE,what=""))
		if(irr-ir>2)
		{
			# ncols <- 0
			# irows <- ir:irr-1
			# for(j in irows)
			# {
			# 	tmp=as.double(scan(fn,skip=j,nlines=1,quiet=TRUE,what=""))
			# 	if(length(tmp)>ncols) ncols <- length(tmp)
			# 	#print(paste(1:ncols))
			# }
			# cname <- paste(1:ncols)
			# dum=as.matrix(read.table(fn,skip=ir,nrow=irr-ir-1,fill=TRUE,col.names=cname))
			# cat("\n ir ",ir," irr ",irr)
			dum=as.matrix(read.table(fn,skip=ir,nrows=irr-ir-1,fill=TRUE,row.names = NULL))
		} 
			

		if(is.numeric(dum))#Logical test to ensure dealing with numbers
		{
			A[[vnam[i]]]=dum
		}
	}
	options(warn=0)
	
	return(A)
}
