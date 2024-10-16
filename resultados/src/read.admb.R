#' Read a set of ADMB outcomes
#' 
#' \code{read.admb} returns the sum of all the values present in its arguments.
#' 
#' This is a generic function: methods can be defined for it directly or via the
#' \code{\link{read.psv}} group generic. For this to work properly, the
#' arguments \code{...} should be unnamed, and dispatch is on the first
#' argument.
#' 
#' @param ifile list of ADMB files after convergence
#' @return list of ADMB objects
#' 
#' @family aggregate functions
#' @seealso \code{\link{prod}} for products, \code{\link{cumsum}} for cumulative
#'   sums, and \code{\link{colSums}}/\code{\link{rowSums}} marginal sums over 
#'   high-dimensional arrays.
#'   
#' @export


read.admb <-
function(ifile)
{	
	ret   <- read.fit(ifile)  # read *.par and *.cor (if last exit)
	fn    <- paste(ifile,'.rep', sep='')
	A     <- read.rep(fn)
	A$fit <- ret
	pfn   <- paste(ifile,'.psv',sep='')
	if(file.exists(pfn))
		A$post.samp <- read.psv(pfn)
	return(A)
}
