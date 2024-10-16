#' Read a ADMB psv file
#' 
#' \code{read.psv} returns the sum of all the values present in its arguments.
#' 
#' This is a generic function: methods can be defined for it directly or via the
#' \code{\link{read.fit}} group generic. For this to work properly, the
#' arguments \code{...} should be unnamed, and dispatch is on the first
#' argument.
#' 
#' @param fn list of ADMB files after convergence
#' @param nsamples total of MCMC samples
#' @return list of ADMB fits
#' 
#' @family aggregate functions
#' @seealso \code{\link{prod}} for products, \code{\link{cumsum}} for cumulative
#'   sums, and \code{\link{colSums}}/\code{\link{rowSums}} marginal sums over 
#'   high-dimensional arrays.


read.psv <-
function(fn, nsamples=10000)
{
	#This function reads the binary output from ADMB
	#-mcsave command line option.
	#fn = paste(ifile,'.psv',sep='')
	filen <- file(fn, "rb")
	nopar <- readBin(filen, what = integer(), n = 1)
	mcmc <- readBin(filen, what = numeric(), n = nopar * nsamples)
	mcmc <- matrix(mcmc, byrow = TRUE, ncol = nopar)
	close(filen)
	return(mcmc)
}
