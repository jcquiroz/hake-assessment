#' Create a palette color
#' 
#' @param colpal A vector
#' @param a A number
#' @return A palette color of length \code{colpal} and alpha \code{a}.
#' @examples
#' colorea(1,0.5)
#' 
#' @export

colorea <- function(colpal,a)
{
  colrgb<-col2rgb(colpal)/255
  coo <- rgb(t(colrgb),alpha=a)
  #coo <- colpal + a
  return(coo)
}
