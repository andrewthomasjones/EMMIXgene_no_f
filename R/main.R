#' Add together two numbers.
#'
#' @param x A number.
#' @param y A number.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' add(1, 1)
#' add(10, 1)
add <- function(x, y) {
  x + y
}

#E-step

tau.estep.wire<-function(dat,pro,mu,sigma,n,m,g)
{


  print(paste0("pro :", pro),na.print="NA")
  print(paste0("mu :", mu),na.print="NA")
  print(paste0("sigma :", sigma),na.print="NA")
  print(paste0("n :", n),na.print="NA")
  print(paste0("m :", m),na.print="NA")
  print(paste0("g :", g),na.print="NA")




  obj <- .Fortran("estepmvn",PACKAGE="EMMIXcontrasts2",
                  as.double(dat),as.integer(n),as.integer(m),as.integer(g),
                  as.double(pro),as.double(mu),as.double(sigma),
                  mtauk=double(n*g),double(g),loglik=double(1),
                  error = integer(1))[8:11]
  if(obj$error) stop("error")
  tau <- array(obj$mtauk,c(n,g))
  list(tau=tau,loglik=obj$loglik,pro=colMeans(tau))
}
