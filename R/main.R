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


select_genes_dat<-function(filename, g=1, k =1, random_starts=5, ll_thresh = 0.001, min_clust_size = 3){
    data<-read.delim(filename, sep=" ", header=F)
    row<-dim(x)[1]
    col<-dim(x)[2]
    res<-select_genes(data, row, col, g, k, random_starts,ll_thresh, min_clust_size)
}

