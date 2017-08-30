#' Selects Genes from data supplied.
#'
#' @param filename Name of file containing gene data.
#' @param g A number.
#' @param k A number.
#' @param random_starts
#' @param ll_thresh
#' @param min_clust_size
#' @return not sure yet
#' @examples
select_genes_dat<-function(filename, g=1, k =1, random_starts=5, ll_thresh = 0.001, min_clust_size = 3){
    data<-read.delim(filename, sep=" ", header=F)
    row<-dim(x)[1]
    col<-dim(x)[2]
    emmix_t(dat, params, 1, 100)
    res<-0#select_genes(data, row, col, g, k, random_starts,ll_thresh, min_clust_size)
}

