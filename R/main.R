
# filname='/home/andrew/projects/EMMIX-GENE/dat/golub_norm.dat'
#/home/andrew/projects/EMMIX-GENE/dat/alon_norm.dat

#' Selects Genes from data supplied.
#'
#' @param filename Name of file containing gene data.
#' @param random_starts
#' @param ll_thresh
#' @param min_clust_size
#' @return matrix of genes
#' @examples
select_genes<-function(filename, random_starts=4, max_it = 400, ll_thresh = 8, min_clust_size = 8){
    data<-as.matrix(read.delim(filename, sep=" ", header=F))
    a<-emmix_gene(data)
    genes <- data[a$g>1,]
    return(genes)
}

#' Clusters genes
#'
#' @param genes matrix of genes
#' @param g number of clusters
#' @return mclust object
#' @examples
cluster_genes<-function(genes, g=NULL){
  clust_genes<-Mclust(genes, G=g, modelNames = "VII")
  return(clust_genes)
}

#' Clusters tissues
#'
#' @param genes matrix of genes
#' @param g number of clusters
#' @return ?
#' @examples
cluster_tissues<-function(genes, method='t'){
  if(method=='t'){
    
  }
  if(method=='mfa')
  return()
}

#' Heat maps
#'
#' @param clust_genes mclust object
#' @return ?
#' @examples
heat_maps<-function(clust_genes){

}








