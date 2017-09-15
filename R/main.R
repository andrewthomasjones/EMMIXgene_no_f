
# filename='/home/andrew/projects/EMMIX-GENE/dat/golub_norm.dat'
# filename='/home/andrew/projects/EMMIX-GENE/dat/alon_norm.dat'

#' Selects Genes from data supplied.
#'
#' @param filename Name of file containing gene data.
#' @param random_starts
#' @param ll_thresh
#' @param min_clust_size
#' @return matrix of genes
#' @examples
select_genes<-function(filename, random_starts=4, max_it = 400, ll_thresh = 8, min_clust_size = 8){
    data<-as.matrix(read.delim(filename, sep=" ", header=F))[,1:72]
    a<-emmix_gene(data)
    genes <- data[a$g>1,]
    return(a)
}

#' Clusters genes
#'
#' @param genes matrix of genes
#' @param g number of clusters
#' @return mclust object
#' @examples
cluster_genes<-function(genes, g=NULL){
  clust_genes<-Mclust((genes), G=g, modelNames = "VII")
  g<-clust_genes$G
  ll_rank_stat<-array(0,g)
  
  for(i in 1:g){
    ll_rank_stat[i]<-each_gene(colMeans(clust_genes$data[clust_genes$classification==i,]))$Ratio
  }
  
  tmp<-factor(clust_genes$classification)
  levels(tmp) <- as.character(order(unlist(ll_rank_stat)))
  clust_genes$classification<-as.numeric(as.character(tmp))
  
  return(clust_genes)
}

#' Clusters tissues
#'
#' @param genes matrix of genes
#' @param g number of clusters
#' @return ?
#' @examples
cluster_tissues<-function(genes, method='t'){
  #clust_genes$data[clust_genes$classification==1,]
  g<-clust_genes$G
  p<-ncol(clust_genes$data)
  clustering<-array(0,c(g,p))
  
  if(method=='t'){
    for(i in 1:g){
        group_means <- colMeans(data[which(a$g>1)[clust_genes$classification==i],])
        t_fit<-emmix_t(group_means, 2, 100)
        clustering[i,]<-as.numeric(xor(t_fit$Clusters, (t_fit$mu[1]>t_fit$mu[2])))
     }
    
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





