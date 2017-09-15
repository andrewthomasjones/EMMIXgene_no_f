
# filename='/home/andrew/projects/EMMIX-GENE/dat/golub_norm.dat'
# filename='/home/andrew/projects/EMMIX-GENE/dat/alon_norm.dat'

#' Selects Genes from data supplied.
#'
#' @param filename Name of file containing gene data. Can be either .csv or space separated .dat. Rows are genes and columns are samples. Must supply one of filename and dat.
#' @param dat A matrix or dataframe containing gene expression data. Rows are genes and columns are samples. Must supply one of filename and dat.
#' @param random_starts The number of random starts used per gene when fitting t-distributions.
#' @param ll_thresh The difference in log-likihood used as a threshold for selecting between g=1 and g=2.
#' @param min_clust_size The minimum number of observations per cluster. 
#' @return emmix-gene object
#' @examples
#' 
#' 
#' @export
select_genes<-function(filename, dat, random_starts=4, max_it = 400, ll_thresh = 8, min_clust_size = 8){
  
    #housekeeping   
  
    if(missing(filename)&missing(dat)){
      stop("Must supply one of filename and dat.")
    }
    
    if(!missing(filename)&!missing(dat)){
      stop("Must supply ONLY one of filename and dat.")
    }  
  
    if(!missing(filename)&missing(dat)){
      if(file_ext(filename)=="dat"){data<-as.matrix(read.delim(filename, sep=" ", header=F))}
      if(file_ext(filename)=="csv"){data<-as.matrix(read.csv(filename, header=F))}
    }
  
    if(missing(filename)&!missing(dat)){
      data<-as.matrix(dat)
    }
    
    
    if(any(complete.cases((data)))|any(complete.cases(t(data)))){
      warning("Incomplete cases removed removed.")
    }
  
    #remove missing data
    data<-data[,complete.cases(t(data))]
    data<-data[complete.cases((data)),]
    
    #actual method
    
    #do emmix_gene C++ routine
    a<-emmix_gene(data)
    
    #add selected genes to result
    a$genes <- data[a$g>1,]
    
    result<-structure(a, class="emmix-gene")
    #return selected genes
    return(result)
}

#' Clusters genes
#'
#' @param gen an emmix-gene object produced by select genes.
#' @param g number of gene clusters . If not specified will be selected automatically on the basis of BIC.
#' @return mclust object
#' @examples
#' 
#' @export
cluster_genes<-function(gen, g=NULL){
  clust_genes<-Mclust((gen$genes), G=g, modelNames = "VII")
  g<-clust_genes$G
  
  ll_rank_stat<-array(0,g)
  
  for(i in 1:g){
    ll_rank_stat[i]<-each_gene(colMeans(clust_genes$data[clust_genes$classification==i,]))$Ratio
  }
  
  #reorder clusters as ranked by mean of ll statistic from emmix-gene fit
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
#' 
#' @export
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
#' 
#' @export
heat_maps<-function(clust_genes){

}





