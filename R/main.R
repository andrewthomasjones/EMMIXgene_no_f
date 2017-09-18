
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
select_genes<-function(dat, filename, random_starts=4, max_it = 400, ll_thresh = 8, min_clust_size = 8){
  
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
    
    
    if(any(!complete.cases((data)))| any(!complete.cases(t(data)))){
      warning("Incomplete cases removed removed.")
    }
  
    #remove missing data
    data<-data[complete.cases((data)),]
    data<-data[,complete.cases(t(data))]
   
    
    #actual method
    
    #do emmix_gene C++ routine
    a<-emmix_gene(data)
    a$selected<-a$g>1
    
    #add selected genes to result
    a$genes <- data[a$g>1,]
    a$all_genes <- data
    
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
#' @param gen emmix-gene object
#' @param clusters mclust object
#' @param method Method for seperating tissue classes. Can be either 't' for a univariate mixture of t-distributions on gene cluster means, or 'mfa' for a mixture of factor analysers. 
#' @param k number of factors if using mfa
#' @return a clustering for each sample (columns) by each group(rows)
#' @examples
#' 
#' @export
cluster_tissues<-function(gen, clusters, method='t', k=6){
  
  g<-clusters$G
  p<-ncol(clusters$data)
  clustering<-array(0,c(g,p))
  clustering2<-array(0,c(g,p))
  
  if(method=='t'){
    for(i in 1:g){
        group_means <- colMeans(gen$genes[clusters$classification==i,])
        t_fit<-emmix_t(group_means, 2)
        clustering[i,]<-as.numeric(xor(t_fit$Clusters, (t_fit$mu[1]>t_fit$mu[2])))
     }
    
  }
  
  if(method=='mfa'){
    for(i in 1:g){
      group <- as.matrix((gen$genes[clusters$classification==i,]))
      #actually mixture of common factor analysers. consider fixing.
      #print(dim(group))
      
        mfa_fit<-mcfa(t(group), 2, k)
        clustering[i,]<- as.numeric(xor(predict_mcfa(mfa_fit, t(group))-1, (diff(mfa_fit$xi[1,])>0))) 
        
    
      
    }
    
    
  }
  
  return(clustering)
}




#' Clusters tissues
#'
#' @param gen emmix-gene object
#' @param n_top number of top genes (as ranked by liklihood) to be selected
#' @param method Method for seperating tissue classes. Can be either 't' for a univariate mixture of t-distributions on gene cluster means, or 'mfa' for a mixture of factor analysers. 
#' @param k number of factors if using mfa
#' @return a clustering for each sample (columns) by each group(rows)
#' @examples
#' 
#' @export
top_genes_cluster_tissues<-function(gen, n_top=100, method='mfa', k=2){
  

  p<-ncol(gen$genes)
  clustering<-array(0,p)
  
  top_genes<-order(test1$stat,decreasing = FALSE)[1:n_top]
  
  if(method=='t'){
    
      group_means <- colMeans(gen$all_genes[top_genes,])
      t_fit<-emmix_t(group_means, 2)
      clustering<-as.numeric(xor(t_fit$Clusters, (t_fit$mu[1]>t_fit$mu[2])))
    
    
  }
  
  if(method=='mfa'){
    
      group <- as.matrix((gen$all_genes[top_genes,]))
      mfa_fit<-mcfa(t(group), 2, k)
      clustering<- as.numeric(xor(predict_mcfa(mfa_fit, t(group))-1, (diff(mfa_fit$xi[1,])>0))) 
      
      
      
    }
    
    
  
  
  return(list(clustering=clustering, top_genes=top_genes, mfa_fit=mfa_fit))
}










#' Heat maps
#'
#' @param clust_genes matrix of genes
#' @return ggplot2 heat map
#' @examples
#' 
#' @export
heat_maps<-function(clust_genes, clustering=NULL){
  colnames(clust_genes) <- NULL
  if(!is.null(clustering)){
    clust_genes<-clust_genes[,order(clustering)]
  }
  
  df_heatmap<-melt(clust_genes)
  names(df_heatmap)<-c("genes", "samples",  "expression_level")
  df_heatmap$genes<-factor(df_heatmap$genes)
  df_heatmap$samples<-factor(df_heatmap$samples)
  
  plot<-ggplot(df_heatmap, aes(samples,genes )) + geom_tile(aes(fill = expression_level),  color = "white") +
    scale_fill_distiller(palette = "RdYlGn") + #, limits=c(-3,3)) + 
    ylab("Genes") +
    xlab("Samples") +
    theme(legend.title = element_text(size = 10),
          legend.text = element_text(size = 12),
          plot.title = element_text(size=16),
          axis.title=element_text(size=14,face="bold"),
          axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y=element_blank(),axis.ticks.x=element_blank()) +
    labs(fill = "Expression level")
  
  return(plot)
}

#' Plot single gene
#'
#' @param dat full matrix of data
#' @param gene_id row number
#' @param random_starts number of random starts for the fit
#' @return ggplot2 histogram with fit
#' @examples
#' 
#' @export
plot_single_gene<-function(dat, gene_id, random_starts=50){ 
  
  
  
  df<-data.frame(x=dat[gene_id,])
  n<-length(df$x)/4
  breaks<-seq(-4,4, length.out=n)
  plot<-ggplot(df, aes(x=x)) + geom_histogram(aes(y=..density..), breaks=breaks, alpha=.5)+theme_bw()
  
  
  df2<-data.frame(x=seq(-4, 4, length.out = 1000))
  res<-each_gene(dat[gene_id,],random_starts)
  for(i in 1:res$components){
    for(j in 1:nrow(df2)){
      df2[[paste0('y',i)]][j]<-res$pi[i]*t_dist(df2$x[j], res$mu[i], res$sigma[i], res$nu[i])
    }
  }
  
  
  plot<-plot+geom_line(data=df2, aes(x=x, y=y1))
  if(res$components>1){
    plot<-plot+geom_line(data=df2, aes(x=x, y=y2))
  }
  
  if(res$components>2){
    plot<-plot+geom_line(data=df2, aes(x=x, y=y3))
  }
  
  return(plot)
}





