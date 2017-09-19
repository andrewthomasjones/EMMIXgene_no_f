library(EMMIXgene)
data(alon_data)
data(golub_data)
#5.1
alon_sel<-select_genes(alon_data)
alon_sel2<-select_genes(dat=alon_data, random_starts=8, max_it = 200, ll_thresh = 8, min_clust_size = 8, tol = 0.0001, start_method ="both", three = T)
n_alon<-sum(alon_sel2$selected)

#5.1.1
alon_top<-top_genes_cluster_tissues(alon_sel, n_alon, k=6)
alon_top_50<-top_genes_cluster_tissues(alon_sel, 50, k=6)

heat_maps(alon_sel$all_genes[alon_top$top_genes,], alon_top$clustering)
heat_maps(alon_sel$all_genes[alon_top_50$top_genes,], alon_top_50$clustering)


#5.1.2
alon_clus_20<-cluster_genes(alon_sel2, 10)
alon_clus_20_t<-cluster_tissues(alon_sel, alon_clus_20, method='t', k=4)

for(j in 1:10){
  plot<-heat_maps(alon_sel$genes[alon_clus_20$classification==j,], alon_clus_20_t[j,]) 
  print(plot)
}

golub_sel<-select_genes(golub_data)
golub_sel2<-select_genes(dat=golub_data, random_starts=8, max_it = 100, ll_thresh = 8, min_clust_size = 8, tol = 0.00001, start_method ="kmeans", three = T)

#5.2
golub_clus_g<-cluster_genes(golub_sel)
n_golub<-sum(golub_sel2$selected)

golub_clus_40<-cluster_genes(golub_sel2, 20)
golub_clus_40_t<-cluster_tissues(golub_sel2, golub_clus_40, method='t', k=8)

golub_top_50<-top_genes_cluster_tissues(golub_sel2, 50, k=6)


heat_maps(golub_sel$all_genes[golub_top_50$top_genes,], golub_top_50$clustering)

for(j in 1:10){
  plot<-heat_maps(golub_sel2$genes[golub_clus_40$classification==j,], golub_clus_40_t[j,]) 
  print(plot)
}


golub_top_50<-top_genes_cluster_tissues(golub_sel2, 50, k=6, g=2)












for(k in 1:g){
  test4<-heat_maps(test1$genes[test2$classification==k,], test3[k,]) 
  plot(test4)
}

test1<-select_genes(alon_data)
test2<-select_genes(alon_data)
test3<-select_genes(alon_data)
test4<-select_genes(alon_data)
test5<-select_genes(alon_data)

load('/home/andrew/Downloads/west_pnas01.RData')

load("/home/andrew/Downloads/chang_lac03.RData")
m<-scale(t(X), center = TRUE, scale = TRUE)
test1<-select_genes(m)
g=10
test2<-cluster_genes(test1, g)
test3<-cluster_tissues(test1, test2, 'mfa')

for(k in 1:g){
  test4<-heat_maps(test1$genes[test2$classification==k,], test3[k,]) 
  plot(test4)
}




test2<-cluster_genes(test1, g)

top50<-(order(test1$stat, decreasing = T))[1:50]


tumour_info<-read.csv('/home/andrew/Documents/test.csv')

heat_maps(m[top50+4000,])
heat_maps(t(X)[top50,order(tumour_info$Lymph.node.status)])
heat_maps(t(X)[top50,])

t(X)[top50,order(tumour_info$Lymph.node.status)]


################################################################################
library(cancerdata)
data(VEER1)
data<-t(as.data.frame(VEER1))
data<-as.matrix(data)
m <- mapply(data, FUN=as.numeric)
m<-matrix(m, dim(data))


test1<-select_genes(dat=m)


g=10
test2<-cluster_genes(test1,g)

test3a<-cluster_tissues(test1, test2)
test3<-cluster_tissues(test1, test2, 'mfa', 2)

for(k in 1:g){
  test4<-heat_maps(test1$genes[test2$classification==k,],test3[k,]) 
  plot(test4)
}



heat_maps(test1$all_genes[test5$top_genes,]) 

rbind(test5$mfa_fit$U[,,1],test5$mfa_fit$U[,,2])

group <- as.matrix((gen$all_genes[top_genes,]))
G1<-factor(predict_mcfa(mfa_fit, t(group)))
qplot(test5$mfa_fit$U[,1,1], test5$mfa_fit$U[,2,1], colour=G1)

