library(EMMIXgene)
filename_b='/home/andrew/projects/EMMIX-GENE/dat/golub_norm.dat'
filename='/home/andrew/projects/EMMIX-GENE/dat/alon_norm.dat'

load("/home/andrew/Downloads/chang_lac03.RData")



test1<-select_genes(filename)
test2<-cluster_genes(test1, g)
test3a<-cluster_tissues(test1, test2)




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


test5<-cluster_tissues_top_genes((test1), 50)
heat_maps(test1$all_genes[test5$top_genes,]) 

qplot(test3$mfa_fit$U[,,1], test3$mfa_fit$U[,,2])

