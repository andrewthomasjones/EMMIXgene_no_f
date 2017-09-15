library(EMMIXgene)
filename_b='/home/andrew/projects/EMMIX-GENE/dat/golub_norm.dat'
filename='/home/andrew/projects/EMMIX-GENE/dat/alon_norm.dat'

g=10

test1<-select_genes(filename)
test2<-cluster_genes(test1, g)
test3a<-cluster_tissues(test1, test2)

test3<-cluster_tissues(test1, test2, 'mfa')


for(k in 1:g){
  test4<-heat_maps(test1$genes[test2$classification==k,], test3[k,]) 
  plot(test4)
}

