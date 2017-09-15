library(EMMIXgene)
filename='/home/andrew/projects/EMMIX-GENE/dat/golub_norm.dat'
filename='/home/andrew/projects/EMMIX-GENE/dat/alon_norm.dat'

test1<-select_genes(filename)
test2<-cluster_genes(test1)
test3<-cluster_tissues(test2)
test4<-heat_maps(test2)


