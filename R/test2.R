

dat1<-iris[,3]



estep(dat1,p1)



p1<-matrix(c(0.3,0.7,0,1.5,30,5,1,1),2,4)

for(i in 1:100){
  e1<-estep(dat1,p1)
  p2<-mstep(dat1, e1[1:2,], e1[3:4,], p1)
  print(p2)
  p1<-p2
}
