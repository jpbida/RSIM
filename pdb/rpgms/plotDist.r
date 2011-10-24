dat3<-read.table(file="ssCCC.dat",sep=",",header=FALSE)
dat4<-read.table(file="ssAAA.dat",sep=",",header=FALSE)

names(dat3)<-c("id","pt","x","y","z")
names(dat4)<-c("id","pt","x","y","z")

dat1<-dat3[dat3$pt=="c1" | dat3$pt=="c3",]
dat2<-dat4[dat4$pt=="c1" | dat4$pt=="c3",]
cls1<-rep(0,length(dat1$pt))
cls1[dat1$pt=="c1"]<-2
cls1[dat1$pt=="c3"]<-4

cls2<-rep(0,length(dat2$pt))
cls2[dat2$pt=="c1"]<-3
cls2[dat2$pt=="c3"]<-5
library("rgl")
open3d()
k<-dat4$pt=="n1" | dat4$pt=="n3"
dat5<-dat4[k,]
#plot3d(dat2$x,dat2$y,dat2$z,col=cls2,size=2,xlim=c(-30,30),ylim=c(-30,30),zlim=c(-30,30))
plot3d(dat2$x,dat2$y,dat2$z,col=cls2,size=2,axes=FALSE,ylab="",xlab="",zlab="",box=FALSE)
#for(j in 1:length(dat2$x)){
#	if(!is.na(dat2$x[j])){
#	segments3d(c(dat2$x[j],dat5$x[j]),c(dat2$y[j],dat5$y[j]),c(dat2$z[j],dat5$z[j]))
#	}}

#points3d(dat1$x,dat1$y,dat1$z,col=cls1,size=2)
