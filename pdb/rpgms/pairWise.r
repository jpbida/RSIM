##### Finds all pair-wise distances between centers ####

### Calculate distance from center ###
dC<-read.table(file="ssAAA.dat", sep=",",header=FALSE)
dC<-dC[!is.na(dC[,4]),]
names(dC)<-c("id","pt","x","y","z")

t<-table(dC$id)
t<-names(t)
t<-as.numeric(t)
dist1<-NULL
dist2<-NULL
w<-1


#### Calculate angles phi and theta ####
theta1<-NULL
ntheta1<-NULL
theta2<-NULL
ntheta2<-NULL
phi1<-NULL
nphi1<-NULL
phi2<-NULL
nphi2<-NULL
dist1<-NULL
ndist1<-NULL
dist2<-NULL
ndist2<-NULL
dist13<-NULL
for(id in t){
d<-dC[dC$id==id,]
pt1<-d[d$pt=="c2",]
pt2<-d[d$pt=="c3",]
pt3<-d[d$pt=="n2",]
pt4<-d[d$pt=="n3",]
pt5<-d[d$pt=="c1",]
dist13[w]<-sqrt((pt2$x-pt5$x)^2+(pt2$z-pt5$z)^2+(pt2$z-pt5$z)^2)
#### Translate pt1 to the origin ###
pt4$x<-pt4$x-pt2$x
pt4$y<-pt4$y-pt2$y
pt4$z<-pt4$z-pt2$z

nphi1[w] = atan2(pt4$y[1],pt4$x[1])
ntheta1[w] = acos(pt4$z[1]/sqrt(pt4$x[1]^2+pt4$y[1]^2+pt4$z[1]^2))
ndist1[w]<-sqrt((pt4$x[1])^2+(pt4$y[1])^2+(pt4$z[1])^2)
phi1[w] = atan2(pt2$y[1],pt2$x[1])
theta1[w] = acos(pt2$z[1]/sqrt(pt2$x[1]^2+pt2$y[1]^2+pt2$z[1]^2))
dist1[w]<-sqrt((pt1$x[1]-pt2$x[1])^2+(pt1$y[1]-pt2$y[1])^2+(pt1$z[1]-pt2$z[1])^2)
pt1<-d[d$pt=="c2",]
pt2<-d[d$pt=="c1",]
pt3<-d[d$pt=="n2",]
pt4<-d[d$pt=="n1",]
pt4$x<-pt4$x-pt2$x
pt4$y<-pt4$y-pt2$y
pt4$z<-pt4$z-pt2$z
nphi2[w] = atan2(pt4$y[1],pt4$x[1])
ntheta2[w] = acos(pt4$z[1]/sqrt(pt4$x[1]^2+pt4$y[1]^2+pt4$z[1]^2))
ndist2[w]<-sqrt((pt4$x[1])^2+(pt4$y[1])^2+(pt4$z[1])^2)
phi2[w] = atan2(pt2$y[1],pt2$x[1])
theta2[w] = acos(pt2$z[1]/sqrt(pt2$x[1]^2+pt2$y[1]^2+pt2$z[1]^2))
dist2[w]<-sqrt((pt1$x[1]-pt2$x[1])^2+(pt1$y[1]-pt2$y[1])^2+(pt1$z[1]-pt2$z[1])^2)
w<-w+1
}


pairs(cbind(phi1,phi2,nphi1,nphi2),ylim=c(-pi,pi),xlim=c(-pi,pi))
library("rgl")
open3d()
k<-dC$pt=="c1"
dg<-dC[k,]
k<-dC$pt=="c3"
dp<-dC[k,]

plot3d(dg$x,dg$y,dg$z,size=3,ylim=c(-10,10),xlim=c(-10,10),zlim=c(-10,10))
points3d(0,0,0,col=4,size=6)
lines3d(c(0,0),c(0,0),c(0,10),col=4)
lines3d(c(0,0),c(0,10),c(0,0))
lines3d(c(0,10),c(0,0),c(0,0))
points3d(dp$x,dp$y,dp$z,size=3,col=5)


#### Based on distributions create a function g2(theta,phi,r) #####
### Fit the distribution to an ellipse ####

ell<-function(theta,phi,a,b){
#### Given theta and phi it outputs R ####
#r<-a*(1-e^2)/(1+e*cos(theta))
r<-a*b/(sqrt(a^2*sin(theta)^2+b^2*cos(theta)^2))
r

}
ell3d<-function(theta,phi,a,b,c,x1,y1,z1){
x<-a*cos(theta)*sin(phi)
y<-b*sin(theta)*sin(phi)
z<-c*cos(phi)
dat<-NULL
dat$x<-x+x1
dat$y<-y+y1
dat$z<-z+z1
dat
}

optFun<-function(p){
a<-p[1]
b<-p[2]
c<-p[3]
x1<-p[4]
y1<-p[5]
z1<-p[6]
ell3d(


}

out<-NULL
j<-1
theta<-seq(0,(2*pi),length=100)
phi<-seq(0,(pi),length=100)

for(t in theta){
	for(p in phi){
out$phi[j]<-p
	out$theta[j]<-t
	j<-j+1	
	}
}


