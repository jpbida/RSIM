#### Features of a network #############
# Degree of node = number of edges coming from the node
# Mean Degree = Expectation of degree
# Variance of node degree
# Degree correlation - Correlation between degree node relative to the neighbors
# 		- Pick to nodes with N edges, make plot of other nodes
# Clustering Coefficient(transitivity)
# 		- Number of triangles/elbows 
#		- Pick three connected nodes A->B->C what is the probability C->A?
# Critical Components 
# 		Plot components with respect to size of network
# Mean shortest path between any two nodes
# Diameter = Largest path between two nodes
# Spacially embbed random network
# 	Protein Interaction Network 
#		Probabilty of interacting based on evolutionary time?  
#               Transcription Factor & Gene regulation
##### Pixelizes the fragment class #####
file1="fRRR.vect"
file2="fRRY.vect"
file3="fRYR.vect"
file4="fRYY.vect"
file5="fYRR.vect"
file6="fYRY.vect"
file7="fYYR.vect"
file8="fYYY.vect"

file1r="rRRR.vect"
file2r="rRRY.vect"
file3r="rRYR.vect"
file4r="rRYY.vect"
file5r="rYRR.vect"
file6r="rYRY.vect"
file7r="rYYR.vect"
file8r="rYYY.vect"
### id,x,y,z,r,gid ###
dat<-read.table(file=file1,sep="\t")
points<-NULL
points$x<-dat[,2]
points$y<-dat[,3]
points$z<-dat[,4]
points$id<-dat[,1]
points<-as.data.frame(points)
points<-points[points$x<30 & points$y < 30 & points$z<30,]
m<-range(c(points$x,points$y,points$z))
r<-diff(m)
x<-m[1]
y<-m[1]
z<-m[1]
names(points)=c("x","y","z","id")


grid<-function(file,num){
dat<-read.table(file=file,sep="\t")
points<-NULL
points$x<-dat[,2]
points$y<-dat[,3]
points$z<-dat[,4]
points$id<-c(0:(length(dat[,4])-1))
print(range(points$id))
points<-as.data.frame(points)
names(points)=c("x","y","z","id")
s<-kmeanSplit(points,num)
#### Get the number of elements in each box #####
write.table(s,file=paste("c",file,sep=""),sep="\t",row.names=FALSE,col.names=FALSE)
}


kmeanSplit<-function(dat,num){
#### use kmeans clustering to identify centers ####
#### Calculate all centers pairwise differences ###
#### Assign radius as minr/2 for each center    ###
#### Determine how many points are ignored      ###
k<-kmeans(cbind(dat$x,dat$y,dat$z),centers=num)
#library("rgl")

rads<-NULL
clus<-k$cluster
dat$gid=clus
ci1<-table(clus)
ci1<-as.numeric(names(ci1))
cen<-k$centers
for(id in c(1:length(ci1))){
##### Loop through each cluster and determine all pairwise distances between members and center ###
tmp<-dat[dat$gid==ci1[id],]
dmax=0
for(j in 1:dim(tmp)[1]){
dist<-(cen[id,1]-tmp[j,1])^2
dist<-dist+(cen[id,2]-tmp[j,2])^2
dist<-dist+(cen[id,3]-tmp[j,3])^2
if(dist>dmax){
dmax=dist
}
}
rads[id]=dmax
}
##### Loop through and determine minimal distance between cluster centers ### 




rads<-cbind(c(1:length(cen[,1])),cen[,1],cen[,2],cen[,3],rads)
rads<-as.data.frame(rads)
names(rads)=c("gid","x1","y1","z1","r")
g<-rads
rgl.spheres(g$x1,g$y1,g$z1,radius=(sqrt(g$r)*1),color="blue",front="line",back="cull",lit=FALSE)
all<-merge(dat,rads,by="gid")
all<-cbind(all$id,all$x1,all$y1,all$z1,sqrt(all$r),all$gid)
all<-as.data.frame(all)
names(all)<-c("id","x","y","z","r","gid")
all
}




splitBox<-function(r,x,y,z,dat,gid,nmax){
### x,y,z is the lower lefti                                        #
### r is the length of an edge                                      #
### dat has columns X1,Y1,Z1,ID are the coordinates of the points   #
### Split into Four Cubes return a dataset with the format          #
### r,x,y,z,id,group                                                #
### where r is new edge length x,y,z is lower left                  # 
print(gid)
c1<-dat$x > (x+r/2) & dat$x <= (x+r) & dat$y > y & dat$y <= (y+r/2) & dat$z > z & dat$z <= (z+r/2)
c2<-dat$x > (x+r/2) & dat$x <= (x+r) & dat$y > y & dat$y <= (y+r/2) & dat$z > (z+r/2) & dat$z <= (z+r)
c3<-dat$x > (x+r/2) & dat$x <= (x+r) & dat$y > (y+r/2) & dat$y <= (y+r) & dat$z > (z+r/2) & dat$z <= (z+r)
c4<-dat$x > (x+r/2) & dat$x <= (x+r) & dat$y > (y+r/2) & dat$y <= (y+r) & dat$z > z & dat$z <= (z+r/2)

c5<-dat$x > (x) & dat$x <= (x+r/2) & dat$y > y & dat$y <= (y+r/2) & dat$z > z & dat$z <= (z+r/2)
c6<-dat$x > (x) & dat$x <= (x+r/2) & dat$y > y & dat$y <= (y+r/2) & dat$z > (z+r/2) & dat$z <= (z+r)
c7<-dat$x > (x) & dat$x <= (x+r/2) & dat$y > (y+r/2) & dat$y <= (y+r) & dat$z > (z+r/2) & dat$z <= (z+r)
c8<-dat$x > (x) & dat$x <= (x+r/2) & dat$y > (y+r/2) & dat$y <= (y+r) & dat$z > z & dat$z <= (z+r/2)
p1=c((x+r/2),y,z)
p2=c((x+r/2),y,(z+r/2))
p3=c((x+r/2),(y+r/2),(z+r/2))
p4=c((x+r/2),(y+r/2),z)

p5=c((x),y,z)
p6=c((x),y,(z+r/2))
p7=c((x),(y+r/2),(z+r/2))
p8=c((x),(y+r/2),z)

g1<-dat[c1,] 
g2<-dat[c2,] 
g3<-dat[c3,] 
g4<-dat[c4,] 
g5<-dat[c5,] 
g6<-dat[c6,] 
g7<-dat[c7,] 
g8<-dat[c8,]
if(dim(g1)[1]==0){g1[1,]<-c(-1,-1,-1,-1)}
if(dim(g2)[1]==0){g2[1,]<-c(-1,-1,-1,-1)}
if(dim(g3)[1]==0){g3[1,]<-c(-1,-1,-1,-1)}
if(dim(g4)[1]==0){g4[1,]<-c(-1,-1,-1,-1)}
if(dim(g5)[1]==0){g5[1,]<-c(-1,-1,-1,-1)}
if(dim(g6)[1]==0){g6[1,]<-c(-1,-1,-1,-1)}
if(dim(g7)[1]==0){g7[1,]<-c(-1,-1,-1,-1)}
if(dim(g8)[1]==0){g8[1,]<-c(-1,-1,-1,-1)}

if(dim(g1)[1]<nmax){
g1$gid=gid
g1$xr=p1[1]+(r/4)
g1$yr=p1[2]+(r/4)
g1$zr=p1[3]+(r/4)
g1$rr=(r/2)
gid=gid+1
}else{
g1<-splitBox((r/2),p1[1],p1[2],p1[3],dat,gid,nmax)
### Get new group_id ######
gid=max(g1$gid)+1
}

if(dim(g2)[1]<nmax){
g2$gid=gid
g2$xr=p2[1]+(r/4)
g2$yr=p2[2]+(r/4)
g2$zr=p2[3]+(r/4)
g2$rr=(r/2)
gid=gid+1
}else{
g2<-splitBox((r/2),p2[1],p2[2],p2[3],dat,gid,nmax)
### Get new group_id ######
gid=max(g2$gid)+1
}

if(dim(g3)[1]<nmax){
g3$gid=gid
g3$xr=p3[1]+(r/4)
g3$yr=p3[2]+(r/4)
g3$zr=p3[3]+(r/4)
g3$rr=(r/2)
gid=gid+1
}else{
g3<-splitBox((r/2),p3[1],p3[2],p3[3],dat,gid,nmax)
### Get new group_id ######
gid=max(g3$gid)+1
}

if(dim(g4)[1]<nmax){
g4$gid=gid
g4$xr=p4[1]+(r/4)
g4$yr=p4[2]+(r/4)
g4$zr=p4[3]+(r/4)
g4$rr=(r/2)
gid=gid+1
}else{
g4<-splitBox((r/2),p4[1],p4[2],p4[3],dat,gid,nmax)
### Get new group_id ######
gid=max(g4$gid)+1
}
 
if(dim(g5)[1]<nmax){
g5$gid=gid
g5$xr=p5[1]+(r/4)
g5$yr=p5[2]+(r/4)
g5$zr=p5[3]+(r/4)
g5$rr=(r/2)
gid=gid+1
}else{
g5<-splitBox((r/2),p5[1],p5[2],p5[3],dat,gid,nmax)
### Get new group_id ######
gid=max(g5$gid)+1
}

if(dim(g6)[1]<nmax){
g6$gid=gid
g6$xr=p6[1]+(r/4)
g6$yr=p6[2]+(r/4)
g6$zr=p6[3]+(r/4)
g6$rr=(r/2)
gid=gid+1
}else{
g6<-splitBox((r/2),p6[1],p6[2],p6[3],dat,gid,nmax)
### Get new group_id ######
gid=max(g6$gid)+1
}

if(dim(g7)[1]<nmax){
g7$gid=gid
g7$xr=p7[1]+(r/4)
g7$yr=p7[2]+(r/4)
g7$zr=p7[3]+(r/4)
g7$rr=(r/2)
gid=gid+1
}else{
g7<-splitBox((r/2),p7[1],p7[2],p7[3],dat,gid,nmax)
### Get new group_id ######
gid=max(g7$gid)+1
}
 
if(dim(g8)[1]<nmax){
g8$gid=gid
g8$xr=p8[1]+(r/4)
g8$yr=p8[2]+(r/4)
g8$zr=p8[3]+(r/4)
g8$rr=(r/2)
gid=gid+1
}else{
g8<-splitBox((r/2),p8[1],p8[2],p8[3],dat,gid,nmax)
### Get new group_id ######
gid=max(g8$gid)+1
}

out<-rbind(g1,g2,g3,g4,g5,g6,g7,g8)
out<-as.data.frame(out)
out<-out[out$id!=-1,]
out
}



grid(file1,300)

#grid(file2,300)
#grid(file3,300)
#grid(file4,300)
#grid(file5,300)
#grid(file6,300)
#grid(file7,300)
#grid(file8,300)
#grid(file1r,300)
#grid(file2r,300)
#grid(file3r,300)
#grid(file4r,300)
#grid(file5r,300)
#grid(file6r,300)
#grid(file7r,300)
#grid(file8r,300)

