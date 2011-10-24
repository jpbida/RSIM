### Read in each table and subset frag_db, frag_nt, frag_at, and frag_id ###
### Put each database in a new directory                                 ###
print("Reading Frag_db")
frag_db=read.table(file="./frag_db",sep="\t")
print(dim(frag_db))
print("Reading Frag_ID")
frag_id=read.table(file="./frag_id",sep="\t")
frag_id[,2]=as.character(frag_id[,2])
frag_id[,3]=as.numeric(as.character(frag_id[,3]))
frag_id[,4]=as.character(frag_id[,4])
frag_id[,5]=as.numeric(as.character(frag_id[,5]))
print(dim(frag_id))
print("Reading Frag_NT")
frag_nt=read.table(file="./frag_nt",sep=",")
print(dim(frag_nt))
print("Reading Frag_AT")
frag_at=read.table(file="./frag_at",sep=",")
print(dim(frag_at))
print("Reading ffrags.vect")
ffrag=read.table(file="./ffrags.vect",sep=",")
print(dim(ffrag))
print("Reading rfrags.vect")
rfrag=read.table(file="./rfrags.vect",sep=",")
print(dim(rfrag))

all_frag_dbs=NULL
frag_num_db=1
### Creates several sets of fragment databases that are based on different distributions of fragment types ###
system("cat ./pdbs/*_combo.txt > all_combos.txt")
combos<-read.table(file="all_combos.txt",sep=",")

### Keep the first model ###

combos_m1=combos[combos$V2==0,]

### Change from nuc to py/pur  ###

combos_m1$V5=gsub("A","R",combos_m1$V5)
combos_m1$V5=gsub("G","R",combos_m1$V5)
combos_m1$V5=gsub("C","Y",combos_m1$V5)
combos_m1$V5=gsub("U","Y",combos_m1$V5)

frag_type=names(table(combos_m1$V5))
### Make fragment databases of various sizes ###
### 200,300, .... 1000 ###
helix=combos_m1[combos_m1$V6==1,]
loop=combos_m1[combos_m1$V6==2,]
bulge=combos_m1[combos_m1$V6==0,]
for(total in seq(200,1000,by=100)){
for(ratio1 in seq(0,1,by=.25)){
for(ratio2 in seq(0,1,by=.25)){
### Randomly Pick Fragments ###
loop_tot=round(ratio1*total)
helix_tot=round((1-ratio1)*total*ratio2)
bulge_tot=round((1-ratio1)*total*(1-ratio2))
print(paste(loop_tot,helix_tot,bulge_tot,sep=":"))
frag_tmp=NULL
start=1
for(f in frag_type){
loop_f=loop[loop$V5==f,]
helix_f=helix[helix$V5==f,]
bulge_f=bulge[bulge$V5==f,]
l_frags=round(runif(loop_tot)*dim(loop_f)[1])
h_frags=round(runif(helix_tot)*dim(helix_f)[1])
b_frags=round(runif(bulge_tot)*dim(bulge_f)[1])
if(start==1){
frag_tmp=rbind(loop_f[l_frags,])
frag_tmp=rbind(frag_tmp,helix_f[h_frags,])
frag_tmp=rbind(frag_tmp,bulge_f[b_frags,])
start=0
}else{

frag_tmp=rbind(frag_tmp,loop_f[l_frags,])
frag_tmp=rbind(frag_tmp,helix_f[h_frags,])
frag_tmp=rbind(frag_tmp,bulge_f[b_frags,])
}
}
### Create Frag Database ###
new_dir=paste("./data/frag_dbs",frag_num_db,sep="")
dir.create(paste("./data/frag_dbs",frag_num_db,sep=""))
write.table(frag_tmp,file=paste("./data/frag_dbs",frag_num_db,"/frag_types",sep=""),row.names=FALSE,col.names=FALSE)
all_frag_dbs$total[frag_num_db]=total
all_frag_dbs$ratio1[frag_num_db]=ratio1
all_frag_dbs$ratio2[frag_num_db]=ratio2
frag_num_db=frag_num_db+1
### Subset frag_db,frag_id,... ####
sub_vect=NULL
w=1
for(j in 1:dim(frag_tmp)[1]){
pdb_id=as.character(frag_tmp[j,1])
model=as.numeric(as.character(frag_tmp[j,2]))
chain=as.character(frag_tmp[j,3])
start=as.numeric(as.character(frag_tmp[j,4]))

v=frag_id[as.character(frag_id[,2])==pdb_id & as.numeric(frag_id[,3])==model & as.character(frag_id[,4])==chain & as.numeric(frag_id[,5])==start,1]+1
if(length(v)==1){
sub_vect[w]=v
w=w+1
}else{
print(paste("Can't find: ",pdb_id," ",model," ",chain," ",start))
}


}

frag_id1=frag_id[sub_vect,]
frag_nt1=frag_nt[sub_vect,]
frag_at1=frag_at[sub_vect,]
frag_db1=frag_db[sub_vect,]
ffrag1  =ffrag[sub_vect,]
rfrag1  =rfrag[sub_vect,]
### Write lib files ###
libs=paste(new_dir,"/libs",sep="")
dir.create(libs)
k1=sub_vect[as.character(frag_tmp[,5])=="RRR"]
k2=sub_vect[as.character(frag_tmp[,5])=="RRY"]
k3=sub_vect[as.character(frag_tmp[,5])=="RYR"]
k4=sub_vect[as.character(frag_tmp[,5])=="RYY"]
k5=sub_vect[as.character(frag_tmp[,5])=="YRR"]
k6=sub_vect[as.character(frag_tmp[,5])=="YRY"]
k7=sub_vect[as.character(frag_tmp[,5])=="YYR"]
k8=sub_vect[as.character(frag_tmp[,5])=="YYY"]

write.table(ffrag[k1,],file=paste(libs,"/fRRR.vect",sep=""),sep=",",row.names=FALSE,col.names=FALSE)
write.table(ffrag[k2,],file=paste(libs,"/fRRY.vect",sep=""),sep=",",row.names=FALSE,col.names=FALSE)
write.table(ffrag[k3,],file=paste(libs,"/fRYR.vect",sep=""),sep=",",row.names=FALSE,col.names=FALSE)
write.table(ffrag[k4,],file=paste(libs,"/fRYY.vect",sep=""),sep=",",row.names=FALSE,col.names=FALSE)
write.table(ffrag[k5,],file=paste(libs,"/fYRR.vect",sep=""),sep=",",row.names=FALSE,col.names=FALSE)
write.table(ffrag[k6,],file=paste(libs,"/fYRY.vect",sep=""),sep=",",row.names=FALSE,col.names=FALSE)
write.table(ffrag[k7,],file=paste(libs,"/fYYR.vect",sep=""),sep=",",row.names=FALSE,col.names=FALSE)
write.table(ffrag[k8,],file=paste(libs,"/fYYY.vect",sep=""),sep=",",row.names=FALSE,col.names=FALSE)

write.table(rfrag[k1,],file=paste(libs,"/rRRR.vect",sep=""),sep=",",row.names=FALSE,col.names=FALSE)
write.table(rfrag[k2,],file=paste(libs,"/rRRY.vect",sep=""),sep=",",row.names=FALSE,col.names=FALSE)
write.table(rfrag[k3,],file=paste(libs,"/rRYR.vect",sep=""),sep=",",row.names=FALSE,col.names=FALSE)
write.table(rfrag[k4,],file=paste(libs,"/rRYY.vect",sep=""),sep=",",row.names=FALSE,col.names=FALSE)
write.table(rfrag[k5,],file=paste(libs,"/rYRR.vect",sep=""),sep=",",row.names=FALSE,col.names=FALSE)
write.table(rfrag[k6,],file=paste(libs,"/rYRY.vect",sep=""),sep=",",row.names=FALSE,col.names=FALSE)
write.table(rfrag[k7,],file=paste(libs,"/rYYR.vect",sep=""),sep=",",row.names=FALSE,col.names=FALSE)
write.table(rfrag[k8,],file=paste(libs,"/rYYY.vect",sep=""),sep=",",row.names=FALSE,col.names=FALSE)

write.table(frag_id1,file=paste(new_dir,"/frag_id",sep=""),row.names=FALSE,col.names=FALSE,sep="\t")
write.table(frag_at1,file=paste(new_dir,"/frag_at",sep=""),row.names=FALSE,col.names=FALSE,sep=",")
write.table(frag_nt1,file=paste(new_dir,"/frag_nt",sep=""),row.names=FALSE,col.names=FALSE,sep=",")
write.table(frag_db1,file=paste(new_dir,"/frag_db",sep=""),row.names=FALSE,col.names=FALSE,sep="\t")


}
}
}


