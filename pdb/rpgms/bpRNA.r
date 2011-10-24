source("functions.r")


#### Takes in a file name ###
files<-read.table(file="files2.txt",sep=",",header=FALSE)
fnames<-as.character(files[,2])
k<-files[,1]>0
files<-fnames[k]

for(f in files){
print(f)
	fo<-""
dat<-read.table(file=paste("bps/",f,sep=""),sep=",",header=FALSE)
names(dat)<-c("pdb","chain1","nt1","pos1","atom1","chain2","nt2","pos2","atom2","dist")
### Keep it simple at first only consider W-C basepairing ###
ids<-paste(dat$chain1,"_",dat$pos1,sep="")

tids<-table(ids)
tids<-names(tids)
### Look for watson-crick basepairing ####
for (j in tids){
ldat<-dat[ids==j,]
#### if NT1 == A  ####
out<-ntAU(ldat)
if(out!=0){fo<-paste(fo,out,"\n",sep="")}
#### NT1 == G #####
out<-ntGC(ldat)
if(out!=0){fo<-paste(fo,out,"\n",sep="")}
}
write(fo,file=paste("bps2/",f,sep=""))
}
