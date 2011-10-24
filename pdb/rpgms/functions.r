#### Given a four points 5' to 3' it calculates the given torsional angle betweent the two middle points ###

torsion<-function(a,b,c,d){
#### b to c is pointing in positive z axis ###
#### a to b is in positive y z axis ###

	
}

##### R file that condenses bp files genereated by bpRNA.pl into structure files ####
###### Write a function that recongnizes each type of base-pairing ####
ntAU<-function(ldat){
out<-0
	ntAU1<-0
	ldat$nt1<-as.character(ldat$nt1)
	ldat$nt2<-as.character(ldat$nt2)
	ldat$atom2<-as.character(ldat$atom2)
	ldat$atom1<-as.character(ldat$atom1)
#### Returns W-C pairing yes/no ####
b1<-ldat$nt1=="  A" & ldat$atom1==" N1 " & ldat$atom2==" N3 " & ldat$nt2=="  U"
b2<-ldat$nt1=="  A" & ldat$atom1==" N6 " & ldat$atom2==" O4 " & ldat$nt2=="  U"
b3<-ldat$nt1=="  U" & ldat$atom1==" N3 " & ldat$atom2==" N1 " & ldat$nt2=="  A"
b4<-ldat$nt1=="  U" & ldat$atom1==" O4 " & ldat$atom2==" N6 " & ldat$nt2=="  A"

#### Finding TRUES ####
#### Both T1 and T2 need to be true ###
t1<-table(b1)
t2<-table(b2)
t1<-names(t1)
t2<-names(t2)
c1<-bpT(t1,t2)
if(c1 ==2){
#### Get the position and chainID of the partner ####
ld1<-ldat[b1,]
ld2<-ldat[b2,]
ld3<-rbind(ld1,ld2)
tld<-table(ld3$pos2)
ps1<-names(tld)[tld==max(tld)]
ps1<-ps1[1]
ret<-ld3[ld3$pos2==ps1,]
ret<-ret[1,]
out<-paste(ret$pdb,ret$nt1,ret$nt2,ret$pos1,ret$pos2,ret$chain1,ret$chain2,sep=",")
ntAU1=1
##### Return NT1, NT2,POS1, POS2,Chain1, Chain2 ######
}else{
# or 
### Both T3 and T4 need to be true ###
t3<-table(b3)
t4<-table(b4)
t3<-names(t3)
t4<-names(t4)
c2<-bpT(t3,t4)
if(c2==2){
ld1<-ldat[b3,]
ld2<-ldat[b4,]
ld3<-rbind(ld1,ld2)
tld<-table(ld3$pos2)
ps1<-names(tld)[tld==max(tld)]
ps1<-ps1[1]
ret<-ld3[ld3$pos2==ps1,]
ret<-ret[1,]
out<-paste(ret$pdb,ret$nt1,ret$nt2,ret$pos1,ret$pos2,ret$chain1,ret$chain2,sep=",")
ntAU1=1
}
	}
out
}

ntGC<-function(ldat){
	ntAU1<-0
out<-0
	ldat$nt1<-as.character(ldat$nt1)
	ldat$nt2<-as.character(ldat$nt2)
	ldat$atom2<-as.character(ldat$atom2)
	ldat$atom1<-as.character(ldat$atom1)
#### Returns W-C pairing yes/no ####
b1<-ldat$nt1=="  G" & ldat$atom1==" N2 " & ldat$atom2==" O2 " & ldat$nt2=="  C"
b2<-ldat$nt1=="  G" & ldat$atom1==" N1 " & ldat$atom2==" N3 " & ldat$nt2=="  C"
b3<-ldat$nt1=="  G" & ldat$atom1==" O6 " & ldat$atom2==" N4 " & ldat$nt2=="  C"
b4<-ldat$nt1=="  C" & ldat$atom1==" O2 " & ldat$atom2==" N2 " & ldat$nt2=="  G"
b5<-ldat$nt1=="  C" & ldat$atom1==" N3 " & ldat$atom2==" N1 " & ldat$nt2=="  G"
b6<-ldat$nt1=="  C" & ldat$atom1==" N4 " & ldat$atom2==" O6 " & ldat$nt2=="  G"
#### Finding TRUES ####
#### Both T1 and T2 need to be true ###
t1<-table(b1)
t2<-table(b2)
t3<-table(b3)

t1<-names(t1)
t2<-names(t2)
t3<-names(t3)

c1<-bpT3(t1,t2,t3)

if(c1 ==3){
#### Get the position and chainID of the partner ####
ld1<-ldat[b1,]
ld2<-ldat[b2,]
ld3<-ldat[b3,]
ld4<-rbind(ld1,ld2,ld3)
tld<-table(ld4$pos2)
ps1<-names(tld)[tld==max(tld)]
if(length(ps1)>1){
ps1<-ps1[1]
}
ret<-ld4[ld4$pos2==ps1,]
ret<-ret[1,]
out<-paste(ret$pdb,ret$nt1,ret$nt2,ret$pos1,ret$pos2,ret$chain1,ret$chain2,sep=",")
ntAU1=1
}else{
# or 
### Both T3 and T4 need to be true ###
t4<-table(b4)
t5<-table(b5)
t6<-table(b6)

t4<-names(t4)
t5<-names(t5)
t6<-names(t6)

c2<-bpT3(t4,t5,t6)
if(c2==3){
ld1<-ldat[b4,]
ld2<-ldat[b5,]
ld3<-ldat[b6,]
ld4<-rbind(ld1,ld2,ld3)
tld<-table(ld4$pos2)
ps1<-names(tld)[tld==max(tld)]
if(length(ps1)>1){
ps1<-ps1[1]
}
ret<-ld4[ld4$pos2==ps1,]
ret<-ret[1,]
out<-paste(ret$pdb,ret$nt1,ret$nt2,ret$pos1,ret$pos2,ret$chain1,ret$chain2,sep=",")
ntAU1=1
}
	}
out
}

bpT3<-function(t1,t2,t3){
### Returns 3 if t1, t2, and t3 are true ###
case1<-0
for(p in 1:length(t1)){
if(t1[p]=="TRUE"){case1=1}
}
if(case1==1){
for(p in 1:length(t2)){
if(t2[p]=="TRUE"){case1=2}
}
if(case1==2){
for(p in 1:length(t3)){
if(t3[p]=="TRUE"){case1=3}
}
}
}else{
case1=0
}
case1
}
bpT<-function(t1,t2){
### Returns 2 if both t1 and t2 are true ###
case1<-0
for(p in 1:length(t1)){
if(t1[p]=="TRUE"){case1=1}
}
if(case1==1){
for(p in 1:length(t2)){
if(t2[p]=="TRUE"){case1=2}
}
}else{
case1=0
}
case1
}

