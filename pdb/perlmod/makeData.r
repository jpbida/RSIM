#### Residues that correspond to RNA nucleotides ####
#55    U  C1' C2  C2' C3' C4  C4' C5  C5' C6  N1  N3  O2  O2' O3' O4  O4'
#71    U  C1' C2  C2' C3' C4  C4' C5  C5' C6  N1  N3  O2  O2' O3' O4  O4' O5'
#8     U  C1' C2  C2' C3' C4  C4' C5  C5' C6  N1  N3  O2  O2' O3' O4  O4' O5' OP1 OP2
#57    U  C1' C2  C2' C3' C4  C4' C5  C5' C6  N1  N3  O2  O2' O3' O4  O4' O5' OP1 OP2 OP3 P  
#43    U  C1' C2  C2' C3' C4  C4' C5  C5' C6  N1  N3  O2  O2' O3' O4  O4' O5' OP1 OP2 P  
#
#12    G  C1' C2  C2' C3' C4  C4' C5  C5' C6  C8  N1  N2  N3  N7  N9  O2' O3' O4' O5' O6 
#9     G  C1' C2  C2' C3' C4  C4' C5  C5' C6  C8  N1  N2  N3  N7  N9  O2' O3' O4' O5' O6  OP1 OP2
#14    G  C1' C2  C2' C3' C4  C4' C5  C5' C6  C8  N1  N2  N3  N7  N9  O2' O3' O4' O5' O6  OP1 OP2 OP3 P  
#15    G  C1' C2  C2' C3' C4  C4' C5  C5' C6  C8  N1  N2  N3  N7  N9  O2' O3' O4' O5' O6  OP1 OP2 P  
#127   G  C1' C2  C2' C3' C4  C4' C5  C5' C6  C8  N1  N2  N3  N7  N9  O2' O3' O4' O5' O6  P  
#13    G  C1' C2  C2' C3' C4  C4' C5  C5' C6  C8  N1  N2  N3  N7  N9  O2' O3' O4' O6 
#
#56    A  C1' C2  C2' C3' C4  C4' C5  C5' C6  C8  N1  N3  N6  N7  N9  O2' O3' O4'
#51    A  C1' C2  C2' C3' C4  C4' C5  C5' C6  C8  N1  N3  N6  N7  N9  O2' O3' O4' O5'
#7     A  C1' C2  C2' C3' C4  C4' C5  C5' C6  C8  N1  N3  N6  N7  N9  O2' O3' O4' O5' OP1 OP2
#42    A  C1' C2  C2' C3' C4  C4' C5  C5' C6  C8  N1  N3  N6  N7  N9  O2' O3' O4' O5' OP1 OP2 P  
#
#6     C  C1' C2  C2' C3' C4  C4' C5  C5' C6  N1  N3  N4  O2  O2' O3' O4'
#11    C  C1' C2  C2' C3' C4  C4' C5  C5' C6  N1  N3  N4  O2  O2' O3' O4' O5'
#10    C  C1' C2  C2' C3' C4  C4' C5  C5' C6  N1  N3  N4  O2  O2' O3' O4' O5' OP1 OP2

if(makenew==0){
dat<-read.table(file="chains.txt",sep=" ")
map<-read.table(file="res.map",sep=" ")
## Good residues ##
rna_res<-c(55,71,8,57,43,12,9,14,15,127,13,56,51,7,42,6,11,10)

bad_residues<-c(97,137,153,86,130,174,172,175,146,184,147,183,187,155,129,173,171,168,98,148,152,132,170,176,169,177,145,144,134,154,151,167,131,178)
map[bad_residues,3]<-"UNK"

chainType<-function(chain,map){
out<-NULL
out$rna=0
out$dna=0
out$protein=0
out$unk=0
for(i in 1:(dim(map)[1]-1)){
j=i+3
if(chain[j]>0){
if(as.character(map[i,3])=="RNA"){out$rna[1]=out$rna[1]+as.numeric(chain[j])}
if(as.character(map[i,3])=="DNA"){out$dna[1]=out$dna[1]+as.numeric(chain[j])}
if(as.character(map[i,3])=="PROTEIN"){out$protein[1]=out$protein[1]+as.numeric(chain[j])}
if(as.character(map[i,3])=="UNK"){out$unk[1]=out$unk[1]+as.numeric(chain[j])}
}
}
out
}
### Determine the length of each chain  ##
all<-NULL
start=0
for(i in 1:dim(dat)[1]){
print(i)
out<-chainType(dat[i,],map)
out<-as.data.frame(out)
out<-cbind(dat[i,1:2],out)
if(start==0){
start=1
all<-out
}else{all<-rbind(all,out)}
}

save(all,file="pdb.stats")
### Determine the type of the chain     ##
### Don't consider hybrid chains  ###
### Foreach PDB label as small-RNA, protein/small-RNA, protein/Large RNA, RNA/DNA
t<-table(all[,1])
pdb_ids<-names(t)

### pdb_id,RNA_Chains,DNA_Chains,Protein_Chains,UNK_Chains,Min_RNA_Chain,Max_RNA_Chain
pdb_stats<-NULL
start=0

for(pdb in pdb_ids){
tmp<-all[all[,1]==pdb,]
stats<-NULL
stats$pdb_id=pdb
stats$hybrid=0
stats$rna=0
stats$dna=0
stats$protein=0
stats$unk=0
stats$min_rna=9999
stats$max_rna=0
stats$tot_chains=dim(tmp)[1]
stats$rna_chain=""
for(j in 1:dim(tmp)[1]){
k<-tmp[j,3:6]!=0
types=rep(1,4)[k]
if(sum(types)>1){
chain_type="HYBRID"
stats$hybrid[1]=stats$hybrid[1]+1
}else{
chain_type=names(tmp)[3:6]
ctype=as.character(chain_type[k])
if(ctype=="dna"){stats$dna[1]=stats$dna[1]+1}
if(ctype=="rna"){
stats$rna_chain=tmp$V2[j]
stats$rna[1]=stats$rna[1]+1
if(stats$min_rna[1] > tmp$rna[j]){stats$min_rna[1]=tmp$rna[j]}
if(stats$max_rna[1] < tmp$rna[j]){stats$max_rna[1]=tmp$rna[j]}
}
if(ctype=="protein"){stats$protein[1]=stats$protein[1]+1}
if(ctype=="unk"){stats$unk[1]=stats$unk[1]+1}
}

}

stats<-as.data.frame(stats)
if(start==0){
start=1
pdb_stats=stats
}else{
pdb_stats=rbind(pdb_stats,stats)
}
print(pdb)
}

save(pdb_stats,file="pdb_stats.stats")
}else{
load("pdb_stats.stats")
}
### Keep only those that don't have hybrid chains ###

k<-pdb_stats$hybrid==0
pstats<-pdb_stats[k,]


### Keep only RNA's > 10nts ###

k<-pstats$min_rna > 10
pstats=pstats[k,]

### Keep only single chain RNA's(no duplex rna)
k<-pstats$rna==1
pstats=pstats[k,]

### Remove any structures containing DNA ##
k=pstats$dna==0
pstats<-pstats[k,]

### Classify all RNA's as small RNA, small RNA + Protein,Ribosome ####
g1=pstats$max_rna < 500
g2=pstats$protein>0

srna=pstats[(g1 & !g2),]
sprna=pstats[(g1 & g2),]
ribosomes=pstats[!g1,]

#### Write out file list path, chain ###
### Get RNA chain, use babel to add hydrogens ####
if(gen==0){

i=1
for(r in as.character(srna[,1])){
newfile=gsub("pdb_db/","pdb_db/small_rna/",r)
print(newfile)
system("rm tmp.pdb")
system(paste("./getChain.pl ",r," ",srna$rna_chain[i]," > tmp.pdb",sep=""))
system(paste("mv tmp.pdb ",newfile,sep=""))
#system(paste("babel -h tmp.pdb ",newfile,sep=""))
i=i+1
}

i=1
for(r in as.character(sprna[,1])){
newfile=gsub("pdb_db/","pdb_db/rna_protein/",r)
print(newfile)
system("rm tmp.pdb")
system(paste("./getChain.pl ",r," ",sprna$rna_chain[i]," > tmp.pdb",sep=""))

system(paste("mv tmp.pdb ",newfile,sep=""))
#system(paste("babel -h tmp.pdb ",newfile,sep=""))
i=i+1
}

#i=1
#for(r in as.character(ribosomes[,1])){
#newfile=gsub("pdb_db/","pdb_db/ribos/",r)
#print(newfile)
#system("rm tmp.pdb")
#system(paste("./getChain.pl ",r," ",ribosomes$rna_chain[i]," > tmp.pdb",sep=""))
#system(paste("cp tmp.pdb ",newfile,sep=""))
#i=i+1
#}
}


### Perl program adds and removes all hydrogens ###









