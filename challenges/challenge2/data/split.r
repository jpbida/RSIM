source("~/projects/rpgms/ct2coord.r")
source("~/projects/rpgms/aptplot.r")

ss="(((...((((....))))..(((...)))..)))"
ct1=makeCt(ss,ss)
ct1=as.data.frame(ct1)
 ct=pseudoKnot(ct1)
 ct=ct[[2]]
 ctc=ct2coord(ct)
 out=NULL
 mx_group=max(ctc$group,na.rm=TRUE)
 w=1
 for(i in 0:mx_group){
 ct1=getComponent(ctc,i)

ct1=ct1[ct1$group==i,]
 for(j in 1:dim(ct1)[1]){
 out$comp[w]=i
 out$pos[w]=ct1$pos[j]
 w=w+1
 }

}




