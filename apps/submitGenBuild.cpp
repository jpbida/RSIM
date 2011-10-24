/* This program tries to replace all fragments in the movement group with one found in the database */
/* This program should become the template program for Monte-carlo simulations */
/* Inputs for this program
-m Master Structure: Used to calculate score of new cluster centers 
-c Cluster Structure: Structure that is moved
-r Min RMSD: Minimum distance moves must be from the original structure 
-n Number of structures in original file
-N Number of Generations: Repeat the process until N structures have been generated  
-p Output file prefix: 
 pre.structs-All centers that have not been moved can be outputted to a file sorted by distance from the master structure
 pre.scores -RMSD distance from master node
-x Max file structure size -- The program will run until 
Clusters each round and keeps best scores from each cluster
Takes an initial structure and final structure and finds a fragment path between the two */
#include <algorithm>
#include <ctime>
#include <boost/lexical_cast.hpp>
#include <fstream>
#include <iostream>
#include "base.h"
#include "algos.h"
#include "fragment.h"
#include "polymere.h"
#include "fragclass.h"
#include "fragatoms.h"
#include "fragall.h"
#include "vclass.h"
#include <string>
#include <cmath>
#include <vector>
#include <istream>
#include <sstream>
#include "matches.h"
#include "graph.h"
#include "element.h"
#include "parse.h"
#include "options.h"
#include "pdb.h"
#include "fabase.h"

struct PairCompareAsc
{
    bool operator() (const std::pair<float,int>  x, const std::pair<float,int>  y) const {
//To Minimize a score       
return x.first > y.first;
//To maximize a score
//return x.score < y.score;
 }
};


struct PairCompareDsc
{
    bool operator() (const std::pair<float,int>  x, const std::pair<float,int>  y) const {
//To maximize a score
return x.first < y.first;
//return x.score < y.score;
 }
};


int main(int argc, char **argv)
{
using namespace std;
//This program takes the inputs from the genCons.r script and generates 
//configuration files for RNA simulations
Options * opt=new Options(argc,argv);
//Initialize buldge contraints
int len_tot=opt->sequence.length();
//printf("writing scripts\n");
int total_scripts=0;
//Read in parameters for building the structure
//sequence
//total number of structures
//round pos1 pos2 seq
std::vector<int> round;
std::vector<int> pos1;  	//Starting Position of movement group
std::vector<int> pos2;		//End position of movement group
std::vector<int> sec_struct; //base-pairing information for the entire molecule



std::ifstream inputb;
inputb.open(opt->gen_confs.c_str());
if(inputb.fail()){
printf("Could not open general constraints file: %s\n",opt->gen_confs.c_str());
abort();
}else{
string line;
string field;
while(getline(inputb, line)) {
//Read in file with each line being x \t y \t z \t group
istringstream iss(line,istringstream::in);	
	getline(iss,field,'\t');
round.push_back(strtoint(field));
	getline(iss,field,'\t');
pos1.push_back(strtoint(field));
	getline(iss,field,'\t');
pos2.push_back(strtoint(field));
}
cout << round.size() << " rounds of constraints loaded from: " << opt->gen_confs.c_str()<<  endl;
}
inputb.close();
//Correct for pos1 to R conversion 
for(int j=0; j<pos1.size(); j++){
pos1[j]=pos1[j]-1;
pos2[j]=pos2[j]-1;
}

sec_struct=opt->sec_structure;
//outputs
//for(int s=0; s< sec_struct.size(); s++){printf("%d - %d\n",s,sec_struct[s]);}

//Vector of Clash Groups
std::vector< std::vector< std::vector<int> > > MasterClash;

//Vector of movement vectors
std::vector< std::vector< std::vector<int> > > MasterMovePos;

//Vector of general constraints
std::vector< std::vector< std::vector<std::pair<std::pair<int,int>,int> > > > MasterGenCons;

//Vectors for holding temporary values of MasterCons and Master Mvoes 
std::vector< std::vector< std::vector<int> > >   cshTemp;
std::vector< std::vector< std::vector<int> > >   moveTemp;
std::vector<std::vector< std::vector< std::pair<std::pair<int,int>,int> > > > consTemp;

for(int r=0; r<round.size(); r++){
//printf("Round: %d Positions: %d - %d\n",r,pos1[r],pos2[r]);
if(round[r]==1){
//If round = 1 just build it
int slen=abs(pos2[r]-pos1[r]);
if(slen < 5){
//Creating a clash vector
std::vector< std::vector<int> > movecsh;
std::vector<int> csh;
for(int c=0; c < len_tot; c++){csh.push_back(0);}
for(int c=0; c<csh.size(); c++){
if(c <=(pos1[r]+5) && c>=(pos1[r]-1)){
csh[c]=1;
}
}
movecsh.push_back(csh);

std::vector< std::vector<int> > movepos;
std::vector< std::vector< std::pair<std::pair<int,int>,int> > >gencons;

//If some kind of constraint exists use predefined loop database
//Use pos1-1 and pos1-2 with even and odd loops 
//with movepos being every other position
//use hard constraints as genCons

int use_loopdb =0;
std::vector<std::pair<std::pair<int,int>,int > > gcons;
for(int s=0; s < slen; s++){
	if(sec_struct[pos1[r]+s] > 0){
	if(sec_struct[pos1[r]+s] > pos1[r]+s){
		use_loopdb=1;
		std::pair<std::pair<int,int>,int> gcon;
		gcon.first.second = sec_struct[pos1[r]+s];
		gcon.first.first = pos1[r]+s;
		gcon.second=-7;
		gcons.push_back(gcon);
}
}
}

if(use_loopdb==1){
gencons.push_back(gcons);
std::vector<int> pos;
pos.push_back(pos1[r]-1);
pos.push_back(pos1[r]+1);
pos.push_back(pos1[r]+3);
movepos.push_back(pos);
}else{
std::vector<int> pos;
pos.push_back(pos1[r]-1);
pos.push_back(pos1[r]+1);
pos.push_back(pos1[r]+3);
movepos.push_back(pos);
std::pair<std::pair<int,int>,int> gcon;
gcon.first.first=pos1[r];
gcon.first.second=pos2[r];
gcon.second=-1; //Near By constraint
gcons.push_back(gcon);
gencons.push_back(gcons);
}
moveTemp.push_back(movepos);
consTemp.push_back(gencons);
cshTemp.push_back(movecsh);
}else{

//Example output
//s1	s2	side1	side2	lp_p	lp_ps1	lp_ps2	lp_typ	minPos	maxPos	metcons	gencons
//6	19	3	5	0	3	5	5	9	14	0	[9 ]-(-2--2) [8 13 ]-(9-14) [7 14 ]-(8-15) (8-16) [6 16 ]-(7-17) (7-18) 
//Pick midpoints that have > 2:1 ratios and abide by the hard constraints
//Identify midpoints that maintain a 0 balance between constraints
//Midpoint Range
int minPos;
int maxPos;
for(int mpt = pos1[r]; mpt<=pos2[r]; mpt++){
//if max-min < 5 don't update and break out
if(sec_struct[mpt] > -1){if(sec_struct[mpt]-mpt >= 5){minPos=mpt; maxPos=sec_struct[mpt];}}
}

//loop through all possible loop types and positions
//printf("Mid Point Range: %d %d\n",minPos,maxPos);
int lp_len = maxPos - minPos;
printf("Round %d\n",round[r]);
printf("s1\ts2\tside1\tside2\tlp_p\tlp_ps1\tlp_ps2\tlp_typ\tminPos\tmaxPos\tmetcons\tgencons\n");
for(int lp_typ=6; lp_typ <=7; lp_typ++){
for(int lp_p = 0; lp_p <= (lp_len - lp_typ); lp_p++){
//Check 2:1 assym
int lp_ps1=(minPos - pos1[r])+lp_p;
int lp_ps2=(pos2[r] - maxPos)+lp_len - lp_p - lp_typ+1;
float asym = (float)lp_ps1/(float)lp_ps2;
if(asym < 2 && asym > 0.5){
//printf("asym: %d %d: %8.3f\n",lp_ps1,lp_ps2,asym);
//Loop brings together pos1+lp_p 
//Generate all asymmetries that are consistent with the hard constraints
int side1 = lp_ps1; 
int side2 = lp_ps2;
//generate all asymettries
int n=0;
std::vector<int> path; 
int tot_len;
if(side1 > side2){tot_len=side2;}else{tot_len=side1;}
//printf("Total Length: %d\n",tot_len);
for(int n=0; n < tot_len; n++){path.push_back(0);}
path[path.size()-1]=-1;
//Copy the current Move and Cons vectors to temporary holdings
while(n==0){
//Each path is a different bulge conformation
//Set reset the current clash group 
//Each helix of round one is built without considering steric clashes from the other 
//helices in round one

std::vector<int> csh_cur;
for(int c=0; c<len_tot; c++){csh_cur.push_back(0);}


std::vector< std::vector<int> > movepos;
std::vector< std::vector< std::pair<std::pair<int,int>,int > > > gencons;
std::vector< std::vector<int> > movecsh;


//for(int p=0; p<path.size(); p++){
//printf("%d ",path[p]);
//}
//printf("\n");
n = nextPath(&path,side1,side2);
//Check if path conforms to all base-pairings
int s1=pos1[r]+side1-1; //Based on Checks [POSSIBLY WRONG] might be side1-1
int s2=pos2[r]-side2+1;
//push loop position onto array
std::vector<int> lpos;
std::vector< std::pair<std::pair<int,int>,int> >  gconsl;
lpos.push_back((pos1[r]+side1));
std::pair<std::pair<int,int>,int> con_loop;
if(lp_typ==6){
con_loop.first.first=-2;
con_loop.first.second=-2;
con_loop.second=-7;
}else{
con_loop.first.first=-3;
con_loop.first.second=-3;
con_loop.second=-7;
}

std::vector<int> csh;
for(int c=0; c < len_tot; c++){csh.push_back(csh_cur[c]);}
for(int c=0; c<csh.size(); c++){
if(c > (pos1[r]+side1-1) && c<(pos1[r]+side1+lp_typ)){
//printf("%d %d\n",c,csh[c]);
csh[c]=1;
csh_cur[c]=1;
}
}

movecsh.push_back(csh);
movepos.push_back(lpos);
gconsl.push_back(con_loop);
gencons.push_back(gconsl);
int n1=0;
int n2=0;
for(int g1=0; g1 < tot_len; g1++){
//Stop once you get to n1 = side1 n2=side2
if(n1!=side1 && n2!=side2){
std::vector<int> mpos;
std::vector< std::pair<std::pair<int,int>,int> > gcons;
std::pair<std::pair<int,int>,int> gcon1;
std::pair<std::pair<int,int>,int> gcon2;
mpos.push_back(s1-1);
mpos.push_back(s2-1);

//Clash Vector
std::vector<int> csh;
for(int c=0; c < len_tot; c++){csh.push_back(csh_cur[c]);}
for(int c=0; c<csh.size(); c++){
if( (c <= (s1+2) && c >= (s1-1)) || ( c< (s2+2) && c>=(s2-1) )   ){
csh[c]=1;
csh_cur[c]=1;
}
}

movecsh.push_back(csh);
movepos.push_back(mpos);
if(path[g1]==0){
gcon1.first.first=s1;
gcon1.first.second=s2;
gcon1.second=-7;
gcons.push_back(gcon1);
s1--;
s2++;
n1++;
n2++;
}

if(path[g1]==1){
gcon2.first.first=s1-1;
gcon2.first.second=s2;
gcon2.second=-9;
gcons.push_back(gcon2);
s2++;
s1--;
s1--;
n1++;
n1++;
n2++;
}

if(path[g1]==2){
gcon2.first.first=s1;
gcon2.first.second=s2+1;
gcon2.second=-9;
gcons.push_back(gcon2);
s2++;
s2++;
s1--;
n1++;
n2++;
n2++;
}

gencons.push_back(gcons);
}
}
//check to make sure gencons contains all bases pairings in the 
//hard constraints
int metcons=0;
for(int ss1=pos1[r]; ss1 <= pos2[r]; ss1++){
int sc2=sec_struct[ss1];
metcons=0;
if(sc2 >=0 && sc2 > ss1){
for(int gc1=0; gc1 < gencons.size(); gc1++){
for(int gc2=0; gc2 < gencons[gc1].size(); gc2++){
//printf("Con Check (%d %d) (%d %d)\n",ss1,sc2,gencons[gc1][gc2].first,gencons[gc1][gc2].second);
if(gencons[gc1][gc2].first.first==ss1 && gencons[gc1][gc2].first.second==sc2){
metcons=1;
break;
}
}
if(metcons==1){break;}
}
if(metcons==0){break;}
}else{
metcons=1;
}
}

if(metcons==1){
//Foreach of the current movepositions push_back the new constraints
if(MasterMovePos.size() > 0){
for(int m=0; m < MasterMovePos.size(); m++){
//Create a new entry that adds to the existing vectors
std::vector< std::vector<int> > tmov=MasterMovePos[m];
for(int mnew=0; mnew < movepos.size(); mnew++){tmov.push_back(movepos[mnew]);}
std::vector< std::vector<int> > tcsh=MasterClash[m];
for(int mnew=0; mnew < movecsh.size(); mnew++){tcsh.push_back(movecsh[mnew]);}
std::vector< std::vector< std::pair<std::pair<int,int>,int> > > tcon=MasterGenCons[m];
for(int mnew=0; mnew < gencons.size(); mnew++){tcon.push_back(gencons[mnew]);}
moveTemp.push_back(tmov);
consTemp.push_back(tcon);
cshTemp.push_back(tcsh);
}
}else{

std::vector< std::vector<int> > tmov;
for(int mnew=0; mnew < movepos.size(); mnew++){tmov.push_back(movepos[mnew]);}
std::vector< std::vector<int> > tcsh;
for(int mnew=0; mnew < movecsh.size(); mnew++){tcsh.push_back(movecsh[mnew]);}
std::vector<std::vector< std::pair<std::pair<int,int>,int> > >  tcon;
for(int mnew=0; mnew < gencons.size(); mnew++){tcon.push_back(gencons[mnew]);}
cshTemp.push_back(tcsh);
moveTemp.push_back(tmov);
consTemp.push_back(tcon);
}

printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t",s1,s2,side1,side2,lp_p,lp_ps1,lp_ps2,lp_typ,minPos,maxPos,metcons);
for(int g=0; g<gencons.size(); g++){
printf("[");
for(int m1=0; m1<movepos[g].size(); m1++){printf("%d ",movepos[g][m1]);}
printf("]-");
for(int g2=0; g2<gencons[g].size(); g2++){
printf("(%d-%d:%d) ",gencons[g][g2].first.first,gencons[g][g2].first.second,gencons[g][g2].second);}

for(int p=0; p<movecsh[g].size(); p++){
printf("%d ",movecsh[g][p]);
}
printf("\n");
}
printf("\n");

}
}
}


}
}


}
 

}else{
//Reset Clash Vector
//csh_all contains all groups that have been moved up until this round

std::vector<int> csh_cur;
for(int c=0; c<len_tot; c++){csh_cur.push_back(0);}

std::vector<int> csh_prv;
for(int c=0; c<len_tot; c++){csh_prv.push_back(0);}

//	make vector of linear positions
std::vector<int> mol_pos;
std::vector<int> mol_pos2;

//	push back everything between pos1 and pos2 excluding things from other groups

for(int p=pos1[r]; p<=pos2[r]; p++){
int keep1=1;
int keep2=1;

for(int r_cur=0; r_cur < round.size(); r_cur++){
if(round[r_cur] < round[r] && p > pos1[r_cur] && p < pos2[r_cur]){keep1=0;}
if(round[r_cur] < round[r] && p >= pos1[r_cur] && p <= pos2[r_cur]){keep2=0;}
}

//    Add previous round positions to csh_cur 


if(keep1==1){mol_pos.push_back(p);}
if(keep1==0){csh_prv[p]=1;}
if(keep2==1){mol_pos2.push_back(p);}



} 
//Need to add pseudo atom foreach base pair to standardize asym generation
//printf("Path2a: ");
//for(int p=0; p<mol_pos.size(); p++){
//printf("%d ",mol_pos[p]);
//}
//printf("\n");
for(int mp=0; mp<(mol_pos.size()-1); mp++){
if((mol_pos[mp+1]-mol_pos[mp])>1){
//printf("Found Break: %d %d\n",mol_pos[mp],mol_pos[mp+1]);
mol_pos.insert(mol_pos.begin()+(mp+1),1,-1);
mp++;
}
}
//printf("\nPath2b: ");
//for(int p=0; p<mol_pos.size(); p++){
//printf("%d ",mol_pos[p]);
//}
//printf("\n");
int minPos;
int maxPos;
int mpt2;


for(int mpt1 = 0; mpt1<mol_pos2.size(); mpt1++){
//loops through mol_pos2, need to map it back to mol_pos
int mpt=mol_pos2[mpt1];
//if max-min < 5 don't update and break out
//only one loop position possible
if(sec_struct[mpt] > -1){
//Determine how far apart the two are in the chain
int dist=0;
int past_it=0;
for(int mlen=0; mlen<mol_pos.size(); mlen++){
if(mol_pos[mlen]==sec_struct[mpt]){break;}else{
if(past_it ==1){dist++;}else{
if(mol_pos[mlen]==mpt){past_it=1;}
}
}
}
if(dist >= 5){
minPos=mpt; maxPos=sec_struct[mpt];
}

}

}

//Map minPos, maxPos to mol_pos vector indices
for(int mp3=0; mp3<mol_pos.size(); mp3++){
if(mol_pos[mp3]==minPos){minPos=mp3;}
if(mol_pos[mp3]==maxPos){maxPos=mp3;}
}
//Determine inner constraints ignoring previous groups
//printf("Inner Constraints: %d - %d\n",minPos,maxPos);
//	loop through all asymmetries 
//for(int p=0; p<mol_pos.size(); p++){
//printf("%d ",mol_pos[p]);
//}
//printf("\n");
//	Identify potential mid-points
//	Generate constraints
//	Remove all constraints containing the base-pair and its pseudo atom
//	Replace with clusterCons constraints
int slen = mol_pos.size()-1;
if(slen < 5){
std::vector< std::vector<int> > movepos;
std::vector< std::vector<int> > movecsh;
std::vector< std::vector< std::pair<std::pair<int,int>,int> > >gencons;
//If some kind of constraint exists use predefined loop database
//Use pos1-1 and pos1-2 with even and odd loops 
//with movepos being every other position
//use hard constraints as genCons
int use_loopdb =0;
std::vector<std::pair<std::pair<int,int>,int> > gcons;
for(int s=0; s < slen; s++){
//It it isn't part of the multi-chain base pairs
if(s > 0){
if(mol_pos[s+1]!=-1 && mol_pos[s-1]!=-1){
if(sec_struct[mol_pos[s]] > 0){
if(sec_struct[mol_pos[s]] > mol_pos[s]){
use_loopdb=1;
std::pair<std::pair<int,int>,int> gcon;
gcon.first.second = sec_struct[mol_pos[s]];
gcon.first.first = mol_pos[s];
gcon.second=-7;
gcons.push_back(gcon);
}
}
}
}else{
if(sec_struct[mol_pos[s]] > 0){
if(sec_struct[mol_pos[s]] > mol_pos[s]){
use_loopdb=1;
std::pair<std::pair<int,int>,int> gcon;
gcon.first.second = sec_struct[mol_pos[s]];
gcon.first.first = mol_pos[s];
gcon.second=-7;
gcons.push_back(gcon);
}
}
}

}

if(use_loopdb==1){
gencons.push_back(gcons);
std::vector<int> pos;
pos.push_back(mol_pos[0]-1);
pos.push_back(mol_pos[0]+2);
pos.push_back(mol_pos[0]+3);
//Set current clash vector and clash vector
movepos.push_back(pos);

std::vector<int> csh;
for(int c=0; c<csh_cur.size(); c++){
csh.push_back(csh_cur[c]);
if((c >=(mol_pos[0]-1) && c<= (mol_pos[0]+5))){
csh[c]=1;
csh_cur[c]=1;
}
}
movecsh.push_back(csh);



}else{
//Use clusterCons with clash minimization
std::vector<int> pos;
pos.push_back(mol_pos[0]-1);
pos.push_back(mol_pos[0]+2);
pos.push_back(mol_pos[0]+3);
movepos.push_back(pos);
std::pair<std::pair<int,int>,int> gcon;
gcon.first.first=mol_pos[0];
gcon.first.second=mol_pos[(mol_pos.size()-1)];
gcon.second=-1; //Near By constraint
gcons.push_back(gcon);
gencons.push_back(gcons);
}

moveTemp.push_back(movepos);
consTemp.push_back(gencons);
}else{
//Pick midpoints that have > 2:1 ratios and abide by the hard constraints
//Identify midpoints that maintain a 0 balance between constraints
//Midpoint Range
//printf("Mid Point Range: %d %d\n",minPos,maxPos);
int lp_len = maxPos - minPos;
printf("Round %d\n",round[r]);
printf("s1\ts2\tside1\tside2\tlp_p\tlp_ps1\tlp_ps2\tlp_typ\tminPos\tmaxPos\tmetcons\tgencons\n");
for(int lp_typ=6; lp_typ <=7; lp_typ++){
for(int lp_p = 0; lp_p <= (lp_len - lp_typ); lp_p++){
//Check 2:1 assym
int lp_ps1=(minPos)+lp_p;  //MIGHT NEED TO BE lp_p +1 [WRONG??]
int lp_ps2=(mol_pos.size() - maxPos)+lp_len - lp_p - lp_typ + 1;//NO + 1??
float asym = (float)lp_ps1/(float)lp_ps2;
if(asym < 2 && asym > 0.5){
//printf("asym: %d %d: %8.3f\n",lp_ps1,lp_ps2,asym);
//Loop brings together pos1+lp_p 
//Generate all asymmetries that are consistent with the hard constraints
int side1 = lp_ps1; 
int side2 = lp_ps2;
//generate all asymettries
int n=0;
std::vector<int> path; 
int tot_len;
if(side1 > side2){tot_len=side2;}else{tot_len=side1;}
//printf("Total Length: %d %d %d %d\n",tot_len,side1,side2,mol_pos.size());
for(int n1=0; n1 < tot_len; n1++){path.push_back(0);}
path[path.size()-1]=-1;
//Copy the current Move and Cons vectors to temporary holdings
n=0;
while(n==0){
std::vector< std::vector<int> > movepos;
std::vector< std::vector< std::pair<std::pair<int,int>,int> > > gencons;

//printf("Path %d %d: ",side1,side2);
//for(int p=0; p<path.size(); p++){
//printf("%d ",path[p]);
//}
//printf(":");
//for(int p=0; p<mol_pos.size(); p++){
//printf("%d ",mol_pos[p]);
//}
//printf("\n");
n = nextPath(&path,side1,side2);
if(n==1){break;}
//Check if path conforms to all base-pairings
int s1=side1-1;
int s2=mol_pos.size()-side2;
//push loop position onto array
std::vector<int> lpos;
std::vector< std::pair<std::pair<int,int>,int> >  gconsl;
std::pair<std::pair<int,int>,int> con_loop;
if(lp_typ==6){
//This is a junction loop 
//Base pairings can occur inside this loop 
//in which case clustered movement groups 
//need to be created 

int bp_pos1=0;
int bp_pos2=0;
int bp_pos3=0;
int bp_pos4=0;
int bp_pos5=0;
int bp_pos6=0;
int bp_pos7=0;
int bp_pos8=0;
if(mol_pos[side1-2]==-1){bp_pos1=1;}
if(mol_pos[side1-1]==-1){bp_pos2=1;}
if(mol_pos[side1+0]==-1){bp_pos3=1;}
if(mol_pos[side1+1]==-1){bp_pos4=1;}
if(mol_pos[side1+2]==-1){bp_pos5=1;}
if(mol_pos[side1+3]==-1){bp_pos6=1;}
if(mol_pos[side1+4]==-1){bp_pos7=1;}
if(mol_pos[side1+5]==-1){bp_pos8=1;}

if(bp_pos1==0 && bp_pos2==0 && bp_pos3==0 && bp_pos4==0 && bp_pos5==0 && bp_pos6==0 && bp_pos7==0 && bp_pos8==0){
//Simple loop with no base-pairings inside 
con_loop.first.first=-2;
con_loop.first.second=-2;
con_loop.second=-7;  //Need to consider if these are base-paired BUG
gconsl.push_back(con_loop);
lpos.push_back(mol_pos[side1-1]);
movepos.push_back(lpos);
gencons.push_back(gconsl);

}else{ 
if(bp_pos1==1 && bp_pos2==0 && bp_pos3==0 && bp_pos4==0 && bp_pos5==0 && bp_pos6==0 && bp_pos7==0 && bp_pos8==0){
//First position is base paired 
con_loop.first.first=mol_pos[side1-3];
con_loop.first.second=mol_pos[side1+5];
con_loop.second=-1;
gconsl.push_back(con_loop);
con_loop.first.first=mol_pos[side1-1];
con_loop.first.second=mol_pos[side1+5];
con_loop.second=-1;
gconsl.push_back(con_loop);

lpos.push_back(mol_pos[side1-1]-1);
lpos.push_back(mol_pos[side1+1]);
lpos.push_back(mol_pos[side1+3]);
lpos.push_back(mol_pos[side1+5]);
movepos.push_back(lpos);
gencons.push_back(gconsl);
}else{ 
if(bp_pos1==0 && bp_pos2==1 && bp_pos3==0 && bp_pos4==0 && bp_pos5==0 && bp_pos6==0 && bp_pos7==0 && bp_pos8==0){
con_loop.first.first=mol_pos[side1-2];
con_loop.first.second=mol_pos[side1+5];
con_loop.second=-1;
gconsl.push_back(con_loop);

lpos.push_back(mol_pos[side1-2]-1);
lpos.push_back(mol_pos[side1]-1);
lpos.push_back(mol_pos[side1+2]);
lpos.push_back(mol_pos[side1+4]);
movepos.push_back(lpos);
gencons.push_back(gconsl);


}else{ 
if(bp_pos1==0 && bp_pos2==0 && bp_pos3==1 && bp_pos4==0 && bp_pos5==0 && bp_pos6==0 && bp_pos7==0 && bp_pos8==0){
con_loop.first.first=mol_pos[side1-2];
con_loop.first.second=mol_pos[side1+5];
con_loop.second=-1;
gconsl.push_back(con_loop);
lpos.push_back(mol_pos[side1-2]);
lpos.push_back(mol_pos[side1-1]-1);
lpos.push_back(mol_pos[side1+1]-1);
lpos.push_back(mol_pos[side1+2]);
lpos.push_back(mol_pos[side1+4]);
movepos.push_back(lpos);
gencons.push_back(gconsl);
}else{ 
///      -2             -1           0               1             2             3           4               5
if(bp_pos1==0 && bp_pos2==0 && bp_pos3==0 && bp_pos4==1 && bp_pos5==0 && bp_pos6==0 && bp_pos7==0 && bp_pos8==0){
con_loop.first.first=mol_pos[side1-2];
con_loop.first.second=mol_pos[side1+5];
con_loop.second=-1;
gconsl.push_back(con_loop);

lpos.push_back(mol_pos[side1-2]);
lpos.push_back(mol_pos[side1]-1);
lpos.push_back(mol_pos[side1+2]-1);
lpos.push_back(mol_pos[side1+3]);
movepos.push_back(lpos);
gencons.push_back(gconsl);


}else{ 
///      -2             -1           0               1             2             3           4               5
if(bp_pos1==0 && bp_pos2==0 && bp_pos3==0 && bp_pos4==0 && bp_pos5==1 && bp_pos6==0 && bp_pos7==0 && bp_pos8==0){

con_loop.first.first=mol_pos[side1-2];
con_loop.first.second=mol_pos[side1+5];
con_loop.second=-1;
gconsl.push_back(con_loop);

lpos.push_back(mol_pos[side1-2]);
lpos.push_back(mol_pos[side1]);
lpos.push_back(mol_pos[side1+3]-1);
lpos.push_back(mol_pos[side1+4]);
movepos.push_back(lpos);
gencons.push_back(gconsl);

}else{ 
///      -2             -1           0               1             2             3           4               5
if(bp_pos1==0 && bp_pos2==0 && bp_pos3==0 && bp_pos4==0 && bp_pos5==0 && bp_pos6==1 && bp_pos7==0 && bp_pos8==0){
con_loop.first.first=mol_pos[side1-2];
con_loop.first.second=mol_pos[side1+5];
con_loop.second=-1;
gconsl.push_back(con_loop);

lpos.push_back(mol_pos[side1-2]);
lpos.push_back(mol_pos[side1+0]);
lpos.push_back(mol_pos[side1+2]);
lpos.push_back(mol_pos[side1+4]-1);
movepos.push_back(lpos);
gencons.push_back(gconsl);

}else{ 
///      -2             -1           0               1             2             3           4               5
if(bp_pos1==0 && bp_pos2==0 && bp_pos3==0 && bp_pos4==0 && bp_pos5==0 && bp_pos6==0 && bp_pos7==1 && bp_pos8==0){

con_loop.first.first=mol_pos[side1-2];
con_loop.first.second=mol_pos[side1+5];
con_loop.second=-1;
gconsl.push_back(con_loop);

lpos.push_back(mol_pos[side1-2]);
lpos.push_back(mol_pos[side1+0]);
lpos.push_back(mol_pos[side1+2]);
lpos.push_back(mol_pos[side1+5]-1);
movepos.push_back(lpos);
gencons.push_back(gconsl);

}else{ 
///      -2             -1           0               1             2             3           4               5
if(bp_pos1==0 && bp_pos2==0 && bp_pos3==0 && bp_pos4==0 && bp_pos5==0 && bp_pos6==0 && bp_pos7==0 && bp_pos8==1){
con_loop.first.first=mol_pos[side1-2];
con_loop.first.second=mol_pos[side1+6];
con_loop.second=-1;
gconsl.push_back(con_loop);
con_loop.first.first=mol_pos[side1-2];
con_loop.first.second=mol_pos[side1+4];
con_loop.second=-1;
gconsl.push_back(con_loop);

lpos.push_back(mol_pos[side1-2]);
lpos.push_back(mol_pos[side1+0]);
lpos.push_back(mol_pos[side1+2]);
lpos.push_back(mol_pos[side1+4]);
movepos.push_back(lpos);
gencons.push_back(gconsl);

}else{ 
///      -2             -1           0               1             2             3           4               5
if(bp_pos1==1 && bp_pos2==0 && bp_pos3==0 && bp_pos4==1 && bp_pos5==0 && bp_pos6==0 && bp_pos7==0 && bp_pos8==0){
con_loop.first.first=mol_pos[side1-3];
con_loop.first.second=mol_pos[side1+5];
con_loop.second=-1;
gconsl.push_back(con_loop);
con_loop.first.first=mol_pos[side1-1];
con_loop.first.second=mol_pos[side1+5];
con_loop.second=-1;
gconsl.push_back(con_loop);

lpos.push_back(mol_pos[side1-3]-1);
lpos.push_back(mol_pos[side1-1]-1);
lpos.push_back(mol_pos[side1+0]);
lpos.push_back(mol_pos[side1+2]-1);
lpos.push_back(mol_pos[side1+3]);
movepos.push_back(lpos);
gencons.push_back(gconsl);

}else{ 
///      -2             -1           0               1             2             3           4               5
if(bp_pos1==1 && bp_pos2==0 && bp_pos3==0 && bp_pos4==0 && bp_pos5==1 && bp_pos6==0 && bp_pos7==0 && bp_pos8==0){
con_loop.first.first=mol_pos[side1-3];
con_loop.first.second=mol_pos[side1+5];
con_loop.second=-1;
gconsl.push_back(con_loop);
con_loop.first.first=mol_pos[side1-1];
con_loop.first.second=mol_pos[side1+5];
con_loop.second=-1;
gconsl.push_back(con_loop);
lpos.push_back(mol_pos[side1-3]-1);
lpos.push_back(mol_pos[side1-1]-1);
lpos.push_back(mol_pos[side1+0]);
lpos.push_back(mol_pos[side1+3]-1);
lpos.push_back(mol_pos[side1+4]);
movepos.push_back(lpos);
gencons.push_back(gconsl);

}else{ 
///      -2             -1           0               1             2             3           4               5
if(bp_pos1==1 && bp_pos2==0 && bp_pos3==0 && bp_pos4==0 && bp_pos5==0 && bp_pos6==1 && bp_pos7==0 && bp_pos8==0){
con_loop.first.first=mol_pos[side1-3];
con_loop.first.second=mol_pos[side1+5];
con_loop.second=-1;
gconsl.push_back(con_loop);
con_loop.first.first=mol_pos[side1-1];
con_loop.first.second=mol_pos[side1+5];
con_loop.second=-1;
gconsl.push_back(con_loop);

lpos.push_back(mol_pos[side1-3]-1);
lpos.push_back(mol_pos[side1-1]-1);
lpos.push_back(mol_pos[side1+0]);
lpos.push_back(mol_pos[side1+2]);
lpos.push_back(mol_pos[side1+4]-1);
movepos.push_back(lpos);
gencons.push_back(gconsl);

}else{ 
///      -2             -1           0               1             2             3           4               5
if(bp_pos1==1 && bp_pos2==0 && bp_pos3==0 && bp_pos4==0 && bp_pos5==0 && bp_pos6==0 && bp_pos7==1 && bp_pos8==0){
con_loop.first.first=mol_pos[side1-3];
con_loop.first.second=mol_pos[side1+5];
con_loop.second=-1;
gconsl.push_back(con_loop);
con_loop.first.first=mol_pos[side1-1];
con_loop.first.second=mol_pos[side1+5];
con_loop.second=-1;
gconsl.push_back(con_loop);

lpos.push_back(mol_pos[side1-3]-1);
lpos.push_back(mol_pos[side1-1]-1);
lpos.push_back(mol_pos[side1+0]);
lpos.push_back(mol_pos[side1+2]);
lpos.push_back(mol_pos[side1+5]-1);
movepos.push_back(lpos);
gencons.push_back(gconsl);

}else{ 
///      -2             -1           0               1             2             3           4               5
if(bp_pos1==1 && bp_pos2==0 && bp_pos3==0 && bp_pos4==0 && bp_pos5==0 && bp_pos6==0 && bp_pos7==0 && bp_pos8==1){
//Try to stack the two base-pairings
con_loop.first.first=mol_pos[side1-3];
con_loop.first.second=mol_pos[side1+6];
con_loop.second=-10;
gconsl.push_back(con_loop);
con_loop.first.first=mol_pos[side1-1];
con_loop.first.second=mol_pos[side1+4];
con_loop.second=-10;
gconsl.push_back(con_loop);

lpos.push_back(mol_pos[side1-3]-1);
lpos.push_back(mol_pos[side1-1]-1);
lpos.push_back(mol_pos[side1+0]);
lpos.push_back(mol_pos[side1+2]);
lpos.push_back(mol_pos[side1+4]-1);
lpos.push_back(mol_pos[side1+6]-1);
movepos.push_back(lpos);
gencons.push_back(gconsl);
}else{ 
///      -2             -1           0               1             2             3           4               5
if(bp_pos1==0 && bp_pos2==1 && bp_pos3==0 && bp_pos4==0 && bp_pos5==1 && bp_pos6==0 && bp_pos7==0 && bp_pos8==0){
con_loop.first.first=mol_pos[side1-2];
con_loop.first.second=mol_pos[side1+5];
con_loop.second=-1;
gconsl.push_back(con_loop);

lpos.push_back(mol_pos[side1-3]-1);
lpos.push_back(mol_pos[side1-1]-1);
lpos.push_back(mol_pos[side1+0]);
lpos.push_back(mol_pos[side1+2]);
lpos.push_back(mol_pos[side1+4]-1);
lpos.push_back(mol_pos[side1+6]-1);
movepos.push_back(lpos);
gencons.push_back(gconsl);

}else{ 
///      -2             -1           0               1             2             3           4               5
if(bp_pos1==0 && bp_pos2==1 && bp_pos3==0 && bp_pos4==0 && bp_pos5==0 && bp_pos6==1 && bp_pos7==0 && bp_pos8==0){
con_loop.first.first=mol_pos[side1-2];
con_loop.first.second=mol_pos[side1+5];
con_loop.second=-1;
gconsl.push_back(con_loop);

lpos.push_back(mol_pos[side1-3]);
lpos.push_back(mol_pos[side1-0]-1);
lpos.push_back(mol_pos[side1+1]);
lpos.push_back(mol_pos[side1+4]-1);
lpos.push_back(mol_pos[side1+5]);
movepos.push_back(lpos);
gencons.push_back(gconsl);
}else{ 
///      -2             -1           0               1             2             3           4               5
if(bp_pos1==0 && bp_pos2==1 && bp_pos3==0 && bp_pos4==0 && bp_pos5==0 && bp_pos6==0 && bp_pos7==1 && bp_pos8==0){
con_loop.first.first=mol_pos[side1-2];
con_loop.first.second=mol_pos[side1+5];
con_loop.second=-1;
gconsl.push_back(con_loop);

lpos.push_back(mol_pos[side1-2]-1);
lpos.push_back(mol_pos[side1-0]-1);
lpos.push_back(mol_pos[side1+1]);
lpos.push_back(mol_pos[side1+3]);
lpos.push_back(mol_pos[side1+5]-1);
movepos.push_back(lpos);
gencons.push_back(gconsl);

}else{ 
///      -2             -1           0               1             2             3           4               5
if(bp_pos1==0 && bp_pos2==1 && bp_pos3==0 && bp_pos4==0 && bp_pos5==0 && bp_pos6==0 && bp_pos7==0 && bp_pos8==1){
con_loop.first.first=mol_pos[side1-2];
con_loop.first.second=mol_pos[side1+6];
con_loop.second=-1;
gconsl.push_back(con_loop);
con_loop.first.first=mol_pos[side1-2];
con_loop.first.second=mol_pos[side1+4];
con_loop.second=-1;
gconsl.push_back(con_loop);

lpos.push_back(mol_pos[side1-2]-1);
lpos.push_back(mol_pos[side1-0]-1);
lpos.push_back(mol_pos[side1+1]);
lpos.push_back(mol_pos[side1+3]);
lpos.push_back(mol_pos[side1+6]-1);
movepos.push_back(lpos);
gencons.push_back(gconsl);
}else{ 
///      -2             -1           0               1             2             3           4               5
if(bp_pos1==0 && bp_pos2==0 && bp_pos3==1 && bp_pos4==0 && bp_pos5==0 && bp_pos6==1 && bp_pos7==0 && bp_pos8==0){
con_loop.first.first=mol_pos[side1-2];
con_loop.first.second=mol_pos[side1+5];
con_loop.second=-1;
gconsl.push_back(con_loop);

lpos.push_back(mol_pos[side1-2]);
lpos.push_back(mol_pos[side1+1]-1);
lpos.push_back(mol_pos[side1+2]);
lpos.push_back(mol_pos[side1+4]-1);
lpos.push_back(mol_pos[side1+5]);
movepos.push_back(lpos);
gencons.push_back(gconsl);

}else{ 
///      -2             -1           0               1             2             3           4               5
if(bp_pos1==0 && bp_pos2==0 && bp_pos3==1 && bp_pos4==0 && bp_pos5==0 && bp_pos6==0 && bp_pos7==1 && bp_pos8==0){
con_loop.first.first=mol_pos[side1-2];
con_loop.first.second=mol_pos[side1+5];
con_loop.second=-1;
gconsl.push_back(con_loop);

lpos.push_back(mol_pos[side1-2]);
lpos.push_back(mol_pos[side1+1]-1);
lpos.push_back(mol_pos[side1+2]);
lpos.push_back(mol_pos[side1+3]);
lpos.push_back(mol_pos[side1+5]-1);
movepos.push_back(lpos);
gencons.push_back(gconsl);

}else{ 
///      -2             -1           0               1             2             3           4               5
if(bp_pos1==0 && bp_pos2==0 && bp_pos3==1 && bp_pos4==0 && bp_pos5==0 && bp_pos6==0 && bp_pos7==0 && bp_pos8==1){
con_loop.first.first=mol_pos[side1-2];
con_loop.first.second=mol_pos[side1+6];
con_loop.second=-1;
gconsl.push_back(con_loop);
con_loop.first.first=mol_pos[side1-2];
con_loop.first.second=mol_pos[side1+4];
con_loop.second=-1;
gconsl.push_back(con_loop);

lpos.push_back(mol_pos[side1-2]);
lpos.push_back(mol_pos[side1+1]-1);
lpos.push_back(mol_pos[side1+2]);
lpos.push_back(mol_pos[side1+4]);
movepos.push_back(lpos);
gencons.push_back(gconsl);
}else{ 
///      -2             -1           0               1             2             3           4               5
if(bp_pos1==0 && bp_pos2==0 && bp_pos3==0 && bp_pos4==1 && bp_pos5==0 && bp_pos6==0 && bp_pos7==1 && bp_pos8==0){

con_loop.first.first=mol_pos[side1-2];
con_loop.first.second=mol_pos[side1+5];
con_loop.second=-1;
gconsl.push_back(con_loop);

lpos.push_back(mol_pos[side1-2]);
lpos.push_back(mol_pos[side1+0]);
lpos.push_back(mol_pos[side1+2]-1);
lpos.push_back(mol_pos[side1+3]);
lpos.push_back(mol_pos[side1+5]-1);
movepos.push_back(lpos);
gencons.push_back(gconsl);

}else{ 
///      -2             -1           0               1             2             3           4               5
if(bp_pos1==0 && bp_pos2==0 && bp_pos3==0 && bp_pos4==1 && bp_pos5==0 && bp_pos6==0 && bp_pos7==0 && bp_pos8==1){
con_loop.first.first=mol_pos[side1-2];
con_loop.first.second=mol_pos[side1+6];
con_loop.second=-1;
gconsl.push_back(con_loop);
con_loop.first.first=mol_pos[side1-2];
con_loop.first.second=mol_pos[side1+4];
con_loop.second=-1;
gconsl.push_back(con_loop);

lpos.push_back(mol_pos[side1-2]);
lpos.push_back(mol_pos[side1+0]);
lpos.push_back(mol_pos[side1+2]-1);
lpos.push_back(mol_pos[side1+3]);
lpos.push_back(mol_pos[side1+4]);
movepos.push_back(lpos);
gencons.push_back(gconsl);


}else{ 
///      -2             -1           0               1             2             3           4               5
if(bp_pos1==0 && bp_pos2==0 && bp_pos3==0 && bp_pos4==0 && bp_pos5==1 && bp_pos6==0 && bp_pos7==0 && bp_pos8==1){
con_loop.first.first=mol_pos[side1-2];
con_loop.first.second=mol_pos[side1+6];
con_loop.second=-1;
gconsl.push_back(con_loop);
con_loop.first.first=mol_pos[side1-2];
con_loop.first.second=mol_pos[side1+4];
con_loop.second=-1;
gconsl.push_back(con_loop);

lpos.push_back(mol_pos[side1-2]);
lpos.push_back(mol_pos[side1+0]);
lpos.push_back(mol_pos[side1+1]);
lpos.push_back(mol_pos[side1+3]-1);
lpos.push_back(mol_pos[side1+4]);
movepos.push_back(lpos);
gencons.push_back(gconsl);

}else{ 
///      -2             -1           0               1             2             3           4               5
if(bp_pos1==1 && bp_pos2==0 && bp_pos3==0 && bp_pos4==1 && bp_pos5==0 && bp_pos6==0 && bp_pos7==1 && bp_pos8==0){
//Try to stack the two base-pairings
con_loop.first.first=mol_pos[side1-1];
con_loop.first.second=mol_pos[side1+5];
con_loop.second=-1;
gconsl.push_back(con_loop);
con_loop.first.first=mol_pos[side1-3];
con_loop.first.second=mol_pos[side1+5];
con_loop.second=-1;
gconsl.push_back(con_loop);

lpos.push_back(mol_pos[side1-3]-1);
lpos.push_back(mol_pos[side1-1]-1);
lpos.push_back(mol_pos[side1+0]);
lpos.push_back(mol_pos[side1+2]-1);
lpos.push_back(mol_pos[side1+3]);
lpos.push_back(mol_pos[side1+5]-1);
movepos.push_back(lpos);
gencons.push_back(gconsl);

}else{ 
///      -2             -1           0               1             2             3           4               5
if(bp_pos1==1 && bp_pos2==0 && bp_pos3==0 && bp_pos4==1 && bp_pos5==0 && bp_pos6==0 && bp_pos7==0 && bp_pos8==1){
//Try to stack the two base-pairings
con_loop.first.first=mol_pos[side1-3];
con_loop.first.second=mol_pos[side1+6];
con_loop.second=-10;
gconsl.push_back(con_loop);
con_loop.first.first=mol_pos[side1-1];
con_loop.first.second=mol_pos[side1+4];
con_loop.second=-10;
gconsl.push_back(con_loop);

lpos.push_back(mol_pos[side1-3]-1);
lpos.push_back(mol_pos[side1-1]-1);
lpos.push_back(mol_pos[side1+0]);
lpos.push_back(mol_pos[side1+2]-1);
lpos.push_back(mol_pos[side1+3]);
lpos.push_back(mol_pos[side1+6]-1);
movepos.push_back(lpos);
gencons.push_back(gconsl);


}else{ 
///      -2             -1           0               1             2             3           4               5
if(bp_pos1==1 && bp_pos2==0 && bp_pos3==0 && bp_pos4==0 && bp_pos5==1 && bp_pos6==0 && bp_pos7==0 && bp_pos8==1){
//Try to stack the two base-pairings
con_loop.first.first=mol_pos[side1-3];
con_loop.first.second=mol_pos[side1+6];
con_loop.second=-10;
gconsl.push_back(con_loop);
con_loop.first.first=mol_pos[side1-1];
con_loop.first.second=mol_pos[side1+4];
con_loop.second=-10;
gconsl.push_back(con_loop);

lpos.push_back(mol_pos[side1-3]-1);
lpos.push_back(mol_pos[side1-1]-1);
lpos.push_back(mol_pos[side1+0]);
lpos.push_back(mol_pos[side1+3]-1);
lpos.push_back(mol_pos[side1+4]);
lpos.push_back(mol_pos[side1+6]-1);
movepos.push_back(lpos);
gencons.push_back(gconsl);

}else{ 
///      -2             -1           0               1             2             3           4               5
if(bp_pos1==0 && bp_pos2==1 && bp_pos3==0 && bp_pos4==0 && bp_pos5==1 && bp_pos6==0 && bp_pos7==0 && bp_pos8==1){
con_loop.first.first=mol_pos[side1-2];
con_loop.first.second=mol_pos[side1+6];
con_loop.second=-1;
gconsl.push_back(con_loop);
con_loop.first.first=mol_pos[side1-2];
con_loop.first.second=mol_pos[side1+4];
con_loop.second=-1;
gconsl.push_back(con_loop);

lpos.push_back(mol_pos[side1-2]-1);
lpos.push_back(mol_pos[side1-0]-1);
lpos.push_back(mol_pos[side1+1]);
lpos.push_back(mol_pos[side1+3]-1);
lpos.push_back(mol_pos[side1+4]);
lpos.push_back(mol_pos[side1+6]-1);
movepos.push_back(lpos);
gencons.push_back(gconsl);


}else{printf("Bad Base-paing constraints\n");} 
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}







}else{
int bp_pos1=0;
int bp_pos2=0;
int bp_pos3=0;
int bp_pos4=0;
int bp_pos5=0;
int bp_pos6=0;
int bp_pos7=0;
if(mol_pos[side1-2]==-1){bp_pos1=1;}
if(mol_pos[side1-1]==-1){bp_pos2=1;}
if(mol_pos[side1+0]==-1){bp_pos3=1;}
if(mol_pos[side1+1]==-1){bp_pos4=1;}
if(mol_pos[side1+2]==-1){bp_pos5=1;}
if(mol_pos[side1+3]==-1){bp_pos6=1;}
if(mol_pos[side1+4]==-1){bp_pos7=1;}
// Need to consider base pairs in odd loops
//UPDATE
///      -2             -1           0               1             2             3           4            
if(bp_pos1==1 && bp_pos2==0 && bp_pos3==0 && bp_pos4==0 && bp_pos5==0 && bp_pos6==0 && bp_pos7==0){
//First position is base paired 
con_loop.first.first=mol_pos[side1-3];
con_loop.first.second=mol_pos[side1+4];
con_loop.second=-1;
gconsl.push_back(con_loop);
con_loop.first.first=mol_pos[side1-1];
con_loop.first.second=mol_pos[side1+4];
con_loop.second=-1;
gconsl.push_back(con_loop);

lpos.push_back(mol_pos[side1-1]-1);
lpos.push_back(mol_pos[side1+0]);
lpos.push_back(mol_pos[side1+2]);
lpos.push_back(mol_pos[side1+4]);
movepos.push_back(lpos);
gencons.push_back(gconsl);

}else{
///      -2             -1           0               1             2             3           4            
if(bp_pos1==1 && bp_pos2==0 && bp_pos3==0 && bp_pos4==1 && bp_pos5==0 && bp_pos6==0 && bp_pos7==0){
//First position is base paired 
con_loop.first.first=mol_pos[side1-3];
con_loop.first.second=mol_pos[side1+4];
con_loop.second=-1;
gconsl.push_back(con_loop);
con_loop.first.first=mol_pos[side1-1];
con_loop.first.second=mol_pos[side1+4];
con_loop.second=-1;
gconsl.push_back(con_loop);

lpos.push_back(mol_pos[side1-1]-1);
lpos.push_back(mol_pos[side1+0]);
lpos.push_back(mol_pos[side1+2]-1);
lpos.push_back(mol_pos[side1+3]);
movepos.push_back(lpos);
gencons.push_back(gconsl);


}else{
///      -2             -1           0               1             2             3           4            
if(bp_pos1==1 && bp_pos2==0 && bp_pos3==0 && bp_pos4==0 && bp_pos5==1 && bp_pos6==0 && bp_pos7==0){
//First position is base paired 
con_loop.first.first=mol_pos[side1-3];
con_loop.first.second=mol_pos[side1+4];
con_loop.second=-1;
gconsl.push_back(con_loop);
con_loop.first.first=mol_pos[side1-1];
con_loop.first.second=mol_pos[side1+4];
con_loop.second=-1;
gconsl.push_back(con_loop);

lpos.push_back(mol_pos[side1-1]-1);
lpos.push_back(mol_pos[side1+0]);
lpos.push_back(mol_pos[side1+3]-1);
lpos.push_back(mol_pos[side1+4]);
movepos.push_back(lpos);
gencons.push_back(gconsl);



}else{
///      -2             -1           0               1             2             3           4            
if(bp_pos1==1 && bp_pos2==0 && bp_pos3==0 && bp_pos4==0 && bp_pos5==0 && bp_pos6==1 && bp_pos7==0){

//First position is base paired 
con_loop.first.first=mol_pos[side1-3];
con_loop.first.second=mol_pos[side1+4];
con_loop.second=-1;
gconsl.push_back(con_loop);
con_loop.first.first=mol_pos[side1-1];
con_loop.first.second=mol_pos[side1+4];
con_loop.second=-1;
gconsl.push_back(con_loop);

lpos.push_back(mol_pos[side1-1]-1);
lpos.push_back(mol_pos[side1+0]);
lpos.push_back(mol_pos[side1+2]);
lpos.push_back(mol_pos[side1+4]-1);
movepos.push_back(lpos);
gencons.push_back(gconsl);


}else{
///      -2             -1           0               1             2             3           4            
if(bp_pos1==1 && bp_pos2==0 && bp_pos3==0 && bp_pos4==0 && bp_pos5==0 && bp_pos6==0 && bp_pos7==1){

//First position is base paired 
con_loop.first.first=mol_pos[side1-3];
con_loop.first.second=mol_pos[side1+5];
con_loop.second=-1;
gconsl.push_back(con_loop);

lpos.push_back(mol_pos[side1-1]-1);
lpos.push_back(mol_pos[side1+0]);
lpos.push_back(mol_pos[side1+2]);
lpos.push_back(mol_pos[side1+5]-1);
movepos.push_back(lpos);
gencons.push_back(gconsl);
}else{
///      -2             -1           0               1             2             3           4            
if(bp_pos1==1 && bp_pos2==0 && bp_pos3==0 && bp_pos4==1 && bp_pos5==0 && bp_pos6==0 && bp_pos7==1){
//First position is base paired 
con_loop.first.first=mol_pos[side1-3];
con_loop.first.second=mol_pos[side1+5];
con_loop.second=-1;
gconsl.push_back(con_loop);
lpos.push_back(mol_pos[side1-1]-1);
lpos.push_back(mol_pos[side1+0]);
lpos.push_back(mol_pos[side1+2]-1);
lpos.push_back(mol_pos[side1+3]);
lpos.push_back(mol_pos[side1+5]-1);
movepos.push_back(lpos);
gencons.push_back(gconsl);

}else{
///      -2             -1           0               1             2             3           4            
if(bp_pos1==0 && bp_pos2==1 && bp_pos3==0 && bp_pos4==0 && bp_pos5==0 && bp_pos6==0 && bp_pos7==0){
con_loop.first.first=mol_pos[side1-2];
con_loop.first.second=mol_pos[side1+4];
con_loop.second=-1;
gconsl.push_back(con_loop);
lpos.push_back(mol_pos[side1-2]-1);
lpos.push_back(mol_pos[side1+0]-1);
lpos.push_back(mol_pos[side1+1]);
lpos.push_back(mol_pos[side1+3]);
movepos.push_back(lpos);
gencons.push_back(gconsl);



}else{
///      -2             -1           0               1             2             3           4            
if(bp_pos1==0 && bp_pos2==1 && bp_pos3==0 && bp_pos4==0 && bp_pos5==1 && bp_pos6==0 && bp_pos7==0){

con_loop.first.first=mol_pos[side1-2];
con_loop.first.second=mol_pos[side1+4];
con_loop.second=-1;
gconsl.push_back(con_loop);
lpos.push_back(mol_pos[side1-2]-1);
lpos.push_back(mol_pos[side1+0]-1);
lpos.push_back(mol_pos[side1+1]);
lpos.push_back(mol_pos[side1+3]-3);
lpos.push_back(mol_pos[side1+4]);
movepos.push_back(lpos);
gencons.push_back(gconsl);


}else{
///      -2             -1           0               1             2             3           4            
if(bp_pos1==0 && bp_pos2==1 && bp_pos3==0 && bp_pos4==0 && bp_pos5==0 && bp_pos6==1 && bp_pos7==0){

con_loop.first.first=mol_pos[side1-2];
con_loop.first.second=mol_pos[side1+4];
con_loop.second=-10; //Stacked 
gconsl.push_back(con_loop);
lpos.push_back(mol_pos[side1-2]-1);
lpos.push_back(mol_pos[side1+0]-1);
lpos.push_back(mol_pos[side1+1]);
lpos.push_back(mol_pos[side1+2]);
lpos.push_back(mol_pos[side1+4]-1);
movepos.push_back(lpos);
gencons.push_back(gconsl);


}else{
///      -2             -1           0               1             2             3           4            
if(bp_pos1==0 && bp_pos2==1 && bp_pos3==0 && bp_pos4==0 && bp_pos5==0 && bp_pos6==0 && bp_pos7==1){
con_loop.first.first=mol_pos[side1-2];
con_loop.first.second=mol_pos[side1+3];
con_loop.second=-1;
gconsl.push_back(con_loop);
con_loop.first.first=mol_pos[side1-2];
con_loop.first.second=mol_pos[side1+5];
con_loop.second=-1;
gconsl.push_back(con_loop);
lpos.push_back(mol_pos[side1-2]-1);
lpos.push_back(mol_pos[side1+0]-1);
lpos.push_back(mol_pos[side1+1]);
lpos.push_back(mol_pos[side1+3]);
lpos.push_back(mol_pos[side1+5]-1);
movepos.push_back(lpos);
gencons.push_back(gconsl);
}else{
///      -2             -1           0               1             2             3           4            
if(bp_pos1==0 && bp_pos2==0 && bp_pos3==1 && bp_pos4==0 && bp_pos5==0 && bp_pos6==0 && bp_pos7==0){
con_loop.first.first=mol_pos[side1-2];
con_loop.first.second=mol_pos[side1+4];
con_loop.second=-1;
gconsl.push_back(con_loop);
lpos.push_back(mol_pos[side1-2]);
lpos.push_back(mol_pos[side1-1]);
lpos.push_back(mol_pos[side1+1]-1);
lpos.push_back(mol_pos[side1+2]);
lpos.push_back(mol_pos[side1+4]);
movepos.push_back(lpos);
gencons.push_back(gconsl);



}else{
///      -2             -1           0               1             2             3           4            
if(bp_pos1==0 && bp_pos2==0 && bp_pos3==1 && bp_pos4==0 && bp_pos5==0 && bp_pos6==1 && bp_pos7==0){

con_loop.first.first=mol_pos[side1-2];
con_loop.first.second=mol_pos[side1+4];
con_loop.second=-1;
gconsl.push_back(con_loop);
lpos.push_back(mol_pos[side1-2]);
lpos.push_back(mol_pos[side1+1]-1);
lpos.push_back(mol_pos[side1+2]);
lpos.push_back(mol_pos[side1+4]-1);
movepos.push_back(lpos);
gencons.push_back(gconsl);



}else{
///      -2             -1           0               1             2             3           4            
if(bp_pos1==0 && bp_pos2==0 && bp_pos3==1 && bp_pos4==0 && bp_pos5==0 && bp_pos6==0 && bp_pos7==1){
con_loop.first.first=mol_pos[side1-2];
con_loop.first.second=mol_pos[side1+5];
con_loop.second=-1;
gconsl.push_back(con_loop);
con_loop.first.first=mol_pos[side1-2];
con_loop.first.second=mol_pos[side1+3];
con_loop.second=-1;
gconsl.push_back(con_loop);
lpos.push_back(mol_pos[side1-2]);
lpos.push_back(mol_pos[side1-1]);
lpos.push_back(mol_pos[side1+1]-1);
lpos.push_back(mol_pos[side1+2]);
lpos.push_back(mol_pos[side1+3]-1);
movepos.push_back(lpos);
gencons.push_back(gconsl);

}else{
///      -2             -1           0               1             2             3           4            
if(bp_pos1==0 && bp_pos2==0 && bp_pos3==0 && bp_pos4==1 && bp_pos5==0 && bp_pos6==0 && bp_pos7==0){
con_loop.first.first=mol_pos[side1-2];
con_loop.first.second=mol_pos[side1+4];
con_loop.second=-1;
gconsl.push_back(con_loop);
lpos.push_back(mol_pos[side1-2]);
lpos.push_back(mol_pos[side1]);
lpos.push_back(mol_pos[side1+2]-1);
lpos.push_back(mol_pos[side1+3]);
movepos.push_back(lpos);
gencons.push_back(gconsl);
}else{
///      -2             -1           0               1             2             3           4            
if(bp_pos1==0 && bp_pos2==0 && bp_pos3==0 && bp_pos4==1 && bp_pos5==0 && bp_pos6==0 && bp_pos7==1){
con_loop.first.first=mol_pos[side1-2];
con_loop.first.second=mol_pos[side1+5];
con_loop.second=-1;
gconsl.push_back(con_loop);
con_loop.first.first=mol_pos[side1-2];
con_loop.first.second=mol_pos[side1+3];
con_loop.second=-1;
gconsl.push_back(con_loop);
lpos.push_back(mol_pos[side1-2]);
lpos.push_back(mol_pos[side1]);
lpos.push_back(mol_pos[side1+2]-1);
lpos.push_back(mol_pos[side1+3]);
lpos.push_back(mol_pos[side1+5]-1);
movepos.push_back(lpos);
gencons.push_back(gconsl);
}else{
///      -2             -1           0               1             2             3           4            
if(bp_pos1==0 && bp_pos2==0 && bp_pos3==0 && bp_pos4==0 && bp_pos5==1 && bp_pos6==0 && bp_pos7==0){
con_loop.first.first=mol_pos[side1-2];
con_loop.first.second=mol_pos[side1+4];
con_loop.second=-1;
gconsl.push_back(con_loop);
lpos.push_back(mol_pos[side1-2]);
lpos.push_back(mol_pos[side1]);
lpos.push_back(mol_pos[side1+1]);
lpos.push_back(mol_pos[side1+3]-1);
lpos.push_back(mol_pos[side1+4]);
movepos.push_back(lpos);
gencons.push_back(gconsl);

}else{
///      -2             -1           0               1             2             3           4            
if(bp_pos1==0 && bp_pos2==0 && bp_pos3==0 && bp_pos4==0 && bp_pos5==0 && bp_pos6==1 && bp_pos7==0){
con_loop.first.first=mol_pos[side1-2];
con_loop.first.second=mol_pos[side1+4];
con_loop.second=-1;
gconsl.push_back(con_loop);
lpos.push_back(mol_pos[side1-2]);
lpos.push_back(mol_pos[side1]);
lpos.push_back(mol_pos[side1+2]);
lpos.push_back(mol_pos[side1+4]-1);
movepos.push_back(lpos);
gencons.push_back(gconsl);

}else{
///      -2             -1           0               1             2             3           4            
if(bp_pos1==0 && bp_pos2==0 && bp_pos3==0 && bp_pos4==0 && bp_pos5==0 && bp_pos6==0 && bp_pos7==1){

con_loop.first.first=mol_pos[side1-2];
con_loop.first.second=mol_pos[side1+3];
con_loop.second=-1;
gconsl.push_back(con_loop);

con_loop.first.first=mol_pos[side1-2];
con_loop.first.second=mol_pos[side1+5];
con_loop.second=-1;
gconsl.push_back(con_loop);

lpos.push_back(mol_pos[side1-2]);
lpos.push_back(mol_pos[side1]);
lpos.push_back(mol_pos[side1+2]);
lpos.push_back(mol_pos[side1+3]);
lpos.push_back(mol_pos[side1+5]-1);
movepos.push_back(lpos);
gencons.push_back(gconsl);

}else{
///      -2             -1           0               1             2             3           4            
if(bp_pos1==0 && bp_pos2==0 && bp_pos3==0 && bp_pos4==0 && bp_pos5==0 && bp_pos6==0 && bp_pos7==0){

//Simple loop with no base-pairings inside 
con_loop.first.first=-3;
con_loop.first.second=-3;
con_loop.second=-7;  //Need to consider if these are base-paired BUG
gconsl.push_back(con_loop);
lpos.push_back(side1-1);
movepos.push_back(lpos);
gencons.push_back(gconsl);


}else{
printf("Bad Constraint File\n");
}

}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}









}


int n1=0;
int n2=0;
for(int g1=0; g1 < tot_len; g1++){
//Stop once you get to n1 = side1 n2=side2
if(n1!=side1 && n2!=side2){
std::vector<int> mpos;
std::vector< std::pair<std::pair<int,int>,int> > gcons;
std::pair<std::pair<int,int>,int> gcon1;
std::pair<std::pair<int,int>,int> gcon2;


int abp1=0;
int abp2=0;
int abp3=0;
int bbp1=0;
int bbp2=0;
int bbp3=0;
int p_type=path[g1];

int s1a=s1-1;
int s1b=s1+1;
if(s1b == mol_pos.size()){s1b=mol_pos.size()-1;}
if(s1a < 0){s1a=0;}

int s2b=s2-1;
int s2a=s2+1;
if(s2a == mol_pos.size()){s2a=mol_pos.size()-1;}
if(s2b < 0){s2b=0;}

if(mol_pos[s1a]==-1){abp1=1;}
if(mol_pos[s1] ==-1){abp2=1;}
if(mol_pos[s1b]==-1){abp3=1;}

if(mol_pos[s2a]==-1){bbp1=1;}
if(mol_pos[s2] ==-1){bbp2=1;}
if(mol_pos[s2b]==-1){bbp3=1;}

if(abp1==1){
mpos.push_back(mol_pos[s1-2]-1);
mpos.push_back(mol_pos[s1]-1);
}else{
if((s1-1) < 0){mpos.push_back(-1);}else{
mpos.push_back(mol_pos[s1-1]);}
}
if(abp2==1){
mpos.push_back(mol_pos[s1-1]-1);
mpos.push_back(mol_pos[s1+1]-1);
}

if(bbp3==1){
mpos.push_back(mol_pos[s2-2]-1);
mpos.push_back(mol_pos[s2]-1);
}else{
mpos.push_back(mol_pos[s2-1]);
}
if(bbp2==1){
mpos.push_back(mol_pos[s2-1]-1);
mpos.push_back(mol_pos[s2+1]-1);
}
movepos.push_back(mpos);


// -1              1            +1 ||       +1         1          -1 //       
if(abp1==1 && abp2==0 && abp3 == 0 && bbp1==0 && bbp2==0 && bbp3==0 && p_type==0){
gcon1.first.first=mol_pos[s1];
gcon1.first.second=mol_pos[s2];
gcon1.second=-1;
gcons.push_back(gcon1);
s1--;s2++;n1++;n2++;
}else{
// -1              1            +1 ||       +1         1          -1 //       
if(abp1==1 && abp2==0 && abp3 == 0 && bbp1==1 && bbp2==0 && bbp3==0 && p_type==0){
gcon1.first.first=mol_pos[s1];
gcon1.first.second=mol_pos[s2];
gcon1.second=-10;//Stacked
gcons.push_back(gcon1);
s1--;s2++;n1++;n2++;
}else{
// -1              1            +1 ||       +1         1          -1 //       
if(abp1==1 && abp2==0 && abp3 == 0 && bbp1==0 && bbp2==1 && bbp3==0 && p_type==0){
gcon1.first.first=mol_pos[s1];
gcon1.first.second=mol_pos[s2+1];
gcon1.second=-10; //stacked
gcons.push_back(gcon1);
s1--;s2++;n1++;n2++;
}else{
// -1              1            +1 ||       +1         1          -1 //       
if(abp1==1 && abp2==0 && abp3 == 0 && bbp1==0 && bbp2==0 && bbp3==1 && p_type==0){
gcon1.first.first=mol_pos[s1];
gcon1.first.second=mol_pos[s2];
gcon1.second=-1;
gcons.push_back(gcon1);
s1--;s2++;n1++;n2++;
}else{
// -1              1            +1 ||       +1         1          -1 //       
if(abp1==0 && abp2==1 && abp3 == 0 && bbp1==0 && bbp2==0 && bbp3==0 && p_type==0){
gcon1.first.first=mol_pos[s1-1];
gcon1.first.second=mol_pos[s2];
gcon1.second=-1;
gcons.push_back(gcon1);
gcon1.first.first=mol_pos[s1+1];
gcon1.first.second=mol_pos[s2];
gcon1.second=-1;
gcons.push_back(gcon1);
s1--;s2++;n1++;n2++;
}else{

// -1              1            +1 ||       +1         1          -1 //       
if(abp1==0 && abp2==1 && abp3 == 0 && bbp1==1 && bbp2==0 && bbp3==0 && p_type==0){
gcon1.first.first=mol_pos[s1-1];
gcon1.first.second=mol_pos[s2];
gcon1.second=-10;//stacked
gcons.push_back(gcon1);
s1--;s2++;n1++;n2++;
}else{
// -1              1            +1 ||       +1         1          -1 //       
if(abp1==0 && abp2==1 && abp3 == 0 && bbp1==0 && bbp2==1 && bbp3==0 && p_type==0){
gcon1.first.first=mol_pos[s1-1];
gcon1.first.second=mol_pos[s2+1];
gcon1.second=-10;//stacked
gcons.push_back(gcon1);
gcon1.first.first=mol_pos[s1+1];
gcon1.first.second=mol_pos[s2-1];
gcon1.second=-10;//stacked
gcons.push_back(gcon1);
s1--;s2++;n1++;n2++;
}else{
// -1              1            +1 ||       +1         1          -1 //       
if(abp1==0 && abp2==1 && abp3 == 0 && bbp1==0 && bbp2==0 && bbp3==1 && p_type==0){
gcon1.first.first=mol_pos[s1+1];
gcon1.first.second=mol_pos[s2];
gcon1.second=-10;//stacked
gcons.push_back(gcon1);
s1--;s2++;n1++;n2++;
}else{

// -1              1            +1 ||       +1         1          -1 //       
if(abp1==0 && abp2==0 && abp3 == 1 && bbp1==0 && bbp2==0 && bbp3==0 && p_type==0){
gcon1.first.first=mol_pos[s1];
gcon1.first.second=mol_pos[s2];
gcon1.second=-1;
gcons.push_back(gcon1);
s1--;s2++;n1++;n2++;
}else{

// -1              1            +1 ||       +1         1          -1 //       
if(abp1==0 && abp2==0 && abp3 == 1 && bbp1==1 && bbp2==0 && bbp3==0 && p_type==0){
gcon1.first.first=mol_pos[s1];
gcon1.first.second=mol_pos[s2];
gcon1.second=-1;
gcons.push_back(gcon1);
s1--;s2++;n1++;n2++;
}else{

// -1              1            +1 ||       +1         1          -1 //       
if(abp1==0 && abp2==0 && abp3 == 1 && bbp1==0 && bbp2==1 && bbp3==0 && p_type==0){
gcon1.first.first=mol_pos[s1];
gcon1.first.second=mol_pos[s2-1];
gcon1.second=-10; //Stacked
gcons.push_back(gcon1);
s1--;s2++;n1++;n2++;
}else{
// -1              1            +1 ||       +1         1          -1 //       
if(abp1==0 && abp2==0 && abp3 == 1 && bbp1==0 && bbp2==0 && bbp3==1 && p_type==0){
gcon1.first.first=mol_pos[s1];
gcon1.first.second=mol_pos[s2];
gcon1.second=-10;//stacked
gcons.push_back(gcon1);
s1--;s2++;n1++;n2++;
}else{
// -1              1            +1 ||       +1         1          -1 //       
if(abp1==0 && abp2==0 && abp3 == 0 && bbp1==0 && bbp2==0 && bbp3==0 && p_type==0){
gcon1.first.first=mol_pos[s1];
gcon1.first.second=mol_pos[s2];
gcon1.second=-9; //planar
gcons.push_back(gcon1);
s1--;s2++;n1++;n2++;
}else{
// -1              1            +1 ||       +1         1          -1 //       
if(abp1==0 && abp2==0 && abp3 == 0 && bbp1==1 && bbp2==0 && bbp3==0 && p_type==0){
gcon1.first.first=mol_pos[s1];
gcon1.first.second=mol_pos[s2];
gcon1.second=-1;
gcons.push_back(gcon1);
s1--;s2++;n1++;n2++;
}else{
// -1              1            +1 ||       +1         1          -1 //       
if(abp1==0 && abp2==0 && abp3 == 0 && bbp1==0 && bbp2==1 && bbp3==0 && p_type==0){
gcon1.first.first=mol_pos[s1];
gcon1.first.second=mol_pos[s2-1];
gcon1.second=-1;
gcons.push_back(gcon1);
gcon1.first.first=mol_pos[s1];
gcon1.first.second=mol_pos[s2+1];
gcon1.second=-1;
gcons.push_back(gcon1);
s1--;s2++;n1++;n2++;
}else{

// -1              1            +1 ||       +1         1          -1 //       
if(abp1==0 && abp2==0 && abp3 == 0 && bbp1==0 && bbp2==0 && bbp3==1 && p_type==0){
gcon1.first.first=mol_pos[s1];
gcon1.first.second=mol_pos[s2];
gcon1.second=-1;
gcons.push_back(gcon1);
s1--;s2++;n1++;n2++;
}else{
////// P_TYPE = 1 implies a buldge and constraints are s2, s1-1
// -1              1            +1 ||       +1         1          -1 //       
if(abp1==1 && abp2==0 && abp3 == 0 && bbp1==0 && bbp2==0 && bbp3==0 && p_type==1){
gcon1.first.first=mol_pos[s1-2];
gcon1.first.second=mol_pos[s2];
gcon1.second=-1;
gcons.push_back(gcon1);
gcon1.first.first=mol_pos[s1];
gcon1.first.second=mol_pos[s2];
gcon1.second=-1;
gcons.push_back(gcon1);
s2++;s1--;s1--;n1++;n1++;n2++;
}else{

// -1              1            +1 ||       +1         1          -1 //       
if(abp1==1 && abp2==0 && abp3 == 0 && bbp1==1 && bbp2==0 && bbp3==0 && p_type==1){
gcon1.first.first=mol_pos[s1-2];
gcon1.first.second=mol_pos[s2];
gcon1.second=-10; //stacked
gcons.push_back(gcon1);
s2++;s1--;s1--;n1++;n1++;n2++;
}else{
// -1              1            +1 ||       +1         1          -1 //       
if(abp1==1 && abp2==0 && abp3 == 0 && bbp1==0 && bbp2==1 && bbp3==0 && p_type==1){

gcon1.first.first=mol_pos[s1-2];
gcon1.first.second=mol_pos[s2+1];
gcon1.second=-10;//stacked
gcons.push_back(gcon1);
s2++;s1--;s1--;n1++;n1++;n2++;

}else{

// -1              1            +1 ||       +1         1          -1 //       
if(abp1==1 && abp2==0 && abp3 == 0 && bbp1==0 && bbp2==0 && bbp3==1 && p_type==1){

gcon1.first.first=mol_pos[s1-2];
gcon1.first.second=mol_pos[s2];
gcon1.second=-1;
gcons.push_back(gcon1);
gcon1.first.first=mol_pos[s1];
gcon1.first.second=mol_pos[s2];
gcon1.second=-1;
gcons.push_back(gcon1);
s2++;s1--;s1--;n1++;n1++;n2++;
}else{
// -1              1            +1 ||       +1         1          -1 //       
if(abp1==0 && abp2==1 && abp3 == 0 && bbp1==0 && bbp2==0 && bbp3==0 && p_type==1){
gcon1.first.first=mol_pos[s1-1];
gcon1.first.second=mol_pos[s2];
gcon1.second=-1;
gcons.push_back(gcon1);
s2++;s1--;s1--;n1++;n1++;n2++;
}else{
// -1              1            +1 ||       +1         1          -1 //       
if(abp1==0 && abp2==1 && abp3 == 0 && bbp1==1 && bbp2==0 && bbp3==0 && p_type==1){
gcon1.first.first=mol_pos[s1-1];
gcon1.first.second=mol_pos[s2];
gcon1.second=-1;
gcons.push_back(gcon1);
s2++;s1--;s1--;n1++;n1++;n2++;
}else{

// -1              1            +1 ||       +1         1          -1 //       
if(abp1==0 && abp2==1 && abp3 == 0 && bbp1==0 && bbp2==1 && bbp3==0 && p_type==1){
gcon1.first.first=mol_pos[s1-1];
gcon1.first.second=mol_pos[s2-1];
gcon1.second=-10;
gcons.push_back(gcon1);
s2++;s1--;s1--;n1++;n1++;n2++;
}else{

// -1              1            +1 ||       +1         1          -1 //       
if(abp1==0 && abp2==1 && abp3 == 0 && bbp1==0 && bbp2==0 && bbp3==1 && p_type==1){
gcon1.first.first=mol_pos[s1-1];
gcon1.first.second=mol_pos[s2];
gcon1.second=-1;
gcons.push_back(gcon1);
s2++;s1--;s1--;n1++;n1++;n2++;
}else{

// -1              1            +1 ||       +1         1          -1 //       
if(abp1==0 && abp2==0 && abp3 == 1 && bbp1==0 && bbp2==0 && bbp3==0 && p_type==1){
gcon1.first.first=mol_pos[s1-1];
gcon1.first.second=mol_pos[s2];
gcon1.second=-9;//planar
gcons.push_back(gcon1);
s2++;s1--;s1--;n1++;n1++;n2++;
}else{

// -1              1            +1 ||       +1         1          -1 //       
if(abp1==0 && abp2==0 && abp3 == 1 && bbp1==1 && bbp2==0 && bbp3==0 && p_type==1){
gcon1.first.first=mol_pos[s1-1];
gcon1.first.second=mol_pos[s2];
gcon1.second=-1;
gcons.push_back(gcon1);
gcon1.first.first=mol_pos[s1-1];
gcon1.first.second=mol_pos[s2+2];
gcon1.second=-1;
gcons.push_back(gcon1);
s2++;s1--;s1--;n1++;n1++;n2++;
}else{
// -1              1            +1 ||       +1         1          -1 //       
if(abp1==0 && abp2==0 && abp3 == 1 && bbp1==0 && bbp2==1 && bbp3==0 && p_type==1){
gcon1.first.first=mol_pos[s1-1];
gcon1.first.second=mol_pos[s2-1];
gcon1.second=-1;
gcons.push_back(gcon1);
gcon1.first.first=mol_pos[s1-1];
gcon1.first.second=mol_pos[s2+1];
gcon1.second=-1;
gcons.push_back(gcon1);
s2++;s1--;s1--;n1++;n1++;n2++;
}else{

// -1              1            +1 ||       +1         1          -1 //       
if(abp1==0 && abp2==0 && abp3 == 1 && bbp1==0 && bbp2==0 && bbp3==1 && p_type==1){
gcon1.first.first=mol_pos[s1-1];
gcon1.first.second=mol_pos[s2];
gcon1.second=-1;
gcons.push_back(gcon1);
s2++;s1--;s1--;n1++;n1++;n2++;
}else{

// -1              1            +1 ||       +1         1          -1 //       
if(abp1==0 && abp2==0 && abp3 == 0 && bbp1==0 && bbp2==0 && bbp3==0 && p_type==1){
gcon1.first.first=mol_pos[s1-1];
gcon1.first.second=mol_pos[s2];
gcon1.second=-9;
gcons.push_back(gcon1);
s2++;s1--;s1--;n1++;n1++;n2++;


}else{

// -1              1            +1 ||       +1         1          -1 //       
if(abp1==0 && abp2==0 && abp3 == 0 && bbp1==1 && bbp2==0 && bbp3==0 && p_type==1){
gcon1.first.first=mol_pos[s1-1];
gcon1.first.second=mol_pos[s2];
gcon1.second=-1;
gcons.push_back(gcon1);
s2++;s1--;s1--;n1++;n1++;n2++;

}else{

// -1              1            +1 ||       +1         1          -1 //       
if(abp1==0 && abp2==0 && abp3 == 0 && bbp1==0 && bbp2==1 && bbp3==0 && p_type==1){
gcon1.first.first=mol_pos[s1-1];
gcon1.first.second=mol_pos[s2-1];
gcon1.second=-1;
gcons.push_back(gcon1);
gcon1.first.first=mol_pos[s1-1];
gcon1.first.second=mol_pos[s2+1];
gcon1.second=-1;
gcons.push_back(gcon1);
s2++;s1--;s1--;n1++;n1++;n2++;
}else{

// -1              1            +1 ||       +1         1          -1 //       
if(abp1==0 && abp2==0 && abp3 == 0 && bbp1==0 && bbp2==0 && bbp3==1 && p_type==1){
gcon1.first.first=mol_pos[s1-1];
gcon1.first.second=mol_pos[s2];
gcon1.second=-1;
gcons.push_back(gcon1);
s2++;s1--;s1--;n1++;n1++;n2++;
}else{
//P_TYPE = 2 implies a buldge with s1 - s2+1
if(abp1==1 && abp2==0 && abp3 == 0 && bbp1==0 && bbp2==0 && bbp3==0 && p_type==2){
gcon1.first.first=mol_pos[s1];
gcon1.first.second=mol_pos[s2+1];
gcon1.second=-1;
gcons.push_back(gcon1);
s2++;s2++;s1--;n1++;n2++;n2++;
}else{

if(abp1==1 && abp2==0 && abp3 == 0 && bbp1==1 && bbp2==0 && bbp3==0 && p_type==2){

gcon1.first.first=mol_pos[s1];
gcon1.first.second=mol_pos[s2+2];
gcon1.second=-10;//planar
gcons.push_back(gcon1);
s2++;s2++;s1--;n1++;n2++;n2++;

}else{

if(abp1==1 && abp2==0 && abp3 == 0 && bbp1==0 && bbp2==1 && bbp3==0 && p_type==2){
gcon1.first.first=mol_pos[s1];
gcon1.first.second=mol_pos[s2+1];
gcon1.second=-1;
gcons.push_back(gcon1);
s2++;s2++;s1--;n1++;n2++;n2++;
}else{

if(abp1==1 && abp2==0 && abp3 == 0 && bbp1==0 && bbp2==0 && bbp3==1 && p_type==2){
gcon1.first.first=mol_pos[s1];
gcon1.first.second=mol_pos[s2+1];
gcon1.second=-1;
gcons.push_back(gcon1);
s2++;s2++;s1--;n1++;n2++;n2++;
}else{

if(abp1==0 && abp2==1 && abp3 == 0 && bbp1==0 && bbp2==0 && bbp3==0 && p_type==2){
gcon1.first.first=mol_pos[s1-1];
gcon1.first.second=mol_pos[s2+1];
gcon1.second=-1;
gcons.push_back(gcon1);
gcon1.first.first=mol_pos[s1+1];
gcon1.first.second=mol_pos[s2+1];
gcon1.second=-1;
gcons.push_back(gcon1);
s2++;s2++;s1--;n1++;n2++;n2++;
}else{

if(abp1==0 && abp2==1 && abp3 == 0 && bbp1==1 && bbp2==0 && bbp3==0 && p_type==2){
gcon1.first.first=mol_pos[s1-1];
gcon1.first.second=mol_pos[s2];
gcon1.second=-10;
gcons.push_back(gcon1);
gcon1.first.first=mol_pos[s1+1];
gcon1.first.second=mol_pos[s2+2];
gcon1.second=-10;
gcons.push_back(gcon1);
s2++;s2++;s1--;n1++;n2++;n2++;
}else{

if(abp1==0 && abp2==1 && abp3 == 0 && bbp1==0 && bbp2==1 && bbp3==0 && p_type==2){
gcon1.first.first=mol_pos[s1+1];
gcon1.first.second=mol_pos[s2+1];
gcon1.second=-10;
gcons.push_back(gcon1);
s2++;s2++;s1--;n1++;n2++;n2++;
}else{

if(abp1==0 && abp2==1 && abp3 == 0 && bbp1==0 && bbp2==0 && bbp3==1 && p_type==2){
gcon1.first.first=mol_pos[s1+1];
gcon1.first.second=mol_pos[s2+1];
gcon1.second=-1;
gcons.push_back(gcon1);
gcon1.first.first=mol_pos[s1-1];
gcon1.first.second=mol_pos[s2+1];
gcon1.second=-1;
gcons.push_back(gcon1);
s2++;s2++;s1--;n1++;n2++;n2++;
}else{

if(abp1==0 && abp2==0 && abp3 == 1 && bbp1==0 && bbp2==0 && bbp3==0 && p_type==2){
gcon1.first.first=mol_pos[s1];
gcon1.first.second=mol_pos[s2+1];
gcon1.second=-1;
gcons.push_back(gcon1);
s2++;s2++;s1--;n1++;n2++;n2++;
}else{

if(abp1==0 && abp2==0 && abp3 == 1 && bbp1==1 && bbp2==0 && bbp3==0 && p_type==2){
gcon1.first.first=mol_pos[s1];
gcon1.first.second=mol_pos[s2];
gcon1.second=-1;
gcons.push_back(gcon1);
gcon1.first.first=mol_pos[s1];
gcon1.first.second=mol_pos[s2+2];
gcon1.second=-1;
gcons.push_back(gcon1);
s2++;s2++;s1--;n1++;n2++;n2++;
}else{

if(abp1==0 && abp2==0 && abp3 == 1 && bbp1==0 && bbp2==1 && bbp3==0 && p_type==2){
gcon1.first.first=mol_pos[s1];
gcon1.first.second=mol_pos[s2+1];
gcon1.second=-1;
gcons.push_back(gcon1);
s2++;s2++;s1--;n1++;n2++;n2++;
}else{

if(abp1==0 && abp2==0 && abp3 == 1 && bbp1==0 && bbp2==0 && bbp3==1 && p_type==2){
gcon1.first.first=mol_pos[s1];
gcon1.first.second=mol_pos[s2+1];
gcon1.second=-1;
gcons.push_back(gcon1);
s2++;s2++;s1--;n1++;n2++;n2++;
}else{

if(abp1==0 && abp2==0 && abp3 == 0 && bbp1==0 && bbp2==0 && bbp3==0 && p_type==2){
gcon1.first.first=mol_pos[s1];
gcon1.first.second=mol_pos[s2+1];
gcon1.second=-9;
gcons.push_back(gcon1);
s2++;s2++;s1--;n1++;n2++;n2++;
}else{

if(abp1==0 && abp2==0 && abp3 == 0 && bbp1==1 && bbp2==0 && bbp3==0 && p_type==2){
gcon1.first.first=mol_pos[s1];
gcon1.first.second=mol_pos[s2];
gcon1.second=-1;
gcons.push_back(gcon1);
gcon1.first.first=mol_pos[s1];
gcon1.first.second=mol_pos[s2+2];
gcon1.second=-1;
gcons.push_back(gcon1);
s2++;s2++;s1--;n1++;n2++;n2++;
}else{

if(abp1==0 && abp2==0 && abp3 == 0 && bbp1==0 && bbp2==1 && bbp3==0 && p_type==2){
gcon1.first.first=mol_pos[s1];
gcon1.first.second=mol_pos[s2+1];
gcon1.second=-1;
gcons.push_back(gcon1);
s2++;s2++;s1--;n1++;n2++;n2++;
}else{

if(abp1==0 && abp2==0 && abp3 == 0 && bbp1==0 && bbp2==0 && bbp3==1 && p_type==2){
gcon1.first.first=mol_pos[s1];
gcon1.first.second=mol_pos[s2+1];
gcon1.second=-9;
gcons.push_back(gcon1);
s2++;s2++;s1--;n1++;n2++;n2++;
}else{
printf("Bad Constraints\n");
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
if(gcon1.first.first==-1 || gcon1.first.second==-1){
printf("%d %d %d, %d-%d,%d %d %d %d %d %d %d\n",mol_pos.size(),s1,s2,gcon1.first.first,gcon1.first.second,abp1,abp2,abp3,bbp1,bbp2,bbp3,p_type);
}
gencons.push_back(gcons);
}
}
//check to make sure gencons contains all bases pairings in the 
//hard constraints

int metcons=0;
for(int ss1=0; ss1 < mol_pos2.size(); ss1++){
int sc1=mol_pos2[ss1];

if(sc1>=0){
int sc2=sec_struct[sc1];
metcons=0;
if(sc2 >=0 && sc2 > sc1){
for(int gc1=0; gc1 < gencons.size(); gc1++){
for(int gc2=0; gc2 < gencons[gc1].size(); gc2++){
//printf("Con Check (%d %d) (%d %d)\n",ss1,sc2,gencons[gc1][gc2].first,gencons[gc1][gc2].second);
if(gencons[gc1][gc2].first.first==sc1 && gencons[gc1][gc2].first.second==sc2){
metcons=1;
break;
}
}
if(metcons==1){break;}
}
if(metcons==0){break;}
}else{
metcons=1;
}
}
}

if(metcons==1){
//Foreach of the current movepositions push_back the new constraints
if(MasterMovePos.size() > 0){
for(int m=0; m < MasterMovePos.size(); m++){
//Create a new entry that adds to the existing vectors
std::vector< std::vector<int> > tmov=MasterMovePos[m];
for(int mnew=0; mnew < movepos.size(); mnew++){tmov.push_back(movepos[mnew]);}
std::vector< std::vector< std::pair< std::pair<int,int>,int > > > tcon=MasterGenCons[m];
for(int mnew=0; mnew < gencons.size(); mnew++){tcon.push_back(gencons[mnew]);}
moveTemp.push_back(tmov);
consTemp.push_back(tcon);
}

}else{
std::vector< std::vector<int> > tmov;

for(int mnew=0; mnew < movepos.size(); mnew++){tmov.push_back(movepos[mnew]);}
std::vector<std::vector< std::pair<std::pair<int,int>,int> > >  tcon;

for(int mnew=0; mnew < gencons.size(); mnew++){tcon.push_back(gencons[mnew]);}
moveTemp.push_back(tmov);
consTemp.push_back(tcon);
}

///Printing out results
////////////////////////////////////////////////////////////
printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t",s1,s2,side1,side2,lp_p,lp_ps1,lp_ps2,lp_typ,minPos,maxPos,metcons);
for(int g=0; g<gencons.size(); g++){
printf("[");
for(int m1=0; m1<movepos[g].size(); m1++){printf("%d ",movepos[g][m1]);}
printf("]-");
for(int g2=0; g2<gencons[g].size(); g2++){
printf("(%d-%d) ",gencons[g][g2].first.first,gencons[g][g2].first.second);}
}
printf("\n");
}
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////




}

}



}
}


}
 
}
//write Temp's to Master

if(moveTemp.size() > 0){
//Filter for all  bp's on sides - combine with next constraint 
//Identify 
//Filter all constraints for bp's internal to loops - replace with cluster cons
for(int gf1=0; gf1<moveTemp.size();gf1++){
//Foreach set of movement positions
for(int gf2=0; gf2<moveTemp[gf1].size(); gf2++){
//Identify blocks that need to be replaced
int bp1=-1;
int bp2=-1;
int gf3=0;

while(gf3<moveTemp[gf1][gf2].size()){

if(bp1==-1){
if(moveTemp[gf1][gf2][gf3]==-99){
//if position = -99 condense the constraints before and after this one
bp1=moveTemp[gf1][gf2][gf3+1];
bp2=moveTemp[gf1][gf2][gf3-1];
gf3=0;


}
}else{



}

gf3++;
}
}
}


}


if(moveTemp.size() > 0){
MasterMovePos.clear();
MasterGenCons.clear();
MasterClash.clear();

//Identify the constraint sets with the minimal number of buldges
int cur_min=0;
for(int mn = 0; mn < consTemp.size(); mn++){
int new_min=0;
for(int p = 0; p < consTemp[mn].size(); p++){
new_min=new_min+consTemp[mn][p].size();
}
if(cur_min==0){
cur_min=new_min;
}else{
if(new_min<cur_min){
cur_min=new_min;
}
}
}

//Keep minimal sized constraint set
for(int g = 0; g< moveTemp.size(); g++){
int keep=1;
int consize=0;
for(int p = 0; p < consTemp[g].size(); p++){
consize=consize+consTemp[g][p].size();
}

if(consize > cur_min+2){keep=0;}
if(keep==1){
MasterMovePos.push_back(moveTemp[g]);
MasterGenCons.push_back(consTemp[g]);
MasterClash.push_back(cshTemp[g]);
}
}
consTemp.clear();
moveTemp.clear();
cshTemp.clear();
}

}



string l_ofile(opt->prefix.c_str());
l_ofile.append("_csh.txt");
string c_ofile(opt->prefix.c_str());
c_ofile.append("_con.txt");
string m_ofile(opt->prefix.c_str());
m_ofile.append("_pos.txt");
FILE * c_oFH = fopen(c_ofile.c_str(),"a");
FILE * m_oFH = fopen(m_ofile.c_str(),"a");
FILE * l_oFH = fopen(l_ofile.c_str(),"a");

for(int id=0; id<MasterClash.size(); id++){
for(int step=0; step<MasterClash[id].size(); step++){
fprintf(l_oFH,"%d %d",id,step);
for(int p=0; p<MasterClash[id][step].size(); p++){
fprintf(l_oFH," %d",MasterClash[id][step][p]);
}
fprintf(l_oFH,"\n");
}
}

for(int id=0; id<MasterMovePos.size(); id++){
for(int step=0; step<MasterMovePos[id].size(); step++){
for(int p=0; p<MasterMovePos[id][step].size(); p++){
fprintf(m_oFH,"%d %d %d\n",id,step,MasterMovePos[id][step][p]);
}
}
}

int d1,d2,d3,d4,d5;
for(int id=0; id<MasterGenCons.size(); id++){
for(int step=0; step<MasterGenCons[id].size(); step++){
for(int p=0; p<MasterGenCons[id][step].size(); p++){
if(step>0 && d4 > 0){
fprintf(c_oFH,"%d %d %d %d %d\n",id,step,d3,d4,d5);
}
fprintf(c_oFH,"%d %d %d %d %d\n",id,step,MasterGenCons[id][step][p].first.first,MasterGenCons[id][step][p].first.second,MasterGenCons[id][step][p].second);
d1=id;
d2=step;
d3=MasterGenCons[id][step][p].first.first;
d4=MasterGenCons[id][step][p].first.second;
d5=MasterGenCons[id][step][p].second;
}
}
}

fflush(m_oFH);
fflush(c_oFH);
fflush(l_oFH);

fclose(m_oFH);
fclose(c_oFH);
fclose(l_oFH);

}//end main
