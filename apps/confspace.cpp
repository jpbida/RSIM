#include <algorithm>
#include <boost/lexical_cast.hpp>
#include <ctime>
#include <fstream>
#include <iostream>
#include "base.h"
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

int main(int argc, char **argv)
{
using namespace std;
Options * opt=new Options(argc,argv);
srand(opt->seed);
//opt->print();

//Reading in Fragments and Scoring Databases
//dat should already load fragments
//cout << "Reading in Fragment Class" << endl;
Fragclass * dat = new Fragclass(opt->db_prefix);
cout << "Reading in master structure" << endl;

Polymere * pol = new Polymere(opt->starting_struct,0);

Polymere * poltmp=new Polymere(opt->master_struct);
pol->copyOver(poltmp,4,0,poltmp->num_mols);
//Reading in doublet parameters
//Wrap this into a function 
cout << "Reading in doublets" << endl;
vector<double> dparams;
loadDoublets(&dparams,opt->doublet_params);
cout << "Number of Doublet Centers: " << dparams.size() << endl;
//Reading in fragments with base atoms included 
vector<Fragatoms *> atoms;
if(opt->atoms_file.c_str()==NULL){
cout << "Not using full atom fragments" << endl;
}else{
std::ifstream input;
input.open(opt->atoms_file.c_str());
string line;
while(getline(input, line)) {atoms.push_back(new Fragatoms(line));}
cout << atoms.size() << " full atom fragments loaded from: " << opt->atoms_file.c_str()<<  endl;
input.close();
}

//Loading all atom fragments
vector<Fragall *> batoms;

//If all atom file exists load it into the simulation

std::ifstream inputb;
inputb.open(opt->batoms_file.c_str());
if(inputb.fail()){
printf("Could not open full atom fragment file: %s\n",opt->batoms_file.c_str());
abort();
}else{
string line;
while(getline(inputb, line)) {batoms.push_back(new Fragall(line));}
cout << batoms.size() << " full atom fragments loaded from: " << opt->batoms_file.c_str()<<  endl;
}
inputb.close();
//Check all atom fragments
vector<Fragment *> frags;
frags.reserve(42000);
//Setting up some test molecules for our polymere 
std::ifstream input;
input.open(opt->fragment_file.c_str());
string line;
while(getline(input, line)){frags.push_back(new Fragment(line));}
input.close();
cout << "Fragments Loaded: " << frags.size() << endl;
//Loop through all fragments
/*
for(int i=0; i<1000; i++){
cen->changeFrag2(4,i,frags);
cen->coordSysNt(4,atoms[i]->atoms);
cen->updateFull();
cen->coordSysAT(4,batoms[i]->atoms);
Pdb * pdbn = new Pdb(cen,'j');
pdbn->print();
}
*/
//Randomizes a structure 
int nmax1=opt->randomize_num;
int imax=opt->iteration_num;
int sim_nums = opt->total_simulations;
//cen->dynamix(nmax1,opt->movegroup,frags,dat,pol,atoms,batoms,0,0);
//Pdb * pdb1 = new Pdb(cen,'J');
//pdb1->write(opt->out_struct);

//cen->dynamix(nmax1,opt->movegroup,frags,dat,pol,atoms,batoms,0,18);
//cen->mixit(nmax1,opt->movegroup,frags,dat,pol,atoms,batoms);
double rmsd_res=2.5;
vector<int> nodes;
vector<int> parent;
vector<int> child;

int out =0;
int pos=0;
//Max structures
Pdb * pdb = new Pdb(opt->starting_struct);
int tots=pdb->models.size();
printf("Total Number of Starting Structures: %d\n",tots);
pdb->clear();
delete pdb;
////////////////

int t=1;
std::string rnaseq;
for(int m=4; m<(pol->num_mols-4); m++){
rnaseq+=boost::lexical_cast<std::string>(pol->mols[m]->seq);
t=2;
if(m<2){t=1;}
if(m > (pol->num_mols-3)){t=1;}
pol->mols[m]->group = t; 
pol->mols[m]->id =opt->sec_structure[m];
}

Polymere * cen = new Polymere(opt->starting_struct,0);

std::vector<Polymere *> centers;
printf("Total Structures Loaded: %d\n",tots);
for(int m=0; m<tots; m++){
printf("loaded %d\n",m);
centers.push_back(new Polymere(opt->starting_struct,m));
nodes.push_back(0);
int t=1;
for(int m=0; m<centers.back()->num_mols; m++){
t=2;
if(m<2){t=1;}
if(m > (centers.back()->num_mols-3)){t=1;}
centers.back()->mols[m]->group = t; 
centers.back()->mols[m]->id =opt->sec_structure[m];
}
}
//Update constraints
t=1;
for(int m=0; m<cen->num_mols; m++){
t=2;
if(m<2){t=1;}
if(m > (cen->num_mols-3)){t=1;}
cen->mols[m]->group = t; 
cen->mols[m]->id =opt->sec_structure[m];
//pol->mols[m]->group = t; 
//pol->mols[m]->id =opt->sec_structure[m];
}

int aclash=cen->groupCheckClashes(1);
int bclash=cen->groupCheckBackBone(1);
printf("Initial Clashes: %d %d\n",aclash,bclash);

printf("RNA Sequnece: %s\n",rnaseq.c_str());
//Generate all cons
std::vector< std::vector< std::pair<int,int> > > allcons;
cen->genAllCons(nmax1,&allcons,rnaseq);
printf("Total Number of Constraints: %zu\n",allcons.size());
/*
printf("Cons: 60\n");
for(int g=0; g<allcons[60].size(); g++){printf("%d-%d\n",allcons[60][g].first,allcons[60][g].second);}
printf("Cons: 57\n");
for(int g=0; g<allcons[57].size(); g++){printf("%d-%d\n",allcons[57][g].first,allcons[57][g].second);}
printf("Cons: 56\n");
for(int g=0; g<allcons[56].size(); g++){printf("%d-%d\n",allcons[56][g].first,allcons[56][g].second);}
printf("Cons: 27\n");
for(int g=0; g<allcons[27].size(); g++){printf("%d-%d\n",allcons[27][g].first,allcons[27][g].second);}
*/
//Determine which constraint class each structure is in

int curmin=0;
int all=0;
//Do a monte-carlo simulation on all the starting structures 
while(pos < tots){
//Change cen to another one of the centers
cen->copy(centers[pos]);
//Update cen's groups so that everything except first and last for are in group
for(int m=0; m<cen->num_mols; m++){
if(m>3 && m < (cen->num_mols-4)){
cen->mols[m]->group=1;}else{
cen->mols[m]->group=0;
}
}
printf("\n\n\nUsing Center: %d\n\n\n",pos);
//HydroMaximixer
cen->score=99999;
int stype=3;
for(int h=0; h<sim_nums; h++){
cen->monteCons(imax,1,frags,dat,pol,atoms,batoms,&dparams,stype,allcons,pos,opt->prefix);
//Calculate Other Scores
//write to overall structures file
Pdb * pdb = new Pdb(cen,'J');
pdb->write("trna_test.pdb");
pdb->clear();
delete pdb;
cen->copy(centers[pos]);
for(int m=0; m<cen->num_mols; m++){
if(m>3 && m < (cen->num_mols-4)){
cen->mols[m]->group=1;}else{
cen->mols[m]->group=0;
}
}
}
pos++;
////////////////
}

//Clean up all the stuff you created 

for(int i=0; i<frags.size(); i++){
delete frags[i];
}
frags.clear();
//Free atom memory 
for(int i=0; i<atoms.size(); i++){
delete atoms[i];
}
atoms.clear();

dat->RRR1.clear();
dat->RRY1.clear();
dat->RYR1.clear();
dat->RYY1.clear();
dat->YRR1.clear();
dat->YRY1.clear();
dat->YYR1.clear();
dat->YYY1.clear();
dat->RRR2.clear();
dat->RRY2.clear();
dat->RYR2.clear();
dat->RYY2.clear();
dat->YRR2.clear();
dat->YRY2.clear();
dat->YYR2.clear();
dat->YYY2.clear();

dat->iRRR1.clear();
dat->iRRY1.clear();
dat->iRYR1.clear();
dat->iRYY1.clear();
dat->iYRR1.clear();
dat->iYRY1.clear();
dat->iYYR1.clear();
dat->iYYY1.clear();
dat->iRRR2.clear();
dat->iRRY2.clear();
dat->iRYR2.clear();
dat->iRYY2.clear();
dat->iYRR2.clear();
dat->iYRY2.clear();
dat->iYYR2.clear();
dat->iYYY2.clear();

for(int c=0; c<centers.size(); c++){
for(int m=0; m<centers[c]->mols.size(); m++){
delete centers[c]->mols[m];
}
centers[c]->full_atoms.clear();
centers[c]->mols.clear();
}
centers.clear();
for(int i = 0; i<pol->mols.size(); i++){
delete pol->mols[i];
}
pol->mols.clear();
delete pol;
pol->mols.clear();
}//end main
