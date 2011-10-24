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
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <algorithm>
#include <ctime>
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

int main(int argc, char **argv)
{
using namespace std;

Options * opt=new Options(argc,argv);
//Set random number
//srand(time(0));
srand(opt->seed);
//Initial starting conformation
float p_o_upper=opt->p_o_upper;
float p_o_lower=opt->p_o_lower;
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////// Loading fragments and simulation parameters /////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
printf("Assemble Structure from general constraints\n");
Fragclass * dat = new Fragclass(opt->db_prefix);
//Wrap this into a function 
printf("Reading in doublets\n");
std::vector<double> dparams;
loadDoublets(&dparams,opt->doublet_params);
printf("Number of Doublet Centers: %d\n",dparams.size());
//Reading in fragments with base atoms included 
vector<Fragatoms *> atoms;
if(opt->atoms_file.c_str()==NULL){
printf("Not using full atom fragments\n");
}else{
std::ifstream input;
input.open(opt->atoms_file.c_str());
string line;
while(getline(input, line)) {atoms.push_back(new Fragatoms(line));}
printf("Full atom fragments loaded from: %s",opt->atoms_file.c_str());
input.close();
}

//Loading all atom fragments
vector<Fragall *> batoms;
//If all atom file exists load it into the simulation

std::ifstream inputb;
inputb.open(opt->batoms_file.c_str());
if(inputb.fail()){
//printf("Could not open full atom fragment file: %s\n",opt->batoms_file);
abort();
}else{
string line;
while(getline(inputb, line)) {batoms.push_back(new Fragall(line));}
cout << batoms.size() << " full atom fragments loaded from: " << opt->batoms_file.c_str()<<  endl;
}
inputb.close();
//batoms.back()->print();
//printf("Size of vector: %d\n",batoms.back()->atoms.size());
//Reading in fragments containing only virtual atoms
vector<Fragment *> frags;
frags.reserve(42000);
//Setting up some test molecules for our polymere 
std::ifstream input;
input.open(opt->fragment_file.c_str());
string line;
while(getline(input, line)){frags.push_back(new Fragment(line));}
input.close();
cout << "Fragments Loaded: " << frags.size() << endl;
std::vector< std::vector<int> > loop_listO;
//Read in loop database
std::ifstream loops_fileO;
string loops_odd=opt->loop_db;
loops_odd.append("O"); //Append E for even library
loops_fileO.open(loops_odd.c_str());
string line_frags;
printf("Reading in odd loops: ");
while(getline(loops_fileO,line_frags)){
string field;
string seq_type;
int frag1;
int frag2;
int frag3;
int frag4;
	istringstream iss(line_frags,istringstream::in);	
   		getline(iss,field,' ');
        seq_type=field;
                getline(iss,field,' ');
        frag1=strtoint(field);
                getline(iss,field,' ');
        frag2=strtoint(field);
                getline(iss,field,' ');
        frag3=strtoint(field);
                getline(iss,field,' ');
        frag4=strtoint(field);
std::vector<int> loop_in;
loop_in.push_back(frag1);
loop_in.push_back(frag2);
loop_in.push_back(frag3);
loop_in.push_back(frag4);
loop_listO.push_back(loop_in);
}
loops_fileO.close();
printf("%d\n",loop_listO.size());

std::vector< std::vector<int> > loop_listE;
//Read in loop database
std::ifstream loops_file;
string loops_even=opt->loop_db;
loops_even.append("E"); //Append E for even library
loops_file.open(loops_even.c_str());
printf("Reading in even loops: ");
while(getline(loops_file,line_frags)){
string field;
string seq_type;
int frag1;
int frag2;
int frag3;
int frag4;
	istringstream iss(line_frags,istringstream::in);	
   		getline(iss,field,' ');
        seq_type=field;
                getline(iss,field,' ');
        frag1=strtoint(field);
                getline(iss,field,' ');
        frag2=strtoint(field);
                getline(iss,field,' ');
        frag3=strtoint(field);
                getline(iss,field,' ');
        frag4=strtoint(field);
std::vector<int> loop_in;
loop_in.push_back(frag1);
loop_in.push_back(frag2);
loop_in.push_back(frag3);
loop_in.push_back(frag4);
loop_listE.push_back(loop_in);
}
printf("%d\n",loop_listE.size());
loops_file.close();

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////  Reading in general constraints ///////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

//This program assumes you are starting from scratch
// movepos - stores vectors of positions that should be moved with the fragment replacement 
// alogorithm at each step 
// gencons - stores vectors of contraints that should be met at each step
// constraint - a constraint is a vector of pairs of position pairs and a constraint type integer

//contype 0 = Any pairing type
//contype 1 = base paired W-C edge
//contype 2 = base paired Hoogstein edge 
//contype 3 = base paired Sugar edge
//contype 4 = stacked 1 above 2 
//contype 5 = stacked 2 above 1
 
std::vector<std::vector<int> > movepos;
std::vector< std::vector<std::pair<std::pair<int,int>,int> > >gencons;

// Read in a particular constraint set from the file //
std::vector<std::vector<int> > clashgroups;
int con_id = opt->constraint_id;
string ofile=opt->prefix;
ofile+=boost::lexical_cast<std::string>(con_id);
ofile.append(".pdb");

string ofile2=opt->prefix;
ofile2+=boost::lexical_cast<std::string>(con_id);
ofile2.append("f.pdb");

string ofile3=opt->prefix;
ofile3+=boost::lexical_cast<std::string>(con_id);
ofile3.append("_steps_");

std::ifstream clash_infile;
string csh_file=opt->input_file;
csh_file.append("_csh.txt");
printf("Reading in clash groups from %s:",csh_file.c_str());
clash_infile.open(csh_file.c_str());
string line_csh;
std::vector<int> tmpcsh;
int rna_size=opt->sequence.size();
int csh_group=0;
while(getline(clash_infile,line_csh)){
string field;
int id;
int step;
std::vector<int> clashGroup;
	istringstream iss(line_csh,istringstream::in);	
   		getline(iss,field,' ');
        id=strtoint(field);
                getline(iss,field,' ');
        step=strtoint(field);
//Read in each clash group
if(id==con_id){
clashGroup.push_back(0);
clashGroup.push_back(0);
clashGroup.push_back(0);
clashGroup.push_back(0);
for(int p=0; p<rna_size; p++){
    int csh;
    getline(iss,field,' ');
    csh=strtoint(field);
    clashGroup.push_back(csh);
}
clashGroup.push_back(0);
clashGroup.push_back(0);
clashGroup.push_back(0);
clashGroup.push_back(0);
clashgroups.push_back(clashGroup);
}
}
printf(" %d\n",clashgroups.size());




std::ifstream m_infile;
string pos_file=opt->input_file;
pos_file.append("_pos.txt");
printf("Reading in movement groups from %s:",pos_file.c_str());
m_infile.open(pos_file.c_str());
string line_pos;
std::vector<int> tmppos;
int start=-1;
while(getline(m_infile,line_pos)){
string field;
int id;
int step;
int pos;
	istringstream iss(line_pos,istringstream::in);	
   		getline(iss,field,' ');
        id=strtoint(field);
                getline(iss,field,' ');
        step=strtoint(field);
                getline(iss,field,' ');
        pos=strtoint(field);
//Add additional NTs

pos=pos+4;

if(id==con_id){
if(start!=step){
if(start!=-1){movepos.push_back(tmppos);}
start=step;
tmppos.clear();
}
tmppos.push_back(pos);
}
}
movepos.push_back(tmppos);
m_infile.close();
printf(" %d\n",movepos.size());
std::ifstream c_infile;
string con_file=opt->input_file;
con_file.append("_con.txt"); 
printf("Reading in general constraints from %s:",con_file.c_str());
c_infile.open(con_file.c_str());
string line_con;
std::vector<std::pair<std::pair<int,int>,int> > tmpcon;
 start=-1;
while(getline(c_infile,line_con)){
string field;
int id;
int step;
int pos1;
int pos2;
int type;
	istringstream iss(line_con,istringstream::in);	
   		getline(iss,field,' ');
        id=strtoint(field);
                getline(iss,field,' ');
        step=strtoint(field);
                getline(iss,field,' ');
        pos1=strtoint(field);
                getline(iss,field,' ');
        pos2=strtoint(field);
                getline(iss,field,' ');
        type=strtoint(field);
//Add additional NTs
if(pos1 >= 0){pos1=pos1+4;}
if(pos2 >= 0){pos2=pos2+4;}

if(id==con_id){
if(start!=step){
if(start!=-1){gencons.push_back(tmpcon);}
start=step;
tmpcon.clear();
}
tmpcon.push_back(std::pair<std::pair<int,int>,int>(std::pair<int,int>(pos1,pos2),type));
}
}
gencons.push_back(tmpcon);
printf(" %d\n",gencons.size());
c_infile.close();

//Build the tertiary structure using the general constraints
for(int tot_out=0; tot_out < 50; tot_out++){

Polymere * outpol = new Polymere(opt->sequence,frags,dat,atoms,batoms,&dparams);
int start_cons=0;
if(opt->starting_struct.compare("")!=0){
printf("Using initial conditions from %s\n",opt->starting_struct.c_str());
Polymere * init_pol = new Polymere(opt->starting_struct,0);
outpol->copyOver(init_pol,0,0,init_pol->num_mols);
start_cons=opt->con_start;
}else{
printf("No initial pdb given: %s\n",opt->starting_struct.c_str());
}
//Print out read in constraints to make sure they are correct


for(int m=0; m<movepos.size(); m++){
printf("step: %d\n",m);
printf(" nucs: ");
for(int n=0; n< outpol->num_mols && n < clashgroups[m].size(); n++){printf("%c",outpol->mols[n]->seq);}
printf("\n");
printf(" cshs: ");
for(int n=0; n< outpol->num_mols && n < clashgroups[m].size(); n++){printf("%d",clashgroups[m][n]);}
printf("\n");
printf(" cons: ");
for(int n=0; n< outpol->num_mols && n < clashgroups[m].size(); n++){
int g1=0;
for(int gc=0; gc < gencons[m].size(); gc++){
if(n==gencons[m][gc].first.second || n==gencons[m][gc].first.first){
printf("c");
g1=1;
break;
}
}
if(g1==0){printf(".");}
}
printf("\n");
printf(" fpos: ");
for(int n=0; n< outpol->num_mols && clashgroups[m].size(); n++){
int f1=0;
for(int m1=0; m1<movepos[m].size(); m1++){
if(n==movepos[m][m1]){
printf("f");
f1=1;
break;
}
}
if(f1==0){
printf(".");
}
}
printf("\n");


for(int p=0; p<movepos[m].size(); p++){
printf("%d ",movepos[m][p]);
}
for(int p=0; p<gencons[m].size(); p++){
printf("(%d %d)-%d",gencons[m][p].first.first,gencons[m][p].first.second,gencons[m][p].second);
}
printf("\n");
}

//Generate a new polymer from the given RNA sequence
for(int c=0; c<clashgroups.size(); c++){
    for(int w=0; w<clashgroups[c].size(); w++){
         printf("-%d",clashgroups[c][w]);
    }
printf("\n");
}

//Skip a certain number of constraints
int wpdb=outpol->assemble(movepos,gencons,clashgroups,frags,dat,atoms,batoms,loop_listO,loop_listE,opt->iteration_num,opt->max_leaves,5,ofile3,opt->write_steps,opt->resolution,opt->max_clust, start_cons,p_o_lower,p_o_upper);
//int wpdb=outpol->consBuild2(movepos,gencons,frags,dat,atoms,batoms,loop_listO,loop_listE,start_cons,opt->write_steps);

if(wpdb==1){
//write pdb file
Pdb * pdb = new Pdb(outpol,'J');
pdb->write(ofile.c_str());
pdb->clear();
delete pdb;
}else{
Pdb * pdb = new Pdb(outpol,'J');
pdb->write(ofile2.c_str());
pdb->clear();
delete pdb;
}

//remove polymere 
for(int i = 0; i<outpol->mols.size(); i++){
delete outpol->mols[i];
}
outpol->mols.clear();
delete outpol;

}
//this->midPointBuild(movegroups,cons,frags,dat,atoms,batoms,dparams);
//Remove leading and lagging sequences
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
}//end main
