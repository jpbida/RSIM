// Montecarlo fractal kinetics simulator //
#include <boost/lexical_cast.hpp>
#include "polymere.h"
#include <string>
#include <istream>
#include <sstream>
#include <vector>
#include <iostream>
#include "parse.h"
#include "model.h"
#include "residue.h"
#include "atom.h"
#include "pdb.h"

string Pdb::getElement(string line,int type){
/*
---------------------------------------------------------------------------
Field |    Column    | FORTRAN |
  No. |     range    | format  | Description
---------------------------------------------------------------------------
   1. |    1 -  6    |   A6    | Record ID (eg ATOM, HETATM)
   2. |    7 - 11    |   I5    | Atom serial number
   -  |   12 - 12    |   1X    | Blank
   3. |   13 - 16    |   A4    | Atom name (eg " CA " , " ND1")
   4. |   17 - 17    |   A1    | Alternative location code (if any)
   5. |   18 - 20    |   A3    | Standard 3-letter amino acid code for residue
   -  |   21 - 21    |   1X    | Blank
   6. |   22 - 22    |   A1    | Chain identifier code
   7. |   23 - 26    |   I4    | Residue sequence number
   8. |   27 - 27    |   A1    | Insertion code (if any)
   -  |   28 - 30    |   3X    | Blank
   9. |   31 - 38    |  F8.3   | Atom's x-coordinate
  10. |   39 - 46    |  F8.3   | Atom's y-coordinate
  11. |   47 - 54    |  F8.3   | Atom's z-coordinate
  12. |   55 - 60    |  F6.2   | Occupancy value for atom
  13. |   61 - 66    |  F6.2   | B-value (thermal factor)
   -  |   67 - 67    |   1X    | Blank
  14. |   68 - 70    |   I3    | Footnote number
---------------------------------------------------------------------------
*/

int start[] = {0,6,11,12,16,17,20,21,22,26,27,30,46,54,60,66,67};
int end[] = {6,5,1,4,1,3,1,1,4,1,3,8,8,8,6,6,1,3};
string out=line.substr(start[type],end[type]);
return out;
}

void Pdb::ModelsToPols(std::vector< Polymere *> * pols){
//Conversts a PDB file with lots of models to a vector of polymeres
for(int model=0; model<models.size(); model++){
int a=0;
Polymere * pol=new Polymere();
int r = pol->mols.size();
while(r < models[model]->chains[0]->residues.size()){
pol->mols.push_back(new Base());
r++;
} 
//Remove all molecules
for(int resnums=0;resnums < models[model]->chains[0]->residues.size(); resnums++){
//Create a new resiude
pol->mols[resnums]-> z    = 0;
pol->mols[resnums]-> x     =0;
pol->mols[resnums]-> y     =0;
pol->mols[resnums]-> group =0;
pol->mols[resnums]-> prob = 0;
pol->mols[resnums]-> seq   ='N';
pol->mols[resnums]-> bx    =0;
pol->mols[resnums]-> by    =0;
pol->mols[resnums]-> bz    =0;
pol->mols[resnums]-> b2x   =0;
pol->mols[resnums]-> b2y   =0;
pol->mols[resnums]-> b2z   =0;
pol->mols[resnums]-> id    =0;
string rname=models[model]->chains[0]->residues[resnums]->res_name;
if(rname.compare("  A")==0 || rname.compare("  G")==0){
vector<float> c1 = models[model]->chains[0]->residues[resnums]->getAtom(" C1'"); 
vector<float> c2 = models[model]->chains[0]->residues[resnums]->getAtom(" N9 "); 
vector<float> c3 = models[model]->chains[0]->residues[resnums]->getAtom(" C4 "); 
//Add coordinates 
pol->mols[resnums]->x=c1[0];
pol->mols[resnums]->y=c1[1];
pol->mols[resnums]->z=c1[2];

pol->mols[resnums]->bx=c2[0];
pol->mols[resnums]->by=c2[1];
pol->mols[resnums]->bz=c2[2];

pol->mols[resnums]->b2x=c3[0];
pol->mols[resnums]->b2y=c3[1];
pol->mols[resnums]->b2z=c3[2];
if(rname.compare("  A")==0){
pol->mols[resnums]->seq='A';
int bpos = a*81;
int ant_i=0;
if(a == 0){ant_i=3;}
for(int atn = ant_i; atn < 26; atn++){
string at_name = models[model]->chains[0]->residues[resnums]->atomnames2[atn];
c1 = models[model]->chains[0]->residues[resnums]->getAtom(at_name);
if(c1.size()==0){
//printf("\nMissing %s in residue %d %s\n",at_name.c_str(),models[model]->chains[0]->residues[resnums]->res_num,models[model]->chains[0]->residues[resnums]->res_name.c_str()); 
}
pol->full_atoms[bpos+ 0]=c1[0];
pol->full_atoms[bpos+ 1]=c1[1];
pol->full_atoms[bpos+ 2]=c1[2];
bpos = bpos+3;
}
}else{
if(rname.compare("  G")==0){
int ant_i=0;
if(a == 0){ant_i=3;}
int bpos = a*81;
pol->mols[resnums]->seq='G';
for(int atn = ant_i; atn < 27; atn++){
if(a == 0 && atn == 0){bpos=9;}
string at_name = models[model]->chains[0]->residues[resnums]->atomnames1[atn];
c1 = models[model]->chains[0]->residues[resnums]->getAtom(at_name);
if(c1.size()==0){
//printf("\nMissing %s in residue %d %s\n",at_name.c_str(),pdb->models[model]->chains[0]->residues[resnums]->res_num,pdb->models[model]->chains[0]->residues[resnums]->res_name.c_str());
//abort();
} 
pol->full_atoms[bpos+ 0]=c1[0];
pol->full_atoms[bpos+ 1]=c1[1];
pol->full_atoms[bpos+ 2]=c1[2];
bpos = bpos+3;
}


}

}


}else{

if(rname.compare("  U")==0 || rname.compare("  C")==0){
vector<float> c1 = models[model]->chains[0]->residues[resnums]->getAtom(" C1'"); 
vector<float> c2 = models[model]->chains[0]->residues[resnums]->getAtom(" N1 "); 
vector<float> c3 = models[model]->chains[0]->residues[resnums]->getAtom(" C2 "); 

pol->mols[resnums]->x=c1[0];
pol->mols[resnums]->y=c1[1];
pol->mols[resnums]->z=c1[2];

pol->mols[resnums]->bx=c2[0];
pol->mols[resnums]->by=c2[1];
pol->mols[resnums]->bz=c2[2];

pol->mols[resnums]->b2x=c3[0];
pol->mols[resnums]->b2y=c3[1];
pol->mols[resnums]->b2z=c3[2];
//Add atoms to full_atoms
if(rname.compare("  U")==0){
pol->mols[resnums]->seq='U';
int bpos=a*81;
int ant_i=0;
if(a == 0){ant_i=3;}
for(int atn = ant_i; atn < 23; atn++){
if(a == 0 && atn == 0){bpos=9;}
string at_name = models[model]->chains[0]->residues[resnums]->atomnames3[atn];
c1 = models[model]->chains[0]->residues[resnums]->getAtom(at_name);
if(c1.size()==0){
//printf("\nMissing %s in residue %d %s\n",at_name.c_str(),pdb->models[model]->chains[0]->residues[resnums]->res_num,pdb->models[model]->chains[0]->residues[resnums]->res_name.c_str());
//abort();
} 
pol->full_atoms[bpos+ 0]=c1[0];
pol->full_atoms[bpos+ 1]=c1[1];
pol->full_atoms[bpos+ 2]=c1[2];
bpos = bpos+3;
}

}else{

if(rname.compare("  C")==0){
pol->mols[resnums]->seq='C';
int bpos=a*81;
int ant_i=0;
if(a == 0){ant_i=3;}
for(int atn = ant_i; atn < 24; atn++){
if(a == 0 && atn == 0){bpos=9;}
string at_name = models[model]->chains[0]->residues[resnums]->atomnames4[atn];
c1 = models[model]->chains[0]->residues[resnums]->getAtom(at_name);
if(c1.size()==0){
//printf("\nMissing %s in residue %d %s\n",at_name.c_str(),pdb->models[model]->chains[0]->residues[resnums]->res_num,pdb->models[model]->chains[0]->residues[resnums]->res_name.c_str());
//abort();
} 
pol->full_atoms[bpos+ 0]=c1[0];
pol->full_atoms[bpos+ 1]=c1[1];
pol->full_atoms[bpos+ 2]=c1[2];
bpos = bpos+3;
}


}

}


}else{
//Failed to load residue
printf("Failed to load PDB: residue names dont' match(A C G U)");
//abort();

}
}
//Get the C1', N9 and C4

//Get the C1', N1 and C2

//Foreach residue in the file extract the C1' atoms coordinates and the base cooridnates and add to mols vector
a++;
}
//Loop through all residues and load all coordinates into the full_atoms vector 
pol->num_mols=a;
pols->push_back(pol);
printf("Model %d\n",model);
}
}

Pdb::Pdb(Polymere * pol,char chain){
//adds the RNA to a model 
models.push_back(new Model(pol,'J'));
}

Pdb::Pdb(string file){
using namespace std;
vector<string> lines;
vector<int> model_start;
vector<int> model_end;
vector<int> chain_start;
vector<int> chain_end;
vector<int> res_start;
vector<int> res_end;
string s_model ( "MODEL ");
string s_endmodel ( "ENDMDL");
string s_atom     ( "ATOM  ");

//Scan through file and find 
std::ifstream input;
input.open(file.c_str());
if(input.fail()){
printf("Reading PDB file failed: %s",file.c_str());
}else{
string line;
int index=0;
string prevchain ("");
string prevres ("");
int prevresnum= -1;
while(getline(input, line)) {
lines.push_back(line);
string rec_id = getElement(line,0);
if(rec_id.compare(s_model)==0){
model_start.push_back(index);

}
if(rec_id.compare(s_endmodel)==0){
model_end.push_back(index);
//End residues,chains
//chain_end.push_back(index);
//res_end.push_back(index);
prevchain="";
prevres="";
}
if(rec_id.compare(s_atom)==0){
string curchain=getElement(line,7);
string curres=getElement(line,5);
int curresnum = strtoint(getElement(line,8));
//Get Chain Name
//printf("Chain: %s, Residue%s\n",curchain.c_str(),curres.c_str());
//printf("Chain: %s, Residue%s\n",prevchain.c_str(),prevres.c_str());
if(curchain.compare(prevchain)!=0){
if(chain_start.size() > 0){
chain_end.push_back(index);
}
chain_start.push_back(index);
prevchain=curchain;
}
//IF res_name changes add to res start
if(curres.compare(prevres)!=0 || curresnum!=prevresnum){
if(res_start.size() > 0){
res_end.push_back(index);
}
res_start.push_back(index);
prevres=curres;
prevresnum=curresnum;
}

}

index++;
}





}

if(model_start.size()==0){model_start.push_back(0);}
if(model_start.size() > model_end.size()){model_end.push_back(lines.size());}
if(chain_start.size() > chain_end.size()){chain_end.push_back(lines.size());}
if(res_start.size() > res_end.size()){res_end.push_back(lines.size());}

//Print out the starting positions
//for(int i =0; i<model_start.size(); i++){printf("model %d: %d - %d",i,model_start[i],model_end[i]);}
//for(int i =0; i<chain_start.size(); i++){printf("chain %d: %d - %d",i,chain_start[i],chain_end[i]);}
//for(int i =0; i<res_start.size(); i++){printf("chain %d: %d - %d",i,res_start[i],res_end[i]);}
//printf("Models:%d\n",model_start.size());
//printf("Chains:%d\n",chain_start.size());
//printf("Residues:%d\n",res_start.size());

//Loop through models and chains to create PDB object
for(int i=0; i < model_start.size(); i++){
//Make a new model
int mstart = model_start[i];
int mend = model_end[i];
models.push_back(new Model(&lines,&chain_start,&chain_end,&res_start,&res_end,mstart,mend));
}

input.close();
}

void Pdb::clear(){
//Delete PDB
for(int m=0; m<models.size(); m++){
for(int c=0; c<models[m]->chains.size(); c++){
for(int r=0; r<models[m]->chains[c]->residues.size(); r++){
for(int a=0; a<models[m]->chains[c]->residues[r]->atoms.size(); a++){
delete models[m]->chains[c]->residues[r]->atoms[a];
}
delete models[m]->chains[c]->residues[r];
}
delete models[m]->chains[c];
}
delete models[m];
} 

}



void Pdb::standard(const char* file){
//converts all atoms in the current pdb file to a common standard used by the simulator
vector<string> current_atom;
vector<string> standard_atom;

vector<string> current_res;
vector<string> standard_res;
//Read through the file and create the mapping vectors
vector<string> lines;
string s_atom     ( "AMAP  ");
string s_res      ( "RMAP  ");
//Scan through file and find 
std::ifstream input;
input.open(file);
if(input.fail()){
printf("Reading Map file failed\n");
}else{
printf("Read in map file\n");
string line;
while(getline(input, line)) {
string map=line.substr(0,6);
if(map.compare(s_atom)==0){
string aold = line.substr(6,4);
string anew = line.substr(10,4);
current_atom.push_back(aold);
standard_atom.push_back(anew);
}else{
if(map.compare(s_res)==0){
string aold = line.substr(6,3);
string anew = line.substr(9,3);
current_res.push_back(aold);
standard_res.push_back(anew);
}

}


}
}

//Loop through all models and chains and rename the residues and atoms using the map
//printf("Models: %d\n",models.size());
for(int m=0; m<models.size(); m++){
//printf("Chains: %d\n",models[m]->chains.size());
for(int c=0; c<models[m]->chains.size(); c++){
//printf("Residues: %d\n",models[m]->chains[c]->residues.size());
for(int r=0; r<models[m]->chains[c]->residues.size(); r++){
//Rename any residues
for(int amaps=0; amaps < current_res.size(); amaps++){
if(models[m]->chains[c]->residues[r]->res_name.compare(current_res[amaps])==0){
models[m]->chains[c]->residues[r]->res_name=standard_res[amaps];
}
}

for(int a=0; a<models[m]->chains[c]->residues[r]->atoms.size(); a++){

for(int amaps=0; amaps < current_atom.size(); amaps++){
if(models[m]->chains[c]->residues[r]->atoms[a]->atom_name.compare(current_atom[amaps])==0){
models[m]->chains[c]->residues[r]->atoms[a]->atom_name=standard_atom[amaps];
}
}


}
}
}
}


}

void Pdb::splits(string file,int tmods){
int start=0;
int sim_id=0;
for(int tm=0; tm<models.size(); tm=tm+tmods){
//Append a number to the end of the file
string fname(file.c_str());
fname += boost::lexical_cast<std::string>(sim_id++);
FILE * sfile = fopen(fname.c_str(),"a");
//Foreach line write the coordinates of each c1' atom of the polymere
//loop through all models and print out the pdb
int end=start+tmods;
printf("File: %s - %d to %d\n",fname.c_str(),start,end);
if(end > models.size()){end=models.size();}
for(int m=start; m<end; m++){
fprintf(sfile,"MODEL  \n");
for(int c=0; c<models[m]->chains.size(); c++){
for(int r=0; r<models[m]->chains[c]->residues.size(); r++){
// Print the residue in PDB format
string rec_id   ="ATOM  ";
string blank=" ";
//Foreach atom in the residue set the values below and print a line
for(int a=0; a < models[m]->chains[c]->residues[r]->atoms.size(); a++){
//Atom Variables
int atom_num		=models[m]->chains[c]->residues[r]->atoms[a]->atom_num;
string atom_name	=models[m]->chains[c]->residues[r]->atoms[a]->atom_name;
string alt_loc		=models[m]->chains[c]->residues[r]->atoms[a]->alt_loc;
string insert_code	=models[m]->chains[c]->residues[r]->atoms[a]->insert_code;
float x			=models[m]->chains[c]->residues[r]->atoms[a]->x;
float y			=models[m]->chains[c]->residues[r]->atoms[a]->y;
float z			=models[m]->chains[c]->residues[r]->atoms[a]->z;
float occ		=models[m]->chains[c]->residues[r]->atoms[a]->occ;
float bval		=models[m]->chains[c]->residues[r]->atoms[a]->bval;
string footnote		=models[m]->chains[c]->residues[r]->atoms[a]->footnote;
fprintf(sfile,"%6s%5d%1s%4s%1s%3s%1s%1s%4d%1s%3s%8.3f%8.3f%8.3f%6.2f%6.2f%1s%3s\n",
rec_id.c_str(),
atom_num,
blank.c_str(),
atom_name.c_str(),
alt_loc.c_str(),
models[m]->chains[c]->residues[r]->res_name.c_str(),
blank.c_str(),
models[m]->chains[c]->residues[r]->chain.c_str(),
models[m]->chains[c]->residues[r]->res_num,
insert_code.c_str(),
blank.c_str(),
x,
y,
z,
occ,
bval,
blank.c_str(),
footnote.c_str()
);
}
}
}
fprintf(sfile,"ENDMDL\n");
}
start=end;
fflush( sfile );
fclose(sfile);
}
}


void Pdb::write(string file){
FILE * sfile = fopen(file.c_str(),"a");
//Foreach line write the coordinates of each c1' atom of the polymere
//loop through all models and print out the pdb
for(int m=0; m<models.size(); m++){
fprintf(sfile,"MODEL  \n");
for(int c=0; c<models[m]->chains.size(); c++){
for(int r=0; r<models[m]->chains[c]->residues.size(); r++){
// Print the residue in PDB format
string rec_id   ="ATOM  ";
string blank=" ";
//Foreach atom in the residue set the values below and print a line
for(int a=0; a < models[m]->chains[c]->residues[r]->atoms.size(); a++){
//Atom Variables
int atom_num		=models[m]->chains[c]->residues[r]->atoms[a]->atom_num;
string atom_name	=models[m]->chains[c]->residues[r]->atoms[a]->atom_name;
string alt_loc		=models[m]->chains[c]->residues[r]->atoms[a]->alt_loc;
string insert_code	=models[m]->chains[c]->residues[r]->atoms[a]->insert_code;
float x			=models[m]->chains[c]->residues[r]->atoms[a]->x;
float y			=models[m]->chains[c]->residues[r]->atoms[a]->y;
float z			=models[m]->chains[c]->residues[r]->atoms[a]->z;
float occ		=models[m]->chains[c]->residues[r]->atoms[a]->occ;
float bval		=models[m]->chains[c]->residues[r]->atoms[a]->bval;
string footnote		=models[m]->chains[c]->residues[r]->atoms[a]->footnote;
fprintf(sfile,"%6s%5d%1s%4s%1s%3s%1s%1s%4d%1s%3s%8.3f%8.3f%8.3f%6.2f%6.2f%1s%3s\n",
rec_id.c_str(),
atom_num,
blank.c_str(),
atom_name.c_str(),
alt_loc.c_str(),
models[m]->chains[c]->residues[r]->res_name.c_str(),
blank.c_str(),
models[m]->chains[c]->residues[r]->chain.c_str(),
models[m]->chains[c]->residues[r]->res_num,
insert_code.c_str(),
blank.c_str(),
x,
y,
z,
occ,
bval,
blank.c_str(),
footnote.c_str()
);
}
}
}
fprintf(sfile,"ENDMDL\n");
}
fflush( sfile );
fclose(sfile);

}

void Pdb::write(char * file){
FILE * sfile = fopen(file,"a");
//Foreach line write the coordinates of each c1' atom of the polymere
//loop through all models and print out the pdb
for(int m=0; m<models.size(); m++){
fprintf(sfile,"MODEL  \n");
for(int c=0; c<models[m]->chains.size(); c++){
for(int r=0; r<models[m]->chains[c]->residues.size(); r++){
// Print the residue in PDB format
string rec_id   ="ATOM  ";
string blank=" ";
//Foreach atom in the residue set the values below and print a line
for(int a=0; a < models[m]->chains[c]->residues[r]->atoms.size(); a++){
//Atom Variables
int atom_num		=models[m]->chains[c]->residues[r]->atoms[a]->atom_num;
string atom_name	=models[m]->chains[c]->residues[r]->atoms[a]->atom_name;
string alt_loc		=models[m]->chains[c]->residues[r]->atoms[a]->alt_loc;
string insert_code	=models[m]->chains[c]->residues[r]->atoms[a]->insert_code;
float x			=models[m]->chains[c]->residues[r]->atoms[a]->x;
float y			=models[m]->chains[c]->residues[r]->atoms[a]->y;
float z			=models[m]->chains[c]->residues[r]->atoms[a]->z;
float occ		=models[m]->chains[c]->residues[r]->atoms[a]->occ;
float bval		=models[m]->chains[c]->residues[r]->atoms[a]->bval;
string footnote		=models[m]->chains[c]->residues[r]->atoms[a]->footnote;
fprintf(sfile,"%6s%5d%1s%4s%1s%3s%1s%1s%4d%1s%3s%8.3f%8.3f%8.3f%6.2f%6.2f%1s%3s\n",
rec_id.c_str(),
atom_num,
blank.c_str(),
atom_name.c_str(),
alt_loc.c_str(),
models[m]->chains[c]->residues[r]->res_name.c_str(),
blank.c_str(),
models[m]->chains[c]->residues[r]->chain.c_str(),
models[m]->chains[c]->residues[r]->res_num,
insert_code.c_str(),
blank.c_str(),
x,
y,
z,
occ,
bval,
blank.c_str(),
footnote.c_str()
);
}
}
}
fprintf(sfile,"ENDMDL\n");
}
fflush( sfile );
fclose(sfile);

}

void Pdb::print(){
//loop through all models and print out the pdb
for(int m=0; m<models.size(); m++){
printf("MODEL  \n");
models[m]->print();
printf("ENDMDL\n");
}
}

