// Montecarlo fractal kinetics simulator //
#include "polymere.h"
#include <string>
#include <istream>
#include <sstream>
#include <vector>
#include <iostream>
#include "parse.h"
#include "residue.h"
#include "atom.h"

//Backbone Names
const char Residue::atomnames1[27][5]={
" P  ", //0
" OP1",
" OP2",
" O5'",
" C5'",
" C4'", //5
" O4'",
" C3'",
" O3'",
" C2'",
" O2'",//10
" C1'",
" N9 ",
" C4 ",
" N3 ",
" C2 ",//15
" N2 ",
" N1 ",
" C6 ",
" O6 ",
" C5 ",//20
" N7 ",
" C8 ",
" H1 ",
" H8 ",
" H21",//25
" H22"
};

const char Residue::atomnames2[26][5]={
" P  ",//0
" OP1",
" OP2",
" O5'",
" C5'",
" C4'",//5
" O4'",
" C3'",
" O3'",
" C2'",
" O2'",//10
" C1'",
" N9 ",
" C4 ",
" N3 ",
" C2 ",//15
" N1 ",
" C6 ",
" N6 ",
" C5 ",
" N7 ",//20
" C8 ",
" H2 ",
" H8 ",
" H61",
" H62"};

//Uracil Atom Names
const char Residue::atomnames3[23][5]=
{
" P  ",//0
" OP1",
" OP2",
" O5'",
" C5'",
" C4'",//5
" O4'",
" C3'",
" O3'",
" C2'",
" O2'",//10
" C1'",
" N1 ",
" C2 ",
" C6 ",
" O2 ",//15
" N3 ",
" C4 ",
" O4 ",
" C5 ",
" H3 ",//20
" H5 ",
" H6 "};

//Cytosine Atom Names
const char Residue::atomnames4[24][5]={
" P  ",//0
" OP1",
" OP2",
" O5'",
" C5'",
" C4'",//5
" O4'",
" C3'",
" O3'",
" C2'",
" O2'",//10
" C1'",
" N1 ",
" C2 ",
" C6 ",
" O2 ",//15
" N3 ",
" C4 ",
" N4 ",
" C5 ",
" H5 ",//20
" H6 ",
" H41",
" H42"};
                      

string Residue::getElement(string line,int type){
int start[] = {0,6,11,12,16,17,20,21,22,26,27,30,46,54,60,66,67};
int end[] = {6,5,1,4,1,3,1,1,4,1,3,8,8,8,6,6,1,3};
string out=line.substr(start[type],end[type]);
return out;
}

Residue::Residue(){

}

Residue::Residue(vector<string> * lines, int start, int end){
using namespace std;
res_name=getElement(lines->at(start),5);
chain=getElement(lines->at(start),7);
res_num=strtoint(getElement(lines->at(start),8));

for(int l=start; l<end; l++){
string line = lines->at(l);
string frontline=line.substr(0,4);
if(frontline=="ATOM"){
res_name=line.substr(17,3);
atoms.push_back(new Atom(line));
}

}
}


int Residue::isBackBone(int a){
int out=0;
for(int j=0; j<12; j++){
if(atoms[a]->atom_name.compare(atomnames1[j])==0){out=1;}
}
return out;
}


vector<float> Residue::getAtom(string aname){
vector<float> out;
out.push_back(0);
out.push_back(0);
out.push_back(0);
int tot=0;
for(int a=0; a<atoms.size(); a++){
if(atoms[a]->atom_name.compare(aname)==0){
out[0]=atoms[a]->x;
out[1]=atoms[a]->y;
out[2]=atoms[a]->z;
tot=1;
return out;
}
}
if(tot==0){
printf("Missing: %s\n",aname.c_str());
}
return out;
}

Residue::Residue(Polymere * pol,char chain1,int resn){
	string rname("   ");
	rname[2]=pol->mols[resn]->seq;
res_name = rname;
	string cname(" ");
	cname[0]=chain1;
chain = cname;
res_num = resn;
//Foreach type of RNA residue loop through and add the appropirate number of atoms
//Get coordinate positions from full_atom vector
int bpos = resn*81;
if(pol->mols[resn]->seq=='G'){
for(int rn=0; rn < 27; rn++){
//create an atom x y z atom_name chain  res_name
int p1 = bpos+rn*3;
int p2 = bpos+rn*3+1;
int p3 = bpos+rn*3+2;
atoms.push_back(new Atom(rn, atomnames1[rn], pol->full_atoms[p1], pol->full_atoms[p2] ,pol->full_atoms[p3],pol->mols[resn]->prob));
}
}else{
if(pol->mols[resn]->seq=='A'){
for(int rn=0; rn < 26; rn++){
//create an atom x y z atom_name chain  res_name
int p1 = bpos+rn*3;
int p2 = bpos+rn*3+1;
int p3 = bpos+rn*3+2;
atoms.push_back(new Atom(rn, atomnames2[rn], pol->full_atoms[p1], pol->full_atoms[p2] ,pol->full_atoms[p3],pol->mols[resn]->prob));
}


}else{
if(pol->mols[resn]->seq=='U'){
for(int rn=0; rn < 23; rn++){
//create an atom x y z atom_name chain  res_name
int p1 = bpos+rn*3;
int p2 = bpos+rn*3+1;
int p3 = bpos+rn*3+2;
atoms.push_back(new Atom(rn, atomnames3[rn], pol->full_atoms[p1], pol->full_atoms[p2] ,pol->full_atoms[p3],pol->mols[resn]->prob));
}


}else{
if(pol->mols[resn]->seq=='C'){
for(int rn=0; rn < 24; rn++){
//create an atom x y z atom_name chain  res_name
int p1 = bpos+rn*3;
int p2 = bpos+rn*3+1;
int p3 = bpos+rn*3+2;
atoms.push_back(new Atom(rn, atomnames4[rn], pol->full_atoms[p1], pol->full_atoms[p2] ,pol->full_atoms[p3],pol->mols[resn]->prob));
}
}else{
}
}
}
}

}


Residue::Residue(const char* file, int seqnum1,char chain1){
using namespace std;
//Scan through file and find 
res_num=seqnum1;
chain = chain1;
std::ifstream input;
input.open(file);
if(input.fail()){
printf("Reading PDB file failed");
}else{
string line;
while(getline(input, line)) {
//If the first line is ATOM parse the line
string frontline=line.substr(0,4);
if(frontline=="ATOM"){
//check seqnum and chain
int seqnum2	= strtoint(line.substr(22,4));
string chain2 	= line.substr(21,1);
if(chain2[0]==chain1 && seqnum1 == seqnum2){
res_name=line.substr(17,3);
atoms.push_back(new Atom(line));
}
}
}
}
input.close();
}

void Residue::print(){
// Print the residue in PDB format
string rec_id   ="ATOM  ";
string blank=" ";
//Foreach atom in the residue set the values below and print a line

for(int a=0; a < atoms.size(); a++){
//Atom Variables
int atom_num		=atoms[a]->atom_num;
string atom_name	=atoms[a]->atom_name;
string alt_loc		=atoms[a]->alt_loc;
string insert_code	=atoms[a]->insert_code;
float x			=atoms[a]->x;
float y			=atoms[a]->y;
float z			=atoms[a]->z;
float occ		=atoms[a]->occ;
float bval		=atoms[a]->bval;
float r			=atoms[a]->r;
string footnote		=atoms[a]->footnote;
printf("%6s%5d%1s%4s%1s%3s%1s%1s%4d%1s%3s%8.3f%8.3f%8.3f%6.2f%6.2f%1s%3s\n",
rec_id.c_str(),
atom_num,
blank.c_str(),
atom_name.c_str(),
alt_loc.c_str(),
res_name.c_str(),
blank.c_str(),
chain.c_str(),
res_num,
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

void Residue::print(string file){
// Print the residue in PDB format
FILE * output=fopen(file.c_str(),"a");
string rec_id   ="ATOM  ";
string blank=" ";
//Foreach atom in the residue set the values below and print a line

for(int a=0; a < atoms.size(); a++){
//Atom Variables
int atom_num		=atoms[a]->atom_num;
string atom_name	=atoms[a]->atom_name;
string alt_loc		=atoms[a]->alt_loc;
string insert_code	=atoms[a]->insert_code;
float x			=atoms[a]->x;
float y			=atoms[a]->y;
float z			=atoms[a]->z;
float occ		=atoms[a]->occ;
float bval		=atoms[a]->bval;
float r			=atoms[a]->r;
string footnote		=atoms[a]->footnote;
fprintf(output,"%6s%5d%1s%4s%1s%3s%1s%1s%4d%1s%3s%8.3f%8.3f%8.3f%6.2f%6.2f%1s%3s\n",
rec_id.c_str(),
atom_num,
blank.c_str(),
atom_name.c_str(),
alt_loc.c_str(),
res_name.c_str(),
blank.c_str(),
chain.c_str(),
res_num,
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

fflush(output);
fclose(output);
}


