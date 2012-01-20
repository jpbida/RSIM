/*This file is part of RSIM.

    RSIM is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    RSIM is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with RSIM.  If not, see <http://www.gnu.org/licenses/>.

*/
// Program determines all planar base-pairings that occur in the RNA structure
// Output File:
// pdb_id, model, chain, Res_id1, Res_id2, Face1, Face2, hydro

#include <fstream>
#include <boost/lexical_cast.hpp>
#include <string>
#include <iostream>
#include <istream>
#include <sstream>
#include <vector>
#include "score.h"
#include "options.h"
#include "residue.h"
#include "polymere.h"
#include "pdb.h"
#include "parse.h"
#include "parser.h"
 

int main(int argc,char **argv) {
using namespace std;

printf("Converting PDB into fuzzy Secondary structure\n");
Options * opt=new Options(argc,argv);
//Reading in Fragments and Scoring Databases
//Identify total number of models and loop through all models in the pdb file //
Pdb * pdb = new Pdb(opt->master_struct);
printf("Read in PDB file with %d models\n",pdb->models.size());
std::vector<Polymere * > pols;
pdb->ModelsToPols(&pols);

//Foreach model write out the statistics //
printf("Read in PDB file with %d models\n",pols.size());
string ofile(opt->prefix.c_str());
FILE * sfile = fopen(ofile.c_str(), "a");

for(int ply=0; ply<pols.size(); ply++){
Polymere * cen = pols[ply];
cen->updateFull();
//Set entire molecule to movement group 
for(int m=0; m<cen->num_mols; m++){cen->mols[m]->group=5;}
//Created secondary structure
string secStructure;
for(int m=0; m<cen->num_mols; m++){
secStructure.append(".");
}
std::vector<std::pair<int,int> > bp_cons;
float total_hydro=0;
printf("Total Nucleotides: %d:%d\n",cen->num_mols,cen->mols.size());
int b1=-1;
int b2=-1;
for(int i=0; i<(cen->num_mols); i++){
//Identify the nucleotide that shares the most hydrogen bonds //
float mx_bond=0;
int mx_pos=-1;
int mx_face1=-1;
int mx_face2=-1;
for(int j=0; j<cen->num_mols; j++){
if(abs(j-i)>2){
std::pair<int,int> faces;
faces.first=-1;
faces.second=-1;
float hydro=cen->hydroType(i,j,&faces);
if(hydro > 0.1){
mx_face1=faces.first;
mx_face2=faces.second;
mx_bond=hydro;
mx_pos=j;
fprintf(sfile,"%s,%d,%d,%d,%d,%d,%8.3f\n",opt->pdb_id.c_str(),ply,i,j,mx_face1,mx_face2,hydro);
fflush(sfile);
}

}
}

}

}
fclose(sfile);
}

