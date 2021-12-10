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

    extractComponent -c starting.pdb -n model_num -k "++++...++++" -p output 

    Takes the string of + & . and creates a pdb file output.pdb
    containing just the molecules at the positions indicated with +'s. If the pdb file 
    contains multiple models then the model_num model is used. If model_num
    is not provided the first model is used by default 

*/
#include <algorithm>
#include <ctime>
#include <fstream>
#include <iostream>
#include "pdb.h"
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
#include "fabase.h"

int main(int argc, char **argv)
{
using namespace std;
Options * opt=new Options(argc,argv);
//Max structures
Pdb * pdb = new Pdb(opt->starting_struct);
printf("Input PDB: %s\n",opt->starting_struct.c_str());

int tots=pdb->models.size();
printf("Total Number of Models: %d\n",tots);

int model = opt->iteration_num;
printf("Model Number: %d\n",model);

string comp = opt->sequence;
printf("Component: %s\n",comp.c_str());

string output=opt->prefix;
output.append(".pdb");
//Creating vector of polymeres
printf("Output: %s\n",output.c_str()); 

std::vector<Polymere * > pols;
pdb->ModelsToPols(&pols);

//Foreach model write out the statistics //
printf("Read in PDB file with %d models\n",pols.size());


vector<char> vcomp(comp.begin(), comp.end());
Polymere * out_pol = new Polymere(); 

//Create new polymer 
for(int n=0; n < vcomp.size(); n++){
  if(vcomp[n]=='+'){
    out_pol->mols.push_back(new Base());
  }
}

for(int ply=0; ply<pols.size(); ply++){
  Polymere * cen = pols[ply];
  cen->updateFull();
  int p = 0;
  int max_mols = vcomp.size();
  if(cen->mols.size() < vcomp.size()){
    max_mols = cen->mols.size();
  }
  for(int m=0; m < max_mols; m++){
    if(vcomp[m]=='+'){
      out_pol->mols[p]-> z     =cen->mols[m]->z;
      out_pol->mols[p]-> x     =cen->mols[m]->x;
      out_pol->mols[p]-> y     =cen->mols[m]->y; 
      out_pol->mols[p]-> seq   =cen->mols[m]->seq;
      out_pol->mols[p]-> group =cen->mols[m]->group;
      out_pol->mols[p]-> bx    =cen->mols[m]->bx;
      out_pol->mols[p]-> by    =cen->mols[m]->by;
      out_pol->mols[p]-> bz    =cen->mols[m]->bz;
      out_pol->mols[p]-> b2x   =cen->mols[m]->b2x;
      out_pol->mols[p]-> b2y   =cen->mols[m]->b2y;
      out_pol->mols[p]-> b2z   =cen->mols[m]->b2z;
      out_pol->mols[p]-> id    =cen->mols[m]->id;
      out_pol->mols[p]-> prob    =cen->mols[m]->prob;
      p++;
    } 
  }

Pdb * out_pdb = new Pdb(out_pol,'J');
out_pdb->write(output);
}
}


