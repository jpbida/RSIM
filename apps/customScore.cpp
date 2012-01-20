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
#include <boost/lexical_cast.hpp>
#include <fstream>
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
//#include "voro/voro++.cc"
#include "parse.h"
#include "parser.h"

//Its only argument is a directory
//The directory contains a bunck of pdb files
//pdbs.txt
// 	pdb_id, file_name,classifier

//Read pdbs.txt and make a list of files to score
//write results to results/file_name.score
 

int main(int argc,char **argv) {
using namespace std;

printf("Incorporation of Packing Probability\n");
Options * opt=new Options(argc,argv);
string output_file("m_");
output_file+=boost::lexical_cast<std::string>(opt->iteration_num);
output_file.append(".json");
FILE * sfile = fopen(output_file.c_str(),"a");
srand(time(0));
std::vector<string> lines;
std::ifstream input;
input.open(opt->starting_struct.c_str());
printf("Model File: %s\n",opt->modelfile.c_str());
Parser * model=new Parser(opt->modelfile);
if(input.fail()){
printf("Reading PDB list file failed");
}else{
string line;
int index=0;
while(getline(input, line)) {
	lines.push_back(line);
}
}
input.close();

std::vector<std::pair<int,int> > edges;
std::ifstream edge_file;
input.open(opt->edge_file.c_str());
printf("Edge File: %s\n",opt->edge_file.c_str());
if(input.fail()){
printf("Reading Edge File Failed");
}else{
string line;
int index=0;
while(getline(input, line)) {
string field;
int node1;
int node2;
	istringstream iss(line,istringstream::in);	
   		getline(iss,field,' ');
        node1=strtoint(field);
                getline(iss,field,' ');
        node2=strtoint(field);
edges.push_back(std::pair<int,int>(node1,node2));
}
}


vector<double> params;
for(int i=0; i<10; i++){
params.push_back(0);
}
vector<double> dparams;
loadDoublets(&dparams,opt->doublet_params);
printf("Scoreing\n");
std::vector<float> scores;
for(int l=0; l<lines.size(); l++){
Polymere * pol=new Polymere(lines[l].c_str());
Score * score = new Score(pol->num_mols);
pol->resscore(score,&dparams);
std::vector<double> resout;
model->score(score,&resout);
float scores1=0;
for(int r=0; r<resout.size(); r++){
scores1=scores1+resout[r];
}
scores.push_back(scores1);
delete score;
}

//calculate all edge differences from the score list
std::vector<float> links;
//Calculate the mean and STD for z-score normalization
float mn=0;
float min=100000;
float n=0;
for(int e1=0; e1<edges.size(); e1++){
float se=abs(scores[edges[e1].first]-scores[edges[e1].second]);
links.push_back(se);
if(min > se){min=se;}
mn=mn+se;
n=n+1;
}
mn=mn/n;

float std=0;
n=0;
for(int e1=0; e1<edges.size(); e1++){
float se=abs(scores[edges[e1].first]-scores[edges[e1].second]);
std=std+pow((se-mn),2);
n=n+1;
}
std=sqrt(std/n);

for(int se=0; se<links.size(); se++){
fprintf(sfile,"%d %d %g\n",edges[se].first,edges[se].second,((links[se]-mn)/std+abs((min-mn)/std)));
}
}
