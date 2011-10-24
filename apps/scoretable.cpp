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
srand(time(0));
std::vector<string> lines;
std::ifstream input;
input.open(opt->starting_struct.c_str());
if(input.fail()){
printf("Reading PDB list file failed");
}else{
string line;
int index=0;
while(getline(input, line)) {
	lines.push_back(line);
}
}

vector<double> params;
for(int i=0; i<10; i++){
params.push_back(0);
}
vector<double> dparams;
loadDoublets(&dparams,opt->doublet_params);
for(int l=0; l<lines.size(); l++){
Polymere * pol=new Polymere(lines[l].c_str());
Score * score = new Score(pol->num_mols);
pol->resscore(score,&dparams);
int badback=pol->checkBackBone();
int totclash=pol->atomClashes();
printf("\nBad Backbone: %d Total Clashes: %d\n",badback,totclash);
float g0=score->gscore(0);
float g1=score->gscore(1);
float g2=score->gscore(2);
float g3=score->gscore(3);
float g4=score->gscore(4);
float g5=score->gscore(5);
float g6=score->gscore(6);
printf("\nG-Score 1: '%8.3f'\n",g1);
printf("\nG-Score 2: '%8.3f'\n",g2);
printf("\nG-Score 3: '%8.3f'\n",g3);
printf("\nG-Score 4: '%8.3f'\n",g4);
printf("\nG-Score 5: '%8.3f'\n",g5);
printf("\nG-Score 6: '%8.3f'\n",g6);
float ogscore1=(g1+g2+g3+g5)/pol->num_mols;
float ogscore2=g6/pol->num_mols;
float ogscore3=g4/pol->num_mols;
float overall_score = (ogscore1) + ((1/(1+exp(-1*(ogscore1-20)))))*(12) + (-1*(1/(1+exp(-1*(ogscore1-20)))-1))*(ogscore2) +  ((-1*(1/(1+exp(-1*(ogscore2-15)))-1))*(-1*(1/(1+exp(-1*(ogscore1-20)))-1)))*(-10*ogscore3);
printf("Overall: %8.3f\n",overall_score);

score->sim_id=pol->num_mols;
score->gscore0=g0;
score->gscore1=g1;
score->gscore2=g2;
score->gscore3=g3;
score->gscore4=g4;
score->gscore5=g5;
score->gscore6=g6;
printf("Scored %s: %d\n\n",lines[l].c_str(),pol->num_mols);
string outfile("scores.txt");
//score->write_globals(outfile);
score->write("scores.txt");
delete score;
}


}
