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
Polymere * pol1=new Polymere(opt->master_struct);
printf("Total number of residues in final structure: %d\n",pol1->mols.size());
Polymere * pol2 = new Polymere(opt->starting_struct);
std::vector<int> helix1;
std::vector<int> helix2;

helix1.push_back(1);
helix1.push_back(2);
helix1.push_back(3);
helix1.push_back(69);
helix1.push_back(70);
helix1.push_back(71);

helix2.push_back(69);
helix2.push_back(70);
helix2.push_back(71);
helix2.push_back(1);
helix2.push_back(2);
helix2.push_back(3);

pol2->alignHelix(pol1,helix1,helix2);

Pdb * pdb1=new Pdb(pol1,'j');
Pdb * pdb2=new Pdb(pol2,'j');
pdb1->write("converted1.pdb");
pdb1->clear();
pdb2->write("converted2.pdb");
pdb2->clear();
delete pdb1;
delete pdb2;


for(int i = 0; i<pol1->mols.size(); i++){delete pol1->mols[i];}
for(int i = 0; i<pol2->mols.size(); i++){delete pol2->mols[i];}
pol1->mols.clear();
pol2->mols.clear();
delete pol1;
delete pol2;
}//end main
