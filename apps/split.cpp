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
Pdb * pdb = new Pdb(opt->starting_struct);
pdb->standard("pdb.map");
int tots=pdb->models.size();
printf("\nTotal Number of Structures: %d\n",tots);
pdb->splits(opt->prefix,opt->iteration_num);
//Read in the current centers
//Delete PDB
pdb->clear();
delete pdb;

}//end main
