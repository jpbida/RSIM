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
srand(time(0));
//srand(opt->seed);
//opt->print();

printf("Creates illustration of closed moves\n");
printf("Loading: %s\n",opt->starting_struct.c_str());

Pdb * in_pdb = new Pdb(opt->starting_struct.c_str());
printf("total number of models: %d",in_pdb->models.size());

int num_models=in_pdb->models.size();

Polymere * cen = new Polymere(opt->starting_struct);
std::vector<Polymere *> all_mods;
if(num_models>1){
for(int i=1; i<num_models; i++){
Polymere * cena = new Polymere(opt->starting_struct);
all_mods.push_back(cena);
}
}

printf("cen num_mols: %d\n",cen->num_mols);
//Find two points that are the maximum distance apart
int pos1=0;
int pos2=cen->num_mols;
float dist = 0;
for(int m2=0; m2<cen->num_mols; m2++){
float ndist = cen->c1distf(pos1,m2);
if(ndist > dist){pos2=m2;dist=ndist;}
}

printf("Pos1: %d Pos2: %d Dist: %8.3f\n",pos1,pos2,dist);
//Translate on of them to the origin 
double tx=cen->mols[0]->x;
double ty=cen->mols[0]->y;
double tz=cen->mols[0]->z;
cen->translate(cen->mols[0]->x,cen->mols[0]->y,cen->mols[0]->z);
//Translate all other models
for(int i=0; i<all_mods.size(); i++){all_mods[i]->translate(tx,ty,tz);}



//Rotate the other to the position Y axis
double ang = atan2(cen->mols[pos2]->z,cen->mols[pos2]->x);
cen->rotateY(ang);
for(int i=0; i<all_mods.size(); i++){all_mods[i]->rotateY(ang);}
double phi=-1*atan2(cen->mols[pos2]->x,cen->mols[pos2]->y);
cen->rotateZ(phi);
for(int i=0; i<all_mods.size(); i++){all_mods[i]->rotateZ(phi);}


//Determine the center of mass for X-axis 
double minx=cen->mols[0]->x;
double maxx=cen->mols[0]->x;
double minz=cen->mols[0]->z;
double maxz=cen->mols[0]->z;
for(int m=0; m<cen->num_mols; m++){
if(cen->mols[m]->x > maxx){maxx=cen->mols[m]->x;}
if(cen->mols[m]->z > maxz){maxz=cen->mols[m]->z;}
if(cen->mols[m]->z < minz){minz=cen->mols[m]->z;}
if(cen->mols[m]->x < minx){minx=cen->mols[m]->x;}
}

double xtrans=minx+(maxx-minx)/2;
double ztrans=minz+(maxz-minz)/2;
cen->translate(xtrans,0,ztrans);
for(int i=0; i<all_mods.size(); i++){all_mods[i]->translate(xtrans,0,ztrans);}

double camera_distZ = (dist/2)/atan((21*PI/180));
double camera_distY = (dist/2);
double camera_distX = 0;
printf("Camera dist: %8.3f,%8.3f,%8.3f\n",camera_distZ,camera_distY,camera_distX);

printf("Pos1: %8.3f,%8.3f,%8.3f\n",cen->mols[pos1]->x,cen->mols[pos1]->y,cen->mols[pos1]->z);
printf("Pos2: %8.3f,%8.3f,%8.3f\n",cen->mols[pos2]->x,cen->mols[pos2]->y,cen->mols[pos2]->z);

//Update full_atoms
cen->updateFull();
for(int i=0; i<all_mods.size(); i++){all_mods[i]->updateFull();}
pos1=pos1*81+3*11;
pos2=pos2*81+3*11;
printf("Pos %d: %8.3f,%8.3f,%8.3f\n",pos1,cen->full_atoms[(pos1)],cen->full_atoms[(pos1+1)],cen->full_atoms[pos1+2]);
printf("Pos %d: %8.3f,%8.3f,%8.3f\n",pos2,cen->full_atoms[(pos2)],cen->full_atoms[(pos2+1)],cen->full_atoms[(pos2+2)]);

std::string ofile("pdb");

std::string outfile;
outfile.append(ofile);
outfile+=boost::lexical_cast<std::string>(0);
outfile.append(".pdb");
cen->topov(outfile,camera_distX,camera_distY,camera_distZ);
for(int i=0; i<all_mods.size(); i++){
std::string outfile;
outfile.append(ofile);
outfile+=boost::lexical_cast<std::string>(i+1);
outfile.append(".pdb");
all_mods[i]->topov(outfile,camera_distX,camera_distY,camera_distZ);
}

//Foreach model make an output file

//std::string povTriangle(std::vector<float> points,std::string color);
//std::string povSphere(std::vector<float> points,float radius,std::string color);
//std::string povCylinder(std::vector<float> points,float radius,std::string color);
//Create spheres for the C4' Atoms
//Create cylinders between C4' and C4' to C1'
//Create polygons for the bases with cylinders for 

for(int i = 0; i<cen->mols.size(); i++){
delete cen->mols[i];
}
cen->mols.clear();
delete cen;
}//end main
