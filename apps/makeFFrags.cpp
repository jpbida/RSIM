#include <ctime>
#include <fstream>
#include <iostream>
#include "base.h"
#include "fragment.h"
#include "polymere.h"
#include <string>
#include <cmath>
#include <vector>
#include <istream>
#include <sstream>
#include "options.h"
const double PI = 3.14159265358979323846;
int main(int argc, char **argv)
{
using namespace std;
Options * opt=new Options(argc,argv);
Polymere * pol = new Polymere(opt->starting_struct);
vector<Fragment *> frags;
frags.reserve(opt->iteration_num);
//Setting up some test molecules for our polymere 
std::ifstream input;
input.open(opt->fragment_file.c_str());
string line;
while(getline(input, line)){frags.push_back(new Fragment(line));}
input.close();
cout << "Fragments Loaded: " << frags.size() << endl;
int tot = frags.size();
pol->alignStruct(1);
//pol->printMols(0);
int w=0;
int f=0;
pol->changeTor2(1,0);
pol->changeBond2(1,PI/2);
pol->alignStruct(1);
pol->changeTor2(4,0);
pol->changeBond2(5,PI/2);
int fid=0;
for(int j=0; j<tot; j++){
pol->changeFrag2(2,j,frags);
pol->alignStruct(1);
fid=frags[j]->frag_type;
//pol->printMols(fid);
//Print only the vector information
//Origin of coordinate system

Base * v1=pol->mols[1];
Base * v5=pol->mols[5];
Base * v6=pol->mols[6];
Base * v4=pol->mols[4];

double ox=v5->x-v1->x;
double oy=v5->y-v1->y;
double oz=v5->z-v1->z;

double vx=v5->x-v4->x;
double vy=v5->y-v4->y;
double vz=v5->z-v4->z;
double dist1 = sqrt(pow(vx,2)+pow(vy,2)+pow(vz,2));
vx=vx/dist1;
vy=vy/dist1;
vz=vz/dist1;
double rx=v6->x-v5->x;
double ry=v6->y-v5->y;
double rz=v6->z-v5->z;
double dist2 = sqrt(pow(rx,2)+pow(ry,2)+pow(rz,2));
rx=rx/dist2;
ry=ry/dist2;
rz=rz/dist2;
cout <<fid << ","<< ox << "," << oy << "," << oz << "," << vx << "," << vy << "," << vz << "," << rx << "," << ry << "," << rz << endl;

}
}//end main
