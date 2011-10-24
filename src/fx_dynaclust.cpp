// Montecarlo fractal kinetics simulator //
#include <string>
#include <cmath>
#include "stdlib.h"
#include "stdio.h"
#include <cstring>
#include <istream>
#include <sstream>
#include <vector>
#include <iostream>
#include "dynaclust.h"

Dynaclust::Dynaclust(int clusts,float res){
using namespace std;
//Sets-up the dynaclust object with a target number of clsuters, clusts, a particular resolution, res and the total number of positions moved numb_pos; 
cur_tot=0;
//Intialize movement positions
Nc = clusts;
resolution = res;

//Initialize frag_set and pos_nc 
float t=0;
for(int c=0; c<=(Nc*3); c++){pos_nc.push_back(t);}
}


//Returns 0 if total number of clusters have been reached 1 if more needed
int Dynaclust::add(float xpos_n, float ypos_n, float zpos_n){
using namespace std;
int outval=0;
int addTo = 0;
if(cur_tot >= Nc){
outval=-1;
}else{
//Check if the coords are within 1A of anything already in the set
for(int i = 0; i < cur_tot; i++){
float xpos_c = pos_nc[(i*3+0)];
float ypos_c = pos_nc[(i*3+1)];
float zpos_c = pos_nc[(i*3+2)];

float dist = sqrt(pow((xpos_c - xpos_n),2) + pow((ypos_c - ypos_n),2) + pow((zpos_c - zpos_n),2));
if(dist < resolution){
//printf("Within Minimal Distance: %8.3f\n",dist);
addTo=1;
break;
}
}
if(addTo==0){
pos_nc[3*cur_tot+0]=xpos_n;
pos_nc[3*cur_tot+1]=ypos_n;
pos_nc[3*cur_tot+2]=zpos_n;
cur_tot++;
outval= 1;
printf("Total Clusts: %d\n",cur_tot);
}else{
outval=0;
}
} 
return outval;
}


