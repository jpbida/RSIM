// Montecarlo fractal kinetics simulator //
#include <boost/lexical_cast.hpp>
#include "algos.h"
#include "ss.h"
#include "ct.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <istream>
#include <sstream>
using namespace std;
const int VECSIZE = 6000;

CT::CT() {
using namespace std;
double z=0;
for (int k=0; k<VECSIZE; k++) {
SS * b = new SS();
b-> z        =z;
b-> x        =z;
b-> y        =z;
b-> group    =0;
b-> seq      =0;
b-> pos      =0;
b-> before   =0;
b-> after    =0;
b-> bp       =0;
nucs.push_back(b);
}
num_mols=0;
}


CT::CT(std::string ss, std::string seq){
//This is were magic happens and the bracket notation is converted into a 2D layout
//Convert secondary structure and sequence into the CT fomrat
double z=0;

//Create two fake base-pairings for the layout algorithm

for(int i=0; i<ss.length(); i++){
SS * b = new SS();
b->x=z;
b->y=z;
b->z=z;
b->pos=i;
b->after=i+1;
b->before=i-1;
b->group=0;
b->seq=seq[i];

int bp=-1;
if(ss.compare(i,1,"(")==0){
bp=forward(ss,i);
}else{
if(ss.compare(i,1,")")==0){
bp=backward(ss,i);
}

}

b->bp=bp;
nucs.push_back(b);
}
num_mols=ss.length();

}

void CT::print(){
for(int w=0; w<nucs.size(); w++){
printf("%8.3f %8.3f %8.3f %d %d %d %d %c %d\n",
nucs[w]->x,
nucs[w]->y,
nucs[w]->z,
nucs[w]->pos,
nucs[w]->after,
nucs[w]->before,
nucs[w]->group,
nucs[w]->seq,
nucs[w]->bp
);
}

}

void CT::layout(){
/*
SS * p1 = new SS();
p1->x=z;
p1->y=z;
p1->z=z;
p1->pos=0;
p1->after=1;
p1->before=-1;
p1->group=0;
p1->seq='N';
p1->bp=ss.length()+1;
nucs->push_back(p1);

SS * p2 = new SS();
p2->x=z;
p2->y=z;
p2->z=z;
p2->pos=1;
p2->after=2;
p2->before=0;
p2->group=0;
p2->seq='N';
p2->bp=ss.length();
nucs->push_back(p1);



SS * p3 = new SS();
p3->x=z;
p3->y=z;
p3->z=z;
p3->pos=ss.length()+2;
p3->after=1;
p3->before=-1;
p3->group=0;
p3->seq='N';
p3->bp=1;
nucs->push_back(p3);

SS * p4 = new SS();
p4->x=z;
p4->y=z;
p4->z=z;
p4->pos=ss.length()+3;
p4->after=2;
p4->before=0;
p4->group=0;
p4->seq='N';
p4->bp=0;
nucs->push_back(p4);
*/
}


int CT::ct2rsim(){



}



