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
// Montecarlo fractal kinetics simulator //
#include <istream>
#include <sstream>
#include <vector>
#include "points.h"
#include "fragatoms.h"
#include <iostream>
#include "parse.h"
#include "fragall.h"

using namespace std;
const int VECSIZE1 = 36; //Number of Fragments
Fragall::Fragall():atoms(VECSIZE1){
using namespace std;
if(atoms.size()<36){
//Need to figure out why this is being called twice or something
for (int k=0; k<VECSIZE1; k++) {
Points pt;
pt.z=0;
pt.x=0;
pt.y=0;
atoms.push_back(pt);
//cout << "k: " << k << endl;
}
}
}

Fragall::Fragall(string line){
//Parses a line from an input file and loads information into atoms class
if(atoms.size()<36){
//Need to figure out why this is being called twice or something
for (int k=0; k<VECSIZE1; k++) {
Points pt;
pt.z=0;
pt.x=0;
pt.y=0;
atoms.push_back(pt);
//cout << "k: " << k << endl;
}
}
string field;
double x,y,z;
	istringstream iss(line,istringstream::in);	
		getline(iss,field,'\t');
int	id2=strtoint(field);
id=id2;
for(int w=0; w<36; w++){
		getline(iss,field,'\t');
	x=strtodouble(field);
		getline(iss,field,'\t');
	y=strtodouble(field);
		getline(iss,field,'\t');
	z=strtodouble(field);
atoms[w].x=x;
atoms[w].y=y;
atoms[w].z=z;
}
}

void Fragall::loadLine(string line){
//Parses a line from an input file and loads information into atoms class
string field;
double x,y,z;
	istringstream iss(line,istringstream::in);	
		getline(iss,field,'\t');
int	id2=strtoint(field);
id=id2;
for(int w=0; w<36; w++){
		getline(iss,field,'\t');
	x=strtodouble(field);
		getline(iss,field,'\t');
	y=strtodouble(field);
		getline(iss,field,'\t');
	z=strtodouble(field);
atoms[w].x=x;
atoms[w].y=y;
atoms[w].z=z;
}
}

int Fragall::check(){
int i=0;
std::vector<float> dists;
dists.push_back(sqrt(pow((atoms[i+0].x -atoms[i+1].x),2) + pow((atoms[i+0].y  - atoms[i+1].y),2)  + pow((atoms[i+0].z - atoms[i+1].z),2)) );
dists.push_back(sqrt(pow((atoms[i+0].x -atoms[i+2].x),2) + pow((atoms[i+0].y  - atoms[i+2].y),2)  + pow((atoms[i+0].z - atoms[i+2].z),2)) );
dists.push_back(sqrt(pow((atoms[i+0].x -atoms[i+3].x),2) + pow((atoms[i+0].y  - atoms[i+3].y),2)  + pow((atoms[i+0].z - atoms[i+3].z),2)) );
dists.push_back(sqrt(pow((atoms[i+3].x -atoms[i+5].x),2) + pow((atoms[i+3].y  - atoms[i+5].y),2)  + pow((atoms[i+3].z - atoms[i+5].z),2)) );
dists.push_back(sqrt(pow((atoms[i+7].x -atoms[i+5].x),2) + pow((atoms[i+7].y  - atoms[i+5].y),2)  + pow((atoms[i+7].z - atoms[i+5].z),2)) );
dists.push_back(sqrt(pow((atoms[i+7].x -atoms[i+9].x),2) + pow((atoms[i+7].y  - atoms[i+9].y),2)  + pow((atoms[i+7].z - atoms[i+9].z),2)) );
dists.push_back(sqrt(pow((atoms[i+9].x -atoms[i+11].x),2) + pow((atoms[i+9].y - atoms[i+11].y),2) + pow((atoms[i+9].z - atoms[i+11].z),2)) );
dists.push_back(sqrt(pow((atoms[i+6].x -atoms[i+11].x),2) + pow((atoms[i+6].y - atoms[i+11].y),2) + pow((atoms[i+6].z - atoms[i+11].z),2)) );
dists.push_back(sqrt(pow((atoms[i+6].x -atoms[i+5].x),2) + pow((atoms[i+6].y  - atoms[i+5].y),2)  + pow((atoms[i+6].z - atoms[i+5].z),2)) );
i=12;
dists.push_back(sqrt(pow((atoms[i+0].x -atoms[i+1].x),2) + pow((atoms[i+0].y - atoms[i+1].y),2) + pow((atoms[i+0].z - atoms[i+1].z),2)) );
dists.push_back(sqrt(pow((atoms[i+0].x -atoms[i+2].x),2) + pow((atoms[i+0].y - atoms[i+2].y),2) + pow((atoms[i+0].z - atoms[i+2].z),2)) );
dists.push_back(sqrt(pow((atoms[i+0].x -atoms[i+3].x),2) + pow((atoms[i+0].y - atoms[i+3].y),2) + pow((atoms[i+0].z - atoms[i+3].z),2)) );
dists.push_back(sqrt(pow((atoms[i+3].x -atoms[i+5].x),2) + pow((atoms[i+3].y - atoms[i+5].y),2) + pow((atoms[i+3].z - atoms[i+5].z),2)) );
dists.push_back(sqrt(pow((atoms[i+7].x -atoms[i+5].x),2) + pow((atoms[i+7].y - atoms[i+5].y),2) + pow((atoms[i+7].z - atoms[i+5].z),2)) );
dists.push_back(sqrt(pow((atoms[i+7].x -atoms[i+9].x),2) + pow((atoms[i+7].y - atoms[i+9].y),2) + pow((atoms[i+7].z - atoms[i+9].z),2)) );
dists.push_back(sqrt(pow((atoms[i+9].x -atoms[i+11].x),2) + pow((atoms[i+9].y -atoms[i+11].y),2) + pow((atoms[i+9].z -atoms[i+11].z),2)) );
dists.push_back(sqrt(pow((atoms[i+6].x -atoms[i+11].x),2) + pow((atoms[i+6].y -atoms[i+11].y),2) + pow((atoms[i+6].z -atoms[i+11].z),2)) );
dists.push_back(sqrt(pow((atoms[i+6].x -atoms[i+5].x),2) + pow((atoms[i+6].y - atoms[i+5].y),2) + pow((atoms[i+6].z - atoms[i+5].z),2)) );
i=24;
dists.push_back(sqrt(pow((atoms[i+0].x -atoms[i+1].x),2) + pow((atoms[i+0].y - atoms[i+1].y),2) + pow((atoms[i+0].z - atoms[i+1].z),2)) );
dists.push_back(sqrt(pow((atoms[i+0].x -atoms[i+2].x),2) + pow((atoms[i+0].y - atoms[i+2].y),2) + pow((atoms[i+0].z - atoms[i+2].z),2)) );
dists.push_back(sqrt(pow((atoms[i+0].x -atoms[i+3].x),2) + pow((atoms[i+0].y - atoms[i+3].y),2) + pow((atoms[i+0].z - atoms[i+3].z),2)) );
dists.push_back(sqrt(pow((atoms[i+3].x -atoms[i+5].x),2) + pow((atoms[i+3].y - atoms[i+5].y),2) + pow((atoms[i+3].z - atoms[i+5].z),2)) );
dists.push_back(sqrt(pow((atoms[i+7].x -atoms[i+5].x),2) + pow((atoms[i+7].y - atoms[i+5].y),2) + pow((atoms[i+7].z - atoms[i+5].z),2)) );
dists.push_back(sqrt(pow((atoms[i+7].x -atoms[i+9].x),2) + pow((atoms[i+7].y - atoms[i+9].y),2) + pow((atoms[i+7].z - atoms[i+9].z),2)) );
dists.push_back(sqrt(pow((atoms[i+9].x -atoms[i+11].x),2) + pow((atoms[i+9].y - atoms[i+11].y),2) + pow((atoms[i+9].z - atoms[i+11].z),2)) );
dists.push_back(sqrt(pow((atoms[i+6].x -atoms[i+11].x),2) + pow((atoms[i+6].y - atoms[i+11].y),2) + pow((atoms[i+6].z - atoms[i+11].z),2)) );
dists.push_back(sqrt(pow((atoms[i+6].x -atoms[i+5].x),2) + pow((atoms[i+6].y - atoms[i+5].y),2) + pow((atoms[i+6].z - atoms[i+5].z),2)) );
float sout=0;
int out;
for(int i=0; i<dists.size(); i++){
sout=sout+dists[i];
}
if(sout > 43.75 || sout < 42.75){out=0;}else{out=1;}
return out;
}


void Fragall::print(){
std::printf("Atom Cordinates for fragment: \%d\nx\ty\tz\n",id);
for(int w=0; w<36; w++){
std::printf("%3.2f\t%3.2f\t%3.2f\n",atoms[w].x,atoms[w].y,atoms[w].z);
}
}
