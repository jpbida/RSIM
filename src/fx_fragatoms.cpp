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

using namespace std;
const int VECSIZE1 = 9; //Number of Fragments
Fragatoms::Fragatoms():atoms(VECSIZE1) {
using namespace std;
if(atoms.size()<9){
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

Fragatoms::Fragatoms(string line){
//Parses a line from an input file and loads information into atoms class
if(atoms.size()<9){
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
for(int w=0; w<9; w++){
		getline(iss,field,'\t');
	x=strtodouble(field);
		getline(iss,field,'\t');
	y=strtodouble(field);
		getline(iss,field,'\t');
	z=strtodouble(field);
atoms[w].x=x;
atoms[w].y=y;
atoms[w].z=z;
if((w+1) % 3 == 0){		
getline(iss,field,'\t');
char	nt1=strtochar(field);
}
}
}

void Fragatoms::loadLine(string line){
//Parses a line from an input file and loads information into atoms class
string field;
double x,y,z;
	istringstream iss(line,istringstream::in);	
		getline(iss,field,'\t');
int	id2=strtoint(field);
id=id2;
for(int w=0; w<9; w++){
		getline(iss,field,'\t');
	x=strtodouble(field);
		getline(iss,field,'\t');
	y=strtodouble(field);
		getline(iss,field,'\t');
	z=strtodouble(field);
atoms[w].x=x;
atoms[w].y=y;
atoms[w++].z=z;
if((w+1) % 3 == 0){		
getline(iss,field,'\t');
char	nt1=strtochar(field);
}
}
}

void Fragatoms::print(){
std::printf("Atom Cordinates for fragment: \%d\nx\ty\tz\n",id);
for(int w=0; w<9; w++){
std::printf("%3.2f\t%3.2f\t%3.2f\n",atoms[w].x,atoms[w].y,atoms[w].z);
}
}
