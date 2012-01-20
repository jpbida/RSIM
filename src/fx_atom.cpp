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
#include <iostream>
#include "parse.h"
#include "atom.h"
#include <string>
/*
Hydrogen: .31
Carbon .730
Nitrogen 0.71
Oxygen: .66
Phosphorus: 1.07
*/
Atom::Atom(int at_num, string at_name, float x1, float y1, float z1,float bval1){
atom_num = at_num;
atom_name = at_name;
x = x1;
y = y1;
z = z1;
//Set other elements to default values
alt_loc=" ";
insert_code = " ";
occ = 0;
bval = bval1;
footnote = "   ";
if (at_name.find("H",0) != string::npos){
r=.31;
}else{
if (at_name.find("C",0) != string::npos){
r=.73;
}else{
if (at_name.find("O",0) != string::npos){
r=.66;
}else{
if (at_name.find("P",0) != string::npos){
r=1.07;
}else{
if (at_name.find("N",0) != string::npos){
r=0.71;
}else{
printf("Can't identify atom type: %s\n",at_name.c_str());
}
}
}
}
}

}
Atom::Atom(int at_num, string at_name, float x1, float y1, float z1){
atom_num = at_num;
atom_name = at_name;
x = x1;
y = y1;
z = z1;
//Set other elements to default values
alt_loc=" ";
insert_code = " ";
occ = 0;
bval = 0;
footnote = "   ";
if (at_name.find("H",0) != string::npos){
r=.31;
}else{
if (at_name.find("C",0) != string::npos){
r=.73;
}else{
if (at_name.find("O",0) != string::npos){
r=.66;
}else{
if (at_name.find("P",0) != string::npos){
r=1.07;
}else{
if (at_name.find("N",0) != string::npos){
r=0.71;
}else{
printf("Can't identify atom type: %s\n",at_name.c_str());
}
}
}
}
}

}

Atom::Atom(string line){
using namespace std;
int length=line.length();
if(length > 11){atom_num		=strtoint(line.substr(6,5).c_str());}
if(length > 15){atom_name	=line.substr(12,4);}
if(length > 16){alt_loc 		=line.substr(16,1);}
if(length > 26){insert_code	=line.substr(26,1);}
if(length > 37){x			=strtofloat(line.substr(30,8).c_str());}
if(length > 45){y			=strtofloat(line.substr(38,8).c_str());}
if(length > 52){z			=strtofloat(line.substr(46,8).c_str());}
if(length > 59){occ		=strtofloat(line.substr(54,6).c_str());}
if(length > 65){bval		=strtofloat(line.substr(60,6).c_str());}
if(length > 68){footnote		=line.substr(67,3);}
if (atom_name.find("H",0) != string::npos){
r=.31;
}else{
if (atom_name.find("C",0) != string::npos){
r=.73;
}else{
if (atom_name.find("O",0) != string::npos){
r=.66;
}else{
if (atom_name.find("P",0) != string::npos){
r=1.07;
}else{
if (atom_name.find("N",0) != string::npos){
r=0.71;
}else{
printf("Can't identify atom type: %s\n",atom_name.c_str());
}
}
}
}
}
}
