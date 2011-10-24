// molecule.h -- molecule class and functions
#ifndef ATOM_H_
#define ATOM_H_
#include <vector>
#include <string>

using namespace std;
class Atom
{
private:
public:
Atom();
Atom(int at_num, string at_name, float x1, float y1, float z1);
Atom(int at_num, string at_name, float x1, float y1, float z1,float bval);
Atom(string line);
//Parse file into the class
int atom_num;
string atom_name;
string alt_loc;
string insert_code;
float x;
float y;
float z;
float r;
float occ;
float bval;
string footnote;


};

#endif
