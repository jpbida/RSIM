// molecule.h -- molecule class and functions
#ifndef RESIDUE_H_
#define RESIDUE_H_
#include <vector>
#include "atom.h"
#include <string>
#include "polymere.h"

using namespace std;
class Residue
{
private:

public:
const static char  atomnames1[27][5];
const static char  atomnames2[26][5];
const static char  atomnames3[23][5];
const static char  atomnames4[24][5];
void print();
void print(std::string file);
Residue();
Residue(Polymere * pol, char chain,int res_num);
Residue(vector<string> * lines, int start,int end);
Residue(const char* file, int seqnum1,char chain1);
string getElement(string,int);
int isBackBone(int a);
vector<float> getAtom(string);
//Parse file into the class
string res_name;
string chain;
int res_num;
vector<Atom *> atoms;
//Vectors storing the atoms names mapping from polymere to PDB

};

#endif
