// molecule.h -- molecule class and functions
#ifndef PDB_H_
#define PDB_H_
#include <vector>
#include "model.h"
#include "polymere.h"
#include <string>
using namespace std;
class Pdb
{
private:
public:
Pdb();
//Parse file into the class
Pdb(Polymere * pol, char chain);
Pdb(string file);
void ModelsToPols(std::vector<Polymere *> * pols);
void splits(string file,int tmods);
void print();
void clear();
void write(char * file);
void write(string file);
//Parses an ATOM line
void standard(const char* file);
string getElement(string,int);
//Data structures
vector <Model *> models;
//	vecotr<Chain *> chains;
//		vector<Residues *> residues;
//			vector<Atom *> atoms;

};

#endif
