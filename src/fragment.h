// molecule.h -- molecule class and functions

#ifndef FRAGMENT_H_
#define FRAGMENT_H_
#include <vector>
#include <string>

using namespace std;
class Fragment
{
private:
public:
Fragment();
Fragment(string line);
void print();
int frag_type;
char nt1;
char nt2;
char nt3;
double len_1;
double tor_1;
double len_2;
double len_3;
double len_4;
double tor_2;
double bond_0;
double bond_1;
double bond_2;
};


#endif



