#ifndef CHAIN_H_
#define CHAIN_H_
#include <vector>
#include <string>
#include "residue.h"
#include "polymere.h"

using namespace std;
class Chain
{
private:
public:
Chain();
Chain(Polymere * pol, char chain);
Chain(std::vector<string> * lines, std::vector<int> * res_start, std::vector<int> * res_end,int start, int end);
void print();
//type RNA, DNA, Protein, Hybrid
string type;
std::vector<Residue *> residues;


};

#endif
