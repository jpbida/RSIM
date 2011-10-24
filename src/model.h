// molecule.h -- molecule class and functions
#ifndef MODEL_H_
#define MODEL_H_
#include <vector>
#include <string>
#include "chain.h"
#include "polymere.h"

using namespace std;
class Model
{
private:
public:
Model();
Model(Polymere * pol,char chain);
Model(vector<string> * lines,vector<int> * chain_start, vector<int> * chain_end, vector<int>* res_start, vector<int> * res_end, int start, int stop);
void print();
//Model(string line);
//Parse file into the class
int model_num;
vector<Chain *> chains;
};

#endif
