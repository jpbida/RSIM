// Montecarlo fractal kinetics simulator //
#include "polymere.h"
#include <string>
#include <istream>
#include <sstream>
#include <vector>
#include <iostream>
#include "parse.h"
#include "residue.h"
#include "atom.h"
#include "model.h"

Model::Model(vector<string> * lines,vector<int> * chain_start, vector<int> * chain_end, vector<int>* res_start, vector<int> * res_end, int start, int stop){
using namespace std;
//Create a vector of chains
for(int c=0; c < chain_start->size(); c++){
if(chain_start->at(c) >= start && chain_start->at(c) < stop){
chains.push_back(new Chain(lines,res_start,res_end,chain_start->at(c),chain_end->at(c)));
}
}

}

Model::Model(Polymere * pol,char chain){
using namespace std;
//Create a vector of chains
chains.push_back(new Chain(pol,chain));
}

void Model::print(){
for(int c=0; c<chains.size(); c++){
chains[c]->print();
}

}
