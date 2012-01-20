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
