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
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "vclass.h"
#include "fragclass.h"
#include <istream>
#include <sstream>
#include "graph.h"
#include <algorithm>
using namespace std;
double Graph::gen_id(Polymere * pol){
double dist=0;
int count = 0;
double bp1 = 0;
double bp2 = 0;
double mind=99999;
double xind=0;
double mdist;
//Check Base-Pair Distances
int j;
//cout << "Checking Base-Pairs" << endl;
for(int i=0; i<pol->mols.size(); i++){
//cout << "id: " << pol->mols[i]->id << endl;
if(pol->mols[i]->id >= 0){
if(pol->mols[i]->id > i){
j=pol->mols[i]->id;
mdist =pow((pol->mols[i]->x - pol->mols[j]->x),2);
mdist =mdist + pow((pol->mols[i]->y - pol->mols[j]->y),2);
mdist =mdist + pow((pol->mols[i]->z - pol->mols[j]->z),2);
//cout << i << ":" << j << ":" << mdist << endl;
if(mdist < mind){
mind=mdist;
}
if(mdist > xind){
xind=mdist;
}

}
}
}

if(mind < 64 || xind > 169){
dist=0;
}else{
mind=99999;
dist=0;
for(int i=0; i<pol->mols.size(); i++){
for(int j=(i+1); j<pol->mols.size(); j++){
count++;
//mol[i]->id contains the base-paired partner of the given molecule
mdist =pow((pol->mols[i]->x - pol->mols[j]->x),2);
mdist =mdist + pow((pol->mols[i]->y - pol->mols[j]->y),2);
mdist =mdist + pow((pol->mols[i]->z - pol->mols[j]->z),2);
dist=mdist+dist;
if(mdist < mind){
if(j>(i+3)){
mind=mdist;
}
} 
}
}
if(mind < 9){
dist=0;
}else{
dist=dist/count;
}
}
return dist;
}

struct PairCompare
{
    bool operator() (const Element * x, const Element * y) const {
        return x->child < y->child;
    }
};


void Graph::sort(){
using namespace std;
std::sort(elements.begin(), elements.end(), PairCompare ());
}

void Graph::print(){
using namespace std;
cout.precision(8);
cout << "Graph g {" << endl;
for(int w=0; w<elements.size();w++){
cout << "n" << elements[w]->parent << " -- n" << elements[w]->child << endl;
}
cout << "}" << endl;
}

void Graph::add(Element * a){
using namespace std;
//// ----- ------ ////
elements.push_back(new Element());
elements[(elements.size()-1)]=a;
}

double Graph::search(double child){
using namespace std;
Element * ele=new Element();
ele->child=child;
vector<Element *>::iterator low=lower_bound(elements.begin(), elements.end(), ele,PairCompare() );
vector<Element *>::iterator high=upper_bound(elements.begin(), elements.end(), ele,PairCompare() );
double parent=-1;
if(low!=high){
//    std::cout << "contains " << child <<  std::endl;
parent=elements[(low - elements.begin())]->parent;
}
  else{
 //   std::cout << "does not contain " << child <<  std::endl;
}
return parent;
}

