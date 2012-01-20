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
#include <algorithm>
#include <ctime>
#include <fstream>
#include <iostream>
#include "pdb.h"
#include "base.h"
#include "fragment.h"
#include "polymere.h"
#include "fragclass.h"
#include "fragatoms.h"
#include "fragall.h"
#include "vclass.h"
#include <string>
#include <cmath>
#include <vector>
#include <istream>
#include <sstream>
#include "matches.h"
#include "graph.h"
#include "element.h"
#include "parse.h"
#include "options.h"
#include "fabase.h"

int main(int argc, char **argv)
{
using namespace std;
Options * opt=new Options(argc,argv);
//Max structures
Pdb * pdb = new Pdb(opt->starting_struct);
int tots=pdb->models.size();
printf("Total Number of Models: %d\n",tots);

string pdb_map=opt->prefix;
//Creating vector of polymeres
pdb->standard(pdb_map.c_str());
pdb->write("converted.pdb");
pdb->clear();
delete pdb;
}//end main
