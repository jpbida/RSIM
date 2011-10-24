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
