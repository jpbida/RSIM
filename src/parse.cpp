#include <istream>
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include <istream>
#include <sstream>
#include "polymere.h"
#include "parse.h"

void loadDoublets(std::vector<double> * dparams, string filedblt){
string field;
	std::ifstream input;
	input.open(filedblt.c_str());
	string line;
double x;
double y;
double z;
double theta2;
double phi2;
double typeRY;
while(getline(input, line)) 
{
//Read in file with each line being x \t y \t z \t group
	istringstream iss(line,istringstream::in);	
		getline(iss,field,'\t'); x=strtofloat(field);
		getline(iss,field,'\t'); y=strtofloat(field);
		getline(iss,field,'\t'); z=strtofloat(field);
		getline(iss,field,'\t'); theta2=strtofloat(field);
		getline(iss,field,'\t'); phi2=strtofloat(field);
		getline(iss,field,'\t'); typeRY=strtofloat(field);
dparams->push_back(x);
dparams->push_back(y);
dparams->push_back(z);
dparams->push_back(theta2);
dparams->push_back(phi2);
dparams->push_back(typeRY);
}
input.close();
}


double strtofloat(std::string what)
{
using namespace std;
	istringstream instr(what);
	float val;
	instr >> val;
	return val;
};
double strtodouble(std::string what)
{
using namespace std;
	istringstream instr(what);
	double val;
	instr >> val;
	return val;
};
int strtoint(std::string what)
{
using namespace std;
	istringstream instr(what);
	int val;
	instr >> val;
	return val;
};
char strtochar(std::string what)
{
using namespace std;
	istringstream instr(what);
	char val;
	instr >> val;
	return val;
};
double round_nplaces(double value, int to)
{
    double places = pow(10.0, to);
    return round(value * places) / places;
};

struct PairCompare
{
    bool operator() (const Polymere *  x, const Polymere * y) const {
//To Minimize a score       
 return x->score > y->score;
//To maximize a score
//return x.score < y.score;
    }
};
