// molecule.h -- molecule class and functions
#ifndef FRAGCLASS_H_
#define FRAGCLASS_H_
#include "vclass.h"
#include "cclass.h"
#include "matches.h"
#include "points.h"
#include <vector>
#include <string>
using namespace std;
class Fragclass
{
private:
public:
Fragclass();
Fragclass(string db_prefix);
vector<Vclass *> RRR1;
vector<Vclass *> RRY1;
vector<Vclass *> RYR1;
vector<Vclass *> RYY1;
vector<Vclass *> YRR1;
vector<Vclass *> YRY1;
vector<Vclass *> YYR1;
vector<Vclass *> YYY1;

vector<Vclass *> RRR2;
vector<Vclass *> RRY2;
vector<Vclass *> RYR2;
vector<Vclass *> RYY2;
vector<Vclass *> YRR2;
vector<Vclass *> YRY2;
vector<Vclass *> YYR2;
vector<Vclass *> YYY2;

vector<Cclass *> iRRR1;
vector<Cclass *> iRRY1;
vector<Cclass *> iRYR1;
vector<Cclass *> iRYY1;
vector<Cclass *> iYRR1;
vector<Cclass *> iYRY1;
vector<Cclass *> iYYR1;
vector<Cclass *> iYYY1;

vector<Cclass *> iRRR2;
vector<Cclass *> iRRY2;
vector<Cclass *> iRYR2;
vector<Cclass *> iRYY2;
vector<Cclass *> iYRR2;
vector<Cclass *> iYRY2;
vector<Cclass *> iYYR2;
vector<Cclass *> iYYY2;

//Return all matched moves for two given positions
int  matchedSingle(int type);
std::vector<int> allSingles(int type1);
vector<Matches> matchedFMoves(int type,double Lmin,double Lmax,Matrix * t,Vclass * o,double xf,double yf, double zf,int fixed,int pos1,int pos2);
vector<Matches> fixedMoves(int type1,int type2,double L,Matrix * t1,Matrix * t2,Vclass * o1,Vclass * o2,int pos1, int pos2);//Returns all pairs of fragments that hold the ends within a particular distance
void printMoves(
		int type1,
	       	int type2,
		double L,
		double r,
		double v,
		double z,
		Matrix * t1,
		Matrix * t2,
		Vclass * o1,
		Vclass * o2,
		double Lx,
		double Ly,
		double Lz,
		int p1,
		int p2
		);
vector<Matches> matchedMoves(
		int type1,
	       	int type2,
		double L,
		double r,
		double v,
		double z,
		Matrix * t1,
		Matrix * t2,
		Vclass * o1,
		Vclass * o2,
		double Lx,
		double Ly,
		double Lz,
		int time,
		int p1,
		int p2
		);
//Load a fragment library 
void loadFrags(string file,int type, std::vector<Vclass *> * oRRR1);
void loadClusters(string file,int type, std::vector<Cclass *> * oRRR1);
int searchfrag(vector<Vclass >,int id);
void sort();
};
#endif



