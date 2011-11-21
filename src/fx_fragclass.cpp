// Montecarlo fractal kinetics simulator //
#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cmath>
#include "vclass.h"
#include "cclass.h"
#include "fragclass.h"
#include <istream>
#include <sstream>


struct PairCompareV
{
    bool operator() (const Vclass * x, const Vclass * y) const {
        return x->id < y->id;
    }
};

struct PairCompareV2
{
    bool operator() (Vclass x, Vclass  y) const {
        return x.id < y.id;
    }
};

Fragclass::Fragclass(string db_prefix){
//Create the vector of Vclasses and Cclasses 
for(int i=0; i < 1800; i++){
RRR1.push_back(new Vclass());
RRR2.push_back(new Vclass());
RRY1.push_back(new Vclass());
RRY2.push_back(new Vclass());
RYR1.push_back(new Vclass());
RYR2.push_back(new Vclass());
RYY1.push_back(new Vclass());
RYY2.push_back(new Vclass());
YRR1.push_back(new Vclass());
YRR2.push_back(new Vclass());
YRY1.push_back(new Vclass());
YRY2.push_back(new Vclass());
YYR1.push_back(new Vclass());
YYR2.push_back(new Vclass());
YYY1.push_back(new Vclass());
YYY2.push_back(new Vclass());
}

for(int i=0; i < 800; i++){
iRRR1.push_back(new Cclass());
iRRR2.push_back(new Cclass());
iRRY1.push_back(new Cclass());
iRRY2.push_back(new Cclass());
iRYR1.push_back(new Cclass());
iRYR2.push_back(new Cclass());
iRYY1.push_back(new Cclass());
iRYY2.push_back(new Cclass());
iYRR1.push_back(new Cclass());
iYRR2.push_back(new Cclass());
iYRY1.push_back(new Cclass());
iYRY2.push_back(new Cclass());
iYYR1.push_back(new Cclass());
iYYR2.push_back(new Cclass());
iYYY1.push_back(new Cclass());
iYYY2.push_back(new Cclass());
}

string file0   = db_prefix;file0.append("fRRR.vect");printf("Reading lib: %s\n" ,file0.c_str() );
string file1   = db_prefix;file1.append("fRRY.vect");printf("Reading lib: %s\n" ,file1.c_str()   );  
string file2   = db_prefix;file2.append("fRYR.vect");printf("Reading lib: %s\n" ,file2.c_str()   );  
string file3   = db_prefix;file3.append("fRYY.vect");printf("Reading lib: %s\n" ,file3.c_str()   );  
string file4   = db_prefix;file4.append("fYRR.vect");printf("Reading lib: %s\n" ,file4.c_str()   );  
string file5   = db_prefix;file5.append("fYRY.vect");printf("Reading lib: %s\n" ,file5.c_str()   );  
string file6   = db_prefix;file6.append("fYYR.vect");printf("Reading lib: %s\n" ,file6.c_str()   );  
string file7   = db_prefix;file7.append("fYYY.vect");printf("Reading lib: %s\n" ,file7.c_str()   );  
string file0r  = db_prefix;file0r.append("rRRR.vect");printf("Reading lib: %s\n" ,file0r.c_str()  );  
string file1r  = db_prefix;file1r.append("rRRY.vect");printf("Reading lib: %s\n" ,file1r.c_str()  );
string file2r  = db_prefix;file2r.append("rRYR.vect");printf("Reading lib: %s\n" ,file2r.c_str()  );
string file3r  = db_prefix;file3r.append("rRYY.vect");printf("Reading lib: %s\n" ,file3r.c_str()  );
string file4r  = db_prefix;file4r.append("rYRR.vect");printf("Reading lib: %s\n" ,file4r.c_str()  );
string file5r  = db_prefix;file5r.append("rYRY.vect");printf("Reading lib: %s\n" ,file5r.c_str()  );
string file6r  = db_prefix;file6r.append("rYYR.vect");printf("Reading lib: %s\n" ,file6r.c_str()  );
string file7r  = db_prefix;file7r.append("rYYY.vect");printf("Reading lib: %s\n" ,file7r.c_str()  );
string file1c  = db_prefix;file1c.append("cfRRR.vect");printf("Reading lib: %s\n",file1c.c_str()  );
string file2c  = db_prefix;file2c.append("cfRRY.vect");printf("Reading lib: %s\n",file2c.c_str()  );
string file3c  = db_prefix;file3c.append("cfRYR.vect");printf("Reading lib: %s\n",file3c.c_str()  );
string file4c  = db_prefix;file4c.append("cfRYY.vect");printf("Reading lib: %s\n",file4c.c_str()  );
string file5c  = db_prefix;file5c.append("cfYRR.vect");printf("Reading lib: %s\n",file5c.c_str()  );
string file6c  = db_prefix;file6c.append("cfYRY.vect");printf("Reading lib: %s\n",file6c.c_str()  );
string file7c  = db_prefix;file7c.append("cfYYR.vect");printf("Reading lib: %s\n",file7c.c_str()  );
string file8c  = db_prefix;file8c.append("cfYYY.vect");printf("Reading lib: %s\n",file8c.c_str()  );
string file1cr = db_prefix;file1cr.append("crRRR.vect");printf("Reading lib: %s\n",file1cr.c_str() );
string file2cr = db_prefix;file2cr.append("crRRY.vect");printf("Reading lib: %s\n",file2cr.c_str() );
string file3cr = db_prefix;file3cr.append("crRYR.vect");printf("Reading lib: %s\n",file3cr.c_str() );
string file4cr = db_prefix;file4cr.append("crRYY.vect");printf("Reading lib: %s\n",file4cr.c_str() );
string file5cr = db_prefix;file5cr.append("crYRR.vect");printf("Reading lib: %s\n",file5cr.c_str() );
string file6cr = db_prefix;file6cr.append("crYRY.vect");printf("Reading lib: %s\n",file6cr.c_str() );
string file7cr = db_prefix;file7cr.append("crYYR.vect");printf("Reading lib: %s\n",file7cr.c_str() );
string file8cr = db_prefix;file8cr.append("crYYY.vect");printf("Reading lib: %s\n",file8cr.c_str() );

loadFrags(file0,0,&RRR1);
loadFrags(file0r,0,&RRR2);
loadFrags(file1,1,&RRY1);
loadFrags(file1r,1,&RRY2);
loadFrags(file2,2,&RYR1);
loadFrags(file2r,2,&RYR2);
loadFrags(file3,3,&RYY1);
loadFrags(file3r,3,&RYY2);
loadFrags(file4,4,&YRR1);
loadFrags(file4r,4,&YRR2);
loadFrags(file5,5,&YRY1);
loadFrags(file5r,5,&YRY2);
loadFrags(file6,6,&YYR1);
loadFrags(file6r,6,&YYR2);
loadFrags(file7,7,&YYY1);
loadFrags(file7r,7,&YYY2);

loadClusters(file1c,7,&iRRR1);
loadClusters(file2c,7,&iRRY1);
loadClusters(file3c,7,&iRYR1);
loadClusters(file4c,7,&iRYY1);
loadClusters(file5c,7,&iYRR1);
loadClusters(file6c,7,&iYRY1);
loadClusters(file7c,7,&iYYR1);
loadClusters(file8c,7,&iYYY1);
loadClusters(file1cr,7,&iRRR2);
loadClusters(file2cr,7,&iRRY2);
loadClusters(file3cr,7,&iRYR2);
loadClusters(file4cr,7,&iRYY2);
loadClusters(file5cr,7,&iYRR2);
loadClusters(file6cr,7,&iYRY2);
loadClusters(file7cr,7,&iYYR2);
loadClusters(file8cr,7,&iYYY2);
}

Fragclass::Fragclass(){
//Create the vector of Vclasses and Cclasses 
for(int i=0; i < 1800; i++){
RRR1.push_back(new Vclass());
RRR2.push_back(new Vclass());
RRY1.push_back(new Vclass());
RRY2.push_back(new Vclass());
RYR1.push_back(new Vclass());
RYR2.push_back(new Vclass());
RYY1.push_back(new Vclass());
RYY2.push_back(new Vclass());
YRR1.push_back(new Vclass());
YRR2.push_back(new Vclass());
YRY1.push_back(new Vclass());
YRY2.push_back(new Vclass());
YYR1.push_back(new Vclass());
YYR2.push_back(new Vclass());
YYY1.push_back(new Vclass());
YYY2.push_back(new Vclass());
}

for(int i=0; i < 800; i++){
iRRR1.push_back(new Cclass());
iRRR2.push_back(new Cclass());
iRRY1.push_back(new Cclass());
iRRY2.push_back(new Cclass());
iRYR1.push_back(new Cclass());
iRYR2.push_back(new Cclass());
iRYY1.push_back(new Cclass());
iRYY2.push_back(new Cclass());
iYRR1.push_back(new Cclass());
iYRR2.push_back(new Cclass());
iYRY1.push_back(new Cclass());
iYRY2.push_back(new Cclass());
iYYR1.push_back(new Cclass());
iYYR2.push_back(new Cclass());
iYYY1.push_back(new Cclass());
iYYY2.push_back(new Cclass());
}


string file0("../libs/fRRR.vect");
string file1("../libs/fRRY.vect");
string file2("../libs/fRYR.vect");
string file3("../libs/fRYY.vect");
string file4("../libs/fYRR.vect");
string file5("../libs/fYRY.vect");
string file6("../libs/fYYR.vect");
string file7("../libs/fYYY.vect");

string file0r(  "../libs/rRRR.vect");
string file1r(  "../libs/rRRY.vect");
string file2r(  "../libs/rRYR.vect");
string file3r(  "../libs/rRYY.vect");
string file4r(  "../libs/rYRR.vect");
string file5r(  "../libs/rYRY.vect");
string file6r(  "../libs/rYYR.vect");
string file7r(  "../libs/rYYY.vect");

string file1c ( "../libs/cfRRR.vect");
string file2c ( "../libs/cfRRY.vect");
string file3c ( "../libs/cfRYR.vect");
string file4c ( "../libs/cfRYY.vect");
string file5c ( "../libs/cfYRR.vect");
string file6c ( "../libs/cfYRY.vect");
string file7c ( "../libs/cfYYR.vect");
string file8c ( "../libs/cfYYY.vect");

string file1cr ( "../libs/crRRR.vect");
string file2cr ( "../libs/crRRY.vect");
string file3cr ( "../libs/crRYR.vect");
string file4cr ( "../libs/crRYY.vect");
string file5cr ( "../libs/crYRR.vect");
string file6cr ( "../libs/crYRY.vect");
string file7cr ( "../libs/crYYR.vect");
string file8cr ( "../libs/crYYY.vect");

loadFrags(file0,0,&RRR1);
loadFrags(file0r,0,&RRR2);
loadFrags(file1,1,&RRY1);
loadFrags(file1r,1,&RRY2);
loadFrags(file2,2,&RYR1);
loadFrags(file2r,2,&RYR2);
loadFrags(file3,3,&RYY1);
loadFrags(file3r,3,&RYY2);
loadFrags(file4,4,&YRR1);
loadFrags(file4r,4,&YRR2);
loadFrags(file5,5,&YRY1);
loadFrags(file5r,5,&YRY2);
loadFrags(file6,6,&YYR1);
loadFrags(file6r,6,&YYR2);
loadFrags(file7,7,&YYY1);
loadFrags(file7r,7,&YYY2);

loadClusters(file1c,7,&iRRR1);
loadClusters(file2c,7,&iRRY1);
loadClusters(file3c,7,&iRYR1);
loadClusters(file4c,7,&iRYY1);
loadClusters(file5c,7,&iYRR1);
loadClusters(file6c,7,&iYRY1);
loadClusters(file7c,7,&iYYR1);
loadClusters(file8c,7,&iYYY1);
loadClusters(file1cr,7,&iRRR2);
loadClusters(file2cr,7,&iRRY2);
loadClusters(file3cr,7,&iRYR2);
loadClusters(file4cr,7,&iRYY2);
loadClusters(file5cr,7,&iYRR2);
loadClusters(file6cr,7,&iYRY2);
loadClusters(file7cr,7,&iYYR2);
loadClusters(file8cr,7,&iYYY2);
}

void Fragclass::sort(){
using namespace std;
std::sort(RRR1.begin(), RRR1.end(), PairCompareV ());
std::sort(RRY1.begin(), RRY1.end(), PairCompareV ());
std::sort(RYR1.begin(), RYR1.end(), PairCompareV ());
std::sort(RYY1.begin(), RYY1.end(), PairCompareV ());
std::sort(YRR1.begin(), YRR1.end(), PairCompareV ());
std::sort(YRY1.begin(), YRY1.end(), PairCompareV ());
std::sort(YYR1.begin(), YYR1.end(), PairCompareV ());
std::sort(YYY1.begin(), YYY1.end(), PairCompareV ());

std::sort(RRR2.begin(), RRR2.end(), PairCompareV ());
std::sort(RRY2.begin(), RRY2.end(), PairCompareV ());
std::sort(RYR2.begin(), RYR2.end(), PairCompareV ());
std::sort(RYY2.begin(), RYY2.end(), PairCompareV ());
std::sort(YRR2.begin(), YRR2.end(), PairCompareV ());
std::sort(YRY2.begin(), YRY2.end(), PairCompareV ());
std::sort(YYR2.begin(), YYR2.end(), PairCompareV ());
std::sort(YYY2.begin(), YYY2.end(), PairCompareV ());
}

int Fragclass::searchfrag(vector<Vclass> frags,int id){
using namespace std;
Vclass * ele=new Vclass();
ele->id=id;
Vclass ele2=*ele;
vector<Vclass>::iterator low=lower_bound(frags.begin(), frags.end(), ele2,PairCompareV2() );
vector<Vclass>::iterator high=upper_bound(frags.begin(), frags.end(), ele2,PairCompareV2() );
int parent=-1;
if(low!=high){
//    std::cout << "contains " << child <<  std::endl;
parent=(low - frags.begin());
}
  else{
 //   std::cout << "does not contain " << child <<  std::endl;
}
return parent;
}

double strtodouble1(std::string what)
{
using namespace std;
	istringstream instr(what);
	double val;
	instr >> val;
	return val;
}


int strtoint1(std::string what)
{
using namespace std;
	istringstream instr(what);
	int val;	
	instr >> val;
	return val;
}

char strtochar1(std::string what)
{
using namespace std;
	istringstream instr(what);
	char val;
	instr >> val;
	return val;
}

using namespace std;
const double PI = 3.14159265358979323846;

//Fragclass has cluster classes 
void Fragclass::loadClusters(string file,int type, std::vector<Cclass*> *oRRR1){
	using namespace std;
double x,y,z,r;
int id;
int gid;
int t0=0;
int newg=0;
int ind=0;
string field;
	std::ifstream input;
	input.open(file.c_str());
	string line;
while(getline(input, line)) 
{
//Read in file with each line being x \t y \t z \t group
istringstream iss(line,istringstream::in);	
	getline(iss,field,'\t');
id=strtoint1(field);
	getline(iss,field,'\t');
x=strtodouble1(field);
	getline(iss,field,'\t');
y=strtodouble1(field);
	getline(iss,field,'\t');
z=strtodouble1(field);
	getline(iss,field,'\t');
r=strtodouble1(field);
	getline(iss,field,'\t');
gid=strtoint1(field);
//If gid alread exists then just add id to members
newg=0;
for(int j=0; j<t0; j++){
if(oRRR1->at(j)->id==gid){
newg=1;
ind=j;
}
}
if(newg==1){
oRRR1->at(ind)->members.push_back(id);
}else{
oRRR1->at(t0)->x=x;
oRRR1->at(t0)->y=y;
oRRR1->at(t0)->z=z;
oRRR1->at(t0)->r=r;
oRRR1->at(t0)->members.push_back(id);
oRRR1->at(t0++)->id=gid;
}
}
cout << "Total Clusters: " << t0 << endl;
}

void Fragclass::loadFrags(string file,int type, vector<Vclass *> *oRRR1){
//Given a file with 2D coordinates generate from secondary structure it loads positions into polymere
using namespace std;
double x,rx,vx;
double y,ry,vy;
double z,rz,vz;
int id;
int gid;
int t0=0;
string field;
	std::ifstream input;
	input.open(file.c_str());
	string line;
while(getline(input, line)) 
{
//Read in file with each line being x \t y \t z \t group
	istringstream iss(line,istringstream::in);	
		getline(iss,field,'\t');
	id=strtoint1(field);
		getline(iss,field,'\t');
	x=strtodouble1(field);
		getline(iss,field,'\t');
	y=strtodouble1(field);
		getline(iss,field,'\t');
	z=strtodouble1(field);
		getline(iss,field,'\t');
	vx=strtodouble1(field);
		getline(iss,field,'\t');
	vy=strtodouble1(field);
		getline(iss,field,'\t');
	vz=strtodouble1(field);
		getline(iss,field,'\t');
	rx=strtodouble1(field);
		getline(iss,field,'\t');
	ry=strtodouble1(field);
		getline(iss,field,'\t');
	rz=strtodouble1(field);
	oRRR1->at(t0)->x=x;
	oRRR1->at(t0)->y=y;
	oRRR1->at(t0)->z=z;
	oRRR1->at(t0)->vx=vx;
	oRRR1->at(t0)->vy=vy;
	oRRR1->at(t0)->vz=vz;
	oRRR1->at(t0)->rx=rx;
	oRRR1->at(t0)->ry=ry;
	oRRR1->at(t0)->rz=rz;
	oRRR1->at(t0)->id=id;
	oRRR1->at(t0++)->group=id;
}
}

vector<Cclass> copyCclass(vector<Cclass*> v){
int tot = v.size();
vector<Cclass> vout;

for(int j=0; j<tot; j++){
Cclass v1 = *v[j];
vout.push_back(v1);
}
return vout;
}

vector<Vclass> copyVclass(vector<Vclass*> v){
int tot = v.size();
vector<Vclass> vout;

for(int j=0; j<tot; j++){
Vclass v1 = *v[j];
vout.push_back(v1);
}
return vout;
}

std::vector<int> Fragclass::allSingles(int type1){
//Returns the id of a fragment with the given type 
vector<Vclass> fgs1;
if(type1==0){fgs1=copyVclass(YYY1);}
if(type1==1){fgs1=copyVclass(RYY1);}
if(type1==2){fgs1=copyVclass(YRY1);}
if(type1==3){fgs1=copyVclass(RRY1);}
if(type1==4){fgs1=copyVclass(YYR1);}
if(type1==5){fgs1=copyVclass(RYR1);}
if(type1==6){fgs1=copyVclass(YRR1);}
if(type1==7){fgs1=copyVclass(RRR1);}
std::vector<int> idlist;
for(int j=0; j<fgs1.size(); j++){
idlist.push_back(fgs1[j].id);   
}
return idlist;
}

int Fragclass::matchedSingle(int type1){
//Returns the id of a fragment with the given type 
vector<Vclass> fgs1;
if(type1==0){fgs1=copyVclass(YYY1);}
if(type1==1){fgs1=copyVclass(RYY1);}
if(type1==2){fgs1=copyVclass(YRY1);}
if(type1==3){fgs1=copyVclass(RRY1);}
if(type1==4){fgs1=copyVclass(YYR1);}
if(type1==5){fgs1=copyVclass(RYR1);}
if(type1==6){fgs1=copyVclass(YRR1);}
if(type1==7){fgs1=copyVclass(RRR1);}

int kr=(std::rand() % fgs1.size());
int id=fgs1[kr].id;   
return id;
}

vector<Matches> Fragclass::matchedFMoves(int type,double Lmin,double Lmax,Matrix * t,Vclass * o,double xf,double yf, double zf,int fixed,int pos1,int pos2){
//Finds a matching fragment for a fixed fragment
//Fixed tells you which position is fixed 1 = position 1 2=position 2
//The fixed position is the base with which you want the 3rd Nt from the
//free fragment to be next to.

//Set the fixed fragment vector to the cluster class with radius = 1
//Create a cluster for the single fragment that is fixed
//You don't need to translate or rotate this cluster
Cclass cc_fixed;
cc_fixed.x=xf;
cc_fixed.y=yf;
cc_fixed.z=zf;
cc_fixed.r=(Lmax-Lmin)/2;
cc_fixed.members.push_back(0); 

Vclass vc_fixed;
vc_fixed.x=xf;
vc_fixed.y=yf;
vc_fixed.z=zf;

//These don't matter because we are only filtering for distance
vc_fixed.vx=0;
vc_fixed.vy=0;
vc_fixed.vz=1;

vc_fixed.rx=1;
vc_fixed.ry=0;
vc_fixed.rz=0;

vc_fixed.group=0;
vc_fixed.id=0;

vector<Matches > matches;
vector<Vclass> fgs1;
vector<Vclass> fgs2;
vector<Cclass> fgi1;
vector<Cclass> fgi2;

if(fixed==1){
fgi1.push_back(cc_fixed);
fgs1.push_back(vc_fixed);
//Cluters
if(type==0){fgi2=copyCclass(iYYY2);}
if(type==1){fgi2=copyCclass(iRYY2);}
if(type==2){fgi2=copyCclass(iYRY2);}
if(type==3){fgi2=copyCclass(iRRY2);}
if(type==4){fgi2=copyCclass(iYYR2);}
if(type==5){fgi2=copyCclass(iRYR2);}
if(type==6){fgi2=copyCclass(iYRR2);}
if(type==7){fgi2=copyCclass(iRRR2);}
//Fragment Types
if(type==0){fgs2=copyVclass(YYY2);}
if(type==1){fgs2=copyVclass(RYY2);}
if(type==2){fgs2=copyVclass(YRY2);}
if(type==3){fgs2=copyVclass(RRY2);}
if(type==4){fgs2=copyVclass(YYR2);}
if(type==5){fgs2=copyVclass(RYR2);}
if(type==6){fgs2=copyVclass(YRR2);}
if(type==7){fgs2=copyVclass(RRR2);}
}else{
if(fixed==2){
//Cluters
fgi2.push_back(cc_fixed);
fgs2.push_back(vc_fixed);
if(type==0){fgi1=copyCclass(iYYY1);}
if(type==1){fgi1=copyCclass(iRYY1);}
if(type==2){fgi1=copyCclass(iYRY1);}
if(type==3){fgi1=copyCclass(iRRY1);}
if(type==4){fgi1=copyCclass(iYYR1);}
if(type==5){fgi1=copyCclass(iRYR1);}
if(type==6){fgi1=copyCclass(iYRR1);}
if(type==7){fgi1=copyCclass(iRRR1);}
//Fragment Types
if(type==0){fgs1=copyVclass(YYY1);}
if(type==1){fgs1=copyVclass(RYY1);}
if(type==2){fgs1=copyVclass(YRY1);}
if(type==3){fgs1=copyVclass(RRY1);}
if(type==4){fgs1=copyVclass(YYR1);}
if(type==5){fgs1=copyVclass(RYR1);}
if(type==6){fgs1=copyVclass(YRR1);}
if(type==7){fgs1=copyVclass(RRR1);}
}
}
//Translate each to there given positions
int i = 0;
int j = 0;
int tot1=fgi1.size();
int tot2=fgi2.size();
int tot3=fgs1.size();
int tot4=fgs2.size();
//Transform to global coordinates
if(fixed==1){
while(j < tot2){
fgi2[j].transform(t);
fgi2[j++].translate(o);
}
while(i < tot4){
fgs2[i].transform(t);
fgs2[i++].translate(o);
}
}else{
if(fixed==2){
while(j < tot1){
fgi1[j].transform(t);
fgi1[j++].translate(o);
}
while(i < tot3){
fgs1[i].transform(t);
fgs1[i++].translate(o);
}
}
}
i=0;
j=0;
double dist,x1,y1,z1,x2,z2,y2,y3,x3,z3,a,dr,dv,dz,xi,yi,zi,xj,yj,zj,angout,x4,y4,z4,r1,r2;
while(i < tot1){
	j=0;
while(j < tot2){
x1=fgi1[i].x;
y1=fgi1[i].y;
z1=fgi1[i].z;

x2=fgi2[j].x;
y2=fgi2[j].y;
z2=fgi2[j].z;
 
x3=x2-x1;
y3=y2-y1;
z3=z2-z1;

dist=pow(x3,2)+pow(y3,2)+pow(z3,2);
//Radius of fragment is stored in rx
double mid=Lmin+(Lmax-Lmin)/2;
a = abs(sqrt(dist)-mid); //L is the desired distance apart
//a = abs(dist - L);
//Look at all positions < 4A difference
//Calculate all pairwise distances for fragments in the two groups
//Set distance threshold 
if(a < (2*(fgi1[i].r+fgi2[j].r))){
//Check pairwise distances of members in the group
//Get fgi1.members and fgi2.members 
int i1,i2;
for(int f1=0; f1<fgi1[i].members.size(); f1++){
i1=fgi1[i].members[f1];
for(int f2=0; f2<fgi2[j].members.size(); f2++){
i2=fgi2[j].members[f2];
x1=fgs1[i1].x;
y1=fgs1[i1].y;
z1=fgs1[i1].z;

x2=fgs2[i2].x;
y2=fgs2[i2].y;
z2=fgs2[i2].z;
 
x3=x2-x1;
y3=y2-y1;
z3=z2-z1;
dist=pow(x3,2)+pow(y3,2)+pow(z3,2);
//Radius of fragment is stored in rx
a = sqrt(dist);
if(a < Lmax && a > Lmin){
//Here we don't care about the angular dependence because we just replace things further sequentially 
//down the two strands
Matches * mat=new Matches();
mat->frag1=fgs1[i1].id;
mat->frag2=fgs2[i2].id;
mat->score = a;
mat->pos1 = pos1;
mat->pos2 = pos2;
matches.push_back(*mat);
delete mat;
}
}	
	}
}
j++;
}
i++;
}
return matches;
}

vector<Matches> Fragclass::fixedMoves(int type1,int type2,double L,Matrix * t1,Matrix * t2,Vclass * o1,Vclass * o2,int pos1, int pos2){
//Transform frag1 and frag2 into same coordinate system with t1 and t2
vector<Matches > matches;
vector<Vclass> fgs1;
vector<Vclass> fgs2;
vector<Cclass> fgi1;
vector<Cclass> fgi2;
//Cluters
if(type1==0){fgi1=copyCclass(iYYY1);}
if(type1==1){fgi1=copyCclass(iRYY1);}
if(type1==2){fgi1=copyCclass(iYRY1);}
if(type1==3){fgi1=copyCclass(iRRY1);}
if(type1==4){fgi1=copyCclass(iYYR1);}
if(type1==5){fgi1=copyCclass(iRYR1);}
if(type1==6){fgi1=copyCclass(iYRR1);}
if(type1==7){fgi1=copyCclass(iRRR1);}
                                
if(type2==0){fgi2=copyCclass(iYYY2);}
if(type2==1){fgi2=copyCclass(iRYY2);}
if(type2==2){fgi2=copyCclass(iYRY2);}
if(type2==3){fgi2=copyCclass(iRRY2);}
if(type2==4){fgi2=copyCclass(iYYR2);}
if(type2==5){fgi2=copyCclass(iRYR2);}
if(type2==6){fgi2=copyCclass(iYRR2);}
if(type2==7){fgi2=copyCclass(iRRR2);}
//Fragment Types
if(type1==0){fgs1=copyVclass(YYY1);}
if(type1==1){fgs1=copyVclass(RYY1);}
if(type1==2){fgs1=copyVclass(YRY1);}
if(type1==3){fgs1=copyVclass(RRY1);}
if(type1==4){fgs1=copyVclass(YYR1);}
if(type1==5){fgs1=copyVclass(RYR1);}
if(type1==6){fgs1=copyVclass(YRR1);}
if(type1==7){fgs1=copyVclass(RRR1);}
                               
if(type2==0){fgs2=copyVclass(YYY2);}
if(type2==1){fgs2=copyVclass(RYY2);}
if(type2==2){fgs2=copyVclass(YRY2);}
if(type2==3){fgs2=copyVclass(RRY2);}
if(type2==4){fgs2=copyVclass(YYR2);}
if(type2==5){fgs2=copyVclass(RYR2);}
if(type2==6){fgs2=copyVclass(YRR2);}
if(type2==7){fgs2=copyVclass(RRR2);}
//Translate each to there given positions
int i = 0;
int j = 0;
int tot1=fgi1.size();
int tot2=fgi2.size();
//Transform to global coordinates
while(i < tot1){fgi1[i++].transform(t1);}
while(j < tot2){fgi2[j++].transform(t2);}
//Translate only x,y,z  
i=0;
j=0;
while(i < tot1){fgi1[i++].translate(o1);}
while(j < tot2){fgi2[j++].translate(o2);}
//Determine cross product of rv,vv
i=0;
j=0;
int tot3=fgs1.size();
int tot4=fgs2.size();
//Transform to global coordinates
while(i < tot3){fgs1[i++].transform(t1);}
while(j < tot4){fgs2[j++].transform(t2);}
//Translate only x,y,z  
i=0;
j=0;
while(i < tot3){fgs1[i++].translate(o1);}
while(j < tot4){fgs2[j++].translate(o2);}
/*
//Print out aligned fragments
i=0;
j=0;
while(i < tot3){cout << fgs1[i].x << "," << fgs1[i].y << "," << fgs1[i++].z << "," << 1 << endl;}
while(j < tot4){cout << fgs2[j].x << "," << fgs2[j].y << "," << fgs2[j++].z << "," << 2 << endl;}
*/
i=0;
j=0;
double dist,x1,y1,z1,x2,z2,y2,y3,x3,z3,a,dr,dv,dz,xi,yi,zi,xj,yj,zj,angout,x4,y4,z4,r1,r2;
//Print Out Positions of fragments
/*
for(int x=0; x<tot3; x++){

}
*/


while(i < tot1){
	j=0; 
while(j < tot2){
x1=fgi1[i].x;
y1=fgi1[i].y;
z1=fgi1[i].z;

x2=fgi2[j].x;
y2=fgi2[j].y;
z2=fgi2[j].z;
 
x3=x2-x1;
y3=y2-y1;
z3=z2-z1;

dist=pow(x3,2)+pow(y3,2)+pow(z3,2);
//Radius of fragment is stored in rx
a = abs(sqrt(dist)-sqrt(L));
//a = abs(dist - L);
//Look at all positions < 4A difference
//Calculate all pairwise distances for fragments in the two groups
//Set distance threshold 
if(a < (2*(fgi1[i].r+fgi2[j].r))){
//Check pairwise distances of members in the group
//Get fgi1.members and fgi2.members 
int i1,i2;
for(int f1=0; f1<fgi1[i].members.size(); f1++){
i1=fgi1[i].members[f1];
for(int f2=0; f2<fgi2[j].members.size(); f2++){
//int i1=searchfrag(fgs1,fgi1[i].members[f1]);
//int i2=searchfrag(fgs2,fgi2[j].members[f2]);
i2=fgi2[j].members[f2];
//cout << "i1: " << i1 << " i2: " << i2 << " type: " << type1 << " type2: " << type2 << " f1: " << f1 << " f2: " << f2 << " i: " << i << " j: " << j << "x: " <<  endl;
x1=fgs1[i1].x;
y1=fgs1[i1].y;
z1=fgs1[i1].z;

x2=fgs2[i2].x;
y2=fgs2[i2].y;
z2=fgs2[i2].z;
 
x3=x2-x1;
y3=y2-y1;
z3=z2-z1;
dist=pow(x3,2)+pow(y3,2)+pow(z3,2);
//Radius of fragment is stored in rx
a = abs(sqrt(dist)-sqrt(L));
if(a < 1){
//convert all coordinate systems to pos1 bottom coordinates
//If the structure is aligned to pos1-1 then the coordinates are already in line 
Matches * mat=new Matches();
mat->frag1=fgs1[i1].id;
mat->frag2=fgs2[i2].id;
mat->score = a;
mat->pos1 = pos1;
mat->pos2 = pos2;
matches.push_back(*mat);
delete mat;
}
}	
	}
}
j++;
}
i++;
}
return matches;
}

void Fragclass::printMoves(int type1,int type2,double L,double odr, double odv, double odz, Matrix * t1,Matrix * t2,Vclass * o1,Vclass * o2,double Lx, double Ly, double Lz,int pos1, int pos2){
//Transform frag1 and frag2 into same coordinate system with t1 and t2
vector<Matches > matches;
vector<Vclass> fgs1;
vector<Vclass> fgs2;
vector<Cclass> fgi1;
vector<Cclass> fgi2;
//Cluters
if(type1==0){fgi1=copyCclass(iYYY1);}
if(type1==1){fgi1=copyCclass(iRYY1);}
if(type1==2){fgi1=copyCclass(iYRY1);}
if(type1==3){fgi1=copyCclass(iRRY1);}
if(type1==4){fgi1=copyCclass(iYYR1);}
if(type1==5){fgi1=copyCclass(iRYR1);}
if(type1==6){fgi1=copyCclass(iYRR1);}
if(type1==7){fgi1=copyCclass(iRRR1);}
                                
if(type2==0){fgi2=copyCclass(iYYY2);}
if(type2==1){fgi2=copyCclass(iRYY2);}
if(type2==2){fgi2=copyCclass(iYRY2);}
if(type2==3){fgi2=copyCclass(iRRY2);}
if(type2==4){fgi2=copyCclass(iYYR2);}
if(type2==5){fgi2=copyCclass(iRYR2);}
if(type2==6){fgi2=copyCclass(iYRR2);}
if(type2==7){fgi2=copyCclass(iRRR2);}
//Fragment Types
if(type1==0){fgs1=copyVclass(YYY1);}
if(type1==1){fgs1=copyVclass(RYY1);}
if(type1==2){fgs1=copyVclass(YRY1);}
if(type1==3){fgs1=copyVclass(RRY1);}
if(type1==4){fgs1=copyVclass(YYR1);}
if(type1==5){fgs1=copyVclass(RYR1);}
if(type1==6){fgs1=copyVclass(YRR1);}
if(type1==7){fgs1=copyVclass(RRR1);}
                               
if(type2==0){fgs2=copyVclass(YYY2);}
if(type2==1){fgs2=copyVclass(RYY2);}
if(type2==2){fgs2=copyVclass(YRY2);}
if(type2==3){fgs2=copyVclass(RRY2);}
if(type2==4){fgs2=copyVclass(YYR2);}
if(type2==5){fgs2=copyVclass(RYR2);}
if(type2==6){fgs2=copyVclass(YRR2);}
if(type2==7){fgs2=copyVclass(RRR2);}
//Translate each to there given positions
int i = 0;
int j = 0;
int tot1=fgi1.size();
int tot2=fgi2.size();
//Transform to global coordinates
while(i < tot1){fgi1[i++].transform(t1);}
while(j < tot2){fgi2[j++].transform(t2);}
//Translate only x,y,z  
i=0;
j=0;
while(i < tot1){fgi1[i++].translate(o1);}
while(j < tot2){fgi2[j++].translate(o2);}
//Determine cross product of rv,vv
i=0;
j=0;
int tot3=fgs1.size();
int tot4=fgs2.size();
//Transform to global coordinates
while(i < tot3){fgs1[i++].transform(t1);}
while(j < tot4){fgs2[j++].transform(t2);}
//Translate only x,y,z  
i=0;
j=0;
while(i < tot3){fgs1[i++].translate(o1);}
while(j < tot4){fgs2[j++].translate(o2);}
//Print out aligned fragments and their clusters
// x,y,z,clus_id,frag_id
// x,y,z,r,clus_id
i=0;
j=0;
FILE * ffile = fopen("aligned_frags.dat","w");
while(i < tot3){
fprintf(ffile,"1,%d,%d,%8.3f,%8.3f,%8.3f\n",fgs1[i].group,fgs1[i].id,fgs1[i].x,fgs1[i].y,fgs1[i].z);
i++;
}
while(j < tot4){
fprintf(ffile,"2,%d,%d,%8.3f,%8.3f,%8.3f\n",fgs2[j].group,fgs2[j].id,fgs2[j].x,fgs2[j].y,fgs2[j].z);
j++;
}
fflush(ffile);
fclose(ffile);


FILE * cfile = fopen("aligned_clusters.dat","w");
i=0;
j=0;
while(i < tot1){
fprintf(cfile,"1,%d,%8.3f,%8.3f,%8.3f,%8.3f\n",fgi1[i].id,fgi1[i].r,fgi1[i].x,fgi1[i].y,fgi1[i].z);
i++;
}
while(j < tot2){
fprintf(cfile,"2,%d,%8.3f,%8.3f,%8.3f,%8.3f\n",fgi2[j].id,fgi2[j].r,fgi2[j].x,fgi2[j].y,fgi2[j].z);
j++;
}
fflush(cfile);
fclose(cfile);

i=0;
j=0;

int ms1=0;
ms1=3240000;
int ms2=0;
double dist,x1,y1,z1,x2,z2,y2,y3,x3,z3,a,dr,dv,dz,xi,yi,zi,xj,yj,zj,angout,x4,y4,z4,r1,r2;
//Print Out Positions of fragments
/*
for(int x=0; x<tot3; x++){

}
*/


while(i < tot1){
	j=0;
while(j < tot2){
x1=fgi1[i].x;
y1=fgi1[i].y;
z1=fgi1[i].z;

x2=fgi2[j].x;
y2=fgi2[j].y;
z2=fgi2[j].z;
 
x3=x2-x1;
y3=y2-y1;
z3=z2-z1;

dist=pow(x3,2)+pow(y3,2)+pow(z3,2);
//Radius of fragment is stored in rx
a = abs(sqrt(dist)-sqrt(L));
//a = abs(dist - L);
//Look at all positions < 4A difference
//Calculate all pairwise distances for fragments in the two groups
//Set distance threshold 
if(a < (2*(fgi1[i].r+fgi2[j].r))){
//Check pairwise distances of members in the group
//Get fgi1.members and fgi2.members 
int i1,i2;
for(int f1=0; f1<fgi1[i].members.size(); f1++){
i1=fgi1[i].members[f1];
if(i1==33890){printf("cluster 1 id: %d",fgi1[i].id);}
for(int f2=0; f2<fgi2[j].members.size(); f2++){
if(i2==16864){printf("cluster 2 id: %d",fgi2[j].id);}
ms2++;
//int i1=searchfrag(fgs1,fgi1[i].members[f1]);
//int i2=searchfrag(fgs2,fgi2[j].members[f2]);
i2=fgi2[j].members[f2];
//cout << "i1: " << i1 << " i2: " << i2 << " type: " << type1 << " type2: " << type2 << " f1: " << f1 << " f2: " << f2 << " i: " << i << " j: " << j << "x: " <<  endl;
x1=fgs1[i1].x;
y1=fgs1[i1].y;
z1=fgs1[i1].z;

x2=fgs2[i2].x;
y2=fgs2[i2].y;
z2=fgs2[i2].z;
 
x3=x2-x1;
y3=y2-y1;
z3=z2-z1;
dist=pow(x3,2)+pow(y3,2)+pow(z3,2);
//Radius of fragment is stored in rx
a = abs(sqrt(dist)-sqrt(L));
if(a > 50){
//cout << fgs1[i1].id << ":" << fgs2[i2].id << ":" << a<< "Type1: " << type1 << "Type2: " << type2 << endl;
}
if(a < 1){
//convert all coordinate systems to pos1 bottom coordinates
//If the structure is aligned to pos1-1 then the coordinates are already in line 
dr = fgs1[i1].dotr(&fgs2[i2]);
dv = fgs1[i1].dotv(&fgs2[i2]);

//fgsi v x r
xi=fgs1[i1].vy*fgs1[i1].rz-fgs1[i1].vz*fgs1[i1].ry;
yi=fgs1[i1].vz*fgs1[i1].rx-fgs1[i1].vx*fgs1[i1].rz;
zi=fgs1[i1].vx*fgs1[i1].ry-fgs1[i1].vy*fgs1[i1].rx;

//Convert to fgs1 coordinate system
x4=fgs1[i1].vx*x3 + fgs1[i1].rx*y3 + xi*z3;
y4=fgs1[i1].vy*x3 + fgs1[i1].ry*y3 + yi*z3;
z4=fgs1[i1].vz*x3 + fgs1[i1].rz*y3 + zi*z3;

//Angle between original vector and new vector
angout = acos((Lx*x4+Ly*y4+Lz*z4)/(sqrt(L)*sqrt(pow(x4,2)+pow(y4,2)+pow(z4,2))+0.00000001));

xj=fgs2[i2].vy*fgs2[i2].rz-fgs2[i2].vz*fgs2[i2].ry;
yj=fgs2[i2].vz*fgs2[i2].rx-fgs2[i2].vx*fgs2[i2].rz;
zj=fgs2[i2].vx*fgs2[i2].ry-fgs2[i2].vy*fgs2[i2].rx;
dz = xi*xj+yi*yj+zi*zj;

dz=abs(acos(dz)-acos(odz));
dr=abs(acos(dr)-acos(odr));
dv=abs(acos(dv)-acos(odv));

//cout << fgs1[i1].id << "," << fgs2[i2].id << "," << a << "," << dr << "," << dv << "," << dz << "," << angout << endl;

if(angout < 0.05){
if(dr < 0.1 & dv < 0.1 & dz < 0.1){
Matches * mat=new Matches();
mat->frag1=fgs1[i1].id;
mat->frag2=fgs2[i2].id;
mat->score = dr+dv+dz+angout+a;
mat->pos1 = pos1;
mat->pos2 = pos2;
matches.push_back(*mat);
delete mat;
//cout << fgs1[i1].id << "," << fgs2[i2].id << "," << a << "," << dr << "," << dv << "," << dz << "," << angout << endl;
}
}
}
//cout << fgs1[i]->id << "," << fgs2[j]->id<<"," << fgs1[i]->x << "," << fgs1[i]->y <<  "," << fgs1[i]->z<<  "," << fgs2[j]->x<<  "," << fgs2[j]->y<<  "," << fgs2[j]->z << endl;  
//cout << fgs1[i]->id << "," << fgs2[j]->id<<"," << fgs1[i]->x << "," << fgs1[i]->y <<  "," << fgs1[i]->z<<  "," << fgs2[j]->x<<  "," << fgs2[j]->y<<  "," << fgs2[j]->z << endl;  
}	
	}
}
j++;
}
i++;
}
//cout << "Eff: "  << "-" <<  ms2 << ":" << ms1 << endl;
//return matches
}

vector<Matches> Fragclass::matchedMoves(int type1,int type2,double L,double odr, double odv, double odz, Matrix * t1,Matrix * t2,Vclass * o1,Vclass * o2,double Lx, double Ly, double Lz,int time,int pos1, int pos2){
//Transform frag1 and frag2 into same coordinate system with t1 and t2
vector<Matches > matches;
vector<Vclass> fgs1;
vector<Vclass> fgs2;
vector<Cclass> fgi1;
vector<Cclass> fgi2;
//Cluters
if(type1==0){fgi1=copyCclass(iYYY1);}
if(type1==1){fgi1=copyCclass(iRYY1);}
if(type1==2){fgi1=copyCclass(iYRY1);}
if(type1==3){fgi1=copyCclass(iRRY1);}
if(type1==4){fgi1=copyCclass(iYYR1);}
if(type1==5){fgi1=copyCclass(iRYR1);}
if(type1==6){fgi1=copyCclass(iYRR1);}
if(type1==7){fgi1=copyCclass(iRRR1);}
                                
if(type2==0){fgi2=copyCclass(iYYY2);}
if(type2==1){fgi2=copyCclass(iRYY2);}
if(type2==2){fgi2=copyCclass(iYRY2);}
if(type2==3){fgi2=copyCclass(iRRY2);}
if(type2==4){fgi2=copyCclass(iYYR2);}
if(type2==5){fgi2=copyCclass(iRYR2);}
if(type2==6){fgi2=copyCclass(iYRR2);}
if(type2==7){fgi2=copyCclass(iRRR2);}
//Fragment Types
if(type1==0){fgs1=copyVclass(YYY1);}
if(type1==1){fgs1=copyVclass(RYY1);}
if(type1==2){fgs1=copyVclass(YRY1);}
if(type1==3){fgs1=copyVclass(RRY1);}
if(type1==4){fgs1=copyVclass(YYR1);}
if(type1==5){fgs1=copyVclass(RYR1);}
if(type1==6){fgs1=copyVclass(YRR1);}
if(type1==7){fgs1=copyVclass(RRR1);}
                               
if(type2==0){fgs2=copyVclass(YYY2);}
if(type2==1){fgs2=copyVclass(RYY2);}
if(type2==2){fgs2=copyVclass(YRY2);}
if(type2==3){fgs2=copyVclass(RRY2);}
if(type2==4){fgs2=copyVclass(YYR2);}
if(type2==5){fgs2=copyVclass(RYR2);}
if(type2==6){fgs2=copyVclass(YRR2);}
if(type2==7){fgs2=copyVclass(RRR2);}
//Translate each to there given positions
int i = 0;
int j = 0;
int tot1=fgi1.size();
int tot2=fgi2.size();
//Transform to global coordinates
while(i < tot1){fgi1[i++].transform(t1);}
while(j < tot2){fgi2[j++].transform(t2);}
//Translate only x,y,z  
i=0;
j=0;
while(i < tot1){fgi1[i++].translate(o1);}
while(j < tot2){fgi2[j++].translate(o2);}
//Determine cross product of rv,vv
i=0;
j=0;
int tot3=fgs1.size();
int tot4=fgs2.size();
//Transform to global coordinates
while(i < tot3){fgs1[i++].transform(t1);}
while(j < tot4){fgs2[j++].transform(t2);}
//Translate only x,y,z  
i=0;
j=0;
while(i < tot3){fgs1[i++].translate(o1);}
while(j < tot4){fgs2[j++].translate(o2);}
/*
//Print out aligned fragments
i=0;
j=0;
while(i < tot3){cout << fgs1[i].x << "," << fgs1[i].y << "," << fgs1[i++].z << "," << 1 << endl;}
while(j < tot4){cout << fgs2[j].x << "," << fgs2[j].y << "," << fgs2[j++].z << "," << 2 << endl;}
*/
i=0;
j=0;
int ms2=0;
double dist,x1,y1,z1,x2,z2,y2,y3,x3,z3,a,dr,dv,dz,xi,yi,zi,xj,yj,zj,angout,x4,y4,z4,r1,r2;
//Print Out Positions of fragments
/*
for(int x=0; x<tot3; x++){

}
*/

if(L > 0){
while(i < tot1){
	j=0;
while(j < tot2){
x1=fgi1[i].x;
y1=fgi1[i].y;
z1=fgi1[i].z;

x2=fgi2[j].x;
y2=fgi2[j].y;
z2=fgi2[j].z;
 
x3=x2-x1;
y3=y2-y1;
z3=z2-z1;

dist=pow(x3,2)+pow(y3,2)+pow(z3,2);
//Radius of fragment is stored in rx
a = abs(sqrt(dist)-sqrt(L));
//a = abs(dist - L);
//Look at all positions < 4A difference
//Calculate all pairwise distances for fragments in the two groups
//Set distance threshold 
if(a < (2*(fgi1[i].r+fgi2[j].r))){
//Check pairwise distances of members in the group
//Get fgi1.members and fgi2.members 
int i1,i2;
for(int f1=0; f1<fgi1[i].members.size(); f1++){
i1=fgi1[i].members[f1];
for(int f2=0; f2<fgi2[j].members.size(); f2++){
ms2++;
//int i1=searchfrag(fgs1,fgi1[i].members[f1]);
//int i2=searchfrag(fgs2,fgi2[j].members[f2]);
i2=fgi2[j].members[f2];
//cout << "i1: " << i1 << " i2: " << i2 << " type: " << type1 << " type2: " << type2 << " f1: " << f1 << " f2: " << f2 << " i: " << i << " j: " << j << "x: " <<  endl;
x1=fgs1[i1].x;
y1=fgs1[i1].y;
z1=fgs1[i1].z;

x2=fgs2[i2].x;
y2=fgs2[i2].y;
z2=fgs2[i2].z;
 
x3=x2-x1;
y3=y2-y1;
z3=z2-z1;
dist=pow(x3,2)+pow(y3,2)+pow(z3,2);
//Radius of fragment is stored in rx
a = abs(sqrt(dist)-sqrt(L));
if(a > 50){
cout << fgs1[i1].id << ":" << fgs2[i2].id << ":" << a<< "Type1: " << type1 << "Type2: " << type2 << endl;
}
if(a < 1){
//convert all coordinate systems to pos1 bottom coordinates
//If the structure is aligned to pos1-1 then the coordinates are already in line 
dr = fgs1[i1].dotr(&fgs2[i2]);
dv = fgs1[i1].dotv(&fgs2[i2]);

//fgsi v x r
xi=fgs1[i1].vy*fgs1[i1].rz-fgs1[i1].vz*fgs1[i1].ry;
yi=fgs1[i1].vz*fgs1[i1].rx-fgs1[i1].vx*fgs1[i1].rz;
zi=fgs1[i1].vx*fgs1[i1].ry-fgs1[i1].vy*fgs1[i1].rx;

//Convert to fgs1 coordinate system
x4=fgs1[i1].vx*x3 + fgs1[i1].rx*y3 + xi*z3;
y4=fgs1[i1].vy*x3 + fgs1[i1].ry*y3 + yi*z3;
z4=fgs1[i1].vz*x3 + fgs1[i1].rz*y3 + zi*z3;

//Angle between original vector and new vector
angout = acos((Lx*x4+Ly*y4+Lz*z4)/(sqrt(L)*sqrt(pow(x4,2)+pow(y4,2)+pow(z4,2))+0.00000001));

xj=fgs2[i2].vy*fgs2[i2].rz-fgs2[i2].vz*fgs2[i2].ry;
yj=fgs2[i2].vz*fgs2[i2].rx-fgs2[i2].vx*fgs2[i2].rz;
zj=fgs2[i2].vx*fgs2[i2].ry-fgs2[i2].vy*fgs2[i2].rx;
dz = xi*xj+yi*yj+zi*zj;

dz=abs(acos(dz)-acos(odz));
dr=abs(acos(dr)-acos(odr));
dv=abs(acos(dv)-acos(odv));

//cout << fgs1[i1].id << "," << fgs2[i2].id << "," << a << "," << dr << "," << dv << "," << dz << "," << angout << endl;

if(angout < 0.05){
if(dr < 0.1 & dv < 0.1 & dz < 0.1){
Matches * mat=new Matches();
mat->frag1=fgs1[i1].id;
mat->frag2=fgs2[i2].id;
mat->score = dr+dv+dz+angout+a;
mat->pos1 = pos1;
mat->pos2 = pos2;
matches.push_back(*mat);
delete mat;
//cout << fgs1[i1].id << "," << fgs2[i2].id << "," << a << "," << dr << "," << dv << "," << dz << "," << angout << endl;
}
}
}
//cout << fgs1[i]->id << "," << fgs2[j]->id<<"," << fgs1[i]->x << "," << fgs1[i]->y <<  "," << fgs1[i]->z<<  "," << fgs2[j]->x<<  "," << fgs2[j]->y<<  "," << fgs2[j]->z << endl;  
//cout << fgs1[i]->id << "," << fgs2[j]->id<<"," << fgs1[i]->x << "," << fgs1[i]->y <<  "," << fgs1[i]->z<<  "," << fgs2[j]->x<<  "," << fgs2[j]->y<<  "," << fgs2[j]->z << endl;  
}	
	}
}
j++;
}
i++;
}
//cout << "Eff: "  << "-" <<  ms2 << ":" << ms1 << endl;
//return matches
}else{
//Filtering for end to end distances within 14
while(i < tot1){
	j=0;
while(j < tot2){
x1=fgi1[i].x;
y1=fgi1[i].y;
z1=fgi1[i].z;

x2=fgi2[j].x;
y2=fgi2[j].y;
z2=fgi2[j].z;
 
x3=x2-x1;
y3=y2-y1;
z3=z2-z1;

dist=pow(x3,2)+pow(y3,2)+pow(z3,2);
if(dist < 144){//set this to 144 to screen in helices
//Check pairwise distances of members in the group
//Get fgi1.members and fgi2.members 
int i1,i2;
for(int f1=0; f1<fgi1[i].members.size(); f1++){
i1=fgi1[i].members[f1];
for(int f2=0; f2<fgi2[j].members.size(); f2++){
ms2++;
//int i1=searchfrag(fgs1,fgi1[i].members[f1]);
//int i2=searchfrag(fgs2,fgi2[j].members[f2]);
i2=fgi2[j].members[f2];
//cout << "i1: " << i1 << " i2: " << i2 << " type: " << type1 << " type2: " << type2 << " f1: " << f1 << " f2: " << f2 << " i: " << i << " j: " << j << "x: " <<  endl;
x1=fgs1[i1].x;
y1=fgs1[i1].y;
z1=fgs1[i1].z;

x2=fgs2[i2].x;
y2=fgs2[i2].y;
z2=fgs2[i2].z;
 
x3=x2-x1;
y3=y2-y1;
z3=z2-z1;
dist=pow(x3,2)+pow(y3,2)+pow(z3,2);

if(dist > 36 && dist < 144){
Matches * mat=new Matches();
mat->frag1=fgs1[i1].id;
mat->frag2=fgs2[i2].id;
mat->score = dr+dv+dz+angout+a;
mat->pos1 = pos1;
mat->pos2 = pos2;
matches.push_back(*mat);
delete mat;
}
//cout << fgs1[i]->id << "," << fgs2[j]->id<<"," << fgs1[i]->x << "," << fgs1[i]->y <<  "," << fgs1[i]->z<<  "," << fgs2[j]->x<<  "," << fgs2[j]->y<<  "," << fgs2[j]->z << endl;  
//cout << fgs1[i]->id << "," << fgs2[j]->id<<"," << fgs1[i]->x << "," << fgs1[i]->y <<  "," << fgs1[i]->z<<  "," << fgs2[j]->x<<  "," << fgs2[j]->y<<  "," << fgs2[j]->z << endl;  
}	
	}
}
j++;
}
i++;
}



}



return matches;
}



