/// polymere.h -- polymere class and functions
#ifndef POLYMERE_H_
#define POLYMERE_H_

#include "dynaclust.h"
#include "score.h"
#include "parser.h"
#include "base.h"
#include "fragment.h"
#include "fragatoms.h"
#include "fragall.h"
#include "doublet.h"
#include "fragclass.h"
#include "matrix.h"
#include <vector>
using namespace std;
class Polymere
{
private:
public:
Polymere();
Polymere(string file);
Polymere(string file,int model);
Polymere(string seq,std::vector< Fragment *>  frags,Fragclass * dat,std::vector< Fragatoms *> atoms,std::vector<Fragall *> batoms,std::vector<double> * dparams);
//Virtual atom storage
vector<Base *> mols;
int num_mols;
double score;
int randtot;
//Full atom storage
//Each residue has the same number of backbone atoms.
vector<double> full_atoms;
//Keep track of which residue is at each position 1=A,2=C,G=3,U=4
vector<int> residues;

//Functions to get particular atoms from the full atom storage
int getResidueStart(int pos);
//Returns index that points to the start position of back bone atoms for a given position. The first 12x3 positions are the coordinates of the backbone atoms.  The following are the atoms for the resiude

//Gets coordinates of a particular atom in residue at position pos
//vector<double> getAtom(char * name, int pos);

void hydroVg(vector<float> * hydroVals,int movegroup);
double rmsAllg(const Polymere * p2,int movegroup);
void alignHelix(Polymere * p2,std::vector<int> helix1, std::vector<int> helix2);
void copyOver(Polymere * pol2,int into_pos,int from_pos,int totlen);
float monteCons(int nmax, int movegroup, const std::vector< Fragment *>  frags,Fragclass * dat,const Polymere * pol,std::vector< Fragatoms *> atoms,std::vector<Fragall *> batoms,std::vector<double> * dparams,int stype,std::vector< std::vector <std::pair<int,int> > > allcons,int sim_id,std::string out_prefix);
std::vector<std::pair<int,int> > confMatch(std::vector< std::vector< std::pair<int,int> > >  confs);
void genAllCons(int nmax,std::vector<std::vector< std::pair<int,int> > > * allcons,std::string rnaseq);
void genCons(std::vector< std::pair<int,int> > * cons,int p1,int lt,int asym_num,std::string rnaseq);
//Writes a PDB file using the all atom model 
void writePDB(string filename);
//void writePDB(char * filename);
void midPointBuild(std::vector<std::pair<int,int> > movegroups,std::vector< std::pair<float,float> > cons,const std::vector< Fragment *>  frags,Fragclass * dat,std::vector< Fragatoms *> atoms,std::vector<Fragall *> batoms,std::vector<double> * dparams);

//Add all atoms to virtual bond model polymere
void coordSysAT(int pos,std::vector<Points> * pts);
float constraints(std::vector< std::pair<int, int> > * constraints,int clash_level, float ctemp);
int hardConstraints(std::vector< std::pair<int,int> > * constraints);
int hardConstraintsL(std::vector< std::pair<int,int> > * cons);
float updateConstraints(std::vector< std::pair<int, int> > * constraints);
vector<Matches> fixed5to3(int pos1, int pos2,double Lmin,double Lmax, Fragclass * frags);
vector<Matches> fixed3to5(int pos1, int pos2,double Lmin,double Lmax, Fragclass * frags);
//Functions for virtual bond model 
void copy(Polymere * pol2);
void readPDB(string);
void updateFull();
void printMols(int);
void topov(std::string outfile,double camx,double camy, double camz);
int c1dist();
float c1distf(int,int);
float rms(Polymere *);
float rms_doublet(std::vector<Doublet> dbl, int i, int j);
float minAA(int i, int j);
void structParams(std::vector<Doublet>, int, int, int);
float minDD(int i, int j);
float acceptorDist(int i,int j);
float donorDist(int i, int j);
float rms_doubletA(std::vector<Doublet> dlt23,int i,int j);
double moldist(int,int);
double hdist(int,int);
float doublet(std::vector<Doublet> db);
float clusterDist(std::vector<Polymere *>);
float cluster(const std::vector<Polymere *> * vect1 , std::vector<int> * vect);
void cluster2(std::vector<Polymere *> * vect1 , std::vector<int> * vect,int max_size,int max_iter);
void readMulti(string, int,int);
double checkbps();
int checkBackBone();
int groupCheckBackBone(int group_id,float p_o_lower=1.96, float p_o_upper=3.81);
int atomClashes();
int groupCheckClashes(int group_id);
double checkhbonds();
double rmsf(Polymere *);
float rmsAllf(Polymere *);
double rmsAll(Polymere *);
double rmsAllc(const Polymere *);
vector<Base *> readStruct(string);
vector<double> readDist(string);
void readStructNew(string);
void writeStruct(string,double);
void writeStructInfo(string,double,double,double);
void writeStructInfo2(string);
void writeDists(string,double);
double scoref();
void resscore(Score * score, vector<double> * dparams);
void group_resscore(Score * score,vector<double> * dparams,int group_id);
void lowres(Score * score,std::vector<float> * resscores);
vector<Fragatoms *> readAtoms(string);
vector<Fragment *> readFrags(string);
vector<Fragment *> readFrags2(string);
Matrix * coordSysF(int pos);
Matrix * coordSysG(int pos);
Matrix * coordSysE(int pos);
Matrix * coordSysR(int pos);
void coordSysNt(int pos,std::vector<Points> * pts);
void baseCoords(int i,int j);
void getParams(int i,int j,vector<double> * p);
double gauss(const vector<double> * p); 
double gaussDbl_ij(const vector<double> * p,int,int,double,double,double,double,double); 
double gaussDbl(const vector<double> * p,double,double,double,double,double,const std::vector<int> * group); 
void doubletV(const vector<double> * p,double,double,double,double,double,std::vector<float> * probs); 
void doubletVg(const vector<double> * p,double,double,double,double,double,std::vector<float> * probs,int movegroup); 
vector<Matches> chooseMoves(int pos1, int pos2,Fragclass * frags,int time);
vector<Matches> chooseMovesPairs(int pos1, int pos2,Fragclass * frags);
void printMoves(int pos1, int pos2,Fragclass * frags);
vector<Matches> fixedMoves(int pos1, int pos2,Fragclass * frags,double L);
int chooseSingle(int pos1,Fragclass * frags);
void chooseAllSingle(int pos1,Fragclass * frags,std::vector<int> * fids);
void cartesian_to_sperical(double x,double y, double z, double *r,double *angy,double * angz);
void rotate3d(double a,double b,double g,double x2,double y2,double z2);
void rotateX(double ang);
void rotateXr(double ang, int);
void rotateXb(double, int);
void rotateY(double ang);
void rotateZ(double ang);
void translate(double x2,double y2,double z2);
void translateb(int pos, double x2,double y2,double z2);
void translater(int pos, double x2,double y2,double z2);
void changeFrag(int, int, std::vector<Fragment *>);
void changeFrag2(int, int, const std::vector<Fragment *>);
void changeFrag1(int, int, const std::vector<Fragment *>);
void changeFragVid(int, int, std::vector<Fragment *>);
void changeFragR(int, int, std::vector<Fragment *>);
void changeTor(int pos,double ang);
void changeTor2(int pos,double ang);
void changeTor2r(int pos,double ang);
float hydro();
void hydroV(vector<float> * hydro);
float hydroG(const std::vector<int> * groups);
float hydroP(int,int);
float hydroType(int i,int j,std::pair<int,int> * faces);
void hydroPairs(std::vector<float> * hydroVals,std::vector<int> * hpair);
float stacked();
float groupStacked(int group_id);
float adjStacked(std::vector<int> * probs);
void genStacked(const std::vector< Fragment *>  frags,Fragclass * dat,std::vector< Fragatoms *> atoms,std::vector<Fragall *> batoms);

int assemble(std::vector<std::vector<int> > movepos,std::vector<std::vector<std::pair<std::pair<int,int>,int> > > gencons,std::vector<std::vector<int> > clashgroups, const std::vector< Fragment *>  frags,Fragclass * dat,std::vector< Fragatoms *> atoms,std::vector<Fragall *> batoms,std::vector<std::vector<int> > loopO,std::vector<std::vector<int> > loopE, int max_iterations,int max_leaves, int group,std::string ofile,int write_step,float min_res, int max_clust,int start_con,float p_o_lower,float p_o_upper);

int consBuild1(int movepos,std::vector<std::pair<std::pair<int,int>,int> > gencons, const std::vector< Fragment *>  frags,Fragclass * dat,std::vector< Fragatoms *> atoms,std::vector<Fragall *> batoms,std::vector<std::vector<int> > loopO,std::vector<int> clashGroup,int max_iterations,int max_leaves,int group,std::vector<float> * con_score,std::vector<int> * constree,std::vector< std::vector<int> > * fragtree,int step,std::string ofile,int write_steps,float p_o_lower,float p_o_upper);
int consBuild2(int movepos,std::vector<std::pair<std::pair<int,int>,int> > gencons, const std::vector< Fragment *>  frags,Fragclass * dat,std::vector< Fragatoms *> atoms,std::vector<Fragall *> batoms,std::vector<std::vector<int> > loopE,std::vector<int> clashGroup,int max_iterations,int max_leaves,int group,std::vector<float> * con_score,std::vector<int> * constree,std::vector< std::vector<int> > * fragtree,int step,std::string ofile,int write_step,float p_o_lower,float p_o_upper);
int consBuild3(int movepos1,int movepos2,std::vector<std::pair<std::pair<int,int>,int> > gencons, const std::vector< Fragment *>  frags,Fragclass * dat,std::vector< Fragatoms *> atoms,std::vector<Fragall *> batoms,std::vector<int> clashGroup,int max_iterations,int max_leaves,int group,std::vector<float> * con_score,std::vector<int> * constree,std::vector< std::vector<int> > * fragtree,int step,std::string ofile,int write_step,float p_o_lower,float p_o_upper);
int consBuild4(std::vector<int> movepos,std::vector<std::pair<std::pair<int,int>,int> > gencons, const std::vector< Fragment *>  frags,Fragclass * dat,std::vector< Fragatoms *> atoms,std::vector<Fragall *> batoms,std::vector<int> clashGroup,int max_iterations,int max_leaves,int group,std::vector<float> * con_score,std::vector<int> * constree,std::vector< std::vector<int> > * fragtree,int step,std::string ofile,int write_step,float min_res, int max_clust,float p_o_lower,float p_o_upper);

int generalCons(std::vector<std::pair<std::pair<int,int> , int> > genCons ,float * total_hydro);
int generalCon(std::pair<int,int>  genCons ,float * hydro);
void buildLoop(int start_pos,int end_pos,int nmax,int iter_max,const std::vector< Fragment *>  frags,Fragclass * dat,std::vector< Fragatoms *> atoms,std::vector<Fragall *> batoms,string gfile,string loopseq);
void buildLoopO(int start_pos,int end_pos,int nmax,int iter_max,const std::vector< Fragment *>  frags,Fragclass * dat,std::vector< Fragatoms *> atoms,std::vector<Fragall *> batoms,string gfile,string loopseq);
int buildBuldge(int s1, int s2,std::vector<int> * small_side,const std::vector< Fragment *>  frags,Fragclass * dat,std::vector< Fragatoms *> atoms,std::vector<Fragall *> batoms,int s_side,std::string outfile);
float zaligned(int i,int j);
int genPairedCons(int clashGroup,std::vector<std::pair<int,int> > * cons, const std::vector< Fragment *>  frags,Fragclass * dat,std::vector< Fragatoms *> atoms,std::vector<Fragall *> batoms,std::vector<std::pair<int,int> > * branch);
void genPairedStacks(int pos1,  int pos2,int clashGroup, const std::vector< Fragment *>  frags,Fragclass * dat,std::vector< Fragatoms *> atoms,std::vector<Fragall *> batoms);
int findLoop(int start_pos,int type,std::vector< std::vector<int> > * frag_list,std::vector<Fragment *> frags, Fragclass * dat,std::vector<Fragatoms* > atoms, std::vector<Fragall *> batoms);
float hydroA();
double coplanar(int i,int j);
int arePlanar(int i,int j);
std::vector<double> planar(int i,int j);
void changeTorR(int pos,double ang);
void changeBaseTor(int pos,double ang);
void changeBond1(int pos,double ang);
void changeBond2(int pos,double ang);
void changeBond2r(int pos,double ang);
void changeBond3(int pos,double ang);
void changeBondLen(int pos,double len);
void changeBondLen2r(int pos,double len);
void changeBondLenR(int pos,double len);
void alignStruct(int pos);
void alignStruct2(int pos);
void align(int pos,double *a, double *b, double *c, double *d,double*e,double*f);
void align3to5(int pos,double *a, double *b, double *c, double *d,double*e,double*f);
void align1(int pos,double *a, double *b, double *c, double *d,double*e,double*f);
void align3(int pos,double *a, double *b, double *c, double *d,double*e,double*f);
void align4(int pos,double *a, double *b, double *c, double *d,double*e,double*f);
void changeBond3a(int pos, double ang);
void rotateZb(double ang, int pos);
//float rmsd(float *v1, float *v2, int N, float *mtx);
void localRMSD(vector<float> * rmsds,const Polymere * master,std::vector<int> * parent,std::vector<int> * child);
void localgraph(std::vector<int> * parent, std::vector<int> * child);
void printVoro();
void clearMols();
void rotateZf(double ang, int pos);
void mixit(int nmax, int movegroup,const std::vector< Fragment *>  frags,Fragclass * dat,const Polymere * pol,std::vector<Fragatoms *> atoms,std::vector<Fragall *> batoms);
float trimix(double rmsd,int nmax, int start_center,int movegroup,const std::vector< Fragment *>  frags,Fragclass * dat,const Polymere * pol,std::vector<Polymere *> *centers,std::vector<Fragatoms *> atoms,std::vector<Fragall *> batoms);
float drimix(double rmsd,int nmax, int start_center,int movegroup,const std::vector< Fragment *>  frags,Fragclass * dat,const Polymere * pol,std::vector<Polymere *> *centers,std::vector<Fragatoms *> atoms,std::vector<Fragall *> batoms);
int modelit(int cycles, int capture,int movegroup,const std::vector< Fragment *>  frags,Fragclass * dat,const Polymere * pol,std::vector<Fragatoms *> atoms,std::vector<Fragall *> batoms,Parser * model,std::vector<double> * dparams,string fileout,std::vector<int> * parent, std::vector<int> * child,string sim_id);
int entropy(int nmax, int movegroup,const std::vector< Fragment *>  frags,Fragclass * dat,const Polymere * pol,std::vector<Fragatoms *> atoms,std::vector<Fragall *> batoms);
void dynamix(int nmax, int movegroup,const std::vector< Fragment *>  frags,Fragclass * dat,const Polymere * pol,std::vector<Fragatoms *> atoms,std::vector<Fragall *> batoms,int nt1,int nt2);
void simulate(int rmax, int nmax, int movegroup,const std::vector< Fragment *>  frags,Fragclass * dat,const Polymere * pol,std::vector<Fragatoms *> atoms,std::vector<Fragall *> batoms,std::vector< std::pair<int,int> > * cons);
void scoreit(int nmax, int movegroup,const std::vector< Fragment *>  frags,std::vector<Fragatoms *> atoms,std::vector<Fragall *> batoms,Fragclass * dat,Polymere * pol,std::vector<double> * p);
float monte(int nmax, int movegroup, const std::vector< Fragment *>  frags,Fragclass * dat,const Polymere * pol,std::vector< Fragatoms *> atoms,std::vector<Fragall *> batoms,std::vector<double> * dparams,int stype);
};

#endif



