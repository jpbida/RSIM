// molecule.h -- molecule class and functions
#ifndef SCORE_H_
#define SCORE_H_
#include <vector>
#include <string>

using namespace std;
class Score
{
private:

public:

//high resolution constants
const static float h_gfreqs[8];
const static float h_afreqs[8];
const static float h_ufreqs[8];
const static float h_cfreqs[8];

//High Resolution - 3-atoms per base exact volumes
/*
FOR some reason i need to initialize the floats inside the class for the UofM machines to compile
static const float h_g1mu=28.5121833333333;
static const float h_g2mu=29.6294623217923;
static const float h_g3mu=29.926917617237;
static const float h_g1sig=7.12826885918243;
static const float h_g2sig=7.13028483216927;
static const float h_g3sig=8.08907525269187;
       
static const float h_a1mu=24.4057940797941;
static const float h_a2mu=29.1543159722222;
static const float h_a3mu=24.5565109947644;
static const float h_a1sig=6.79689351941335;
static const float h_a2sig=7.50295316102331;
static const float h_a3sig=8.2462165085006;
       
static const float h_c1mu=28.8828419047619;
static const float h_c2mu=31.3192229508197;
static const float h_c3mu=27.4985533063428;
static const float h_c1sig=6.39935886576906;
static const float h_c2sig=7.36920424291617;
static const float h_c3sig=6.68391267867415;
       
static const float h_u1mu=28.2866161048689;
static const float h_u2mu=30.4679118942731;
static const float h_u3mu=28.4437012987013;
static const float h_u1sig=7.58531846171727;
static const float h_u2sig=8.2005337254098;
static const float h_u3sig=6.59684401632079;
*/

const float static h_g1mu;
const float static h_g2mu;
const float static h_g3mu;
const float static h_g1sig;
const float static h_g2sig;
const float static h_g3sig;

const float static h_a1mu;
const float static h_a2mu;
const float static h_a3mu;
const float static h_a1sig;
const float static h_a2sig;
const float static h_a3sig;

const float static h_u1mu;
const float static h_u2mu;
const float static h_u3mu;
const float static h_u1sig;
const float static h_u2sig;
const float static h_u3sig;

const float static h_c1mu;
const float static h_c2mu;
const float static h_c3mu;
const float static h_c1sig;
const float static h_c2sig;
const float static h_c3sig;


const static float freqs1[12];
const static float mu1[12];
const static float sig1[12];

const static float freqs2[12];
const static float    mu2[12];
const static float   sig2[12];

const static float freqs3[12];
const static float    mu3[12];
const static float   sig3[12];

const static float freqs4[12];
const static float    mu4[12];
const static float   sig4[12];

const static float freqs5[12];
const static float    mu5[12];
const static float   sig5[12];

const static float freqs6[12];
const static float    mu6[12];
const static float   sig6[12];

const static float freqs7[12];
const static float    mu7[12];
const static float   sig7[12];

const static float freqs8[12];
const static float    mu8[12];
const static float   sig8[12];

const static float freqs9[12];
const static float    mu9[12];
const static float   sig9[12];

const static float freqs10[12];
const static float    mu10[12];
const static float   sig10[12];
const static float graph_order[15];

Score(int num_mols);
void write(string file);
void write_globals(string file);
float gscore(int gtype);
float gaussianF(float mu,float sig,float x);
int sim_id;
float gscore0;
float gscore1;
float gscore2;
float gscore3;
float gscore4;
float gscore5;
float gscore6;
vector<float> local_rmsd;
  vector<int>  s0;// score->s0[vind]=gorder[vind];
  vector<int>  s1;// score->s1[vind]=gtype1[vind];
  vector<int>  s2;// score->s2[vind]=gtype2[vind];
  vector<int>  s3;// score->s3[vind]=gtype3[vind];
  vector<int>  s4;// score->s4[vind]=gtype4[vind];
  vector<int>  s5;// score->s5[vind]=rtype[vind];
vector<float>  s6;//  score->s6[vind]=volumes[vind];
  vector<int>  s7;//  score->s7[vind]=unpacked[vind];
vector<float>  s8;//  score->s8[vind]=surface[vind];
vector<float>  s9;//  score->s9[vind]=bvolumes[vind];
  vector<int> s10;//  score->s10[vind]=bunpacked[vind];
vector<float> s11;//  score->s11[vind]=bsurface[vind];
vector<float> s12;//  score->s12[vind]=surfG[vind];
vector<float> s13;//  score->s13[vind]=surfA[vind];
vector<float> s14;//  score->s14[vind]=surfU[vind];
vector<float> s15;//  score->s15[vind]=surfC[vind] ;
vector<float> s16;//  score->s16[vind]=bsurfG[vind];
vector<float> s17;//  score->s17[vind]=bsurfA[vind];
vector<float> s18;//  score->s18[vind]=bsurfU[vind];
vector<float> s19;//  score->s19[vind]=bsurfC[vind];
   vector<int> s20;// score->s20[vind]=orderG[vind];
   vector<int> s21;// score->s21[vind]=orderA[vind];
   vector<int> s22;// score->s22[vind]=orderU[vind];
   vector<int> s23;// score->s23[vind]=orderC[vind];
   vector<int> s24;// score->s24[vind]=borderG[vind];
   vector<int> s25;// score->s25[vind]=borderA[vind];
   vector<int> s26;// score->s26[vind]=borderU[vind];
   vector<int> s27;// score->s27[vind]=borderC[vind];
vector<float> s28;//  score->s28[vind]=vatom0[vind];
vector<float> s29;//  score->s29[vind]=vatom1[vind];
vector<float> s30;//  score->s30[vind]=vatom2[vind];
vector<float> s31;//  score->s31[vind]=vatom3[vind];
vector<float> s32;//  score->s32[vind]=vatom4[vind];
vector<float> s33;//  score->s33[vind]=vatom5[vind];
vector<float> s34;//  score->s34[vind]=vatom6[vind];
vector<float> s35;//  score->s35[vind]=vatom7[vind];
vector<float> s36;//  score->s36[vind]=vatom8[vind];
vector<float> s37;//  score->s37[vind]=vatom9[vind];
vector<float> s38;//  score->s38[vind]=vatom10[vind];
vector<float> s39;//  score->s39[vind]=vatom11[vind];
vector<float> s40;//  score->s40[vind]=vatom12[vind];
vector<float> s41;//  score->s41[vind]=vatom13[vind];
vector<float> s42;//  score->s42[vind]=vatom14[vind];
vector<float> s43;//  score->s43[vind]=vatom15[vind];
vector<float> s44;//  score->s44[vind]=vatom16[vind];
vector<float> s45;//  score->s45[vind]=vatom17[vind];
vector<float> s46;//  score->s46[vind]=vatom18[vind];
vector<float> s47;//  score->s47[vind]=vatom19[vind];
vector<float> s48;//  score->s48[vind]=vatom20[vind];
vector<float> s49;//  score->s49[vind]=vatom21[vind];
vector<float> s50;//  score->s50[vind]=vatom22[vind];
vector<float> s51;//  score->s51[vind]=vatom23[vind];
vector<float> s52;//  score->s52[vind]=vatom24[vind];
vector<float> s53;//  score->s53[vind]=vatom25[vind];
vector<float> s54;//  score->s54[vind]=vatom26[vind];
  vector<int> s55;//  score->s55[vind]=atom0[vind];
  vector<int> s56;//  score->s56[vind]=atom1[vind];
  vector<int> s57;//  score->s57[vind]=atom2[vind];
  vector<int> s58;//  score->s58[vind]=atom3[vind];
  vector<int> s59;//  score->s59[vind]=atom4[vind];
  vector<int> s60;//  score->s60[vind]=atom5[vind];
  vector<int> s61;//  score->s61[vind]=atom6[vind];
  vector<int> s62;//  score->s62[vind]=atom7[vind];
  vector<int> s63;//  score->s63[vind]=atom8[vind];
  vector<int> s64;//  score->s64[vind]=atom9[vind];
  vector<int> s65;//  score->s65[vind]=atom10[vind];
  vector<int> s66;//  score->s66[vind]=atom11[vind];
  vector<int> s67;//  score->s67[vind]=atom12[vind];
  vector<int> s68;//  score->s68[vind]=atom13[vind];
  vector<int> s69;//  score->s69[vind]=atom14[vind];
  vector<int> s70;//  score->s70[vind]=atom15[vind];
  vector<int> s71;//  score->s71[vind]=atom16[vind];
  vector<int> s72;//  score->s72[vind]=atom17[vind];
  vector<int> s73;//  score->s73[vind]=atom18[vind];
  vector<int> s74;//  score->s74[vind]=atom19[vind];
  vector<int> s75;//  score->s75[vind]=atom20[vind];
  vector<int> s76;//  score->s76[vind]=atom21[vind];
  vector<int> s77;//  score->s77[vind]=atom22[vind];
  vector<int> s78;//  score->s78[vind]=atom23[vind];
vector<int> s79;//   score->s79[vind]=atom24[vind];
vector<int> s80;//   score->s80[vind]=atom25[vind];
vector<int> s81;//   score->s81[vind]=atom26[vind];
vector<float> s82;// score->s82[vind]=hydros[vind];
vector<float> s83;// score->s83[vind]=doublets[vind];
vector<int> s84;//   type of resiude that is contributing the most hydrogen bonds  
vector<float> s85;// volume probability - based on small_rna'
vector<float> s86;// Total Residue Volume
vector<float> s87;// Total Unpacked Atoms 
vector<float> s88;
vector<float> s89;
//Summary Scores 
//packing prob - for the residue
//volume prob 

};                                       


#endif                                     

