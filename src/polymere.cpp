#include <boost/lexical_cast.hpp>
#include "algos.h"
#include <algorithm>
#include "pdb.h"
#include "fabase.h"
#include "doublet.h"
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "polymere.h"
#include "base.h"
#include "fragclass.h"
#include "matrix.h"
#include "fragatoms.h"
#include "fragall.h"
#include "fragment.h"
#include <istream>
#include <sstream>
#include "parse.h"
#include "parser.h"
#include "score.h"
#include "voro/voro++.cc"
/* RMSD.c Program http://www.personal.leeds.ac.uk/~bgy1mm/Bioinformatics/rmsd.html*/
struct PairCompareAscPol
{
    bool operator() (const std::pair<float,int>  x, const std::pair<float,int>  y) const {
//To Minimize a score       
return x.first > y.first;
//To maximize a score
//return x.score < y.score;
 }
};
struct PairComparePol
{
    bool operator() (const Polymere *  x, const Polymere * y) const {
//To Minimize a score       
 return x->score > y->score;
//To maximize a score
//return x.score < y.score;
    }
};

typedef struct
{
  float m[4][4];
} MATRIX;

#define vdiff2(a,b) ( ((a)[0]-(b)[0]) * ((a)[0]-(b)[0]) +	\
  ((a)[1]-(b)[1]) * ((a)[1]-(b)[1]) + \
  ((a)[2]-(b)[2]) * ((a)[2]-(b)[2]) )

static float alignedrmsd(float *v1, float *v2, int N);
static void centroid(float *ret, float *v, int N);
static int getalignmtx(float *v1, float *v2, int N, MATRIX *mtx);
static void crossproduct(float *ans, float *pt1, float *pt2);
static void mtx_root(MATRIX *mtx);
static int almostequal(MATRIX *a, MATRIX *b);
static void mulpt(MATRIX *mtx, float *pt);
static void mtx_mul(MATRIX *ans, MATRIX *x, MATRIX *y);
static void mtx_identity(MATRIX *mtx);
static void mtx_trans(MATRIX *mtx, float x, float y, float z);
static int mtx_invert(float *mtx, int N);
static float absmaxv(float *v, int N);
static float gaussFreq(float x,int type);




/*
  calculate rmsd between two structures
  Params: v1 - first set of points
          v2 - second set of points
          N - number of points
          mtx - return for transfrom matrix used to align structures
  Returns: rmsd score
  Notes: mtx can be null. Transform will be rigid. Inputs must
         be previously aligned for sequence alignment
 */



float rmsd(float *v1, float *v2, int N, float *mtx)
{
  float cent1[3];
  float cent2[3];
  MATRIX tmtx;
  MATRIX tempmtx;
  MATRIX move1;
  MATRIX move2;
  int i;
  float answer;
  float *temp1 = 0;
  float *temp2 = 0;
  int err;

  assert(N > 3);

  temp1 =(float*) malloc(N * 3 * sizeof(float));
  temp2 =(float*) malloc(N * 3 * sizeof(float));
  if(!temp1 || !temp2)
    goto error_exit;

  centroid(cent1, v1, N);
  centroid(cent2, v2, N);
  for(i=0;i<N;i++)
  {
    temp1[i*3+0] = v1[i*3+0] - cent1[0];
    temp1[i*3+1] = v1[i*3+1] - cent1[1];
    temp1[i*3+2] = v1[i*3+2] - cent1[2];
    
    temp2[i*3+0] = v2[i*3+0] - cent2[0];
    temp2[i*3+1] = v2[i*3+1] - cent2[1];
    temp2[i*3+2] = v2[i*3+2] - cent2[2];
  }

  err = getalignmtx(temp1, temp2, N, &tmtx);
  if(err == -1)
    goto error_exit;
 
  mtx_trans(&move1, -cent2[0], -cent2[1], -cent2[2]);
  mtx_mul(&tempmtx, &move1, &tmtx);
  mtx_trans(&move2, cent1[0], cent1[1], cent1[2]);
  mtx_mul(&tmtx, &tempmtx, &move2);
  memcpy(temp2, v2, N * sizeof(float) * 3);
  for(i=0;i<N;i++)
    mulpt(&tmtx, temp2 + i * 3);
//set v2 to temp2
  memcpy(v2, temp2, N * sizeof(float) * 3);

  answer = alignedrmsd(v1, temp2, N);
  free(temp1);
  free(temp2);
  if(mtx)
    memcpy(mtx, &tmtx.m, 16 * sizeof(float));

  return answer;
 error_exit:
  free(temp1);
  free(temp2);
  if(mtx)
  {
    for(i=0;i<16;i++)
      mtx[i] = 0;
  }
  return sqrt(-1.0);
}

/*
  calculate rmsd between two aligned structures (trivial)
  Params: v1 - first structure
          v2 - second structure
          N - number of points
  Returns: rmsd
 */
static float alignedrmsd(float *v1, float *v2, int N)
{
  float answer =0;
  int i;

  for(i=0;i<N;i++)
    answer += vdiff2(v1 + i *3, v2 + i * 3);
  return sqrt(answer/N);
}

/*
  compute the centroid
 */
static void centroid(float *ret, float *v, int N)
{
  int i;

  ret[0] = 0;
  ret[1] = 0;
  ret[2] = 0;
  for(i=0;i<N;i++)
  {
    ret[0] += v[i*3+0];
    ret[1] += v[i*3+1];
    ret[2] += v[i*3+2];
  }
  ret[0] /= N;
  ret[1] /= N;
  ret[2] /= N;
}

/*
  get the matrix needed to align two structures
  Params: v1 - reference structure
          v2 - structure to align
          N - number of points
          mtx - return for rigid body alignment matrix
  Notes: only calculates rotation part of matrix.
         assumes input has been aligned to centroids 
 */
static int getalignmtx(float *v1, float *v2, int N, MATRIX *mtx)
{
  MATRIX A = { {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,1}} };
  MATRIX At;
  MATRIX Ainv;
  MATRIX temp;
  float tv[3];
  float tw[3];
  float tv2[3];
  float tw2[3];
  int k, i, j;
  int flag = 0;
  float correction;

  correction = absmaxv(v1, N * 3) * absmaxv(v2, N * 3);
  
  for(k=0;k<N;k++)
    for(i=0;i<3;i++)
      for(j=0;j<3;j++)
	A.m[i][j] += (v1[k*3+i] * v2[k*3+j])/correction;
  
  while(flag < 3)
  {
    for(i=0;i<4;i++)
      for(j=0;j<4;j++)
        At.m[i][j] = A.m[j][i];

    memcpy(&Ainv, &A, sizeof(MATRIX));
    /* this will happen if all points are in a plane */
    if( mtx_invert((float *) &Ainv, 4) == -1)
    {
      if(flag == 0)
      {
        crossproduct(tv, v1, v1+3);
        crossproduct(tw, v2, v2+3);
      }
      else
      {
        crossproduct(tv2, tv, v1);
        crossproduct(tw2, tw, v2);
        memcpy(tv, tv2, 3 * sizeof(float));
        memcpy(tw, tw2, 3 * sizeof(float));
      }
      for(i=0;i<3;i++)
        for(j=0;j<3;j++)
	  A.m[i][j] += tv[i] * tw[j];
      
      flag++;
    }
    else
      flag = 5;
  }
  if(flag != 5)
    return -1;

  mtx_mul(&temp, &At, &A);
  mtx_root(&temp);
  mtx_mul(mtx, &temp, &Ainv); 
  return 0;
}

/*
  get the crossproduct of two vectors.
  Params: ans - return pinter for answer.
          pt1 - first vector
		  pt2 - second vector.
  Notes: crossproduct is at right angles to the two vectors.
*/
static void crossproduct(float *ans, float *pt1, float *pt2)
{
  ans[0] = pt1[1] * pt2[2] - pt1[2] * pt2[1];
  ans[1] = pt1[0] * pt2[2] - pt1[2] * pt2[0];
  ans[2] = pt1[0] * pt2[1] - pt1[1] * pt2[0];
}

/*
  Denman-Beavers square root iteration
 */
static void mtx_root(MATRIX *mtx)
{
  MATRIX Y = *mtx;
  MATRIX Z;
  MATRIX Y1;
  MATRIX Z1;
  MATRIX invY;
  MATRIX invZ;
  MATRIX Y2;
  int iter = 0;
  int i, ii;

  mtx_identity(&Z);

  do
  {
    invY = Y;
    invZ = Z;
    if( mtx_invert((float *) &invY, 4) == -1)
      return;
    if( mtx_invert((float *) &invZ, 4) == -1)
      return;
    for(i=0;i<4;i++)
      for(ii=0;ii<4;ii++)
      {
        Y1.m[i][ii] = 0.5 * (Y.m[i][ii] + invZ.m[i][ii]);
	Z1.m[i][ii] = 0.5 * (Z.m[i][ii] + invY.m[i][ii]);
      }
    Y = Y1;
    Z = Z1;

    mtx_mul(&Y2, &Y, &Y);
  }
  while(!almostequal(&Y2, mtx) && iter++ < 20 );

  *mtx = Y;
}

/*
  Check two matrices for near-enough equality
  Params: a - first matrix
          b - second matrix
  Returns: 1 if almost equal, else 0, epsilon 0.0001f.
 */
static int almostequal(MATRIX *a, MATRIX *b)
{
  int i, ii;
  float epsilon = 0.001f;

  for(i=0;i<4;i++)
    for(ii=0;ii<4;ii++)
      if(fabs(a->m[i][ii] - b->m[i][ii]) > epsilon) 
        return 0;
  return 1;
}  

/*
  multiply a point by a matrix.
  Params: mtx - matrix
          pt - the point (transformed)
*/
static void mulpt(MATRIX *mtx, float *pt)
{
  float ans[4] = {0};
  int i;
  int ii;

  for(i=0;i<4;i++)
  {
    for(ii=0;ii<3;ii++)
	{
	  ans[i] += pt[ii] * mtx->m[ii][i];
	}
	ans[i] += mtx->m[3][i];
  }
  pt[0] = ans[0];
  pt[1] = ans[1];
  pt[2] = ans[2];
} 

/*
  multiply two matrices.
  Params: ans - return pointer for answer.
          x - first matrix
		  y - second matrix.
  Notes: ans may not be equal to x or y.
*/
static void mtx_mul(MATRIX *ans, MATRIX *x, MATRIX *y)
{
  int i;
  int ii;
  int iii;

  for(i=0;i<4;i++)
    for(ii=0;ii<4;ii++)
	{
	  ans->m[i][ii] = 0;
	  for(iii=0;iii<4;iii++)
	    ans->m[i][ii] += x->m[i][iii] * y->m[iii][ii];
    }
}


/*
  create an identity matrix.
  Params: mtx - return pointer.
*/
static void mtx_identity(MATRIX *mtx)
{
  int i;
  int ii;

  for(i=0;i<4;i++)
    for(ii=0;ii<4;ii++)
	{
	  if(i==ii)
	    mtx->m[i][ii] = 1.0f;
	  else
	    mtx->m[i][ii] = 0;
	}
}

/*
  create a translation matrix.
  Params: mtx - return pointer for matrix.
          x - x translation.
		  y - y translation.
		  z - z translation
*/
static void mtx_trans(MATRIX *mtx, float x, float y, float z)
{
  mtx->m[0][0] = 1;
  mtx->m[0][1] = 0;
  mtx->m[0][2] = 0;
  mtx->m[0][3] = 0;

  mtx->m[1][0] = 0;
  mtx->m[1][1] = 1;
  mtx->m[1][2] = 0;
  mtx->m[1][3] = 0;
  
  mtx->m[2][0] = 0;
  mtx->m[2][1] = 0;
  mtx->m[2][2] = 1;
  mtx->m[2][3] = 0;
  
  mtx->m[3][0] = x;
  mtx->m[3][1] = y;
  mtx->m[3][2] = z;
  mtx->m[3][3] = 1; 
}

/*
   matrix invert routine
  Params: mtx - the matrix in raw format, in/out
          N - width and height
  Returns: 0 on success, -1 on fail
 */
static int mtx_invert(float *mtx, int N)
{
  int indxc[100]; /* these 100s are the only restriction on matrix size */
  int indxr[100];
  int ipiv[100];
  int i, j, k;
  int irow, icol;
  double big;
  double pinv;
  int l, ll;
  double dum;
  double temp;
  
  assert(N <= 100);

  for(i=0;i<N;i++)
    ipiv[i] = 0;

  for(i=0;i<N;i++)
  {
    big = 0.0;

    /* find biggest element */
    for(j=0;j<N;j++)
      if(ipiv[j] != 1)
        for(k=0;k<N;k++)
          if(ipiv[k] == 0)
            if(fabs(mtx[j*N+k]) >= big)
	    {
	       big = fabs(mtx[j*N+k]);
               irow = j;
               icol = k;
	    }       
	      
    ipiv[icol]=1;

    if(irow != icol)
      for(l=0;l<N;l++)
      {
        temp = mtx[irow * N + l];
        mtx[irow * N + l] = mtx[icol * N + l];
	mtx[icol * N + l] = temp;
      }

    indxr[i] = irow;
    indxc[i] = icol;

       
    /* if biggest element is zero matrix is singular, bail */
    if(mtx[icol* N + icol] == 0)
       goto error_exit;
		  
    pinv = 1.0/mtx[icol * N + icol];
           
    mtx[icol * N + icol] = 1.0;
 
    for(l=0;l<N;l++)
      mtx[icol * N + l] *= pinv;
                   
    for(ll=0;ll<N;ll++)
      if(ll != icol)
      {
        dum = mtx[ll * N + icol];
        mtx[ll * N + icol] = 0.0;
        for(l=0;l<N;l++) 
          mtx[ll * N + l] -= mtx[icol * N + l]*dum;
      }
  }                
   

  /* unscramble matrix */
  for (l=N-1;l>=0;l--) 
  {
    if (indxr[l] != indxc[l])
    for (k=0;k<N;k++)
    {
      temp = mtx[k * N + indxr[l]];
      mtx[k * N + indxr[l]] = mtx[k * N + indxc[l]];
      mtx[k * N + indxc[l]] = temp;
    }
  } 

  return 0;

 error_exit:
  return -1;
}

/*
  get the asolute maximum of an array
 */
static float absmaxv(float *v, int N)
{
  float answer;
  int i;

  for(i=0;i<N;i++)
    if(answer < fabs(v[i]))
      answer = fabs(v[i]);
  return answer;
}

#include <stdio.h>

/*
  debug utlitiy
 */

static void printmtx(FILE *fp, MATRIX *mtx)
{
  int i, ii;

  for(i=0;i<4;i++)
  {
    for(ii=0;ii<4;ii++)
      fprintf(fp, "%f, ", mtx->m[i][ii]);
    fprintf(fp, "\n");
  }
}

/*
int rmsdmain(void)
{
  float one[4*3] = {0,0,0, 1,0,0, 2,1,0, 0,3,1};
  float two[4*3] = {0,0,0, 0,1,0, 1,2,0, 3,0,1};
  MATRIX mtx;
  double diff;
  int i;

  diff = rmsd(one, two, 4, (float *) &mtx.m);
  printf("%f\n", diff);
  printmtx(stdout, &mtx);
  for(i=0;i<4;i++)
  {
    mulpt(&mtx, two + i * 3);
    printf("%f %f %f\n", two[i*3], two[i*3+1], two[i*3+2]);
  }
  return 0;
}
*/
using namespace std;
const int VECSIZE = 6000;

Polymere::Polymere() {
using namespace std;
double z=0;
for (int k=0; k<VECSIZE; k++) {
Base * b = new Base();
b-> z     =z;
b-> x     =z;
b-> y     =z;
b-> group =0;
b-> prob = 0;
b-> seq   ='N';
b-> bx    =z;
b-> by    =z;
b-> bz    =z;
b-> b2x   =z;
b-> b2y   =z;
b-> b2z   =z;
b-> id    =0;
mols.push_back(b);
}
score=0;
num_mols=0;
randtot=0;
double x=0;
for(int i=0; i < (VECSIZE*81); i++){
full_atoms.push_back(x);
}
}


Polymere::Polymere(string seq,const std::vector< Fragment *>  frags,Fragclass * dat,std::vector< Fragatoms *> atoms,std::vector<Fragall *> batoms,std::vector<double> * dparams){

//Initializes a polymere to a random single strand 
int seq_size = seq.size();
int failed=0;
for(int resnum=0; resnum<seq_size+8; resnum++){
//Generate random coordinates
double zc = rand()/((double)RAND_MAX)+5*(double)resnum;
double xc = rand()/((double)RAND_MAX)*2;
double yc = rand()/((double)RAND_MAX)*2;
Base * b = new Base();
b-> z     =zc;
b-> x     =xc;
b-> y     =yc;
b-> group =0;
b-> prob = 0;
b-> seq   ='A';
b-> bx    =0;
b-> by    =0;
b-> bz    =0;
b-> b2x   =0;
b-> b2y   =0;
b-> b2z   =0;
b-> id    =0;
//add 78 elements to the full_atoms
for(int psw=0; psw<81; psw++){
full_atoms.push_back(0);
} 
if(resnum > 3 && resnum < (seq_size + 4)){b-> seq   =seq[(resnum-4)];}
mols.push_back(b);
}
//printf("Sequence: \n"); 
//for(int m=0; m<mols.size(); m++){
//printf("%c ",mols[m]->seq);
//}
//printf("\n"); 
num_mols=seq_size+8;
//Set the initial group
for(int m=0; m<mols.size(); m++){mols[m]->group=0;}
mols[1]->group=5;
mols[2]->group=5;
mols[3]->group=5;
int i=1;
printf("Smoothing Backbone\n");
Polymere * tmp=new Polymere();
Polymere * outpol=new Polymere();
for(int att=0; att<100; att++){
int iter=0;
outpol->copy(this);
i=1;
while(i < (seq_size+6)){
tmp->copy(outpol);
int fid=chooseSingle(i,dat);
if(batoms[fid]->check()==1 ){
tmp->changeFrag2(i,fid,frags);
tmp->coordSysNt(i,&atoms[fid]->atoms);
tmp->updateFull();
tmp->coordSysAT(i,&batoms[fid]->atoms);
if(tmp->groupCheckBackBone(5) < 1 && tmp->groupCheckClashes(5) < 1){
//printf("%d %c -- updated\n",i,tmp->mols[i]->seq);
i++;
if(i < (seq_size+6)){tmp->mols[(i+2)]->group=5;}
outpol->copy(tmp);
}
}
iter++;
if(iter > 5000){
failed=1;
break; }
}
if(iter < 5000){this->copy(outpol); break;}else{
//printf("Smoothing again: %d\n",att);
}
}
//Identify midpoint
this->genStacked(frags,dat,atoms,batoms);
}



Polymere::Polymere(string file,int model) {
using namespace std;
score=0;
num_mols=0;
randtot=0;
Pdb * pdb = new Pdb(file);
//pdb->print();
//Convert to a the pdb format used in the simulator
//printf("Standarizing Model: %d\n",model);
pdb->standard("./pdb.map");
//printf("Loading Model: %d\n",model);
//pdb->print();
//printf("Number of models: %d\n",pdb->models.size());
//Get the model, chain, and residues
int a=0;
for(int resnums=0;resnums < pdb->models[model]->chains[0]->residues.size(); resnums++){
//printf("Residue: %d %d\n",resnums,pdb->models[model]->chains.size());
//Create a new resiude
Base * b = new Base();
b-> z    = 0;
b-> x     =0;
b-> y     =0;
b-> group =0;
b-> prob = 0;
b-> seq   ='N';
b-> bx    =0;
b-> by    =0;
b-> bz    =0;
b-> b2x   =0;
b-> b2y   =0;
b-> b2z   =0;
b-> id    =0;
//add 78 elements to the full_atoms
for(int psw=0; psw<81; psw++){
full_atoms.push_back(0);
} 
//printf("Created Base: %d\n",resnums);
//Loop through all the resiudes 
string rname=pdb->models[model]->chains[0]->residues[resnums]->res_name;
//printf("Resname: '%s'",rname.c_str());
if(rname.compare("  A")==0 || rname.compare("  G")==0){
vector<float> c1 = pdb->models[model]->chains[0]->residues[resnums]->getAtom(" C1'"); 
vector<float> c2 = pdb->models[model]->chains[0]->residues[resnums]->getAtom(" N9 "); 
vector<float> c3 = pdb->models[model]->chains[0]->residues[resnums]->getAtom(" C4 "); 
//Add coordinates 
b->x=c1[0];
b->y=c1[1];
b->z=c1[2];

b->bx=c2[0];
b->by=c2[1];
b->bz=c2[2];

b->b2x=c3[0];
b->b2y=c3[1];
b->b2z=c3[2];
if(rname.compare("  A")==0){
b->seq='A';
int bpos = a*81;
for(int atn = 0; atn < 26; atn++){
string at_name = pdb->models[model]->chains[0]->residues[resnums]->atomnames2[atn];
c1 = pdb->models[model]->chains[0]->residues[resnums]->getAtom(at_name);
full_atoms[bpos+ 0]=c1[0];
full_atoms[bpos+ 1]=c1[1];
full_atoms[bpos+ 2]=c1[2];
bpos = bpos+3;
}
}else{
if(rname.compare("  G")==0){
int bpos = a*81;
b->seq='G';
for(int atn = 0; atn < 27; atn++){
string at_name = pdb->models[model]->chains[0]->residues[resnums]->atomnames1[atn];
c1 = pdb->models[model]->chains[0]->residues[resnums]->getAtom(at_name);
full_atoms[bpos+ 0]=c1[0];
full_atoms[bpos+ 1]=c1[1];
full_atoms[bpos+ 2]=c1[2];
bpos = bpos+3;
}


}

}


}else{

if(rname.compare("  U")==0 || rname.compare("  C")==0){
vector<float> c1 = pdb->models[model]->chains[0]->residues[resnums]->getAtom(" C1'"); 
vector<float> c2 = pdb->models[model]->chains[0]->residues[resnums]->getAtom(" N1 "); 
vector<float> c3 = pdb->models[model]->chains[0]->residues[resnums]->getAtom(" C2 "); 

b->x=c1[0];
b->y=c1[1];
b->z=c1[2];

b->bx=c2[0];
b->by=c2[1];
b->bz=c2[2];

b->b2x=c3[0];
b->b2y=c3[1];
b->b2z=c3[2];
//Add atoms to full_atoms
if(rname.compare("  U")==0){
b->seq='U';
int bpos=a*81;
for(int atn = 0; atn < 23; atn++){
string at_name = pdb->models[model]->chains[0]->residues[resnums]->atomnames3[atn];
c1 = pdb->models[model]->chains[0]->residues[resnums]->getAtom(at_name);
full_atoms[bpos+ 0]=c1[0];
full_atoms[bpos+ 1]=c1[1];
full_atoms[bpos+ 2]=c1[2];
bpos = bpos+3;
}

}else{

if(rname.compare("  C")==0){
b->seq='C';
int bpos=a*81;
for(int atn = 0; atn < 24; atn++){
string at_name = pdb->models[model]->chains[0]->residues[resnums]->atomnames4[atn];
c1 = pdb->models[model]->chains[0]->residues[resnums]->getAtom(at_name);
full_atoms[bpos+ 0]=c1[0];
full_atoms[bpos+ 1]=c1[1];
full_atoms[bpos+ 2]=c1[2];
bpos = bpos+3;
}


}

}


}else{
//Failed to load residue
printf("Failed to load PDB: residue names dont' match(A C G U)");
//abort();

}
}
//Get the C1', N9 and C4

//Get the C1', N1 and C2

//Foreach residue in the file extract the C1' atoms coordinates and the base cooridnates and add to mols vector
a++;
mols.push_back(b);
}
//Loop through all residues and load all coordinates into the full_atoms vector 
num_mols=a;
while(full_atoms.size()<(VECSIZE*81)){
full_atoms.push_back(0);
}

//Delete PDB
for(int m=0; m<pdb->models.size(); m++){
for(int c=0; c<pdb->models[m]->chains.size(); c++){
for(int r=0; r<pdb->models[m]->chains[c]->residues.size(); r++){
for(int a=0; a<pdb->models[m]->chains[c]->residues[r]->atoms.size(); a++){
delete pdb->models[m]->chains[c]->residues[r]->atoms[a];
}
delete pdb->models[m]->chains[c]->residues[r];
}
delete pdb->models[m]->chains[c];
}
delete pdb->models[m];
} 
////////////////
delete pdb;

}

Polymere::Polymere(string file) {
using namespace std;
score=0;
num_mols=0;
randtot=0;
Pdb * pdb = new Pdb(file);
printf("read into pdb..\n");
//Convert to a the pdb format used in the simulator
pdb->standard("./pdb.map");
printf("standardized %d residues\n",pdb->models[0]->chains[0]->residues.size());

//Get the model, chain, and residues
int a=0;
for(int resnums=0;resnums < pdb->models[0]->chains[0]->residues.size(); resnums++){
//Create a new resiude
Base * b = new Base();
b-> z    = 0;
b-> x     =0;
b-> y     =0;
b-> group =0;
b-> prob = 0;
b-> seq   ='N';
b-> bx    =0;
b-> by    =0;
b-> bz    =0;
b-> b2x   =0;
b-> b2y   =0;
b-> b2z   =0;
b-> id    =0;
//add 78 elements to the full_atoms
for(int psw=0; psw<81; psw++){
full_atoms.push_back(0);
} 

//Loop through all the resiudes 
string rname=pdb->models[0]->chains[0]->residues[resnums]->res_name;
printf("Resname: '%s'",rname.c_str());
if(rname.compare("  A")==0 || rname.compare("  G")==0){
vector<float> c1 = pdb->models[0]->chains[0]->residues[resnums]->getAtom(" C1'"); 
vector<float> c2 = pdb->models[0]->chains[0]->residues[resnums]->getAtom(" N9 "); 
vector<float> c3 = pdb->models[0]->chains[0]->residues[resnums]->getAtom(" C4 "); 
//Add coordinates 
b->x=c1[0];
b->y=c1[1];
b->z=c1[2];

b->bx=c2[0];
b->by=c2[1];
b->bz=c2[2];

b->b2x=c3[0];
b->b2y=c3[1];
b->b2z=c3[2];
if(rname.compare("  A")==0){
b->seq='A';
int bpos = a*81;
for(int atn = 0; atn < 26; atn++){
string at_name = pdb->models[0]->chains[0]->residues[resnums]->atomnames2[atn];
c1 = pdb->models[0]->chains[0]->residues[resnums]->getAtom(at_name);
full_atoms[bpos+ 0]=c1[0];
full_atoms[bpos+ 1]=c1[1];
full_atoms[bpos+ 2]=c1[2];
bpos = bpos+3;
}
}else{
if(rname.compare("  G")==0){
int bpos = a*81;
b->seq='G';
for(int atn = 0; atn < 27; atn++){
string at_name = pdb->models[0]->chains[0]->residues[resnums]->atomnames1[atn];
c1 = pdb->models[0]->chains[0]->residues[resnums]->getAtom(at_name);
full_atoms[bpos+ 0]=c1[0];
full_atoms[bpos+ 1]=c1[1];
full_atoms[bpos+ 2]=c1[2];
bpos = bpos+3;
}
}
}


}else{

if(rname.compare("  U")==0 || rname.compare("  C")==0){
vector<float> c1 = pdb->models[0]->chains[0]->residues[resnums]->getAtom(" C1'"); 
vector<float> c2 = pdb->models[0]->chains[0]->residues[resnums]->getAtom(" N1 "); 
vector<float> c3 = pdb->models[0]->chains[0]->residues[resnums]->getAtom(" C2 "); 

b->x=c1[0];
b->y=c1[1];
b->z=c1[2];

b->bx=c2[0];
b->by=c2[1];
b->bz=c2[2];

b->b2x=c3[0];
b->b2y=c3[1];
b->b2z=c3[2];
//Add atoms to full_atoms
if(rname.compare("  U")==0){
b->seq='U';
int bpos=a*81;
for(int atn = 0; atn < 23; atn++){
string at_name = pdb->models[0]->chains[0]->residues[resnums]->atomnames3[atn];
c1 = pdb->models[0]->chains[0]->residues[resnums]->getAtom(at_name);
full_atoms[bpos+ 0]=c1[0];
full_atoms[bpos+ 1]=c1[1];
full_atoms[bpos+ 2]=c1[2];
bpos = bpos+3;
}
}else{

if(rname.compare("  C")==0){
b->seq='C';
int bpos=a*81;
for(int atn = 0; atn < 24; atn++){
string at_name = pdb->models[0]->chains[0]->residues[resnums]->atomnames4[atn];
c1 = pdb->models[0]->chains[0]->residues[resnums]->getAtom(at_name);
full_atoms[bpos+ 0]=c1[0];
full_atoms[bpos+ 1]=c1[1];
full_atoms[bpos+ 2]=c1[2];
bpos = bpos+3;
}


}

}


}else{
//Failed to load residue
printf("Failed to load PDB: residue names dont' match(A C G U)");
//abort();

}
}
//Get the C1', N9 and C4

//Get the C1', N1 and C2

//Foreach residue in the file extract the C1' atoms coordinates and the base cooridnates and add to mols vector
a++;
mols.push_back(b);
}
//Loop through all residues and load all coordinates into the full_atoms vector 
num_mols=a;
while(full_atoms.size()<(VECSIZE*81)){
full_atoms.push_back(0);
}

//Delete PDB
for(int m=0; m<pdb->models.size(); m++){
for(int c=0; c<pdb->models[m]->chains.size(); c++){
for(int r=0; r<pdb->models[m]->chains[c]->residues.size(); r++){
for(int a=0; a<pdb->models[m]->chains[c]->residues[r]->atoms.size(); a++){
delete pdb->models[m]->chains[c]->residues[r]->atoms[a];
}
delete pdb->models[m]->chains[c]->residues[r];
}
delete pdb->models[m]->chains[c];
}
delete pdb->models[m];
} 
////////////////
delete pdb;
}

void Polymere::alignHelix(Polymere* p2,std::vector<int> helix1, std::vector<int> helix2){
//helix1 matches the Polymere object calling alignHelix
//helix2 matches the p2 Polymere object

//Check the helix1 is a helix
//Check that the size of helix1 is the same as the size of helix2

/*
The object of alignResidues is find the tranformation matrix 
that minimizes the root mean square distance between the coordinates 
(mols[helix1[i]]->x,mols[helix1[i]]->y,mols[helix1[i]]->z) and (p2->mols[helix1[i]]->x,p2->mols[helix1[i]]->y,p2->mols[helix1[i]]->z)
after identifying the transformation matrix the full_atom vector should be 
transformed to refect the alingment.  An example is given below that 
uses the rmsd function to align the coordinates
*/

float d;
float *v1,*v2;
float *mtx;
v1=NULL;
v2=NULL;
mtx=NULL;
mtx=(float*) malloc(16*sizeof(float));
  v1 =(float*) malloc(num_mols * 9 * sizeof(float));
  v2 =(float*) malloc(num_mols * 9 * sizeof(float));
if(v1==NULL || v2==NULL || mtx==NULL){
d=-1;
cout << "Allocation Failed" << endl;
}else{
int i=0;
for(int j=0; j<num_mols; j++){
//if(p2->mols[j]->group==1 || p2->mols[j]->group==2){
float x1 = (float)mols[i]->x;
float y1 = (float)mols[i]->y;
float z1 = (float)mols[i]->z;
float xb1 = (float)mols[i]->bx;
float yb1 = (float)mols[i]->by;
float zb1 = (float)mols[i]->bz;
float xb21 = (float)mols[i]->b2x;
float yb21 = (float)mols[i]->b2y;
float zb21 = (float)mols[i]->b2z;
v1[i*9+0]=x1;
v1[i*9+1]=y1;
v1[i*9+2]=z1;
v1[i*9+3]=xb1;
v1[i*9+4]=yb1;
v1[i*9+5]=zb1;
v1[i*9+6]=xb21;
v1[i*9+7]=yb21;
v1[i*9+8]=zb21;

float x2 =(float) p2->mols[i]->x;
float y2 =(float) p2->mols[i]->y;
float z2 =(float) p2->mols[i]->z;
float xb2 = (float)p2->mols[i]->bx;
float yb2 = (float)p2->mols[i]->by;
float zb2 = (float)p2->mols[i]->bz;
float xb22 = (float)p2->mols[i]->b2x;
float yb22 = (float)p2->mols[i]->b2y;
float zb22 = (float)p2->mols[i]->b2z;
v2[i*9+0]=x2;
v2[i*9+1]=y2;
v2[i*9+2]=z2;
v2[i*9+3]=xb2;
v2[i*9+4]=yb2;
v2[i*9+5]=zb2;
v2[i*9+6]=xb22;
v2[i*9+7]=yb22;
v2[i*9+8]=zb22;
i++;
//}
}
//cout << "Groups Numbers: " << i << endl;
d=rmsd(v1,v2,i,mtx);
printf("%8.3f\n",d);
//update vectors
for(int j=0; j<num_mols; j++){
p2->mols[j]->x=v2[j*9+0];
p2->mols[j]->y=v2[j*9+1];
p2->mols[j]->z=v2[j*9+2];
p2->mols[j]->bx=v2[j*9+3];
p2->mols[j]->by=v2[j*9+4];
p2->mols[j]->bz=v2[j*9+5];
p2->mols[j]->b2x=v2[j*9+6];
p2->mols[j]->b2y=v2[j*9+7];
p2->mols[j]->b2z=v2[j*9+8];

mols[j]->x=v1[j*9+0];
mols[j]->y=v1[j*9+1];
mols[j]->z=v1[j*9+2];
mols[j]->bx=v1[j*9+3];
mols[j]->by=v1[j*9+4];
mols[j]->bz=v1[j*9+5];
mols[j]->b2x=v1[j*9+6];
mols[j]->b2y=v1[j*9+7];
mols[j]->b2z=v1[j*9+8];

//printf("%8.3f %8.3f %8.3f, %8.3f %8.3f %8.3f\n",mols[j]->x, mols[j]->y, mols[j]->z, p2->mols[j]->x, p2->mols[j]->y, p2->mols[j]->z);
}
this->updateFull();
p2->updateFull();  
free(v1);
  free(v2);
free(mtx);
}
}







void Polymere::clearMols(){
//Deletes all pointers to the base clase molecule
for(int m=0; m<mols.size(); m++){delete mols[m];}
mols.clear();
full_atoms.clear();
}

void Polymere::printMols(int time){
using namespace std;

for(int j = 0; j < num_mols; j++){
cout<<time<<","<<mols[j]->x<< "," << mols[j]->y << "," << mols[j]->z <<","<<mols[j]->bx<< "," << mols[j]->by << "," << mols[j]->bz  <<","<<mols[j]->b2x<< "," << mols[j]->b2y << "," << mols[j]->b2z  <<"," << mols[j]->seq<< ","<< mols[j]->group<<endl;
}
}
void Polymere::topov(std::string outf,double camera_distX,double  camera_distY,double camera_distZ){
int bpos=0;
std::vector<std::vector<double> > C4_spheres;
for(int m=0; m<this->num_mols; m++){
bpos=m*81;
std::vector<double> points;
points.push_back(this->full_atoms[bpos+3*5+0]);
points.push_back(this->full_atoms[bpos+3*5+1]);
points.push_back(this->full_atoms[bpos+3*5+2]);
C4_spheres.push_back(points);
}

std::vector<std::vector<double> > C1_spheres;
for(int m=0; m<this->num_mols; m++){
bpos=m*81;
std::vector<double> points;
points.push_back(this->full_atoms[bpos+3*11+0]);
points.push_back(this->full_atoms[bpos+3*11+1]);
points.push_back(this->full_atoms[bpos+3*11+2]);
C1_spheres.push_back(points);
}

std::vector<std::vector<double> > Base_polygons;
for(int m=0; m<this->num_mols; m++){
bpos=m*81;
if(this->mols[m]->seq=='G'){
std::vector<double> points;
//C4
points.push_back(this->full_atoms[bpos+3*13+0]);
points.push_back(this->full_atoms[bpos+3*13+1]);
points.push_back(this->full_atoms[bpos+3*13+2]);
//N3
points.push_back(this->full_atoms[bpos+3*14+0]);
points.push_back(this->full_atoms[bpos+3*14+1]);
points.push_back(this->full_atoms[bpos+3*14+2]);
//C2
points.push_back(this->full_atoms[bpos+3*15+0]);
points.push_back(this->full_atoms[bpos+3*15+1]);
points.push_back(this->full_atoms[bpos+3*15+2]);
//N1
points.push_back(this->full_atoms[bpos+3*17+0]);
points.push_back(this->full_atoms[bpos+3*17+1]);
points.push_back(this->full_atoms[bpos+3*17+2]);
//C6
points.push_back(this->full_atoms[bpos+3*18+0]);
points.push_back(this->full_atoms[bpos+3*18+1]);
points.push_back(this->full_atoms[bpos+3*18+2]);

//C5
points.push_back(this->full_atoms[bpos+3*20+0]);
points.push_back(this->full_atoms[bpos+3*20+1]);
points.push_back(this->full_atoms[bpos+3*20+2]);

//N7
points.push_back(this->full_atoms[bpos+3*21+0]);
points.push_back(this->full_atoms[bpos+3*21+1]);
points.push_back(this->full_atoms[bpos+3*21+2]);

//C8
points.push_back(this->full_atoms[bpos+3*22+0]);
points.push_back(this->full_atoms[bpos+3*22+1]);
points.push_back(this->full_atoms[bpos+3*22+2]);

//N9
points.push_back(this->full_atoms[bpos+3*12+0]);
points.push_back(this->full_atoms[bpos+3*12+1]);
points.push_back(this->full_atoms[bpos+3*12+2]);
Base_polygons.push_back(points);
points.clear();
}
if(this->mols[m]->seq=='A'){
std::vector<double> points;
//C4
points.push_back(this->full_atoms[bpos+3*13+0]);
points.push_back(this->full_atoms[bpos+3*13+1]);
points.push_back(this->full_atoms[bpos+3*13+2]);
//N3
points.push_back(this->full_atoms[bpos+3*14+0]);
points.push_back(this->full_atoms[bpos+3*14+1]);
points.push_back(this->full_atoms[bpos+3*14+2]);
//C2
points.push_back(this->full_atoms[bpos+3*15+0]);
points.push_back(this->full_atoms[bpos+3*15+1]);
points.push_back(this->full_atoms[bpos+3*15+2]);
//N1
points.push_back(this->full_atoms[bpos+3*16+0]);
points.push_back(this->full_atoms[bpos+3*16+1]);
points.push_back(this->full_atoms[bpos+3*16+2]);
//C6
points.push_back(this->full_atoms[bpos+3*17+0]);
points.push_back(this->full_atoms[bpos+3*17+1]);
points.push_back(this->full_atoms[bpos+3*17+2]);
//C5
points.push_back(this->full_atoms[bpos+3*19+0]);
points.push_back(this->full_atoms[bpos+3*19+1]);
points.push_back(this->full_atoms[bpos+3*19+2]);
//N7
points.push_back(this->full_atoms[bpos+3*20+0]);
points.push_back(this->full_atoms[bpos+3*20+1]);
points.push_back(this->full_atoms[bpos+3*20+2]);
//C8
points.push_back(this->full_atoms[bpos+3*21+0]);
points.push_back(this->full_atoms[bpos+3*21+1]);
points.push_back(this->full_atoms[bpos+3*21+2]);
//N9
points.push_back(this->full_atoms[bpos+3*12+0]);
points.push_back(this->full_atoms[bpos+3*12+1]);
points.push_back(this->full_atoms[bpos+3*12+2]);
Base_polygons.push_back(points);
points.clear();
}

if(this->mols[m]->seq=='U'){
std::vector<double> points;
//N1
points.push_back(this->full_atoms[bpos+3*12+0]);
points.push_back(this->full_atoms[bpos+3*12+1]);
points.push_back(this->full_atoms[bpos+3*12+2]);
//C2
points.push_back(this->full_atoms[bpos+3*13+0]);
points.push_back(this->full_atoms[bpos+3*13+1]);
points.push_back(this->full_atoms[bpos+3*13+2]);
//N3
points.push_back(this->full_atoms[bpos+3*16+0]);
points.push_back(this->full_atoms[bpos+3*16+1]);
points.push_back(this->full_atoms[bpos+3*16+2]);
//C4
points.push_back(this->full_atoms[bpos+3*17+0]);
points.push_back(this->full_atoms[bpos+3*17+1]);
points.push_back(this->full_atoms[bpos+3*17+2]);
//C5
points.push_back(this->full_atoms[bpos+3*19+0]);
points.push_back(this->full_atoms[bpos+3*19+1]);
points.push_back(this->full_atoms[bpos+3*19+2]);

//C6
points.push_back(this->full_atoms[bpos+3*14+0]);
points.push_back(this->full_atoms[bpos+3*14+1]);
points.push_back(this->full_atoms[bpos+3*14+2]);
Base_polygons.push_back(points);
points.clear();
}

if(this->mols[m]->seq=='C'){
std::vector<double> points;
//N1
points.push_back(this->full_atoms[bpos+3*12+0]);
points.push_back(this->full_atoms[bpos+3*12+1]);
points.push_back(this->full_atoms[bpos+3*12+2]);
//C2
points.push_back(this->full_atoms[bpos+3*13+0]);
points.push_back(this->full_atoms[bpos+3*13+1]);
points.push_back(this->full_atoms[bpos+3*13+2]);
//N3
points.push_back(this->full_atoms[bpos+3*16+0]);
points.push_back(this->full_atoms[bpos+3*16+1]);
points.push_back(this->full_atoms[bpos+3*16+2]);
//C4
points.push_back(this->full_atoms[bpos+3*17+0]);
points.push_back(this->full_atoms[bpos+3*17+1]);
points.push_back(this->full_atoms[bpos+3*17+2]);
//C5
points.push_back(this->full_atoms[bpos+3*19+0]);
points.push_back(this->full_atoms[bpos+3*19+1]);
points.push_back(this->full_atoms[bpos+3*19+2]);
//C6
points.push_back(this->full_atoms[bpos+3*14+0]);
points.push_back(this->full_atoms[bpos+3*14+1]);
points.push_back(this->full_atoms[bpos+3*14+2]);
Base_polygons.push_back(points);
points.clear();
}
}

FILE * povout=fopen(outf.c_str(),"w");
std::string pheader=povHeader();
std::string pcam   =povCamera(camera_distX,camera_distY,camera_distZ,0,camera_distY,0);
fprintf(povout,"%s%s",pheader.c_str(),pcam.c_str());
for(int s=0; s<C4_spheres.size(); s++){
std::string scol;
if(this->mols[s]->seq=='G'){scol.append("G_Texture");}
if(this->mols[s]->seq=='A'){scol.append("A_Texture");}
if(this->mols[s]->seq=='U'){scol.append("U_Texture");}
if(this->mols[s]->seq=='C'){scol.append("C_Texture");}
std::string sphere_pov=povSphere(C4_spheres[s],1,scol);
fprintf(povout,"%s",sphere_pov.c_str());
}

for(int s=0; s<C1_spheres.size(); s++){
std::string scol;
int start=0;
if(this->mols[s]->seq=='G'){scol.append("G_Texture");start=8;}
if(this->mols[s]->seq=='A'){scol.append("A_Texture");start=8;}
if(this->mols[s]->seq=='U'){scol.append("U_Texture");}
if(this->mols[s]->seq=='C'){scol.append("C_Texture");}
std::string sphere_pov=povSphere(C1_spheres[s],0.5,scol);
std::string scol2("Cylinder_Texture");
std::vector<double> ends;
for(int n=0; n<3; n++){ends.push_back(C1_spheres[s][n]);}
for(int n=0; n<3; n++){ends.push_back(Base_polygons[s][(start*3)+n]);}
std::string cyl_pov=povCylinder(ends,0.25,scol2);
fprintf(povout,"%s%s",sphere_pov.c_str(),cyl_pov.c_str());

}

for(int s=0; s<Base_polygons.size(); s++){
std::string scol;
if(this->mols[s]->seq=='G'){scol.append("G_Texture");std::string sphere_pov=povBaseAG(Base_polygons[s],scol);fprintf(povout,"%s",sphere_pov.c_str());}
if(this->mols[s]->seq=='A'){scol.append("A_Texture");std::string sphere_pov=povBaseAG(Base_polygons[s],scol);fprintf(povout,"%s",sphere_pov.c_str());}
if(this->mols[s]->seq=='U'){scol.append("U_Texture");std::string sphere_pov=povBaseCU(Base_polygons[s],scol);fprintf(povout,"%s",sphere_pov.c_str());}
if(this->mols[s]->seq=='C'){scol.append("C_Texture");std::string sphere_pov=povBaseCU(Base_polygons[s],scol);fprintf(povout,"%s",sphere_pov.c_str());}
}

for(int s=0; s<(C4_spheres.size()-1); s++){
std::vector<double> ends;
for(int n=0; n<3; n++){ends.push_back(C4_spheres[s][n]);}
for(int n=0; n<3; n++){ends.push_back(C4_spheres[s+1][n]);}
std::string scol("Cylinder_Texture");
std::string cyl_pov=povCylinder(ends,0.5,scol);
//std::string povCylinder(std::vector<float> points,float radius,std::string color);
fprintf(povout,"%s",cyl_pov.c_str());
}

for(int s=0; s<(C4_spheres.size()); s++){
std::vector<double> ends;
for(int n=0; n<3; n++){ends.push_back(C4_spheres[s][n]);}
for(int n=0; n<3; n++){ends.push_back(C1_spheres[s][n]);}
std::string scol("Cylinder_Texture");
std::string cyl_pov=povCylinder(ends,0.25,scol);
//std::string povCylinder(std::vector<float> points,float radius,std::string color);
fprintf(povout,"%s",cyl_pov.c_str());
}

fflush(povout);
fclose(povout);
}

void Polymere::genCons(std::vector< std::pair<int,int> > * cons,int p1,int lt,int asym_num,std::string rnaseq){
//Check that lenght small_side.size + small_size==1 = e2-s2
int lcheck=0;
int success;
int m=0;
std::vector<int> small_side;
int s_side=makeSmallSide(&small_side,rnaseq,p1,lt,asym_num);
//Side 1 Ads for to adjust for the extra sequence added to the ends 
int s1 = p1+4;
int incrLoop=7+lt;
//Side 2
int s2 = s1+incrLoop;
if(s_side==1){//Small side is position 1
for(int j=0; j<(small_side.size()); j++){
std::vector< std::pair<int,int> > bpcons;
int k=0;
int icr=1;	
int icr1=1;	
	if(small_side.at(j)==0){
		std::pair<int,int> scon((s1-j+0),(s2+m+0));	
		cons->push_back(scon);
		icr=1;
		icr1=1;
		}else{
		std::pair<int,int> scon1((s1-j+0),(s2+m+0));	
		cons->push_back(scon1);
		std::pair<int,int> scon2((s1-j+0),(s2+m+1));	
		cons->push_back(scon2);
		icr=2;
		icr1=2;
		}

m=m+icr1;
}
}else{
if(s_side==2){
for(int j=0; j<(small_side.size()); j++){
std::vector< std::pair<int,int> > bpcons;
int k=0;
int icr=1;	
int icr1=1;	
	if(small_side.at(j)==0){
		std::pair<int,int> scon((s1-m+0),(s2+j+0));	
		cons->push_back(scon);
		icr=1;
		icr1=1;
		}else{
		std::pair<int,int> scon1((s1-m+0),(s2+j+0));	
		cons->push_back(scon1);
		std::pair<int,int> scon2((s1-m-1),(s2+j+0));	
		cons->push_back(scon2);
		icr=2;
		icr1=2;
		}
m=m+icr1;
}


}else{
printf("Small side must be either 1 or 2 for buildBuldge\n");
}
}
//Print out constraints for debugging purpose
for(int j=0; j<cons->size(); j++){
	printf("pos %d\n",j);
		printf("\t%d-%d\n",cons->at(j).first,cons->at(j).second);
}

}

void Polymere::genAllCons(int nmax,std::vector<std::vector< std::pair<int,int> > > * allcons,string rnaseq){
std::vector<int> small_side;
//Write the top nmax structures
int len_tot = rnaseq.length();
for(int lt=0; lt<2; lt++){
for(int p1=0; p1<len_tot; p1++){
int s_side=1;
int asym_num=0;

std::vector<std::pair<float,int> > scores;
std::vector<std::vector<int> > params;
int asym_i=0;
//Foreach position find the top nmax scoring structures
//limit to 5k loops
int tot=0;
while(s_side >0 && tot < 5000){
tot++;
s_side=makeSmallSide(&small_side,rnaseq,p1,lt,asym_num);
if(s_side >0){
float score=scoreSmallSide(&small_side,rnaseq,p1,lt,s_side);
std::vector<int> pms;
pms.push_back(s_side);
pms.push_back(p1);
pms.push_back(lt);
pms.push_back(rnaseq.length());
pms.push_back(asym_num);
//printf("%8.3f small_side: %d pos1: %d loop_type: %d: seqlength: %d sslen %d\n",score,s_side,p1,lt,opt->sequence.length(),small_side.size());
std::pair<float,int> plsc(score,asym_i++);
scores.push_back(plsc);
params.push_back(pms);
}
asym_num++;
}
std::sort(scores.begin(), scores.end(), PairCompareAscPol ());
///////////////////////////////////////////////////////////
//                Get first n_max                        //
///////////////////////////////////////////////////////////
int mnum=0;
if(nmax > scores.size()){mnum=scores.size();}else{mnum=nmax;}
for(int n=0; n<mnum; n++){
std::vector< std::pair<int,int> > ncon;
this->genCons(&ncon,params[scores[n].second][1],params[scores[n].second][2],params[scores[n].second][4],rnaseq);
allcons->push_back(ncon);
}
}
}
}


std::vector<std::pair<int,int> > Polymere::confMatch(std::vector< std::vector< std::pair<int,int> > >  confs){
std::vector< std::pair<int,int> > lastcons1;
std::vector< std::pair<int,int> > lastcons2;
//Returns the integer position of the conformation that contains the most matches 
//To a vector of std::pair in confs
std::vector<std::pair<int,int> > outvect;
std::vector<int> cons_scores;
//Generate the current conf
int cur_score,min_score,lastc;
lastc=5;
for(int g=0; g<confs.size(); g++){
int mdiff = abs((confs[g][0].first-confs[g].back().first) - (confs[g].back().second-confs[g][0].second));
//printf("mdiff: %d\n",mdiff);
lastcons1.push_back(confs[g].back());
lastcons2.push_back(confs[g][(confs[g].size()-2)]);
int lastc1=this->hardConstraints(&lastcons1);
int lastc2=this->hardConstraints(&lastcons2);
if(lastc1<2 || lastc2<2){lastc=0;}
lastcons1.clear();
lastcons2.clear();

cur_score=this->hardConstraints(&confs[g]);

cur_score = cur_score - mdiff;
//printf("cur_score: %d\n",cur_score);
if(g==0){min_score=cur_score;}
cons_scores.push_back(cur_score);
if(min_score > cur_score){min_score=cur_score;}
}
//Identify the minimal score and return it
for(int g=0; g<cons_scores.size(); g++){
if(cons_scores[g]==min_score){
std::pair<int,int> pcon(g,lastc);
outvect.push_back(pcon);
}
}
//printf("min_score: %d\n",min_score);
return outvect;
}


void Polymere::cartesian_to_sperical(double x,double y, double z, double * r,double * angy,double * angz){
//returns values in r,angy, angz
using namespace std;
*r = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
if(*r==0){
*angy=0;
*angz=0;
}else{
*angy = atan2(y,x);
*angz = acos(z/(*r));
} 
}

/*
void Polymere::pivot(float anchor_x,float anchor_y,float anchor_z,float x1,float y1,float z1,float x1p,float y1p,float z1p){
//Rotates the point x1,y1,z1 into the point x1p,y1p,z1p
//while keeping the anchor point fixed
//If it can't be rotated into the position it will lie along the vector between 
//the anchor and the new point






}

*/




void Polymere::rotate3d(double a,double b,double g,double x2,double y2,double z2)
{
using namespace std;
//For each molecule in the polymere rotate it around x2,y2,z2
for(int j = 0; j < num_mols; j++){
mols[j]->rotate(a,b,g,x2,y2,z2);
}
}

void Polymere::rotateX(double ang)
{
using namespace std;
//For each molecule in the polymere rotate it around x2,y2,z2
for(int j = 0; j < num_mols; j++){
mols[j]->rotateX(ang);
}
}

void Polymere::rotateXr(double ang,int pos)
{
using namespace std;
//For each molecule in the polymere rotate it around x2,y2,z2
for(int j = pos; j >= 0; j=j-1){
mols[j]->rotateX(ang);
}
}
void Polymere::rotateXb(double ang,int pos)
{
using namespace std;
//For each molecule in the polymere rotate it around x2,y2,z2
for(int j = pos; j < num_mols; j++){
mols[j]->rotateX(ang);
}
}

void Polymere::rotateY(double ang)
{
using namespace std;
//For each molecule in the polymere rotate it around x2,y2,z2
for(int j = 0; j < num_mols; j++){
mols[j]->rotateY(ang);
}
}

void Polymere::rotateZb(double ang,int pos)
{
using namespace std;
//For each molecule in the polymere rotate it around x2,y2,z2
for(int j = pos; j >=0; j=j-1){
mols[j]->rotateZ(ang);
}
}
void Polymere::rotateZf(double ang,int pos)
{
using namespace std;
//For each molecule in the polymere rotate it around x2,y2,z2
for(int j = pos; j <num_mols; j++){
mols[j]->rotateZ(ang);
}
}
void Polymere::rotateZ(double ang)
{
using namespace std;
//For each molecule in the polymere rotate it around x2,y2,z2
for(int j = 0; j < num_mols; j++){
mols[j]->rotateZ(ang);
}
}

void Polymere::translater(int pos,double x2,double y2,double z2)
{
	using namespace std;
//For each molecule in the polymere rotate it around x2,y2,z2
for(int j = pos; j >= 0; j=j-1){
mols[j]->translate(x2,y2,z2);
}
}
void Polymere::translateb(int pos,double x2,double y2,double z2)
{
	using namespace std;
//For each molecule in the polymere rotate it around x2,y2,z2
for(int j = pos; j < num_mols; j++){
mols[j]->translate(x2,y2,z2);
}
}

void Polymere::translate(double x2,double y2,double z2)
{
	using namespace std;
//For each molecule in the polymere rotate it around x2,y2,z2
for(int j = 0; j < num_mols; j++){
mols[j]->translate(x2,y2,z2);
}
}

void Polymere::changeBondLenR(int pos, double length){
using namespace std;
//Changes the lenght between pos and pos+1
//Generate the unit vector in the given direction
double x1 = mols[pos]->x;
double y1 = mols[pos]->y;
double z1 = mols[pos]->z;
int pos2=pos-1;
double x2 = mols[pos2]->x;
double y2 = mols[pos2]->y;
double z2 = mols[pos2]->z;

double x3=(x2-x1);
double y3=(y2-y1);
double z3=(z2-z1);
double r = sqrt(pow(x3,2)+pow(y3,2)+pow(z3,2));
length=length-r;
x3=-1*x3/r*length;
y3=-1*y3/r*length;
z3=-1*z3/r*length;
translater(pos2,x3,y3,z3);
}
void Polymere::changeBondLen2r(int pos, double length){
using namespace std;
//Changes the lenght between pos and pos+1

//Generate the unit vector in the given direction
int pos2=pos+1;
int pos1=pos;

double x1 = mols[pos2]->x;
double y1 = mols[pos2]->y;
double z1 = mols[pos2]->z;

double x2 = mols[pos1]->x;
double y2 = mols[pos1]->y;
double z2 = mols[pos1]->z;

double x3=(x2-x1);
double y3=(y2-y1);
double z3=(z2-z1);
double r = sqrt(pow(x3,2)+pow(y3,2)+pow(z3,2));
length=length-r;
x3=-1*x3/r*length;
y3=-1*y3/r*length;
z3=-1*z3/r*length;
translater(pos1,x3,y3,z3);
}
void Polymere::changeBondLen(int pos, double length){
using namespace std;
//Changes the lenght between pos and pos+1

//Generate the unit vector in the given direction
double x1 = mols[pos]->x;
double y1 = mols[pos]->y;
double z1 = mols[pos]->z;
int pos2=pos+1;
double x2 = mols[pos2]->x;
double y2 = mols[pos2]->y;
double z2 = mols[pos2]->z;

double x3=(x2-x1);
double y3=(y2-y1);
double z3=(z2-z1);
double r = sqrt(pow(x3,2)+pow(y3,2)+pow(z3,2));
length=length-r;
x3=-1*x3/r*length;
y3=-1*y3/r*length;
z3=-1*z3/r*length;
translateb(pos2,x3,y3,z3);
}
void Polymere::changeBond3a(int pos,double ang){
using namespace std;
//determine angle around x-axis to rotate bond
double unx1,uny1,unz1,rotx,roty,rotz; 
align4(pos,&rotx,&roty,&rotz,&unx1,&uny1,&unz1);
double Bx,By,Bz,Ax,Ay,Az;
int p3=pos+1;
Bx=mols[p3]->x;
By=mols[p3]->y;
Bz=mols[p3]->z;
p3=pos;
Ax=mols[p3]->bx;
Ay=mols[p3]->by;
Az=mols[p3]->bz;
double AA=sqrt(pow(Ax,2)+pow(Ay,2)+pow(Az,2));
double BB=sqrt(pow(Bx,2)+pow(By,2)+pow(Bz,2));
//double val =(AA*BB*cos(ang)-Bx*Ax)/By;
//val=(val+sqrt(pow(Ay,2)+pow(Az,2)))/(sqrt(pow(Az,2)+pow(Ay,2))*2);
//			cout<< "Ax: " << Ax << "Ay: " << Ay << "Az: " << Az << endl;
//			cout <<"Bx: "<< Bx << "By: "<<  By << "Bz: " << Bz << endl;
//			if(val > 0){
//val=acos(-1*sqrt(val));
//			}else{val=0;}
//Brute Force
//       		        cout <<	"Angle of rotation " << val <<  endl;
//Rotate Bond around X by acos(val)
//Determine what the current angle is
double angy=atan2(Az,Ay);
double bang=0;
double ax=0;
double mang=9999;
for(int a=0; a<=200; a++){
ax=a*2*PI/2000;
p3=pos;
mols[p3]->rotateX(2*PI/200);
p3=pos+1;
Bx=mols[p3]->x;
By=mols[p3]->y;
Bz=mols[p3]->z;
p3=pos;
Ax=mols[p3]->bx;
Ay=mols[p3]->by;
Az=mols[p3]->bz;
//cout<< "Ax: " << Ax << "Ay: " << Ay << "Az: " << Az << endl;
//cout <<"Bx: "<< Bx << "By: "<<  By << "Bz: " << Bz << endl;
//cout <<	"Angle of rotation " << val <<  endl;
//Calculate angle between the two vectors
AA=sqrt(pow(Ax,2)+pow(Ay,2)+pow(Az,2));
BB=sqrt(pow(Bx,2)+pow(By,2)+pow(Bz,2));
double angout = acos((Ax*Bx+Ay*By+Az*Bz)/(AA*BB));
angout=abs(angout-ang);
if(angout < mang){
	mang=angout;
bang=ax;
}
}
mols[p3]->rotateX(-1*bang);
p3=pos+1;
Bx=mols[p3]->x;
By=mols[p3]->y;
Bz=mols[p3]->z;
p3=pos;
Ax=mols[p3]->bx;
Ay=mols[p3]->by;
Az=mols[p3]->bz;
//cout<< "Ax: " << Ax << "Ay: " << Ay << "Az: " << Az << endl;
//cout <<"Bx: "<< Bx << "By: "<<  By << "Bz: " << Bz << endl;
//cout <<	"Angle of rotation " << val <<  endl;
//Calculate angle between the two vectors
AA=sqrt(pow(Ax,2)+pow(Ay,2)+pow(Az,2));
BB=sqrt(pow(Bx,2)+pow(By,2)+pow(Bz,2));
double angout = acos((Ax*Bx+Ay*By+Az*Bz)/(AA*BB));
//cout << "Desired: " << angout << " Actual: " << ang << endl;
}
void Polymere::changeBond3(int pos,double ang)
{
using namespace std;
//Bond angle should not eqaul zero or PI Otherwise it throws off the torsion angle 
//align first two along positive x-axis with initial point in the negative x positive y plane
if(ang==0){ang=0.00000000001;}
if(ang==PI){ang=PI-0.00000000001;}
double unx1,uny1,unz1,rotx,roty,rotz; 
align3(pos,&rotx,&roty,&rotz,&unx1,&uny1,&unz1);
int p2 = pos+1;
double x1=mols[p2]->x;
double y1=mols[p2]->y;
double z1=mols[p2]->z;
double angy = atan2(y1,x1);
ang=angy-ang;
rotateZf(ang,p2);
//cout << "bond3: " << ang << endl;
//Move back to original position 
//rotateX(rotx);
//rotateY(roty);
//rotateZ(rotz);
//translate(unx1,uny1,unz1);
}

void Polymere::changeBond1(int pos,double ang)
{
using namespace std;
//Bond angle should not eqaul zero or PI Otherwise it throws off the torsion angle 
//align first two along positive x-axis with initial point in the negative x positive y plane
if(ang==0){ang=0.00000000001;}
if(ang==PI){ang=PI-0.00000000001;}
double unx1,uny1,unz1,rotx,roty,rotz; 
align1(pos,&rotx,&roty,&rotz,&unx1,&uny1,&unz1);
int p2 = pos;
double x1=mols[p2]->bx;
double y1=mols[p2]->by;
double z1=mols[p2]->bz;
double angy = atan2(y1,x1);
//cout << "ang: " << ang << " angy: " << angy;
ang=angy-ang;
//cout << " ang2: " << ang<< endl;
mols[p2]->rotateZ(ang);
//cout << "bond1: " << ang << endl;
//Move back to original position 
//rotateX(rotx);
//rotateY(roty);
//rotateZ(rotz);
//translate(unx1,uny1,unz1);
}

void Polymere::changeBond2r(int pos,double ang)
{
using namespace std;
//Bond angle should not eqaul zero or PI Otherwise it throws off the torsion angle 
//align first two along positive x-axis with initial point in the negative x positive y plane
if(ang==0){ang=0.0000000001;}
if(ang==PI){ang=PI-0.00000000001;}
double unx1,uny1,unz1,rotx,roty,rotz; 
align3to5(pos,&rotx,&roty,&rotz,&unx1,&uny1,&unz1);
int p2 = pos-1;
double x1=mols[p2]->x;
double y1=mols[p2]->y;
double z1=mols[p2]->z;
double angy = atan2(y1,x1);
ang=angy-ang;
rotateZb(ang,p2);
}

void Polymere::changeBond2(int pos,double ang)
{
using namespace std;
//Bond angle should not eqaul zero or PI Otherwise it throws off the torsion angle 
//align first two along positive x-axis with initial point in the negative x positive y plane
if(ang==0){ang=0.0000000001;}
if(ang==PI){ang=PI-0.00000000001;}
double unx1,uny1,unz1,rotx,roty,rotz; 
align(pos,&rotx,&roty,&rotz,&unx1,&uny1,&unz1);
int p2 = pos+1;
double x1=mols[p2]->x;
double y1=mols[p2]->y;
double z1=mols[p2]->z;
double angy = atan2(y1,x1);
ang=angy-ang;
rotateZf(ang,p2);
}



void Polymere::changeBaseTor(int pos,double ang)
{
//Align entire molecule to the first three points
//Rotate all points before pos1 around the x-axis
using namespace std;
//Align Structure 
//Translate first position to origin
double x1=mols[pos]->x;
double y1=mols[pos]->y;
double z1=mols[pos]->z;

double unx1=-1*x1;
double uny1=-1*y1;
double unz1=-1*z1;

translate(x1,y1,z1);
//rotate second position to x-axis
double r;
double angz;
double angy;
int p2 = pos;
x1=mols[p2]->bx;
y1=mols[p2]->by;
z1=mols[p2]->bz;
//find angles between pos2 and x-axis
cartesian_to_sperical(x1,y1,z1,&r,&angy,&angz);
//rotate structure so that x1,y1,z1 is on the x-axis 
double unrotz = -1*angy;
rotateZ(angy);
//Need to set position 2 along the positive x-xis 
angz=PI/2-angz;

double unroty = -1*angz;
rotateY(angz);
//Rotate Position -1 into positive Y plane
p2 = pos-1;
y1=mols[p2]->y;
x1=mols[p2]->x;
z1=mols[p2]->z;
angz = atan2(y1,z1);
angz=angz-PI/2;
double unrotx = -1*angz;
rotateX(angz);

//Rotate to the new torsion angle
p2=pos;
x1=mols[p2]->b2x;
y1=mols[p2]->b2y;
z1=mols[p2]->b2z;
angz =atan2(y1,z1);
angz=angz-PI/2+ang;
mols[p2]->rotateX(angz);
//Undo rotations 
//rotateX(unrotx);
//rotateY(unroty);
//rotateZ(unrotz);
//translate(unx1,uny1,unz1);
}
void Polymere::changeTorR(int pos,double ang)
{
//Align entire molecule to the first three points
//Rotate all points before pos1 around the x-axis
using namespace std;
//Align Structure 
//Translate first position to origin
double x1=mols[pos]->x;
double y1=mols[pos]->y;
double z1=mols[pos]->z;

double unx1=-1*x1;
double uny1=-1*y1;
double unz1=-1*z1;

translate(x1,y1,z1);
//rotate second position to x-axis
double r;
double angz;
double angy;
int p2 = pos-1;
x1=mols[p2]->x;
y1=mols[p2]->y;
z1=mols[p2]->z;
//find angles between pos2 and x-axis
cartesian_to_sperical(x1,y1,z1,&r,&angy,&angz);
//rotate structure so that x1,y1,z1 is on the x-axis 
double unrotz = -1*angy;
rotateZ(angy);
//Need to set position 2 along the positive x-xis 
angz=PI/2-angz;

double unroty = -1*angz;
rotateY(angz);
//Rotate Position -1 into positive Y plane
p2 = pos+1;
y1=mols[p2]->y;
x1=mols[p2]->x;
z1=mols[p2]->z;
angz = atan2(y1,z1);
angz=angz-PI/2;
double unrotx = -1*angz;
rotateX(angz);

//Rotate to the new torsion angle
p2=pos-2;
x1=mols[p2]->x;
y1=mols[p2]->y;
z1=mols[p2]->z;
angz =atan2(y1,z1);
angz=angz-PI/2+ang;
rotateXr(angz,p2);
}
void Polymere::changeTor2r(int pos,double ang)
{
//Align entire molecule to the first three points
//Rotate all points before pos1 around the x-axis
using namespace std;
//Align Structure 
//Translate first position to origin
double x1=mols[pos]->x;
double y1=mols[pos]->y;
double z1=mols[pos]->z;

double unx1=-1*x1;
double uny1=-1*y1;
double unz1=-1*z1;

translate(x1,y1,z1);
//rotate second position to x-axis
double r;
double angz;
double angy;
int p2 = pos-1;
x1=mols[p2]->x;
y1=mols[p2]->y;
z1=mols[p2]->z;
//find angles between pos2 and x-axis
cartesian_to_sperical(x1,y1,z1,&r,&angy,&angz);
//rotate structure so that x1,y1,z1 is on the x-axis 
double unrotz = -1*angy;
rotateZ(angy);
//Need to set position 2 along the positive x-xis 
angz=PI/2-angz;

double unroty = -1*angz;
rotateY(angz);
//Rotate Position -1 into positive Y plane
p2 = pos+1;
y1=mols[p2]->y;
x1=mols[p2]->x;
z1=mols[p2]->z;
angz = atan2(y1,z1);
angz=angz-PI/2;
double unrotx = -1*angz;
rotateX(angz);

//Rotate to the new torsion angle
p2=pos-2;
x1=mols[p2]->x;
y1=mols[p2]->y;
z1=mols[p2]->z;
angz =atan2(y1,z1);
angz=(angz-PI/2+ang);
rotateXr(angz,p2);
}

void Polymere::changeTor2(int pos,double ang)
{
//Align entire molecule to the first three points
//Rotate all points before pos1 around the x-axis
using namespace std;
//Align Structure 
//Translate first position to origin
double x1=mols[pos]->x;
double y1=mols[pos]->y;
double z1=mols[pos]->z;

double unx1=-1*x1;
double uny1=-1*y1;
double unz1=-1*z1;

translate(x1,y1,z1);
//rotate second position to x-axis
double r;
double angz;
double angy;
int p2 = pos+1;
x1=mols[p2]->x;
y1=mols[p2]->y;
z1=mols[p2]->z;
//find angles between pos2 and x-axis
cartesian_to_sperical(x1,y1,z1,&r,&angy,&angz);
//rotate structure so that x1,y1,z1 is on the x-axis 
double unrotz = -1*angy;
rotateZ(angy);
//Need to set position 2 along the positive x-xis 
angz=PI/2-angz;

double unroty = -1*angz;
rotateY(angz);
//Rotate Position -1 into positive Y plane
p2 = pos-1;
y1=mols[p2]->y;
x1=mols[p2]->x;
z1=mols[p2]->z;
angz = atan2(y1,z1);
angz=angz-PI/2;
double unrotx = -1*angz;
rotateX(angz);

//Rotate to the new torsion angle
p2=pos+2;
x1=mols[p2]->x;
y1=mols[p2]->y;
z1=mols[p2]->z;
angz =atan2(y1,z1);
angz=angz-PI/2+ang;
rotateXb(angz,p2);
}

void Polymere::changeTor(int pos,double ang)
{
//Align entire molecule to the first three points
//Rotate all points before pos1 around the x-axis
using namespace std;
//Align Structure 
//Translate first position to origin
double x1=mols[pos]->x;
double y1=mols[pos]->y;
double z1=mols[pos]->z;

double unx1=-1*x1;
double uny1=-1*y1;
double unz1=-1*z1;

translate(x1,y1,z1);
//rotate second position to x-axis
double r;
double angz;
double angy;
int p2 = pos+1;
x1=mols[p2]->x;
y1=mols[p2]->y;
z1=mols[p2]->z;
//find angles between pos2 and x-axis
cartesian_to_sperical(x1,y1,z1,&r,&angy,&angz);
//rotate structure so that x1,y1,z1 is on the x-axis 
double unrotz = -1*angy;
rotateZ(angy);
//Need to set position 2 along the positive x-xis 
angz=PI/2-angz;

double unroty = -1*angz;
rotateY(angz);
//Rotate Position -1 into positive Y plane
p2 = pos;
y1=mols[p2]->by;
x1=mols[p2]->bx;
z1=mols[p2]->bz;
angz = atan2(y1,z1);
angz=angz-PI/2;
double unrotx = -1*angz;
rotateX(angz);

//Rotate to the new torsion angle
p2=pos+1;
x1=mols[p2]->bx;
y1=mols[p2]->by;
z1=mols[p2]->bz;
angz =atan2(y1,z1);
angz=angz-PI/2+ang;
rotateXb(angz,p2);
//ang=-1*ang;
//rotateXb(ang,p2);
//Undo rotations 
//rotateX(unrotx);
//rotateY(unroty);
//rotateZ(unrotz);
//translate(unx1,uny1,unz1);
}

void Polymere::align1(int pos,double * unrotx, double * unroty, double * unrotz, double * unx1, double * uny1, double * unz1){
//align structures and returns angles and translations to move back to original position after modification 
//Align entire molecule to the first three points
//Rotate all points before pos1 around the x-axis
using namespace std;
//Align Structure 
//Translate first position to origin
double x1=mols[pos]->x;
double y1=mols[pos]->y;
double z1=mols[pos]->z;

*unx1=-1*x1;
*uny1=-1*y1;
*unz1=-1*z1;

translate(x1,y1,z1);
//rotate second position to x-axis
double r;
double angz;
double angy;
int p2 = pos-1;
x1=mols[p2]->x;
y1=mols[p2]->y;
z1=mols[p2]->z;
//find angles between pos2 and x-axis
cartesian_to_sperical(x1,y1,z1,&r,&angy,&angz);
//rotate structure so that x1,y1,z1 is on the x-axis 
*unrotz = -1*angy;
rotateZ(angy);
//Need to set position 2 along the positive x-xis 
angz=PI/2-angz;
*unroty = -1*angz;
rotateY(angz);
//Rotate Position -1 into positive Y plane
p2 = pos;
y1=mols[p2]->by;
x1=mols[p2]->bx;
z1=mols[p2]->bz;
angz = atan2(y1,z1);
angz=angz-PI/2;
*unrotx = -1*angz;
rotateX(angz);
}
void Polymere::align4(int pos,double * unrotx, double * unroty, double * unrotz, double * unx1, double * uny1, double * unz1){
//align structures and returns angles and translations to move back to original position after modification 
//Align entire molecule to the first three points
//Rotate all points before pos1 around the x-axis
using namespace std;
//Align Structure 
//Translate first position to origin
double x1=mols[pos]->x;
double y1=mols[pos]->y;
double z1=mols[pos]->z;

*unx1=-1*x1;
*uny1=-1*y1;
*unz1=-1*z1;

translate(x1,y1,z1);
//rotate second position to x-axis
double r;
double angz;
double angy;
int p2 = pos-1;
x1=mols[p2]->x;
y1=mols[p2]->y;
z1=mols[p2]->z;
//find angles between pos2 and x-axis
cartesian_to_sperical(x1,y1,z1,&r,&angy,&angz);
//rotate structure so that x1,y1,z1 is on the x-axis 
*unrotz = -1*angy;
rotateZ(angy);
//Need to set position 2 along the positive x-xis 
angz=PI/2-angz;
*unroty = -1*angz;
rotateY(angz);
//Rotate Position -1 into positive Y plane
p2 = pos+1;
y1=mols[p2]->y;
x1=mols[p2]->x;
z1=mols[p2]->z;
angz = atan2(y1,z1);
angz=angz-PI/2;
*unrotx = -1*angz;
rotateX(angz);
}
void Polymere::align3(int pos,double * unrotx, double * unroty, double * unrotz, double * unx1, double * uny1, double * unz1){
//align structures and returns angles and translations to move back to original position after modification 
//Align entire molecule to the first three points
//Rotate all points before pos1 around the x-axis
using namespace std;
//Align Structure 
//Translate first position to origin
double x1=mols[pos]->x;
double y1=mols[pos]->y;
double z1=mols[pos]->z;

*unx1=-1*x1;
*uny1=-1*y1;
*unz1=-1*z1;

translate(x1,y1,z1);
//rotate second position to x-axis
double r;
double angz;
double angy;
int p2 = pos;
x1=mols[p2]->bx;
y1=mols[p2]->by;
z1=mols[p2]->bz;
//find angles between pos2 and x-axis
cartesian_to_sperical(x1,y1,z1,&r,&angy,&angz);
//rotate structure so that x1,y1,z1 is on the x-axis 
*unrotz = -1*angy;
rotateZ(angy);
//Need to set position 2 along the positive x-xis 
angz=PI/2-angz;
*unroty = -1*angz;
rotateY(angz);
//Rotate Position -1 into positive Y plane
p2 = pos+1;
y1=mols[p2]->y;
x1=mols[p2]->x;
z1=mols[p2]->z;
angz = atan2(y1,z1);
angz=angz-PI/2;
*unrotx = -1*angz;
rotateX(angz);
}

void Polymere::align3to5(int pos,double * unrotx, double * unroty, double * unrotz, double * unx1, double * uny1, double * unz1){
//align structures and returns angles and translations to move back to original position after modification 
//Align entire molecule to the first three points
//Rotate all points before pos1 around the x-axis
using namespace std;
//Align Structure 
//Translate first position to origin
double x1=mols[pos]->x;
double y1=mols[pos]->y;
double z1=mols[pos]->z;

*unx1=-1*x1;
*uny1=-1*y1;
*unz1=-1*z1;

translate(x1,y1,z1);
//rotate second position to x-axis
double r;
double angz;
double angy;
int p2 = pos+1;
x1=mols[p2]->x;
y1=mols[p2]->y;
z1=mols[p2]->z;
//find angles between pos2 and x-axis
cartesian_to_sperical(x1,y1,z1,&r,&angy,&angz);
//rotate structure so that x1,y1,z1 is on the x-axis 
*unrotz = -1*angy;
rotateZ(angy);
//Need to set position 2 along the positive x-xis 
angz=PI/2-angz;
*unroty = -1*angz;
rotateY(angz);
//Rotate Position -1 into positive Y plane
p2 = pos-1;
y1=mols[p2]->y;
x1=mols[p2]->x;
z1=mols[p2]->z;
angz = atan2(y1,z1);
angz=angz-PI/2;
*unrotx = -1*angz;
rotateX(angz);
}

void Polymere::align(int pos,double * unrotx, double * unroty, double * unrotz, double * unx1, double * uny1, double * unz1){
//align structures and returns angles and translations to move back to original position after modification 
//Align entire molecule to the first three points
//Rotate all points before pos1 around the x-axis
using namespace std;
//Align Structure 
//Translate first position to origin
double x1=mols[pos]->x;
double y1=mols[pos]->y;
double z1=mols[pos]->z;

*unx1=-1*x1;
*uny1=-1*y1;
*unz1=-1*z1;

translate(x1,y1,z1);
//rotate second position to x-axis
double r;
double angz;
double angy;
int p2 = pos-1;
x1=mols[p2]->x;
y1=mols[p2]->y;
z1=mols[p2]->z;
//find angles between pos2 and x-axis
cartesian_to_sperical(x1,y1,z1,&r,&angy,&angz);
//rotate structure so that x1,y1,z1 is on the x-axis 
*unrotz = -1*angy;
rotateZ(angy);
//Need to set position 2 along the positive x-xis 
angz=PI/2-angz;
*unroty = -1*angz;
rotateY(angz);
//Rotate Position -1 into positive Y plane
p2 = pos+1;
y1=mols[p2]->y;
x1=mols[p2]->x;
z1=mols[p2]->z;
angz = atan2(y1,z1);
angz=angz-PI/2;
*unrotx = -1*angz;
rotateX(angz);
}

void Polymere::alignStruct2(int pos){
using namespace std;
//Translate first position to origin
double x1=mols[pos]->x;
double y1=mols[pos]->y;
double z1=mols[pos]->z;
translate(x1,y1,z1);
//rotate second position to x-axis
double r;
double angz;
double angy;
int p2 = pos-1;
x1=mols[p2]->x;
y1=mols[p2]->y;
z1=mols[p2]->z;
//find angles between pos2 and x-axis
cartesian_to_sperical(x1,y1,z1,&r,&angy,&angz);
//rotate structure so that x1,y1,z1 is on the x-axis 
rotateZ(angy);
//Need to set position 2 along the positive x-xis 
angz=PI/2-angz;
rotateY(angz);
//Rotate Position -1 into positive XY plane
p2 = pos+1;
y1=mols[p2]->y;
x1=mols[p2]->x;
z1=mols[p2]->z;
angz = atan2(y1,z1);
angz=angz-PI/2;
rotateX(angz);
}


void Polymere::updateFull(){
//Aligns each full atom backbone to the current virtual bond model

for(int m=0; m < num_mols; m++){
double ax,ay,az,bx,by,bz,cx,cy,cz;
ax = mols[m]->x;
ay = mols[m]->y;
az = mols[m]->z;

bx = mols[m]->bx;
by = mols[m]->by;
bz = mols[m]->bz;

cx = mols[m]->b2x;
cy = mols[m]->b2y;
cz = mols[m]->b2z;

//Create a forth point using the cross product of a-b X c-b + b
double A1=ax-bx;
double A2=ay-by;
double A3=az-bz;

double B1=cx-bx;
double B2=cy-by;
double B3=cz-bz;

double Bd = sqrt(pow(B1,2)+pow(B2,2)+pow(B3,2));
if(mols[m]->seq=='A' || mols[m]->seq=='G' ){
B1=B1/Bd * 1.3745;
B2=B2/Bd * 1.3745;
B3=B3/Bd * 1.3745;
mols[m]->b2x = mols[m]->bx + B1;
mols[m]->b2y = mols[m]->by + B2;
mols[m]->b2z = mols[m]->bz + B3;
}else{
if(mols[m]->seq=='C' || mols[m]->seq=='U'){
B1=B1/Bd * 1.3965;
B2=B2/Bd * 1.3965;
B3=B3/Bd * 1.3965;
//Adjust positions in pseudo atoms
mols[m]->b2x = mols[m]->bx + B1;
mols[m]->b2y = mols[m]->by + B2;
mols[m]->b2z = mols[m]->bz + B3;
}
}

double Z1=B2*A3-B3*A2;
double Z2=B3*A1-B1*A3;
double Z3=B1*A2-B2*A1;

//Forth point is cross product of A and B
double dx=Z1+bx;
double dy=Z2+by;
double dz=Z3+bz;

//Determine transformation Matrix 

Matrix * Q = new Matrix();
Q->a11=bx-ax;
Q->a21=cx-ax;
Q->a31=dx-ax;
Q->a12=by-ay;
Q->a22=cy-ay;
Q->a32=dy-ay;
Q->a13=bz-az;
Q->a23=cz-az;
Q->a33=dz-az;

double ax2,ay2,az2,bx2,by2,bz2,cx2,cy2,cz2;

int bpos=getResidueStart(m);
ax2=full_atoms[bpos+33];
ay2=full_atoms[bpos+34];
az2=full_atoms[bpos+35];

bx2=full_atoms[bpos+36];
by2=full_atoms[bpos+37];
bz2=full_atoms[bpos+38];

cx2=full_atoms[bpos+39];
cy2=full_atoms[bpos+40];
cz2=full_atoms[bpos+41];

//create the  point d
//Create a forth point using the cross product of a-b X c-b + b
 A1=ax2-bx2;
 A2=ay2-by2;
 A3=az2-bz2;

 B1=cx2-bx2;
 B2=cy2-by2;
 B3=cz2-bz2;

 Z1=B2*A3-B3*A2;
 Z2=B3*A1-B1*A3;
 Z3=B1*A2-B2*A1;

//Forth point is cross product of A and B
double dx2=Z1+bx2;
double dy2=Z2+by2;
double dz2=Z3+bz2;

 
Matrix * P = new Matrix();
P->a11=bx2-ax2;
P->a21=cx2-ax2;
P->a31=dx2-ax2;
P->a12=by2-ay2;
P->a22=cy2-ay2;
P->a32=dy2-ay2;
P->a13=bz2-az2;
P->a23=cz2-az2;
P->a33=dz2-az2;
Matrix Ip = P->inverse();
Matrix * T = new Matrix();

T->a11=Ip.a11*Q->a11+Ip.a12*Q->a21+Ip.a13*Q->a31;
T->a21=Ip.a21*Q->a11+Ip.a22*Q->a21+Ip.a23*Q->a31;
T->a31=Ip.a31*Q->a11+Ip.a32*Q->a21+Ip.a33*Q->a31;

T->a12=Ip.a11*Q->a12+Ip.a12*Q->a22+Ip.a13*Q->a32;
T->a22=Ip.a21*Q->a12+Ip.a22*Q->a22+Ip.a23*Q->a32;
T->a32=Ip.a31*Q->a12+Ip.a32*Q->a22+Ip.a33*Q->a32;

T->a13=Ip.a11*Q->a13+Ip.a12*Q->a23+Ip.a13*Q->a33;
T->a23=Ip.a21*Q->a13+Ip.a22*Q->a23+Ip.a23*Q->a33;
T->a33=Ip.a31*Q->a13+Ip.a32*Q->a23+Ip.a33*Q->a33;

double ax3= ax-(ax2*T->a11+ay2*T->a21+az2*T->a31);
double ay3= ay-(ax2*T->a12+ay2*T->a22+az2*T->a32);
double az3= az-(ax2*T->a13+ay2*T->a23+az2*T->a33);

double x,y,z,w1,w2,w3;

for(int f=bpos; f<(bpos+81); f=f+3){
x=full_atoms[f];y=full_atoms[f+1];z=full_atoms[f+2];

w1=x*T->a11+y*T->a21+z*T->a31+ax3;
w2=x*T->a12+y*T->a22+z*T->a32+ay3;
w3=x*T->a13+y*T->a23+z*T->a33+az3;

full_atoms[f]=w1;
full_atoms[f+1]=w2;
full_atoms[f+2]=w3;

}

delete P;
delete T; 
delete Q; 

}



}


void Polymere::alignStruct(int pos){
using namespace std;
//Translate first position to origin
double x1=mols[pos]->x;
double y1=mols[pos]->y;
double z1=mols[pos]->z;
translate(x1,y1,z1);
//rotate second position to x-axis
double r;
double angz;
double angy;
int p2 = pos+1;
x1=mols[p2]->x;
y1=mols[p2]->y;
z1=mols[p2]->z;
//find angles between pos2 and x-axis
cartesian_to_sperical(x1,y1,z1,&r,&angy,&angz);
//rotate structure so that x1,y1,z1 is on the x-axis 
rotateZ(angy);
//Need to set position 2 along the positive x-xis 
angz=PI/2-angz;
rotateY(angz);
//Rotate Position -1 into positive Y plane
p2 = pos-1;
y1=mols[p2]->y;
x1=mols[p2]->x;
z1=mols[p2]->z;
angz = atan2(y1,z1);
angz=angz-PI/2;
rotateX(angz);
}

void Polymere::changeFragVid(int pos,int a,std::vector<Fragment *> frags){
//Pos = First NT in the fragment
using namespace std;
double bond_0=frags[a]->bond_0;
double len_1=frags[a]->len_1;
double tor_1=frags[a]->tor_1;
double bond_1=frags[a]->bond_1;
double len_2=frags[a]->len_2;
double tor_2=frags[a]->tor_2;
double bond_2=frags[a]->bond_2;
double len_3=frags[a]->len_3;
double len_4=frags[a]->len_4;
//Update Bond angles
//cout<< nt1 << nt2 << nt3 << len_1<< na3_1 << tor_1 << len_2 << btr_2 << na1_2 << na2_2 << na3_2 << tor_2 << btr_3 << na1_3 << endl;
//if(pos>0){changeBond3a(pos,na3_1);}
int p0 = pos -1;
int p2 = pos + 1;
int p3 = pos + 2;
int p4 = pos + 3;
//Bond between C1-C2-C3
//write out structure
writeStruct("vid_structs",0);
changeBond2(p2,bond_1);
writeStruct("vid_structs",1);
changeBond2(pos,bond_0);
writeStruct("vid_structs",2);
changeBond2(p3,bond_2);
writeStruct("vid_structs",3);
changeTor2(pos,tor_1);
writeStruct("vid_structs",4);
changeTor2(p2,tor_2);
writeStruct("vid_structs",5);
//Update Bond lengths;
changeBondLen(p0,len_1);
writeStruct("vid_structs",6);
changeBondLen(pos,len_2);
writeStruct("vid_structs",7);
changeBondLen(p2,len_3);
writeStruct("vid_structs",8);
changeBondLen(p3,len_4);
writeStruct("vid_structs",9);
}

void Polymere::changeFrag1(int pos,int a,const std::vector<Fragment *> frags){
//Pos = First NT in the fragment
using namespace std;
double bond_0=frags[a]->bond_0;
double len_1=frags[a]->len_1;
double tor_1=frags[a]->tor_1;
double bond_1=frags[a]->bond_1;
double len_2=frags[a]->len_2;
double tor_2=frags[a]->tor_2;
double bond_2=frags[a]->bond_2;
double len_3=frags[a]->len_3;
double len_4=frags[a]->len_4;
//Update Bond angles
//cout<< nt1 << nt2 << nt3 << len_1<< na3_1 << tor_1 << len_2 << btr_2 << na1_2 << na2_2 << na3_2 << tor_2 << btr_3 << na1_3 << endl;
//if(pos>0){changeBond3a(pos,na3_1);}
int p0 = pos -1;
int p2 = pos + 1;
int p3 = pos + 2;
int p4 = pos + 3;
//Bond between C1-C2-C3
changeBond2r(p2,bond_1);
changeBond2r(pos,bond_0);
changeBond2r(p3,bond_2);
changeTor2r(p2,tor_2);
changeTor2r(pos,tor_1);
//Update Bond lengths;
changeBondLen2r(p0,len_1);
changeBondLen2r(pos,len_2);
changeBondLen2r(p2,len_3);
changeBondLen2r(p3,len_4);
}

void Polymere::changeFrag2(int pos,int a,const std::vector<Fragment *> frags){
//Pos = First NT in the fragment
using namespace std;
double bond_0=frags[a]->bond_0;
double len_1=frags[a]->len_1;
double tor_1=frags[a]->tor_1;
double bond_1=frags[a]->bond_1;
double len_2=frags[a]->len_2;
double tor_2=frags[a]->tor_2;
double bond_2=frags[a]->bond_2;
double len_3=frags[a]->len_3;
double len_4=frags[a]->len_4;
//Update Bond angles
//cout<< nt1 << nt2 << nt3 << len_1<< na3_1 << tor_1 << len_2 << btr_2 << na1_2 << na2_2 << na3_2 << tor_2 << btr_3 << na1_3 << endl;
//if(pos>0){changeBond3a(pos,na3_1);}
int p0 = pos -1;
int p2 = pos + 1;
int p3 = pos + 2;
int p4 = pos + 3;
//Bond between C1-C2-C3
changeBond2(p2,bond_1);
changeBond2(pos,bond_0);
changeBond2(p3,bond_2);
changeTor2(pos,tor_1);
changeTor2(p2,tor_2);
//Update Bond lengths;
changeBondLen(p0,len_1);
changeBondLen(pos,len_2);
changeBondLen(p2,len_3);
changeBondLen(p3,len_4);
}


void Polymere::changeFragR(int pos,int a,std::vector<Fragment *> frags){
//Pos = First NT in the fragment
using namespace std;
double bond_0=frags[a]->bond_0;
double len_1=frags[a]->len_1;
double tor_1=frags[a]->tor_1;
double bond_1=frags[a]->bond_1;
double len_2=frags[a]->len_2;
double tor_2=frags[a]->tor_2;
double bond_2=frags[a]->bond_2;
double len_3=frags[a]->len_3;
double len_4=frags[a]->len_4;
//Update Bond angles
//cout<< nt1 << nt2 << nt3 << len_1<< na3_1 << tor_1 << len_2 << btr_2 << na1_2 << na2_2 << na3_2 << tor_2 << btr_3 << na1_3 << endl;
//if(pos>0){changeBond3a(pos,na3_1);}
int p0 = pos -1;
int p2 = pos + 1;
int p3 = pos + 2;
int p4 = pos + 3;
//Bond between C1-C2-C3
changeBond2(p2,bond_1);
changeBond2(pos,bond_0);
changeBond2(p3,bond_2);
changeTorR(p2,tor_1);
changeTorR(p3,tor_2);
//Update Bond lengths;
changeBondLenR(pos,len_1);
changeBondLenR(p2,len_2);
changeBondLenR(p3,len_3);
changeBondLenR(p4,len_4);
}
float Polymere::c1distf(int i, int j){
//If any C1dist < 3.3 return 0
//Caclulate rho
double Dx,Dy,Dz,rho;
Dx=mols[i]->x - mols[j]->x;
Dy=mols[i]->y - mols[j]->y;
Dz=mols[i]->z - mols[j]->z;
rho = sqrt(pow(Dx,2)+pow(Dy,2)+pow(Dz,2));
return rho;
}
int Polymere::c1dist(){
//If any C1dist < 3.3 return 0
//Caclulate rho
int i=0;
int j=0;
int out=1;
double Dx,Dy,Dz,rho;
while(i < num_mols){
j=i+1;
while(j < num_mols){ 
Dx=mols[i]->x - mols[j]->x;
Dy=mols[i]->y - mols[j]->y;
Dz=mols[i]->z - mols[j]->z;
rho = sqrt(pow(Dx,2)+pow(Dy,2)+pow(Dz,2));
if(rho < 3.3){
//cout << "rho: " << i << ":" << j << ":" << rho << endl; 
i=num_mols;
j=num_mols;
out=0;

}
j++;
}
i++;
}
return out;
}


std::vector<double> Polymere::planar(int i, int j){
//Given the positions of two nucleotides it calculates the rho,theta,phi.. coordinates
double Ax1=mols[i]->x - mols[i]->bx;
double Ay1=mols[i]->y - mols[i]->by;
double Az1=mols[i]->z - mols[i]->bz;

double Bx1=mols[i]->b2x - mols[i]->bx;
double By1=mols[i]->b2y - mols[i]->by;
double Bz1=mols[i]->b2z - mols[i]->bz;

double Ax2=mols[j]->x - mols[j]->bx;
double Ay2=mols[j]->y - mols[j]->by;
double Az2=mols[j]->z - mols[j]->bz;

double Bx2=mols[j]->b2x - mols[j]->bx;
double By2=mols[j]->b2y - mols[j]->by;
double Bz2=mols[j]->b2z - mols[j]->bz;

double Cx1=Ay1*Bz1 - Az1*By1;
double Cy1=Az1*Bx1 - Ax1*Bz1;
double Cz1=Ax1*By1 - Ay1*Bx1;

//Making orthonormal set
double Ex1=Cy1*Az1 - Cz1*Ay1;
double Ey1=Cz1*Ax1 - Cx1*Az1;
double Ez1=Cx1*Ay1 - Cy1*Ax1;

double Cx2 =Ay2*Bz2 - Az2*By2;
double Cy2 =Az2*Bx2 - Ax2*Bz2;
double Cz2 =Ax2*By2 - Ay2*Bx2;

double Ex2=Cy2*Az2 - Cz2*Ay2;
double Ey2=Cz2*Ax2 - Cx2*Az2;
double Ez2=Cx2*Ay2 - Cy2*Ax2;

//A,E,C is a coordinate frame
//Determine theta1,phi1,rho,theta2,phi2 
//plug into energy function
//Caclulate angle between normal vectors
double CC1=sqrt(pow(Cx1,2)+pow(Cy1,2)+pow(Cz1,2));
double CC2=sqrt(pow(Cx2,2)+pow(Cy2,2)+pow(Cz2,2));
double AA1=sqrt(pow(Ax1,2)+pow(Ay1,2)+pow(Az1,2));
double AA2=sqrt(pow(Ax2,2)+pow(Ay2,2)+pow(Az2,2));
double EE1=sqrt(pow(Ex1,2)+pow(Ey1,2)+pow(Ez1,2));
double EE2=sqrt(pow(Ex2,2)+pow(Ey2,2)+pow(Ez2,2));

//Create orthonormal set
Cx1=Cx1/CC1;
Cy1=Cy1/CC1;
Cz1=Cz1/CC1;

Ax1=Ax1/AA1;
Ay1=Ay1/AA1;
Az1=Az1/AA1;

Ex1=Ex1/EE1;
Ey1=Ey1/EE1;
Ez1=Ez1/EE1;

Cx2=Cx2/CC2;
Cy2=Cy2/CC2;
Cz2=Cz2/CC2;

Ax2=Ax2/AA2;
Ay2=Ay2/AA2;
Az2=Az2/AA2;
       
Ex2=Ex2/EE2;
Ey2=Ey2/EE2;
Ez2=Ez2/EE2;

//Caclulate rho 
double Dx=mols[j]->bx - mols[i]->bx;
double Dy=mols[j]->by - mols[i]->by;
double Dz=mols[j]->bz - mols[i]->bz;

double rho = sqrt(pow(Dx,2)+pow(Dy,2)+pow(Dz,2));
//Transform D into A,B,C coordinate system
double dx_i=(Dx*Ax1+Dy*Ay1+Dz*Az1);
double dy_i=(Dx*Ex1+Dy*Ey1+Dz*Ez1);
double dz_i=(Dx*Cx1+Dy*Cy1+Dz*Cz1);

//Calculate theta
double xout=dx_i;
double yout=dy_i;
double zout=dz_i;
std::vector<double> out;
out.push_back(xout);
out.push_back(yout);
out.push_back(zout);
return out;
}

int Polymere::arePlanar(int i,int j){
//returns 1 is the two are in the same plane and 0 is they are not
int areP=0;
std::vector<double> params1;
std::vector<double> params2;
for(int w=0; w<11; w++){params1.push_back(0); params2.push_back(0);}
this->getParams(i,j,&params1);
this->getParams(j,i,&params2);
if(params1[8] < 2 && params1[8] > -2 && params2[8] < 2 && params2[8] > -2){
areP=1;
} 
return areP;
}

double Polymere::coplanar(int i, int j){
//Given the positions of two nucleotides it calculates the rho,theta,phi.. coordinates
double Ax1=mols[i]->x - mols[i]->bx;
double Ay1=mols[i]->y - mols[i]->by;
double Az1=mols[i]->z - mols[i]->bz;

double Bx1=mols[i]->b2x - mols[i]->bx;
double By1=mols[i]->b2y - mols[i]->by;
double Bz1=mols[i]->b2z - mols[i]->bz;

double Ax2=mols[j]->x - mols[j]->bx;
double Ay2=mols[j]->y - mols[j]->by;
double Az2=mols[j]->z - mols[j]->bz;

double Bx2=mols[j]->b2x - mols[j]->bx;
double By2=mols[j]->b2y - mols[j]->by;
double Bz2=mols[j]->b2z - mols[j]->bz;

double Cx1=Ay1*Bz1 - Az1*By1;
double Cy1=Az1*Bx1 - Ax1*Bz1;
double Cz1=Ax1*By1 - Ay1*Bx1;

//Making orthonormal set
double Ex1=Cy1*Az1 - Cz1*Ay1;
double Ey1=Cz1*Ax1 - Cx1*Az1;
double Ez1=Cx1*Ay1 - Cy1*Ax1;

double Cx2 =Ay2*Bz2 - Az2*By2;
double Cy2 =Az2*Bx2 - Ax2*Bz2;
double Cz2 =Ax2*By2 - Ay2*Bx2;

double Ex2=Cy2*Az2 - Cz2*Ay2;
double Ey2=Cz2*Ax2 - Cx2*Az2;
double Ez2=Cx2*Ay2 - Cy2*Ax2;

//A,E,C is a coordinate frame
//Determine theta1,phi1,rho,theta2,phi2 
//plug into energy function
//Caclulate angle between normal vectors
double CC1=sqrt(pow(Cx1,2)+pow(Cy1,2)+pow(Cz1,2));
double CC2=sqrt(pow(Cx2,2)+pow(Cy2,2)+pow(Cz2,2));
double AA1=sqrt(pow(Ax1,2)+pow(Ay1,2)+pow(Az1,2));
double AA2=sqrt(pow(Ax2,2)+pow(Ay2,2)+pow(Az2,2));
double EE1=sqrt(pow(Ex1,2)+pow(Ey1,2)+pow(Ez1,2));
double EE2=sqrt(pow(Ex2,2)+pow(Ey2,2)+pow(Ez2,2));

//Create orthonormal set
Cx1=Cx1/CC1;
Cy1=Cy1/CC1;
Cz1=Cz1/CC1;

Ax1=Ax1/AA1;
Ay1=Ay1/AA1;
Az1=Az1/AA1;

Ex1=Ex1/EE1;
Ey1=Ey1/EE1;
Ez1=Ez1/EE1;

Cx2=Cx2/CC2;
Cy2=Cy2/CC2;
Cz2=Cz2/CC2;

Ax2=Ax2/AA2;
Ay2=Ay2/AA2;
Az2=Az2/AA2;
       
Ex2=Ex2/EE2;
Ey2=Ey2/EE2;
Ez2=Ez2/EE2;

//Caclulate rho 
double Dx=mols[i]->bx - mols[j]->bx;
double Dy=mols[i]->by - mols[j]->by;
double Dz=mols[i]->bz - mols[j]->bz;

double rho = sqrt(pow(Dx,2)+pow(Dy,2)+pow(Dz,2));
//Transform D into A,B,C coordinate system
double dx_i=(Dx*Ax1+Dy*Ay1+Dz*Az1);
double dy_i=(Dx*Ex1+Dy*Ey1+Dz*Ez1);
double dz_i=(Dx*Cx1+Dy*Cy1+Dz*Cz1);

double dx_j=-1*(Dx*Ax2+Dy*Ay2+Dz*Az2);
double dy_j=-1*(Dx*Ex2+Dy*Ey2+Dz*Ez2);
double dz_j=-1*(Dx*Cx2+Dy*Cy2+Dz*Cz2);

//Calculate theta
double zout=abs(dz_i);
return zout;
}

void Polymere::getParams(int i, int j,vector<double> * params){
//Given the positions of two nucleotides it calculates the rho,theta,phi.. coordinates
double Ax1=mols[i]->x - mols[i]->bx;
double Ay1=mols[i]->y - mols[i]->by;
double Az1=mols[i]->z - mols[i]->bz;

double Bx1=mols[i]->b2x - mols[i]->bx;
double By1=mols[i]->b2y - mols[i]->by;
double Bz1=mols[i]->b2z - mols[i]->bz;

double Ax2=mols[j]->x - mols[j]->bx;
double Ay2=mols[j]->y - mols[j]->by;
double Az2=mols[j]->z - mols[j]->bz;

double Bx2=mols[j]->b2x - mols[j]->bx;
double By2=mols[j]->b2y - mols[j]->by;
double Bz2=mols[j]->b2z - mols[j]->bz;

double Cx1=Ay1*Bz1 - Az1*By1;
double Cy1=Az1*Bx1 - Ax1*Bz1;
double Cz1=Ax1*By1 - Ay1*Bx1;

//Making orthonormal set
double Ex1=Cy1*Az1 - Cz1*Ay1;
double Ey1=Cz1*Ax1 - Cx1*Az1;
double Ez1=Cx1*Ay1 - Cy1*Ax1;

double Cx2 =Ay2*Bz2 - Az2*By2;
double Cy2 =Az2*Bx2 - Ax2*Bz2;
double Cz2 =Ax2*By2 - Ay2*Bx2;

double Ex2=Cy2*Az2 - Cz2*Ay2;
double Ey2=Cz2*Ax2 - Cx2*Az2;
double Ez2=Cx2*Ay2 - Cy2*Ax2;

//A,E,C is a coordinate frame
//Determine theta1,phi1,rho,theta2,phi2 
//plug into energy function
//Caclulate angle between normal vectors
double CC1=sqrt(pow(Cx1,2)+pow(Cy1,2)+pow(Cz1,2));
double CC2=sqrt(pow(Cx2,2)+pow(Cy2,2)+pow(Cz2,2));
double AA1=sqrt(pow(Ax1,2)+pow(Ay1,2)+pow(Az1,2));
double AA2=sqrt(pow(Ax2,2)+pow(Ay2,2)+pow(Az2,2));
double EE1=sqrt(pow(Ex1,2)+pow(Ey1,2)+pow(Ez1,2));
double EE2=sqrt(pow(Ex2,2)+pow(Ey2,2)+pow(Ez2,2));

//Create orthonormal set
Cx1=Cx1/CC1;
Cy1=Cy1/CC1;
Cz1=Cz1/CC1;

Ax1=Ax1/AA1;
Ay1=Ay1/AA1;
Az1=Az1/AA1;

Ex1=Ex1/EE1;
Ey1=Ey1/EE1;
Ez1=Ez1/EE1;

Cx2=Cx2/CC2;
Cy2=Cy2/CC2;
Cz2=Cz2/CC2;

Ax2=Ax2/AA2;
Ay2=Ay2/AA2;
Az2=Az2/AA2;
       
Ex2=Ex2/EE2;
Ey2=Ey2/EE2;
Ez2=Ez2/EE2;

//Calculate Angle between Z-axis
double zzang=acos((Cx1*Cx2)+(Cy1*Cy2)+(Cz1*Cz2));

//Caclulate rho 
double Dx=mols[i]->bx - mols[j]->bx;
double Dy=mols[i]->by - mols[j]->by;
double Dz=mols[i]->bz - mols[j]->bz;

double rho = sqrt(pow(Dx,2)+pow(Dy,2)+pow(Dz,2));
//Transform D into A,B,C coordinate system
double dx_i=(Dx*Ax1+Dy*Ay1+Dz*Az1);
double dy_i=(Dx*Ex1+Dy*Ey1+Dz*Ez1);
double dz_i=(Dx*Cx1+Dy*Cy1+Dz*Cz1);

double dx_j=-1*(Dx*Ax2+Dy*Ay2+Dz*Az2);
double dy_j=-1*(Dx*Ex2+Dy*Ey2+Dz*Ez2);
double dz_j=-1*(Dx*Cx2+Dy*Cy2+Dz*Cz2);

//Calculate theta
double theta1=acos(dz_i/rho);
//Calculate phi
double s = sqrt(pow(dx_i,2)+pow(dy_i,2));
double phi1 = atan2(dy_i,dx_i);

//Calculate theta
double theta2=acos(dz_j/rho);
//Calculate phi
double s2 = sqrt(pow(dx_j,2)+pow(dy_j,2));
double phi2 = atan2(dy_j,dx_j);
int type=0;
if(mols[i]->seq=='A' && mols[j]->seq=='U'){type=0;}
if(mols[i]->seq=='A' && mols[j]->seq=='C'){type=1;}
if(mols[i]->seq=='A' && mols[j]->seq=='G'){type=2;}
if(mols[i]->seq=='A' && mols[j]->seq=='A'){type=3;}
if(mols[i]->seq=='U' && mols[j]->seq=='C'){type=4;}
if(mols[i]->seq=='U' && mols[j]->seq=='G'){type=5;}
if(mols[i]->seq=='U' && mols[j]->seq=='U'){type=6;}
if(mols[i]->seq=='U' && mols[j]->seq=='A'){type=7;}
if(mols[i]->seq=='C' && mols[j]->seq=='U'){type=8;}
if(mols[i]->seq=='C' && mols[j]->seq=='A'){type=9;}
if(mols[i]->seq=='C' && mols[j]->seq=='G'){type=10;}
if(mols[i]->seq=='C' && mols[j]->seq=='C'){type=11;}
if(mols[i]->seq=='G' && mols[j]->seq=='G'){type=12;}
if(mols[i]->seq=='G' && mols[j]->seq=='A'){type=13;}
if(mols[i]->seq=='G' && mols[j]->seq=='C'){type=14;}
if(mols[i]->seq=='G' && mols[j]->seq=='U'){type=15;}
double typeRY;
//RY
if(type==0 || type==1 || type==14 || type==15){typeRY=0;} 
//RR
if(type==2 || type==3 || type==12 || type==13){typeRY=1;} 
//YY
if(type==4 || type==6 || type==8 || type==11){typeRY=2;} 
//YR
if(type==5 || type==7 || type==9|| type==10){typeRY=3;} 
params->at(0)=type;
params->at(1)=rho;
params->at(2)=theta1;
params->at(3)=phi1;
params->at(4)=phi2;
params->at(5)=theta2;
params->at(6)=dx_i;
params->at(7)=dy_i;
params->at(8)=dz_i;
params->at(9)=typeRY;
//Z-angle was added as an after thought and might not be in parts of the program
if(params->size() > 10){params->at(10)=zzang;}
//cout << dx_i << "\t" << dy_i << "\t" << dz_i << "\t" << theta2 << "\t" << phi2 << "\t" << typeRY << endl;
}


void Polymere::baseCoords(int i, int j){
//Given the positions of two nucleotides it calculates the rho,theta,phi.. coordinates
double Ax1=mols[i]->x - mols[i]->bx;
double Ay1=mols[i]->y - mols[i]->by;
double Az1=mols[i]->z - mols[i]->bz;

double Bx1=mols[i]->b2x - mols[i]->bx;
double By1=mols[i]->b2y - mols[i]->by;
double Bz1=mols[i]->b2z - mols[i]->bz;

double Ax2=mols[j]->x - mols[j]->bx;
double Ay2=mols[j]->y - mols[j]->by;
double Az2=mols[j]->z - mols[j]->bz;

double Bx2=mols[j]->b2x - mols[j]->bx;
double By2=mols[j]->b2y - mols[j]->by;
double Bz2=mols[j]->b2z - mols[j]->bz;

double Cx1=Ay1*Bz1 - Az1*By1;
double Cy1=Az1*Bx1 - Ax1*Bz1;
double Cz1=Ax1*By1 - Ay1*Bx1;

//Making orthonormal set
double Ex1=Cy1*Az1 - Cz1*Ay1;
double Ey1=Cz1*Ax1 - Cx1*Az1;
double Ez1=Cx1*Ay1 - Cy1*Ax1;

double Cx2 =Ay2*Bz2 - Az2*By2;
double Cy2 =Az2*Bx2 - Ax2*Bz2;
double Cz2 =Ax2*By2 - Ay2*Bx2;

double Ex2=Cy2*Az2 - Cz2*Ay2;
double Ey2=Cz2*Ax2 - Cx2*Az2;
double Ez2=Cx2*Ay2 - Cy2*Ax2;

//A,E,C is a coordinate frame
//Determine theta1,phi1,rho,theta2,phi2 
//plug into energy function
//Caclulate angle between normal vectors
double CC1=sqrt(pow(Cx1,2)+pow(Cy1,2)+pow(Cz1,2));
double CC2=sqrt(pow(Cx2,2)+pow(Cy2,2)+pow(Cz2,2));
double AA1=sqrt(pow(Ax1,2)+pow(Ay1,2)+pow(Az1,2));
double AA2=sqrt(pow(Ax2,2)+pow(Ay2,2)+pow(Az2,2));
double EE1=sqrt(pow(Ex1,2)+pow(Ey1,2)+pow(Ez1,2));
double EE2=sqrt(pow(Ex2,2)+pow(Ey2,2)+pow(Ez2,2));

//Create orthonormal set
Cx1=Cx1/CC1;
Cy1=Cy1/CC1;
Cz1=Cz1/CC1;

Ax1=Ax1/AA1;
Ay1=Ay1/AA1;
Az1=Az1/AA1;

Ex1=Ex1/EE1;
Ey1=Ey1/EE1;
Ez1=Ez1/EE1;

Cx2=Cx2/CC2;
Cy2=Cy2/CC2;
Cz2=Cz2/CC2;

Ax2=Ax2/AA2;
Ay2=Ay2/AA2;
Az2=Az2/AA2;
       
Ex2=Ex2/EE2;
Ey2=Ey2/EE2;
Ez2=Ez2/EE2;

//Caclulate rho 
double Dx=mols[i]->bx - mols[j]->bx;
double Dy=mols[i]->by - mols[j]->by;
double Dz=mols[i]->bz - mols[j]->bz;

double rho = sqrt(pow(Dx,2)+pow(Dy,2)+pow(Dz,2));
//Transform D into A,B,C coordinate system
double dx_i=(Dx*Ax1+Dy*Ay1+Dz*Az1);
double dy_i=(Dx*Ex1+Dy*Ey1+Dz*Ez1);
double dz_i=(Dx*Cx1+Dy*Cy1+Dz*Cz1);

double dx_j=-1*(Dx*Ax2+Dy*Ay2+Dz*Az2);
double dy_j=-1*(Dx*Ex2+Dy*Ey2+Dz*Ez2);
double dz_j=-1*(Dx*Cx2+Dy*Cy2+Dz*Cz2);

//Calculate theta
double theta1=acos(dz_i/rho);
//Calculate phi
double s = sqrt(pow(dx_i,2)+pow(dy_i,2));
double phi1 = atan2(dy_i,dx_i);

//Calculate theta
double theta2=acos(dz_j/rho);
//Calculate phi
double s2 = sqrt(pow(dx_j,2)+pow(dy_j,2));
double phi2 = atan2(dy_j,dx_j);
//cout << mols[i]->seq << "\t" << mols[j]->seq << "\t" << rho << "\t" << theta1 << "\t" << phi1 << "\t" << theta2 << "\t" << phi2 << endl; 
}
void Polymere::hydroPairs(std::vector<float> * hydroVals,std::vector<int> * hpair){
for(int i=0; i<num_mols; i++){
hydroVals->push_back(0);
hpair->push_back(0);
float mallprob=0;
int cur_pair=0;
for(int j=0; j<num_mols; j++){
float allprob=0;
if(abs(i-j) > 1){
//Only Consider Pairs that are in the movement groups
//cout << "Molsi: " << i << ":" << mols[i]->group << " Molsj: " << j << ":" << mols[j]->group << endl; 
double mdist=100;
mdist =pow((mols[i]->x - mols[j]->x),2);
mdist =mdist + pow((mols[i]->y -mols[j]->y),2);
mdist =mdist + pow((mols[i]->z -mols[j]->z),2);
mdist = sqrt(mdist);
//cout << "Type: " << params[0] << ":" << i << "," << j << " i: " << mols[i]->seq << " j: " <<  mols[j]->seq << " mdist: " << mdist << endl;
//cout << "mdist: " << mdist << endl;
double zplane=coplanar(i,j);
if(zplane < 1.5){
if(mdist< 20){
vector <Points> p;
vector <Points> p2;
Points c1,b1,b2;
int type1,type2;
//Returns a score for hydrogen bonding
c1.x=mols[i]->x;
c1.y=mols[i]->y;
c1.z=mols[i]->z;
b1.x=mols[i]->bx;
b1.y=mols[i]->by;
b1.z=mols[i]->bz;
b2.x=mols[i]->b2x;
b2.y=mols[i]->b2y;
b2.z=mols[i]->b2z;
p.push_back(c1);
p.push_back(b1);
p.push_back(b2);

c1.x=mols[j]->x;
c1.y=mols[j]->y;
c1.z=mols[j]->z;
b1.x=mols[j]->bx;
b1.y=mols[j]->by;
b1.z=mols[j]->bz;
b2.x=mols[j]->b2x;
b2.y=mols[j]->b2y;
b2.z=mols[j]->b2z;
p2.push_back(c1);
p2.push_back(b1);
p2.push_back(b2);

if(mols[i]->seq=='G'){type1=1;}else{
if(mols[i]->seq=='A'){type1=2;}else{
if(mols[i]->seq=='U'){type1=3;}else{
if(mols[i]->seq=='C'){type1=4;}else{
}
}
}
}
if(mols[j]->seq=='G'){type2=1;}else{
if(mols[j]->seq=='A'){type2=2;}else{
if(mols[j]->seq=='U'){type2=3;}else{
if(mols[j]->seq=='C'){type2=4;}else{
}
}
}
}
FABase * base1 = new FABase(type1);  
FABase * base2 = new FABase(type2);  
//Move to coordinate system
base1->alignBase(p);
base2->alignBase(p2);
//base1->print();
//base2->print();
//cout << "Types: " << type1 << ":" << type2 << endl; 
//Identify all donors-acceptor pairs

vector<Points> donors1;base1->donors(&donors1);
vector<Points> accept1;base1->acceptors(&accept1);
vector<Points> donors2;base2->donors(&donors2);
vector<Points> accept2;base2->acceptors(&accept2);

//Caclulate distance and two angles
//Acceptors have hydrogen
//donors have LP
//cout << "acceptors i -> donors j" << endl;
//cout << "donors i -> acceptors j" << endl;
for(int d=0; d < donors1.size(); d=d+2){
for(int a=0; a < accept2.size(); a=a+2){
//angle between H and D
double dx = donors1[d].x  -accept2[a].x;
double dy = donors1[d].y  -accept2[a].y;
double dz = donors1[d].z  -accept2[a].z;
double dd = sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));

double px = accept2[a+1].x-donors1[d+1].x;
double py = accept2[a+1].y-donors1[d+1].y;
double pz = accept2[a+1].z-donors1[d+1].z;
double pp =  sqrt(pow(px,2)+pow(py,2)+pow(pz,2));

double hx = accept2[a+1].x-accept2[a].x;
double hy = accept2[a+1].y-accept2[a].y;
double hz = accept2[a+1].z-accept2[a].z;
double hh =  sqrt(pow(hx,2)+pow(hy,2)+pow(hz,2));

double lx = donors1[d+1].x-donors1[d].x;
double ly = donors1[d+1].y-donors1[d].y;
double lz = donors1[d+1].z-donors1[d].z;
double ll =  sqrt(pow(lx,2)+pow(ly,2)+pow(lz,2));

//angle between H and D
float ang1=(float)atanh(cos(acos((dx*hx+dy*hy+dz*hz)/(dd*hh))));
float ang2=(float)atanh(cos(acos(-1*(dx*lx+dy*ly+dz*lz)/(dd*ll))));
float dd1=(float)log(pow(pp,3));
float prob=base1->Eval(dd1,ang1,ang2,1);
allprob=allprob+prob;
//cout << dd1 << "," << ang1 << "," << ang2 << "," << prob << endl;
}
}
delete base1;
delete base2;
}
}
}
if(allprob > mallprob){
mallprob=allprob;
if(mols[j]->seq=='G'){cur_pair=1;}
if(mols[j]->seq=='A'){cur_pair=2;}
if(mols[j]->seq=='U'){cur_pair=3;}
if(mols[j]->seq=='C'){cur_pair=4;}
}

}
hydroVals->at(i)=mallprob;
hpair->at(i)=cur_pair;
printf("cur_pair: %d\n",cur_pair);
cur_pair=0;

}
}

void Polymere::hydroVg(vector<float> * hydroVals,int movegroup){
int hpos25=0;
hydroVals->clear();
for(int i=0; i<(num_mols); i++){
if(mols[i]->group==movegroup){
hydroVals->push_back(0);
float allprob=0;
for(int j=0; j<(num_mols); j++){
float mallprob=0;
if(abs(i-j) > 1 && mols[i]->group==movegroup && mols[j]->group==movegroup){
//Only Consider Pairs that are in the movement groups
//cout << "Molsi: " << i << ":" << mols[i]->group << " Molsj: " << j << ":" << mols[j]->group << endl; 
double mdist=100;
mdist =pow((mols[i]->x - mols[j]->x),2);
mdist =mdist + pow((mols[i]->y -mols[j]->y),2);
mdist =mdist + pow((mols[i]->z -mols[j]->z),2);
mdist = sqrt(mdist);
//cout << "Type: " << params[0] << ":" << i << "," << j << " i: " << mols[i]->seq << " j: " <<  mols[j]->seq << " mdist: " << mdist << endl;
//cout << "mdist: " << mdist << endl;
double zplane=coplanar(i,j);
if(zplane < 1.5){
if(mdist< 20){
vector <Points> p;
vector <Points> p2;
Points c1,b1,b2;
int type1,type2;
//Returns a score for hydrogen bonding
c1.x=mols[i]->x;
c1.y=mols[i]->y;
c1.z=mols[i]->z;
b1.x=mols[i]->bx;
b1.y=mols[i]->by;
b1.z=mols[i]->bz;
b2.x=mols[i]->b2x;
b2.y=mols[i]->b2y;
b2.z=mols[i]->b2z;
p.push_back(c1);
p.push_back(b1);
p.push_back(b2);

c1.x=mols[j]->x;
c1.y=mols[j]->y;
c1.z=mols[j]->z;
b1.x=mols[j]->bx;
b1.y=mols[j]->by;
b1.z=mols[j]->bz;
b2.x=mols[j]->b2x;
b2.y=mols[j]->b2y;
b2.z=mols[j]->b2z;
p2.push_back(c1);
p2.push_back(b1);
p2.push_back(b2);

if(mols[i]->seq=='G'){type1=1;}else{
if(mols[i]->seq=='A'){type1=2;}else{
if(mols[i]->seq=='U'){type1=3;}else{
if(mols[i]->seq=='C'){type1=4;}else{
}
}
}
}
if(mols[j]->seq=='G'){type2=1;}else{
if(mols[j]->seq=='A'){type2=2;}else{
if(mols[j]->seq=='U'){type2=3;}else{
if(mols[j]->seq=='C'){type2=4;}else{
}
}
}
}
FABase * base1 = new FABase(type1);  
FABase * base2 = new FABase(type2);  
//Move to coordinate system
base1->alignBase(p);
base2->alignBase(p2);
//base1->print();
//base2->print();
//cout << "Types: " << type1 << ":" << type2 << endl; 
//Identify all donors-acceptor pairs

vector<Points> donors1;base1->donors(&donors1);
vector<Points> accept1;base1->acceptors(&accept1);
vector<Points> donors2;base2->donors(&donors2);
vector<Points> accept2;base2->acceptors(&accept2);

//Caclulate distance and two angles
//Acceptors have hydrogen
//donors have LP
//cout << "acceptors i -> donors j" << endl;
//cout << "donors i -> acceptors j" << endl;
for(int d=0; d < donors1.size(); d=d+2){
for(int a=0; a < accept2.size(); a=a+2){
//angle between H and D
double dx = donors1[d].x  -accept2[a].x;
double dy = donors1[d].y  -accept2[a].y;
double dz = donors1[d].z  -accept2[a].z;
double dd = sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));

double px = accept2[a+1].x-donors1[d+1].x;
double py = accept2[a+1].y-donors1[d+1].y;
double pz = accept2[a+1].z-donors1[d+1].z;
double pp =  sqrt(pow(px,2)+pow(py,2)+pow(pz,2));

double hx = accept2[a+1].x-accept2[a].x;
double hy = accept2[a+1].y-accept2[a].y;
double hz = accept2[a+1].z-accept2[a].z;
double hh =  sqrt(pow(hx,2)+pow(hy,2)+pow(hz,2));

double lx = donors1[d+1].x-donors1[d].x;
double ly = donors1[d+1].y-donors1[d].y;
double lz = donors1[d+1].z-donors1[d].z;
double ll =  sqrt(pow(lx,2)+pow(ly,2)+pow(lz,2));

//angle between H and D
float ang1=(float)atanh(cos(acos((dx*hx+dy*hy+dz*hz)/(dd*hh))));
float ang2=(float)atanh(cos(acos(-1*(dx*lx+dy*ly+dz*lz)/(dd*ll))));
float dd1=(float)log(pow(pp,3));
float prob=base1->Eval(dd1,ang1,ang2,1);
mallprob=mallprob+prob;
//cout << dd1 << "," << ang1 << "," << ang2 << "," << prob << endl;
}
}
if(mallprob > allprob){allprob=mallprob;}
delete base1;
delete base2;
}
}
}
}
hydroVals->at(hpos25++)=allprob;
allprob=0;
}

}
}
void Polymere::hydroV(vector<float> * hydroVals){
for(int i=0; i<(num_mols); i++){
hydroVals->push_back(0);
float allprob=0;
for(int j=0; j<(num_mols); j++){
if(abs(i-j) > 1){
//Only Consider Pairs that are in the movement groups
//cout << "Molsi: " << i << ":" << mols[i]->group << " Molsj: " << j << ":" << mols[j]->group << endl; 
double mdist=100;
mdist =pow((mols[i]->x - mols[j]->x),2);
mdist =mdist + pow((mols[i]->y -mols[j]->y),2);
mdist =mdist + pow((mols[i]->z -mols[j]->z),2);
mdist = sqrt(mdist);
//cout << "Type: " << params[0] << ":" << i << "," << j << " i: " << mols[i]->seq << " j: " <<  mols[j]->seq << " mdist: " << mdist << endl;
//cout << "mdist: " << mdist << endl;
double zplane=coplanar(i,j);
if(zplane < 1.5){
if(mdist< 20){
vector <Points> p;
vector <Points> p2;
Points c1,b1,b2;
int type1,type2;
//Returns a score for hydrogen bonding
c1.x=mols[i]->x;
c1.y=mols[i]->y;
c1.z=mols[i]->z;
b1.x=mols[i]->bx;
b1.y=mols[i]->by;
b1.z=mols[i]->bz;
b2.x=mols[i]->b2x;
b2.y=mols[i]->b2y;
b2.z=mols[i]->b2z;
p.push_back(c1);
p.push_back(b1);
p.push_back(b2);

c1.x=mols[j]->x;
c1.y=mols[j]->y;
c1.z=mols[j]->z;
b1.x=mols[j]->bx;
b1.y=mols[j]->by;
b1.z=mols[j]->bz;
b2.x=mols[j]->b2x;
b2.y=mols[j]->b2y;
b2.z=mols[j]->b2z;
p2.push_back(c1);
p2.push_back(b1);
p2.push_back(b2);

if(mols[i]->seq=='G'){type1=1;}else{
if(mols[i]->seq=='A'){type1=2;}else{
if(mols[i]->seq=='U'){type1=3;}else{
if(mols[i]->seq=='C'){type1=4;}else{
}
}
}
}
if(mols[j]->seq=='G'){type2=1;}else{
if(mols[j]->seq=='A'){type2=2;}else{
if(mols[j]->seq=='U'){type2=3;}else{
if(mols[j]->seq=='C'){type2=4;}else{
}
}
}
}
FABase * base1 = new FABase(type1);  
FABase * base2 = new FABase(type2);  
//Move to coordinate system
base1->alignBase(p);
base2->alignBase(p2);
//base1->print();
//base2->print();
//cout << "Types: " << type1 << ":" << type2 << endl; 
//Identify all donors-acceptor pairs

vector<Points> donors1;base1->donors(&donors1);
vector<Points> accept1;base1->acceptors(&accept1);
vector<Points> donors2;base2->donors(&donors2);
vector<Points> accept2;base2->acceptors(&accept2);

//Caclulate distance and two angles
//Acceptors have hydrogen
//donors have LP
//cout << "acceptors i -> donors j" << endl;
//cout << "donors i -> acceptors j" << endl;
for(int d=0; d < donors1.size(); d=d+2){
for(int a=0; a < accept2.size(); a=a+2){
//angle between H and D
double dx = donors1[d].x  -accept2[a].x;
double dy = donors1[d].y  -accept2[a].y;
double dz = donors1[d].z  -accept2[a].z;
double dd = sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));

double px = accept2[a+1].x-donors1[d+1].x;
double py = accept2[a+1].y-donors1[d+1].y;
double pz = accept2[a+1].z-donors1[d+1].z;
double pp =  sqrt(pow(px,2)+pow(py,2)+pow(pz,2));

double hx = accept2[a+1].x-accept2[a].x;
double hy = accept2[a+1].y-accept2[a].y;
double hz = accept2[a+1].z-accept2[a].z;
double hh =  sqrt(pow(hx,2)+pow(hy,2)+pow(hz,2));

double lx = donors1[d+1].x-donors1[d].x;
double ly = donors1[d+1].y-donors1[d].y;
double lz = donors1[d+1].z-donors1[d].z;
double ll =  sqrt(pow(lx,2)+pow(ly,2)+pow(lz,2));

//angle between H and D
float ang1=(float)atanh(cos(acos((dx*hx+dy*hy+dz*hz)/(dd*hh))));
float ang2=(float)atanh(cos(acos(-1*(dx*lx+dy*ly+dz*lz)/(dd*ll))));
float dd1=(float)log(pow(pp,3));
float prob=base1->Eval(dd1,ang1,ang2,1);
allprob=allprob+prob;
//cout << dd1 << "," << ang1 << "," << ang2 << "," << prob << endl;
}
}
delete base1;
delete base2;
}
}
}
}
hydroVals->at(i)=allprob;
allprob=0;

}
}



float Polymere::hydroA(){
float allProb=0;
for(int i=0; i<num_mols; i++){
for(int j=(i+3); j<num_mols; j++){
//Only Consider Pairs that are in the movement groups
//cout << "Molsi: " << i << ":" << mols[i]->group << " Molsj: " << j << ":" << mols[j]->group << endl; 
double mdist=100;
mdist =pow((mols[i]->x - mols[j]->x),2);
mdist =mdist + pow((mols[i]->y -mols[j]->y),2);
mdist =mdist + pow((mols[i]->z -mols[j]->z),2);
mdist = sqrt(mdist);
//cout << "Type: " << params[0] << ":" << i << "," << j << " i: " << mols[i]->seq << " j: " <<  mols[j]->seq << " mdist: " << mdist << endl;
//cout << "mdist: " << mdist << endl;
double zplane=coplanar(i,j);
if(zplane < 1.5){
if(mdist< 20){
vector <Points> p;
vector <Points> p2;
Points c1,b1,b2;
int type1,type2;
//Returns a score for hydrogen bonding
c1.x=mols[i]->x;
c1.y=mols[i]->y;
c1.z=mols[i]->z;
b1.x=mols[i]->bx;
b1.y=mols[i]->by;
b1.z=mols[i]->bz;
b2.x=mols[i]->b2x;
b2.y=mols[i]->b2y;
b2.z=mols[i]->b2z;
p.push_back(c1);
p.push_back(b1);
p.push_back(b2);

c1.x=mols[j]->x;
c1.y=mols[j]->y;
c1.z=mols[j]->z;
b1.x=mols[j]->bx;
b1.y=mols[j]->by;
b1.z=mols[j]->bz;
b2.x=mols[j]->b2x;
b2.y=mols[j]->b2y;
b2.z=mols[j]->b2z;
p2.push_back(c1);
p2.push_back(b1);
p2.push_back(b2);

if(mols[i]->seq=='G'){type1=1;}else{
if(mols[i]->seq=='A'){type1=2;}else{
if(mols[i]->seq=='U'){type1=3;}else{
if(mols[i]->seq=='C'){type1=4;}else{
}
}
}
}
if(mols[j]->seq=='G'){type2=1;}else{
if(mols[j]->seq=='A'){type2=2;}else{
if(mols[j]->seq=='U'){type2=3;}else{
if(mols[j]->seq=='C'){type2=4;}else{
}
}
}
}
FABase * base1 = new FABase(type1);  
FABase * base2 = new FABase(type2);  
//Move to coordinate system
base1->alignBase(p);
base2->alignBase(p2);
//base1->print();
//base2->print();
//cout << "Types: " << type1 << ":" << type2 << endl; 
//Identify all donors-acceptor pairs

vector<Points> donors1;base1->donors(&donors1);
vector<Points> accept1;base1->acceptors(&accept1);
vector<Points> donors2;base2->donors(&donors2);
vector<Points> accept2;base2->acceptors(&accept2);

//Caclulate distance and two angles
//Acceptors have hydrogen
//donors have LP
//cout << "acceptors i -> donors j" << endl;
for(int a=0; a < accept1.size(); a=a+2){
for(int d=0; d < donors2.size(); d=d+2){
//angle between H and D
double dx = donors2[d].x  -accept1[a].x;
double dy = donors2[d].y  -accept1[a].y;
double dz = donors2[d].z  -accept1[a].z;
double dd = sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));

double px = accept1[a+1].x-donors2[d+1].x;
double py = accept1[a+1].y-donors2[d+1].y;
double pz = accept1[a+1].z-donors2[d+1].z;
double pp =  sqrt(pow(px,2)+pow(py,2)+pow(pz,2));

double hx = accept1[a+1].x-accept1[a].x;
double hy = accept1[a+1].y-accept1[a].y;
double hz = accept1[a+1].z-accept1[a].z;
double hh =  sqrt(pow(hx,2)+pow(hy,2)+pow(hz,2));

double lx = donors2[d+1].x-donors2[d].x;
double ly = donors2[d+1].y-donors2[d].y;
double lz = donors2[d+1].z-donors2[d].z;
double ll =  sqrt(pow(lx,2)+pow(ly,2)+pow(lz,2));

//angle between H and D
float ang1=(float)atanh((dx*hx+dy*hy+dz*hz)/(dd*hh));
float ang2=(float)atanh(-1*(dx*lx+dy*ly+dz*lz)/(dd*ll));
float dd1=(float)log(pow(pp,3));
float prob=base1->Eval(dd1,ang1,ang2,1);
allProb=allProb+prob;
//cout << dd1 << "," << ang1 << "," << ang2 << "," << prob << endl;
}
}
//cout << "donors i -> acceptors j" << endl;
for(int d=0; d < donors1.size(); d=d+2){
for(int a=0; a < accept2.size(); a=a+2){
//angle between H and D
double dx = donors1[d].x  -accept2[a].x;
double dy = donors1[d].y  -accept2[a].y;
double dz = donors1[d].z  -accept2[a].z;
double dd = sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));

double px = accept2[a+1].x-donors1[d+1].x;
double py = accept2[a+1].y-donors1[d+1].y;
double pz = accept2[a+1].z-donors1[d+1].z;
double pp =  sqrt(pow(px,2)+pow(py,2)+pow(pz,2));

double hx = accept2[a+1].x-accept2[a].x;
double hy = accept2[a+1].y-accept2[a].y;
double hz = accept2[a+1].z-accept2[a].z;
double hh =  sqrt(pow(hx,2)+pow(hy,2)+pow(hz,2));

double lx = donors1[d+1].x-donors1[d].x;
double ly = donors1[d+1].y-donors1[d].y;
double lz = donors1[d+1].z-donors1[d].z;
double ll =  sqrt(pow(lx,2)+pow(ly,2)+pow(lz,2));

//angle between H and D
float ang1=(float)atanh(cos(acos((dx*hx+dy*hy+dz*hz)/(dd*hh))));
float ang2=(float)atanh(cos(acos(-1*(dx*lx+dy*ly+dz*lz)/(dd*ll))));
float dd1=(float)log(pow(pp,3));
float prob=base1->Eval(dd1,ang1,ang2,1);

allProb=allProb+prob;
//cout << dd1 << "," << ang1 << "," << ang2 << "," << prob << endl;
}
}
delete base1;
delete base2;
}
}
}
}
return allProb;
}

float Polymere::groupStacked(int group_id){
//Counts the number of bases that have another base between 3 and 6.2 angstrums in the positive Z axis
//Plus is within 4 angstrums in the x,y plane
float allProb=0;
for(int i=0; i<num_mols; i++){
if(mols[i]->group==group_id){
for(int j=i+1; j<num_mols; j++){
if(mols[j]->group==group_id){
double mdist=100;
mdist =pow((mols[i]->x - mols[j]->x),2);
mdist =mdist + pow((mols[i]->y -mols[j]->y),2);
mdist =mdist + pow((mols[i]->z -mols[j]->z),2);
mdist = sqrt(mdist);
if(mdist< 10){
std::vector<double> plane=planar(i,j);
//cout << "i: " << i << " j: " << j << " : " << plane[0] << ":" << plane[1] << ":" << plane[2] << endl;
if(plane[2] < 6.2 && plane[2] > 3){
double rad=sqrt(pow(plane[0],2)+pow(plane[1],2));
if(rad < 4){
allProb=allProb+1;
}
}
}
}
}
}
}
return allProb;
}

float Polymere::adjStacked(std::vector<int> * stackprobs){
//Counts the number of bases that have another base between 3 and 6.2 angstrums in the positive Z axis
//Plus is within 4 angstrums in the x,y plane
//Adjusts the movement probability
float allProb=0;
for(int i=0; i<(num_mols-1); i++){
stackprobs->push_back(0);
int j=i+1;
double mdist=100;
mdist =pow((mols[i]->x - mols[j]->x),2);
mdist =mdist + pow((mols[i]->y -mols[j]->y),2);
mdist =mdist + pow((mols[i]->z -mols[j]->z),2);
mdist = sqrt(mdist);
if(mdist< 10){
std::vector<double> plane=planar(i,j);
//cout << "i: " << i << " j: " << j << " : " << plane[0] << ":" << plane[1] << ":" << plane[2] << endl;
if(plane[2] < 6.2 && plane[2] > 3){
double rad=sqrt(pow(plane[0],2)+pow(plane[1],2));
if(rad < 4){
stackprobs->at(i)=1;
allProb=allProb+1;
}
}
}
}
return allProb;
}


float Polymere::zaligned(int i,int j){
//Aligns the bases along the z-axis
double a1_1=mols[i]->x-mols[i]->bx;
double a2_1=mols[i]->y-mols[i]->by;
double a3_1=mols[i]->z-mols[i]->bz;

double b1_1=mols[i]->b2x-mols[i]->bx;
double b2_1=mols[i]->b2y-mols[i]->by;
double b3_1=mols[i]->b2z-mols[i]->bz;

double a1_2=mols[j]->x-mols[j]->bx;
double a2_2=mols[j]->y-mols[j]->by;
double a3_2=mols[j]->z-mols[j]->bz;

double b1_2=mols[j]->b2x-mols[j]->bx;
double b2_2=mols[j]->b2y-mols[j]->by;
double b3_2=mols[j]->b2z-mols[j]->bz;

double AxB1_x = (a2_1*b3_1 - a3_1*b2_1);
double AxB1_y = (a3_1*b1_1 - a1_1*b3_1);
double AxB1_z = (a1_1*b2_1 - a2_1*b1_1);

double AxB2_x = (a2_2*b3_2 - a3_2*b2_2);
double AxB2_y = (a3_2*b1_2 - a1_2*b3_2);
double AxB2_z = (a1_2*b2_2 - a2_2*b1_2);

float mag1=sqrt(pow(AxB1_x,2)+pow(AxB1_y,2)+pow(AxB1_z,2));
float mag2=sqrt(pow(AxB2_x,2)+pow(AxB2_y,2)+pow(AxB2_z,2));
float dprod=AxB1_x*AxB2_x+AxB1_y*AxB2_y+AxB1_z*AxB2_z;
float ang = acos(dprod/(mag1*mag2)); 
return ang;
}

float Polymere::stacked(){
//Counts the number of bases that have another base between 3 and 6.2 angstrums in the positive Z axis
//Plus is within 4 angstrums in the x,y plane
float allProb=0;
for(int i=0; i<num_mols; i++){
for(int j=i+1; j<num_mols; j++){
double mdist=100;
mdist =pow((mols[i]->x - mols[j]->x),2);
mdist =mdist + pow((mols[i]->y -mols[j]->y),2);
mdist =mdist + pow((mols[i]->z -mols[j]->z),2);
mdist = sqrt(mdist);
if(mdist< 10){
std::vector<double> plane=planar(i,j);
//cout << "i: " << i << " j: " << j << " : " << plane[0] << ":" << plane[1] << ":" << plane[2] << endl;
if(plane[2] < 6.2 && plane[2] > 3){
double rad=sqrt(pow(plane[0],2)+pow(plane[1],2));
if(rad < 4){
allProb=allProb+1;
}
}
}
}
}
return allProb;
}

float Polymere::hydroType(int i,int j,std::pair<int,int> * faces){
float allProb=0;

vector <Points> p;
vector <Points> p2;
Points c1,b1,b2;
int type1,type2;
//Returns a score for hydrogen bonding
c1.x=mols[i]->x;
c1.y=mols[i]->y;
c1.z=mols[i]->z;
b1.x=mols[i]->bx;
b1.y=mols[i]->by;
b1.z=mols[i]->bz;
b2.x=mols[i]->b2x;
b2.y=mols[i]->b2y;
b2.z=mols[i]->b2z;
p.push_back(c1);
p.push_back(b1);
p.push_back(b2);

c1.x=mols[j]->x;
c1.y=mols[j]->y;
c1.z=mols[j]->z;
b1.x=mols[j]->bx;
b1.y=mols[j]->by;
b1.z=mols[j]->bz;
b2.x=mols[j]->b2x;
b2.y=mols[j]->b2y;
b2.z=mols[j]->b2z;
p2.push_back(c1);
p2.push_back(b1);
p2.push_back(b2);

if(mols[i]->seq=='G'){type1=1;}else{
if(mols[i]->seq=='A'){type1=2;}else{
if(mols[i]->seq=='U'){type1=3;}else{
if(mols[i]->seq=='C'){type1=4;}else{
}
}
}
}
if(mols[j]->seq=='G'){type2=1;}else{
if(mols[j]->seq=='A'){type2=2;}else{
if(mols[j]->seq=='U'){type2=3;}else{
if(mols[j]->seq=='C'){type2=4;}else{
}
}
}
}
FABase * base1 = new FABase(type1);  
FABase * base2 = new FABase(type2);  
//Move to coordinate system
base1->alignBase(p);
base2->alignBase(p2);
//base1->print();
//base2->print();
//cout << "Types: " << type1 << ":" << type2 << endl; 
//Identify all donors-acceptor pairs

std::vector<Points> donors1;base1->donors(&donors1);
std::vector<Points> accept1;base1->acceptors(&accept1);
std::vector<Points> donors2;base2->donors(&donors2);
std::vector<Points> accept2;base2->acceptors(&accept2);

std::vector<int> aface_type1;base1->acceptorsFace(&aface_type1);
std::vector<int> dface_type1;base1->donorsFace(&dface_type1);

std::vector<int> aface_type2;base2->acceptorsFace(&aface_type2);
std::vector<int> dface_type2;base2->donorsFace(&dface_type2);

std::vector<int> faces1;
std::vector<int> faces2;

for(int a=0; a < accept1.size(); a=a+2){
float a_max=0;
int aface_type=-1;
int dface_type=-1;
if(aface_type1[a]>=0){
for(int d=0; d < donors2.size(); d=d+2){
//angle between H and D
if(dface_type2[d]>=0){
float dx =(float) donors2[d].x  -accept1[a].x;
float dy =(float) donors2[d].y  -accept1[a].y;
float dz =(float) donors2[d].z  -accept1[a].z;
float dd =(float) sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));
          
float px =(float) accept1[a+1].x-donors2[d+1].x;
float py =(float) accept1[a+1].y-donors2[d+1].y;
float pz =(float) accept1[a+1].z-donors2[d+1].z;
float pp =(float)  sqrt(pow(px,2)+pow(py,2)+pow(pz,2));
         
float hx =(float) accept1[a+1].x-accept1[a].x;
float hy =(float) accept1[a+1].y-accept1[a].y;
float hz =(float) accept1[a+1].z-accept1[a].z;
float hh =(float)  sqrt(pow(hx,2)+pow(hy,2)+pow(hz,2));
        
float lx =(float) donors2[d+1].x-donors2[d].x;
float ly =(float) donors2[d+1].y-donors2[d].y;
float lz =(float) donors2[d+1].z-donors2[d].z;
float ll =(float)  sqrt(pow(lx,2)+pow(ly,2)+pow(lz,2));

//angle between H and D
float ang1=(float)atanh((dx*hx+dy*hy+dz*hz)/(dd*hh));
float ang2=(float)atanh(-1*(dx*lx+dy*ly+dz*lz)/(dd*ll));
float dd1=(float)log(pow(pp,3));
float prob=base1->Eval(dd1,ang1,ang2,1);
if(prob > a_max){
a_max=prob;
aface_type=aface_type1[a];
dface_type=dface_type2[d];
}
//cout << dd1 << "," << ang1 << "," << ang2 << "," << prob << "," << dd << endl;
}
}
}
allProb=allProb+a_max;
if(a_max > 0.25){
faces1.push_back(aface_type);
faces2.push_back(dface_type);
}
}

for(int d=0; d < donors1.size(); d=d+2){
float a_max=0;
int aface_type=-1;
int dface_type=-1;
if(dface_type1[d]>=0){
for(int a=0; a < accept2.size(); a=a+2){
//angle between H and D
if(dface_type2[a]>=0){
float dx =(float)  donors1[d].x  -accept2[a].x;
float dy =(float)  donors1[d].y  -accept2[a].y;
float dz =(float)  donors1[d].z  -accept2[a].z;
float dd =(float)  sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));
           
float px =(float)  accept2[a+1].x-donors1[d+1].x;
float py =(float)  accept2[a+1].y-donors1[d+1].y;
float pz =(float)  accept2[a+1].z-donors1[d+1].z;
float pp =(float)   sqrt(pow(px,2)+pow(py,2)+pow(pz,2));
                   
float hx =(float)  accept2[a+1].x-accept2[a].x;
float hy =(float)  accept2[a+1].y-accept2[a].y;
float hz =(float)  accept2[a+1].z-accept2[a].z;
float hh =(float)   sqrt(pow(hx,2)+pow(hy,2)+pow(hz,2));
                   
float lx =(float)  donors1[d+1].x-donors1[d].x;
float ly =(float)  donors1[d+1].y-donors1[d].y;
float lz =(float)  donors1[d+1].z-donors1[d].z;
float ll =(float)   sqrt(pow(lx,2)+pow(ly,2)+pow(lz,2));

//angle between H and D
float ang1=(float)atanh(cos(acos((dx*hx+dy*hy+dz*hz)/(dd*hh))));
float ang2=(float)atanh(cos(acos(-1*(dx*lx+dy*ly+dz*lz)/(dd*ll))));
float dd1=(float)log(pow(pp,3));
float prob=base1->Eval(dd1,ang1,ang2,1);
if(prob > a_max){
a_max=prob;
aface_type=aface_type2[a];
dface_type=dface_type1[d];
}
}

}
}
allProb=allProb+a_max;
if(a_max > 0.25){
faces1.push_back(dface_type);
faces2.push_back(aface_type);
}
}

//If faces1 and faces2 are uniformly on side or the other 
if(faces1.size() > 0 && faces2.size() > 0){
int ftype1=faces1[0];
int ftype2=faces2[0];
for(int i=0; i<faces1.size(); i++){
//printf("faces1: %d\n",faces1[i]); 
if(faces1[i]!=ftype1){ftype1=-1;} }
for(int i=0; i<faces2.size(); i++){
//printf("faces2: %d\n",faces2[i]);
 if(faces1[i]!=ftype1){ftype1=-1;} }

faces->first=ftype1;
faces->second=ftype2;
}
delete base1;
delete base2;
return allProb;
}


float Polymere::hydroP(int i,int j){
float allProb=0;
//Only Consider Pairs that are in the movement groups
//cout << "Molsi: " << i << ":" << mols[i]->group << " Molsj: " << j << ":" << mols[j]->group << endl; 
vector <Points> p;
vector <Points> p2;
Points c1,b1,b2;
int type1,type2;
//Returns a score for hydrogen bonding
c1.x=mols[i]->x;
c1.y=mols[i]->y;
c1.z=mols[i]->z;
b1.x=mols[i]->bx;
b1.y=mols[i]->by;
b1.z=mols[i]->bz;
b2.x=mols[i]->b2x;
b2.y=mols[i]->b2y;
b2.z=mols[i]->b2z;
p.push_back(c1);
p.push_back(b1);
p.push_back(b2);

c1.x=mols[j]->x;
c1.y=mols[j]->y;
c1.z=mols[j]->z;
b1.x=mols[j]->bx;
b1.y=mols[j]->by;
b1.z=mols[j]->bz;
b2.x=mols[j]->b2x;
b2.y=mols[j]->b2y;
b2.z=mols[j]->b2z;
p2.push_back(c1);
p2.push_back(b1);
p2.push_back(b2);

if(mols[i]->seq=='G'){type1=1;}else{
if(mols[i]->seq=='A'){type1=2;}else{
if(mols[i]->seq=='U'){type1=3;}else{
if(mols[i]->seq=='C'){type1=4;}else{
}
}
}
}
if(mols[j]->seq=='G'){type2=1;}else{
if(mols[j]->seq=='A'){type2=2;}else{
if(mols[j]->seq=='U'){type2=3;}else{
if(mols[j]->seq=='C'){type2=4;}else{
}
}
}
}
FABase * base1 = new FABase(type1);  
FABase * base2 = new FABase(type2);  
//Move to coordinate system
base1->alignBase(p);
base2->alignBase(p2);
//base1->print();
//base2->print();
//cout << "Types: " << type1 << ":" << type2 << endl; 
//Identify all donors-acceptor pairs

vector<Points> donors1;base1->donors(&donors1);
vector<Points> accept1;base1->acceptors(&accept1);
vector<Points> donors2;base2->donors(&donors2);
vector<Points> accept2;base2->acceptors(&accept2);

//Caclulate distance and two angles
//Acceptors have hydrogen
//donors have LP
//cout << "acceptors i -> donors j" << endl;
//acceptor can only match one donor. Need to update this in other scores as well
for(int a=0; a < accept1.size(); a=a+2){
float a_max=0;
for(int d=0; d < donors2.size(); d=d+2){
//angle between H and D
float dx =(float) donors2[d].x  -accept1[a].x;
float dy =(float) donors2[d].y  -accept1[a].y;
float dz =(float) donors2[d].z  -accept1[a].z;
float dd =(float) sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));
          
float px =(float) accept1[a+1].x-donors2[d+1].x;
float py =(float) accept1[a+1].y-donors2[d+1].y;
float pz =(float) accept1[a+1].z-donors2[d+1].z;
float pp =(float)  sqrt(pow(px,2)+pow(py,2)+pow(pz,2));
         
float hx =(float) accept1[a+1].x-accept1[a].x;
float hy =(float) accept1[a+1].y-accept1[a].y;
float hz =(float) accept1[a+1].z-accept1[a].z;
float hh =(float)  sqrt(pow(hx,2)+pow(hy,2)+pow(hz,2));
        
float lx =(float) donors2[d+1].x-donors2[d].x;
float ly =(float) donors2[d+1].y-donors2[d].y;
float lz =(float) donors2[d+1].z-donors2[d].z;
float ll =(float)  sqrt(pow(lx,2)+pow(ly,2)+pow(lz,2));

//angle between H and D
float ang1=(float)atanh((dx*hx+dy*hy+dz*hz)/(dd*hh));
float ang2=(float)atanh(-1*(dx*lx+dy*ly+dz*lz)/(dd*ll));
float dd1=(float)log(pow(pp,3));
float prob=base1->Eval(dd1,ang1,ang2,1);
if(prob > a_max){
a_max=prob;
}
//cout << dd1 << "," << ang1 << "," << ang2 << "," << prob << "," << dd << endl;
}
allProb=allProb+a_max;
}

for(int d=0; d < donors1.size(); d=d+2){
float a_max=0;
for(int a=0; a < accept2.size(); a=a+2){
//angle between H and D
float dx =(float)  donors1[d].x  -accept2[a].x;
float dy =(float)  donors1[d].y  -accept2[a].y;
float dz =(float)  donors1[d].z  -accept2[a].z;
float dd =(float)  sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));
           
float px =(float)  accept2[a+1].x-donors1[d+1].x;
float py =(float)  accept2[a+1].y-donors1[d+1].y;
float pz =(float)  accept2[a+1].z-donors1[d+1].z;
float pp =(float)   sqrt(pow(px,2)+pow(py,2)+pow(pz,2));
                   
float hx =(float)  accept2[a+1].x-accept2[a].x;
float hy =(float)  accept2[a+1].y-accept2[a].y;
float hz =(float)  accept2[a+1].z-accept2[a].z;
float hh =(float)   sqrt(pow(hx,2)+pow(hy,2)+pow(hz,2));
                   
float lx =(float)  donors1[d+1].x-donors1[d].x;
float ly =(float)  donors1[d+1].y-donors1[d].y;
float lz =(float)  donors1[d+1].z-donors1[d].z;
float ll =(float)   sqrt(pow(lx,2)+pow(ly,2)+pow(lz,2));

//angle between H and D
float ang1=(float)atanh(cos(acos((dx*hx+dy*hy+dz*hz)/(dd*hh))));
float ang2=(float)atanh(cos(acos(-1*(dx*lx+dy*ly+dz*lz)/(dd*ll))));
float dd1=(float)log(pow(pp,3));
float prob=base1->Eval(dd1,ang1,ang2,1);
if(prob > a_max){
a_max=prob;
}
//cout << dd1 << "," << ang1 << "," << ang2 << "," << prob << endl;
}
allProb=allProb+a_max;
}
delete base1;
delete base2;
return allProb;
}

float Polymere::minDD(int i, int j){
float prob = 100;
vector <Points> p;
vector <Points> p2;
Points c1,b1,b2;
int type1,type2;
//Returns a score for hydrogen bonding
c1.x=mols[i]->x;
c1.y=mols[i]->y;
c1.z=mols[i]->z;
b1.x=mols[i]->bx;
b1.y=mols[i]->by;
b1.z=mols[i]->bz;
b2.x=mols[i]->b2x;
b2.y=mols[i]->b2y;
b2.z=mols[i]->b2z;
p.push_back(c1);
p.push_back(b1);
p.push_back(b2);

c1.x=mols[j]->x;
c1.y=mols[j]->y;
c1.z=mols[j]->z;
b1.x=mols[j]->bx;
b1.y=mols[j]->by;
b1.z=mols[j]->bz;
b2.x=mols[j]->b2x;
b2.y=mols[j]->b2y;
b2.z=mols[j]->b2z;
p2.push_back(c1);
p2.push_back(b1);
p2.push_back(b2);

if(mols[i]->seq=='G'){type1=1;}else{
if(mols[i]->seq=='A'){type1=2;}else{
if(mols[i]->seq=='U'){type1=3;}else{
if(mols[i]->seq=='C'){type1=4;}else{
}
}
}
}
if(mols[j]->seq=='G'){type2=1;}else{
if(mols[j]->seq=='A'){type2=2;}else{
if(mols[j]->seq=='U'){type2=3;}else{
if(mols[j]->seq=='C'){type2=4;}else{
}
}
}
}
FABase * base1 = new FABase(type1);  
FABase * base2 = new FABase(type2);  
//Move to coordinate system
base1->alignBase(p);
base2->alignBase(p2);
//base1->print();
//base2->print();
//cout << "Types: " << type1 << ":" << type2 << endl; 
//Identify all donors-acceptor pairs
vector<Points> donors1;base1->donors(&donors1);
vector<Points> donors2;base2->donors(&donors2);
for(int a=0; a < donors1.size(); a=a+2){
for(int d=0; d < donors2.size(); d=d+2){
//angle between H and D
double dx = donors2[d].x  -donors1[a].x;
double dy = donors2[d].y  -donors1[a].y;
double dz = donors2[d].z  -donors1[a].z;
double dd = sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));
if(dd < prob){
prob=dd;
}
}
}
delete base1;
delete base2;
return prob;
}

float Polymere::minAA(int i,int j){
float prob = 100;
vector <Points> p;
vector <Points> p2;
Points c1,b1,b2;
int type1,type2;
//Returns a score for hydrogen bonding
c1.x=mols[i]->x;
c1.y=mols[i]->y;
c1.z=mols[i]->z;
b1.x=mols[i]->bx;
b1.y=mols[i]->by;
b1.z=mols[i]->bz;
b2.x=mols[i]->b2x;
b2.y=mols[i]->b2y;
b2.z=mols[i]->b2z;
p.push_back(c1);
p.push_back(b1);
p.push_back(b2);

c1.x=mols[j]->x;
c1.y=mols[j]->y;
c1.z=mols[j]->z;
b1.x=mols[j]->bx;
b1.y=mols[j]->by;
b1.z=mols[j]->bz;
b2.x=mols[j]->b2x;
b2.y=mols[j]->b2y;
b2.z=mols[j]->b2z;
p2.push_back(c1);
p2.push_back(b1);
p2.push_back(b2);

if(mols[i]->seq=='G'){type1=1;}else{
if(mols[i]->seq=='A'){type1=2;}else{
if(mols[i]->seq=='U'){type1=3;}else{
if(mols[i]->seq=='C'){type1=4;}else{
}
}
}
}
if(mols[j]->seq=='G'){type2=1;}else{
if(mols[j]->seq=='A'){type2=2;}else{
if(mols[j]->seq=='U'){type2=3;}else{
if(mols[j]->seq=='C'){type2=4;}else{
}
}
}
}
FABase * base1 = new FABase(type1);  
FABase * base2 = new FABase(type2);  
//Move to coordinate system
base1->alignBase(p);
base2->alignBase(p2);
//base1->print();
//base2->print();
//cout << "Types: " << type1 << ":" << type2 << endl; 
//Identify all donors-acceptor pairs
vector<Points> accept1;base1->acceptors(&accept1);
vector<Points> accept2;base2->acceptors(&accept2);
for(int a=0; a < accept1.size(); a=a+2){
for(int d=0; d < accept2.size(); d=d+2){
//angle between H and D
double dx = accept2[d].x  -accept1[a].x;
double dy = accept2[d].y  -accept1[a].y;
double dz = accept2[d].z  -accept1[a].z;
double dd = sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));
if(dd < prob){
prob=dd;
}
}
}
delete base1;
delete base2;
return prob;
}


float Polymere::hydroG(const std::vector<int> * group){
//Calculates Hydrogen Bonding for a particular movment group
float allProb=0;
for(int t1=0; t1 < group->size(); t1++){
int i = group->at(t1);
for(int t2=0; t2< group->size(); t2++){
int j = group->at(t2);
// Only consider bases that are separated by at least 2 nucleotides
if(j > i+2){
double mdist=100;
mdist =pow((mols[i]->x - mols[j]->x),2);
mdist =mdist + pow((mols[i]->y -mols[j]->y),2);
mdist =mdist + pow((mols[i]->z -mols[j]->z),2);
mdist = sqrt(mdist);
//cout << "Type: " << params[0] << ":" << i << "," << j << " i: " << mols[i]->seq << " j: " <<  mols[j]->seq << " mdist: " << mdist << endl;
//cout << "mdist: " << mdist << endl;
double zplane=coplanar(i,j);
if(zplane < 1){
if(mdist< 20){
vector <Points> p;
vector <Points> p2;
Points c1,b1,b2;
int type1,type2;
//Returns a score for hydrogen bonding
c1.x=mols[i]->x;
c1.y=mols[i]->y;
c1.z=mols[i]->z;
b1.x=mols[i]->bx;
b1.y=mols[i]->by;
b1.z=mols[i]->bz;
b2.x=mols[i]->b2x;
b2.y=mols[i]->b2y;
b2.z=mols[i]->b2z;
p.push_back(c1);
p.push_back(b1);
p.push_back(b2);

c1.x=mols[j]->x;
c1.y=mols[j]->y;
c1.z=mols[j]->z;
b1.x=mols[j]->bx;
b1.y=mols[j]->by;
b1.z=mols[j]->bz;
b2.x=mols[j]->b2x;
b2.y=mols[j]->b2y;
b2.z=mols[j]->b2z;
p2.push_back(c1);
p2.push_back(b1);
p2.push_back(b2);

if(mols[i]->seq=='G'){type1=1;}else{
if(mols[i]->seq=='A'){type1=2;}else{
if(mols[i]->seq=='U'){type1=3;}else{
if(mols[i]->seq=='C'){type1=4;}else{
}
}
}
}
if(mols[j]->seq=='G'){type2=1;}else{
if(mols[j]->seq=='A'){type2=2;}else{
if(mols[j]->seq=='U'){type2=3;}else{
if(mols[j]->seq=='C'){type2=4;}else{
}
}
}
}
FABase * base1 = new FABase(type1);  
FABase * base2 = new FABase(type2);  
//Move to coordinate system
base1->alignBase(p);
base2->alignBase(p2);
//base1->print();
//base2->print();
//cout << "Types: " << type1 << ":" << type2 << endl; 
//Identify all donors-acceptor pairs

vector<Points> donors1;base1->donors(&donors1);
vector<Points> accept1;base1->acceptors(&accept1);
vector<Points> donors2;base2->donors(&donors2);
vector<Points> accept2;base2->acceptors(&accept2);

//Caclulate distance and two angles
//Acceptors have hydrogen
//donors have LP
//cout << "acceptors i -> donors j" << endl;
for(int a=0; a < accept1.size(); a=a+2){
for(int d=0; d < donors2.size(); d=d+2){
//angle between H and D
double dx = donors2[d].x  -accept1[a].x;
double dy = donors2[d].y  -accept1[a].y;
double dz = donors2[d].z  -accept1[a].z;
double dd = sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));

double px = accept1[a+1].x-donors2[d+1].x;
double py = accept1[a+1].y-donors2[d+1].y;
double pz = accept1[a+1].z-donors2[d+1].z;
double pp =  sqrt(pow(px,2)+pow(py,2)+pow(pz,2));

double hx = accept1[a+1].x-accept1[a].x;
double hy = accept1[a+1].y-accept1[a].y;
double hz = accept1[a+1].z-accept1[a].z;
double hh =  sqrt(pow(hx,2)+pow(hy,2)+pow(hz,2));

double lx = donors2[d+1].x-donors2[d].x;
double ly = donors2[d+1].y-donors2[d].y;
double lz = donors2[d+1].z-donors2[d].z;
double ll =  sqrt(pow(lx,2)+pow(ly,2)+pow(lz,2));

//angle between H and D
float ang1=(float)atanh((dx*hx+dy*hy+dz*hz)/(dd*hh));
float ang2=(float)atanh(-1*(dx*lx+dy*ly+dz*lz)/(dd*ll));
float dd1=(float)log(pow(pp,3));
float prob=base1->Eval(dd1,ang1,ang2,1);
allProb=allProb+prob;
//cout << dd1 << "," << ang1 << "," << ang2 << "," << prob << endl;
}
}
//Penalize for acceptor acceptor distance
for(int a=0; a < donors1.size(); a=a+2){
for(int d=0; d < donors2.size(); d=d+2){
//angle between H and D
double dx = donors2[d].x  -donors1[a].x;
double dy = donors2[d].y  -donors1[a].y;
double dz = donors2[d].z  -donors1[a].z;
double dd = sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));
float prob = 0;
if(dd < 2.2){prob = -10;}
allProb=allProb+prob;
//cout << dd1 << "," << ang1 << "," << ang2 << "," << prob << endl;
}
}
for(int a=0; a < accept1.size(); a=a+2){
for(int d=0; d < accept2.size(); d=d+2){
//angle between H and D
double dx = accept2[d].x  -accept1[a].x;
double dy = accept2[d].y  -accept1[a].y;
double dz = accept2[d].z  -accept1[a].z;
double dd = sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));
float prob = 0;
if(dd < 2.2){prob = -10;}
allProb=allProb+prob;
//cout << dd1 << "," << ang1 << "," << ang2 << "," << prob << endl;
}
}
//cout << "donors i -> acceptors j" << endl;
for(int d=0; d < donors1.size(); d=d+2){
for(int a=0; a < accept2.size(); a=a+2){
//angle between H and D
double dx = donors1[d].x  -accept2[a].x;
double dy = donors1[d].y  -accept2[a].y;
double dz = donors1[d].z  -accept2[a].z;
double dd = sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));

double px = accept2[a+1].x-donors1[d+1].x;
double py = accept2[a+1].y-donors1[d+1].y;
double pz = accept2[a+1].z-donors1[d+1].z;
double pp =  sqrt(pow(px,2)+pow(py,2)+pow(pz,2));

double hx = accept2[a+1].x-accept2[a].x;
double hy = accept2[a+1].y-accept2[a].y;
double hz = accept2[a+1].z-accept2[a].z;
double hh =  sqrt(pow(hx,2)+pow(hy,2)+pow(hz,2));

double lx = donors1[d+1].x-donors1[d].x;
double ly = donors1[d+1].y-donors1[d].y;
double lz = donors1[d+1].z-donors1[d].z;
double ll =  sqrt(pow(lx,2)+pow(ly,2)+pow(lz,2));

//angle between H and D
float ang1=(float)atanh(cos(acos((dx*hx+dy*hy+dz*hz)/(dd*hh))));
float ang2=(float)atanh(cos(acos(-1*(dx*lx+dy*ly+dz*lz)/(dd*ll))));
float dd1=(float)log(pow(pp,3));
float prob=base1->Eval(dd1,ang1,ang2,1);

allProb=allProb+prob;
//cout << dd1 << "," << ang1 << "," << ang2 << "," << prob << endl;
}
}
delete base1;
delete base2;
}
}
}
}
}
return allProb;
}
float Polymere::hydro(){
float allProb=0;
for(int i=0; i<num_mols; i++){
for(int j=(i+3); j<num_mols; j++){
//Only Consider Pairs that are in the movement groups
//cout << "Molsi: " << i << ":" << mols[i]->group << " Molsj: " << j << ":" << mols[j]->group << endl; 
double mdist=100;
if((mols[i]->group==1 || mols[i]->group==2) && (mols[j]->group==1 || mols[j]->group==2) && mols[i]->id < 0 && mols[j]->id < 0){ 
mdist =pow((mols[i]->x - mols[j]->x),2);
mdist =mdist + pow((mols[i]->y -mols[j]->y),2);
mdist =mdist + pow((mols[i]->z -mols[j]->z),2);
mdist = sqrt(mdist);
}else{mdist=100;}
//cout << "Type: " << params[0] << ":" << i << "," << j << " i: " << mols[i]->seq << " j: " <<  mols[j]->seq << " mdist: " << mdist << endl;
//cout << "mdist: " << mdist << endl;
double zplane=coplanar(i,j);
if(zplane < 1){
if(mdist< 20){
vector <Points> p;
vector <Points> p2;
Points c1,b1,b2;
int type1,type2;
//Returns a score for hydrogen bonding
c1.x=mols[i]->x;
c1.y=mols[i]->y;
c1.z=mols[i]->z;
b1.x=mols[i]->bx;
b1.y=mols[i]->by;
b1.z=mols[i]->bz;
b2.x=mols[i]->b2x;
b2.y=mols[i]->b2y;
b2.z=mols[i]->b2z;
p.push_back(c1);
p.push_back(b1);
p.push_back(b2);

c1.x=mols[j]->x;
c1.y=mols[j]->y;
c1.z=mols[j]->z;
b1.x=mols[j]->bx;
b1.y=mols[j]->by;
b1.z=mols[j]->bz;
b2.x=mols[j]->b2x;
b2.y=mols[j]->b2y;
b2.z=mols[j]->b2z;
p2.push_back(c1);
p2.push_back(b1);
p2.push_back(b2);

if(mols[i]->seq=='G'){type1=1;}else{
if(mols[i]->seq=='A'){type1=2;}else{
if(mols[i]->seq=='U'){type1=3;}else{
if(mols[i]->seq=='C'){type1=4;}else{
}
}
}
}
if(mols[j]->seq=='G'){type2=1;}else{
if(mols[j]->seq=='A'){type2=2;}else{
if(mols[j]->seq=='U'){type2=3;}else{
if(mols[j]->seq=='C'){type2=4;}else{
}
}
}
}
FABase * base1 = new FABase(type1);  
FABase * base2 = new FABase(type2);  
//Move to coordinate system
base1->alignBase(p);
base2->alignBase(p2);
//base1->print();
//base2->print();
//cout << "Types: " << type1 << ":" << type2 << endl; 
//Identify all donors-acceptor pairs

vector<Points> donors1;base1->donors(&donors1);
vector<Points> accept1;base1->acceptors(&accept1);
vector<Points> donors2;base2->donors(&donors2);
vector<Points> accept2;base2->acceptors(&accept2);

//Caclulate distance and two angles
//Acceptors have hydrogen
//donors have LP
//cout << "acceptors i -> donors j" << endl;
for(int a=0; a < accept1.size(); a=a+2){
for(int d=0; d < donors2.size(); d=d+2){
//angle between H and D
double dx = donors2[d].x  -accept1[a].x;
double dy = donors2[d].y  -accept1[a].y;
double dz = donors2[d].z  -accept1[a].z;
double dd = sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));

double px = accept1[a+1].x-donors2[d+1].x;
double py = accept1[a+1].y-donors2[d+1].y;
double pz = accept1[a+1].z-donors2[d+1].z;
double pp =  sqrt(pow(px,2)+pow(py,2)+pow(pz,2));

double hx = accept1[a+1].x-accept1[a].x;
double hy = accept1[a+1].y-accept1[a].y;
double hz = accept1[a+1].z-accept1[a].z;
double hh =  sqrt(pow(hx,2)+pow(hy,2)+pow(hz,2));

double lx = donors2[d+1].x-donors2[d].x;
double ly = donors2[d+1].y-donors2[d].y;
double lz = donors2[d+1].z-donors2[d].z;
double ll =  sqrt(pow(lx,2)+pow(ly,2)+pow(lz,2));

//angle between H and D
float ang1=(float)atanh((dx*hx+dy*hy+dz*hz)/(dd*hh));
float ang2=(float)atanh(-1*(dx*lx+dy*ly+dz*lz)/(dd*ll));
float dd1=(float)log(pow(pp,3));
float prob=base1->Eval(dd1,ang1,ang2,1);
allProb=allProb+prob;
//cout << dd1 << "," << ang1 << "," << ang2 << "," << prob << endl;
}
}
//Penalize for acceptor acceptor distance
for(int a=0; a < donors1.size(); a=a+2){
for(int d=0; d < donors2.size(); d=d+2){
//angle between H and D
double dx = donors2[d].x  -donors1[a].x;
double dy = donors2[d].y  -donors1[a].y;
double dz = donors2[d].z  -donors1[a].z;
double dd = sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));
float prob = 0;
if(dd < 2.2){prob = -10;}
allProb=allProb+prob;
//cout << dd1 << "," << ang1 << "," << ang2 << "," << prob << endl;
}
}
for(int a=0; a < accept1.size(); a=a+2){
for(int d=0; d < accept2.size(); d=d+2){
//angle between H and D
double dx = accept2[d].x  -accept1[a].x;
double dy = accept2[d].y  -accept1[a].y;
double dz = accept2[d].z  -accept1[a].z;
double dd = sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));
float prob = 0;
if(dd < 2.2){prob = -10;}
allProb=allProb+prob;
//cout << dd1 << "," << ang1 << "," << ang2 << "," << prob << endl;
}
}
//cout << "donors i -> acceptors j" << endl;
for(int d=0; d < donors1.size(); d=d+2){
for(int a=0; a < accept2.size(); a=a+2){
//angle between H and D
double dx = donors1[d].x  -accept2[a].x;
double dy = donors1[d].y  -accept2[a].y;
double dz = donors1[d].z  -accept2[a].z;
double dd = sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));

double px = accept2[a+1].x-donors1[d+1].x;
double py = accept2[a+1].y-donors1[d+1].y;
double pz = accept2[a+1].z-donors1[d+1].z;
double pp =  sqrt(pow(px,2)+pow(py,2)+pow(pz,2));

double hx = accept2[a+1].x-accept2[a].x;
double hy = accept2[a+1].y-accept2[a].y;
double hz = accept2[a+1].z-accept2[a].z;
double hh =  sqrt(pow(hx,2)+pow(hy,2)+pow(hz,2));

double lx = donors1[d+1].x-donors1[d].x;
double ly = donors1[d+1].y-donors1[d].y;
double lz = donors1[d+1].z-donors1[d].z;
double ll =  sqrt(pow(lx,2)+pow(ly,2)+pow(lz,2));

//angle between H and D
float ang1=(float)atanh(cos(acos((dx*hx+dy*hy+dz*hz)/(dd*hh))));
float ang2=(float)atanh(cos(acos(-1*(dx*lx+dy*ly+dz*lz)/(dd*ll))));
float dd1=(float)log(pow(pp,3));
float prob=base1->Eval(dd1,ang1,ang2,1);

allProb=allProb+prob;
//cout << dd1 << "," << ang1 << "," << ang2 << "," << prob << endl;
}
}
delete base1;
delete base2;
}
}
}
}
return allProb;
}
void Polymere::doubletV(const vector<double> * vec,double a1,double a2,double a3,double a4,double a5,vector<float> * probs){

float dist;
float score2=0;
probs->clear();
//Initialize Probs vectors
for(int i = 0; i<(num_mols); i++){probs->push_back(0);}
for(int i = 0; i<(num_mols); i++){
for(int j = 0; j<(num_mols); j++){
if(abs(i-j)>1){
dist=c1distf(i,j);
if(dist<20){
score2=(float)gaussDbl_ij(vec,i,j,a1,a2,a3,a4,a5);
probs->at(i)=probs->at(i)+score2;
probs->at(j)=probs->at(j)+score2;
}
}
}
}

}

void Polymere::doubletVg(const vector<double> * vec,double a1,double a2,double a3,double a4,double a5,vector<float> * probs,int movegroup){

float dist;
float score2=0;
probs->clear();
//Initialize Probs vectors
std::vector<int> gmap;
int t_pos=0;
for(int i = 0; i<(num_mols); i++){
if(mols[i]->group==movegroup){
probs->push_back(0);
gmap.push_back(t_pos);
t_pos++;
}else{
gmap.push_back(-1);
}

}
for(int i = 0; i<(num_mols); i++){
for(int j = 0; j<(num_mols); j++){
if(abs(i-j)>1 && mols[i]->group==movegroup && mols[j]->group==movegroup){
dist=c1distf(i,j);
if(dist<20){
score2=(float)gaussDbl_ij(vec,i,j,a1,a2,a3,a4,a5);
probs->at(gmap[i])=probs->at(gmap[i])+score2;
probs->at(gmap[j])=probs->at(gmap[j])+score2;
}
}
}
}

}

double Polymere::gaussDbl(const vector<double> * vec,double a1,double a2,double a3,double a4,double a5,const vector<int> * group){
//Adjust Group Values to reflect local score
//Fragments affect +0 +1 +2 positions
//Sum these scores to get a local score
std::vector<float> probs;
for(int i=0; i<group->size(); i++){
probs.push_back(0);
}

double score1=0;
float dist;
double score2;
for(int t1 = 0; t1<group->size(); t1++){
int i = group->at(t1);
for(int t2 = 0; t2<group->size(); t2++){
int j = group->at(t2);
if(i!=j){
dist=c1distf(i,j);
if(dist<20){
score2=gaussDbl_ij(vec,i,j,a1,a2,a3,a4,a5);
probs[t1]=probs[t1]+(float)score2;
//probs[(i+1)]=probs[(i+1)]+(float)score2;
//probs[(i+2)]=probs[(i+2)]+(float)score2;
probs[t2]=probs[t2]+(float)score2;
//probs[(j+1)]=probs[(j+1)]+(float)score2;
//probs[(j+2)]=probs[(j+2)]+(float)score2;
score1 = score1+score2;
}
}
}
}
//Adjust groups to reflect percentage of score
for(int i=0; i<probs.size(); i++){
//The score is for 3 nucleotide positions
probs[i]=probs[i]/(float)score1 * ((float)probs.size());
}
for(int i=0; i<group->size(); i++){
mols[group->at(i)]->prob=probs[i];
}
//Adjusted score that rewards structures that have even distributions of goodness
float avescore = score1/group->size();
float newscore = 0;
for(int i=0; i<probs.size(); i++){
if(probs[i] > avescore){
newscore=newscore+avescore;
}else{
newscore=newscore+probs[i];
}
}

score1=(double)newscore;
return score1;
}

void Polymere::cluster2(std::vector<Polymere *> * polymeres, vector<int>  * centers,int max_size,int max_iter){
//max_size = maximum size of any cluster
//This tries to cluster the vector of polymeres by RMSD
//Doesn't store the pair-wise rmsd's

//Create ave weights vectors;
std::vector<float> aves(centers->size(),0);
//1) Identify which cluster each polymere belongs too
//array of all pairwise rmsds
int size1 = polymeres->size();
int w=0;
float minrmsd=9999;   //stores the current minimal rmsd to the current center
int minclus;          //stores the index of the cluster that is closest
std::vector< std::vector<int> > clusters; //Vector of vectors containing indices of polymers in each cluster
for(int c=0; c<centers->size(); c++){
clusters.push_back(std::vector<int>());
}
float cur_ave=0;
printf("Total Clusters: %d\n",centers->size());
for(int iter=0; iter < max_iter; iter++){
float prev_ave=0;
for(int a=0; a<aves.size(); a++){
prev_ave=(aves[a]/(1+clusters[a].size()))+prev_ave;
}

printf("Iteration: %d - %8.3f\n",iter,prev_ave);

for(int c=0; c <  centers->size(); c++){
clusters[c].clear();
aves[c]=0;
} 

//Loop through all structures and assign the to the cluster with the nearest center
for(int i = 0; i<size1; i++){
minrmsd=9999;
float currmsd=99999;
minclus=0;
for(int j = 0; j<centers->size(); j++){
currmsd=polymeres->at(i)->rmsAllc(polymeres->at(centers->at(j)));
if(currmsd<minrmsd){
minrmsd=currmsd;
minclus=j;
}
}
clusters[minclus].push_back(i);
aves[minclus]=aves[minclus]+currmsd;
}

for(int a=0; a<aves.size(); a++){
cur_ave=aves[a]/(1+clusters[a].size())+cur_ave;
}


//Print out the current cluster assignments
for(int x=0; x < clusters.size(); x++){
printf("Cluster: %d N: %d Ave: %8.3f \n",centers->at(x),clusters[x].size(),(aves[x]/(clusters[x].size()+1)));
}
printf("Compared: %8.3f - %8.3f\n",cur_ave,prev_ave);
if(cur_ave==prev_ave){
//Write out id pairs and edge lengths to edge.dat file 
FILE * sfile = fopen("edge.dat","w");
cout.precision(8);
//Foreach line write the coordinates of each c1' atom of the polymere

for(int j = 0; j<centers->size(); j++){
for(int w = 0; w<clusters[j].size(); w++){
float edge=polymeres->at(clusters[j][w])->rmsAllc(polymeres->at(centers->at(j)));
fprintf(sfile,"%d\t%d\t%8.3f\n",clusters[j][w],centers->at(j),edge);
}
}

for(int i = 0; i<centers->size(); i++){
for(int j = 0; j<centers->size(); j++){
float edge=polymeres->at(centers->at(i))->rmsAllc(polymeres->at(centers->at(j)));
fprintf(sfile,"%d\t%d\t%8.3f\n",centers->at(i),centers->at(j),edge);
}
}
fflush( sfile );
fclose(sfile);

break;
} 
cur_ave=0;
prev_ave=0;
//2) Move the cluster center to the center of mass of all assigned to the cluster
//For each cluster find the polymere with the total minimum distance to all others
//This could also be used to maximize the number within 3A
for(int j=0; j<clusters.size(); j++){
int msize=0;
if(clusters[j].size() > max_size){
msize=max_size;
}else{
msize=clusters[j].size();
}
std::vector<float> weights(msize,0);
std::vector<int> ind(msize,0);
//Pick random elements from each group to determine next center
if(clusters[j].size() > max_size){
for(int c1=0; c1<msize; c1++){
ind[c1]=(std::rand() % clusters[j].size());
//printf("Rands %d\n",ind[c1]);
}
}else{
for(int c1=0; c1<msize; c1++){
ind[c1]=c1;
}
}
	for(int i1=0; i1<ind.size(); i1++){
        for(int i2=(i1+1); i2<ind.size(); i2++){
	float radd=polymeres->at(clusters[j][ind[i1]])->rmsAllc(polymeres->at(clusters[j][ind[i2]]));
	weights[i2]=weights[i2]+radd;
	weights[i1]=weights[i1]+radd;
	}
}
float minw=weights[0];
for(int w=0; w<weights.size(); w++){
if(minw > weights[w]){
minw=weights[w];
centers->at(j)=clusters[j][ind[w]];
}
}
}	
}


}

float Polymere::cluster(const std::vector<Polymere *> * polymeres, vector<int>  * centers){
//Returns the index of the centers of kmeans clustering
//Use this for less than 1,000 structures 
//Number of clusters is taken from size of centers
//OpenCL tutorial: sdg0919@gmail.com


//Iterate through the process until all centers converge 
//or you hit the maximum interation limit

//1) Identify which cluster each polymere belongs too
//array of all pairwise rmsds
int size1 = polymeres->size();
if(size1 > 100){size1=100;}
int tot = size1*size1;
vector<float> pwrmsd;
pwrmsd.reserve(10000);
for(int i = 0; i < tot; i++){
pwrmsd.push_back(0);
}
int size2 = (*centers).size();
for(int j = 0; j<size1; j++){
for(int i = 0; i<size1; i++){
float r=polymeres->at(i)->rmsAll(polymeres->at(j));
if(r > 100){
cout << "Bad RMSD: " << i << ":" << j << endl;
//polymeres[i]->printMols(0);
//polymeres[j]->printMols(0);
}
pwrmsd[(j*size1+i)]=r;
}
}

//Determine which clusters the structures are in
int w=0;
float minrmsd=100;
int minclus;
vector<int> clus;
clus.reserve(size1);
//Identify structures in each cluster
int attempts=0;
int diff = 1;

while(diff!=0 && attempts < 100){
attempts++;
for(int i = 0; i < size1; i++){
minrmsd=100;
for(int j=0; j<size2; j++){
if(minrmsd > pwrmsd[((i*size1)+(*centers)[j])]){

minrmsd = pwrmsd[(i*size1+(*centers)[j])];
minclus=j;

}
}
if(attempts == 1){clus.push_back(minclus);}else{
clus[i]=minclus;
}
}

//2) Caclulate new centers by finding center of each cluster
//Find the structure with the minimal average pairwise distance in the cluster center
//Average RMSD
//Determine minimal rmsd between any of the clusters
vector<float> clusave;
for(int c=0; c<size2; c++){clusave.push_back(9999);}
vector <int> clusind;
for(int i=0; i < size2; i++){clusind.push_back(0);}
for(int c=0; c<size1; c++){
float tot=0;
float csize1=0;
//Determine the element in the cluster with the minimal average rmsd to all others in the cluster
//loop through clus
for(int j=0; j<size1; j++){
if(clus[j]==clus[c]){
tot=tot+pwrmsd[c*size1+j];
csize1++;
}
}
//Average PW distance for c in its cluster
tot=tot/csize1;
//cout << "tot: " << tot << endl;
//If its ave is the minimum set clusave and clusind
if(tot < clusave[clus[c]]){ clusave[clus[c]]=tot;clusind[clus[c]]=c; 
//cout << "New Ave: " << clus[c] << ":" << c << endl; 
}
//Calculate the minimum 
}
//Center has the minimum average pairwise rmsd
//Compare new clusind to old clusind 
diff=0;
for(int jt=0; jt< size2; jt++){if(clusind[jt]!=(*centers)[jt]){diff=1;
}}
cout << "Centers: ";
for(int jt=0; jt<size2; jt++){
(*centers)[jt]=clusind[jt];
cout << " " << clusind[jt];
}
cout << endl;
}

float mrmsd1=9999;
for(int j=0; j< size2; j++){
for(int i=(j+1); i<size2; i++){
int ind1=(*centers)[i]*size1 + (*centers)[j];
float mrmsd=pwrmsd[ind1];
if(mrmsd<mrmsd1){mrmsd1=mrmsd;}
//cout << "Cluster: " << clus[j] << endl;
}
} 
pwrmsd.clear();
return mrmsd1;
}

double Polymere::gaussDbl_ij(const vector<double> * vec,int i,int j,double a1,double a2,double a3,double a4,double a5){
//Dataset
//x,y,z,theta2,phi2,type
//Given a dataset of centers and a,b,c for the distrubtion it fits the given pair
//a = center taken from dataset
//b = width  
//c = height high/med/low
double scoreO = 0;
vector<double> params;
int pts=vec->size();
double x,y,z,phi2,theta2,score1,score2,score3,score4,score5;
for(int pi=0; pi<10; pi++){params.push_back(0);}
getParams(i,j,&params);
for(int g=0; g < pts;g=g+6){
//Check to make sure the types are the same
//cout << "Params: " << params[9] << " Vector: " << vec->at((g+5)) << endl;
if(params[9]==vec->at((g+5))){
y=params[6];
x=params[7];
z=params[8];
theta2=params[5];
phi2=params[4];
score1=pow((x-vec->at(g)),2);
score2=pow((y-vec->at(g+1)),2);
score3=pow((z-vec->at(g+2)),2);
score4=pow((theta2-vec->at(g+3)),2);
score5=pow((phi2-vec->at(g+4)),2);
//cout << scoreO << ":" << score1 << "," << score2 << "," << score3 << "," << score4 << "," << score5 << endl;
//cout << x << "," << y << "," << z << "," << theta2 << "," << phi2 << endl;
scoreO=5*exp(-1*(score1/a1+score2/a2+score4/a4+score3/a3+score5/a5))+scoreO;
}
}
scoreO = -1*scoreO;
return scoreO;
}

double Polymere::gauss(const vector<double> * vec){
//Dataset
//Nt1(A/G/U/C=1/2/3/4),Nt2,rho,theta1,phi1,theta2,phi2
//Since theta1==theta2 for almost all cases reduce to one variable
//Given a dataset of centers and a,b,c for the distrubtion it fits the given pair
//a = center taken from dataset
//b = width  
//c = height high/med/low
vector<double> params;
int pts=vec->size();
double score = 0;
double score1=0;
double score2=0;
double score3=0;
double score4=0;
double score5=0;
double z=10;
double mdist=50;
double x,y,x1,y1;
for(int i=0; i<6; i++){params.push_back(0);}
//For each base-pair
for(int i=0; i<num_mols; i++){
for(int j=(i+3); j<num_mols; j++){
//Only Consider Pairs that are in the movement groups
//cout << "Molsi: " << i << ":" << mols[i]->group << " Molsj: " << j << ":" << mols[j]->group << endl; 
if((mols[i]->group==1 || mols[i]->group==2) && (mols[j]->group==1 || mols[j]->group==2)){ 
mdist =pow((mols[i]->x - mols[j]->x),2);
mdist =mdist + pow((mols[i]->y -mols[j]->y),2);
mdist =mdist + pow((mols[i]->z -mols[j]->z),2);
mdist = sqrt(mdist);
//cout << "Type: " << params[0] << ":" << i << "," << j << " i: " << mols[i]->seq << " j: " <<  mols[j]->seq << " mdist: " << mdist << endl;
//cout << "mdist: " << mdist << endl;
if(mdist< 20){
getParams(i,j,&params);
//cout << params[0] << ":" << params[1] << ":" << params[2] << ":" << params[3] << ":" << params[4] << endl;
//Calculate score
for(int g=0; g < pts;g=g+6){
//cout << "g: " << g << "-V:" << vec->at(g) << ":" << params[0] << endl;
if(params[0]==vec->at(g)){
//cout << "rho: " << params[1] << " rscore: " << 5*exp(-1*(pow((vec->at(g+1)-params[1]),2)/0.2)) << endl;
//cout << params[2] << ":" << vec->at(g+2) <<endl;
/*
score1=5*exp(-1*(pow((vec->at(g+1)-params[1]),2)/0.5));
score2=5*exp(-1*(pow((vec->at(g+3)-params[3]),2)/0.005));
score3=5*exp(-1*(pow((vec->at(g+4)-params[4]),2)/0.005));
score4=5*exp(-1*(pow(params[2],2)/0.02));
*/
//type,rho,theta,phi1,phi2
y=sin(params[3])*params[1];
x=cos(params[3])*params[1];
y1=sin(vec->at(g+3))*vec->at(g+1);
x1=cos(vec->at(g+3))*vec->at(g+1);
/*
score1=pow((vec->at(g+1)-params[1]),2);
score2=pow((vec->at(g+3)-params[3]),2);
score3=pow((vec->at(g+4)-params[4]),2);
*/
score1=pow((x-x1),2);
score2=pow((y-y1),2);
score3=pow((vec->at(g+4)-params[4]),2);
score4=pow((params[2]-1.57),2);
score5=pow((params[5]-1.57),2);
z=abs(params[1]*cos((params[2])));
//cout << "Type: " << params[0] << ":" << i << "," << j << ":" << vec->at(g) << " i: " << mols[i]->seq << " j: " << mols[j]->seq << " x: " << x << " y: " << y << " x1: " << x1 << " y1: " << y1 << endl; 
//cout << "p1: " << params[0] << " p2: " << params[1] << " p3: " << params[2] << " p4: " << params[3] << " p5: " << params[4] << endl;
//cout << "v1: " << vec->at(g) << " v2: " << vec->at(g+1) << " v3: " << vec->at(g+2) << " v4: " << vec->at(g+3) << " v5: " << vec->at(g+4) << endl;
//cout << "s2: " << score1 << " s3: " << score4 << " s4: " << score2 << " s5: " << score3 << endl;
if(z<1.5){
//cout << "i: " << i << " j: " << j << " nt: " << mols[i]->seq << " nt: " << mols[j]->seq << endl;
score=5*exp(-1*(score1/0.2+score2/0.2+score3/0.02+score5/0.02))+score;}
}
}

   }

}

}
}
score = -1*score;
//cout<<"Scores: " << score1 << ":" << score2 << ":" << score3 << ":" << score4 << endl;
return score;
}

double Polymere::scoref(){
//Writes Distances used to generate scores 
double distAA=0;
double distAC=0;
double distAU=0;
double distAG=0;
double distCC=0;
double distCU=0;
double distCG=0;
double distGG=0;
double distGU=0;
double distUU=0;
int numAA=0;
int numAC=0;
int numAU=0;
int numAG=0;
int numCC=0;
int numCU=0;
int numCG=0;
int numGG=0;
int numGU=0;
int numUU=0;
double mdist=0;
double score = 0;
for(int i=0; i<num_mols; i++){
for(int j=(i+3); j<num_mols; j++){
//Only Consider Pairs that are in the movement groups
if((mols[i]->group==1 || mols[i]->group==2) && (mols[j]->group==1 || mols[j]->group==2)){ 
mdist =pow((mols[i]->x - mols[j]->x),2);
mdist =mdist + pow((mols[i]->y -mols[j]->y),2);
mdist =mdist + pow((mols[i]->z -mols[j]->z),2);
mdist = sqrt(mdist);
if(mdist>20){
}else{
//Cacluate angle between planes

double Ax1=mols[i]->x - mols[i]->bx;
double Ay1=mols[i]->y - mols[i]->by;
double Az1=mols[i]->z - mols[i]->bz;

double Bx1=mols[i]->b2x - mols[i]->bx;
double By1=mols[i]->b2y - mols[i]->by;
double Bz1=mols[i]->b2z - mols[i]->bz;

double Ax2=mols[j]->x - mols[j]->bx;
double Ay2=mols[j]->y - mols[j]->by;
double Az2=mols[j]->z - mols[j]->bz;

double Bx2=mols[j]->b2x - mols[j]->bx;
double By2=mols[j]->b2y - mols[j]->by;
double Bz2=mols[j]->b2z - mols[j]->bz;

double Cx1=Ay1*Bz1 - Az1*By1;
double Cy1=Az1*Bx1 - Ax1*Bz1;
double Cz1=Ax1*By1 - Ay1*Bx1;

//Making orthonormal set
double Ex1=Cy1*Az1 - Cz1*Ay1;
double Ey1=Cz1*Ax1 - Cx1*Az1;
double Ez1=Cx1*Ay1 - Cy1*Ax1;

double Cx2 =Ay2*Bz2 - Az2*By2;
double Cy2 =Az2*Bx2 - Ax2*Bz2;
double Cz2 =Ax2*By2 - Ay2*Bx2;

double Ex2=Cy2*Az2 - Cz2*Ay2;
double Ey2=Cz2*Ax2 - Cx2*Az2;
double Ez2=Cx2*Ay2 - Cy2*Ax2;

//A,E,C is a coordinate frame
//Determine theta1,phi1,rho,theta2,phi2 
//plug into energy function
//Caclulate angle between normal vectors
double CC1=sqrt(pow(Cx1,2)+pow(Cy1,2)+pow(Cz1,2));
double CC2=sqrt(pow(Cx2,2)+pow(Cy2,2)+pow(Cz2,2));
double AA1=sqrt(pow(Ax1,2)+pow(Ay1,2)+pow(Az1,2));
double AA2=sqrt(pow(Ax2,2)+pow(Ay2,2)+pow(Az2,2));
double EE1=sqrt(pow(Ex1,2)+pow(Ey1,2)+pow(Ez1,2));
double EE2=sqrt(pow(Ex2,2)+pow(Ey2,2)+pow(Ez2,2));

//Create orthonormal set
Cx1=Cx1/CC1;
Cy1=Cy1/CC1;
Cz1=Cz1/CC1;

Ax1=Ax1/AA1;
Ay1=Ay1/AA1;
Az1=Az1/AA1;

Ex1=Ex1/EE1;
Ey1=Ey1/EE1;
Ez1=Ez1/EE1;

Cx2=Cx2/CC2;
Cy2=Cy2/CC2;
Cz2=Cz2/CC2;

Ax2=Ax2/AA2;
Ay2=Ay2/AA2;
Az2=Az2/AA2;
       
Ex2=Ex2/EE2;
Ey2=Ey2/EE2;
Ez2=Ez2/EE2;

//Caclulate rho 
double Dx=mols[i]->bx - mols[j]->bx;
double Dy=mols[i]->by - mols[j]->by;
double Dz=mols[i]->bz - mols[j]->bz;

double rho = sqrt(pow(Dx,2)+pow(Dy,2)+pow(Dz,2));
//Transform D into A,B,C coordinate system
double dx_i=(Dx*Ax1+Dy*Ay1+Dz*Az1);
double dy_i=(Dx*Ex1+Dy*Ey1+Dz*Ez1);
double dz_i=(Dx*Cx1+Dy*Cy1+Dz*Cz1);

double dx_j=-1*(Dx*Ax2+Dy*Ay2+Dz*Az2);
double dy_j=-1*(Dx*Ex2+Dy*Ey2+Dz*Ez2);
double dz_j=-1*(Dx*Cx2+Dy*Cy2+Dz*Cz2);

//Calculate theta
double theta1=acos(dz_i/rho);
//Calculate phi
double s = sqrt(pow(dx_i,2)+pow(dy_i,2));
double phi1 = atan2(dy_i,dx_i);

//Calculate theta
double theta2=acos(dz_i/rho);
//Calculate phi
double s2 = sqrt(pow(dx_j,2)+pow(dy_j,2));
double phi2 = atan2(dy_j,dx_j);
//cout << "rho: " << rho << " theta1: " << theta1 << " phi1: " << phi1 << " theta2: " << theta2 << " phi2: " << phi2 << endl; 

//Use parameters as inputs for scoring
double dprod=(Cx1*Cx2+Cy1*Cy2+Cz2*Cz1);
double ang_score = -1*pow(abs(dprod),20);
//cout << "ang_score: " << ang_score << endl;
//Projection of D onto C1
double proj=abs((Dx*Cx1+Dy*Cy1+Dz*Cz1));
double proj_score = -1/(1+pow(proj,2));  
//cout << "projdist: " << proj_score << endl;
double bp_score=0;
if((mols[i]->seq == 'A' && mols[j]->seq=='A') || (mols[i]->seq == 'A' && mols[j]->seq=='A')  ){  distAA=mdist+distAA;numAA++;}else{
if((mols[i]->seq == 'A' && mols[j]->seq=='U') || (mols[i]->seq == 'U' && mols[j]->seq=='A')  ){  distAU=mdist+distAU;numAU++;bp_score=-2;}else{
if((mols[i]->seq == 'A' && mols[j]->seq=='C') || (mols[i]->seq == 'C' && mols[j]->seq=='A')  ){  distAC=mdist+distAC;numAC++;}else{
if((mols[i]->seq == 'A' && mols[j]->seq=='G') || (mols[i]->seq == 'G' && mols[j]->seq=='A')  ){  distAG=mdist+distAG;numAG++;}else{
if((mols[i]->seq == 'U' && mols[j]->seq=='U') || (mols[i]->seq == 'U' && mols[j]->seq=='U')  ){  distUU=mdist+distUU;numUU++;}else{
if((mols[i]->seq == 'U' && mols[j]->seq=='C') || (mols[i]->seq == 'C' && mols[j]->seq=='U')  ){  distCU=mdist+distCU;numCU++;}else{
if((mols[i]->seq == 'U' && mols[j]->seq=='G') || (mols[i]->seq == 'G' && mols[j]->seq=='U')  ){  distGU=mdist+distGU;numGU++;bp_score=-1;}else{
if((mols[i]->seq == 'C' && mols[j]->seq=='C') || (mols[i]->seq == 'C' && mols[j]->seq=='C')  ){  distCC=mdist+distCC;numCC++;}else{
if((mols[i]->seq == 'C' && mols[j]->seq=='G') || (mols[i]->seq == 'C' && mols[j]->seq=='G')  ){  distCG=mdist+distCG;numCG++;bp_score=-2;}else{
if((mols[i]->seq == 'G' && mols[j]->seq=='G') || (mols[i]->seq == 'G' && mols[j]->seq=='G')  ){  distGG=mdist+distGG;numGG++;}else{
}
}
}
}
}
}
}
}
}
}
score=score+2*ang_score+2*proj_score+bp_score;
//score=score+2*proj_score;
}

}

}
}
//Foreach line write the coordinates of each c1' atom of the polymere
return score;
}
void Polymere::writeDists(string file, double id){
//Writes Distances used to generate scores 
ofstream SaveFile(file.c_str(),ios::app);
cout.precision(8);
double distAA=0;
double distAC=0;
double distAU=0;
double distAG=0;
double distCC=0;
double distCU=0;
double distCG=0;
double distGG=0;
double distGU=0;
double distUU=0;
int numAA=0;
int numAC=0;
int numAU=0;
int numAG=0;
int numCC=0;
int numCU=0;
int numCG=0;
int numGG=0;
int numGU=0;
int numUU=0;
double mdist=0;
double score = 0;
for(int i=0; i<num_mols; i++){
for(int j=(i+1); j<num_mols; j++){
//Only Consider Pairs that are in the movement groups
if((mols[i]->group==1 || mols[i]->group==2) && (mols[j]->group==1 || mols[j]->group==2)){ 
mdist =pow((mols[i]->x - mols[j]->x),2);
mdist =mdist + pow((mols[i]->y -mols[j]->y),2);
mdist =mdist + pow((mols[i]->z -mols[j]->z),2);
mdist = sqrt(mdist);
score=score+abs(mdist-10);
if((mols[i]->seq == 'A' && mols[j]->seq=='A') || (mols[i]->seq == 'A' && mols[j]->seq=='A')  ){  distAA=mdist+distAA;numAA++;}else{
if((mols[i]->seq == 'A' && mols[j]->seq=='U') || (mols[i]->seq == 'U' && mols[j]->seq=='A')  ){  distAU=mdist+distAU;numAU++;}else{
if((mols[i]->seq == 'A' && mols[j]->seq=='C') || (mols[i]->seq == 'C' && mols[j]->seq=='A')  ){  distAC=mdist+distAC;numAC++;}else{
if((mols[i]->seq == 'A' && mols[j]->seq=='G') || (mols[i]->seq == 'G' && mols[j]->seq=='A')  ){  distAG=mdist+distAG;numAG++;}else{
if((mols[i]->seq == 'U' && mols[j]->seq=='U') || (mols[i]->seq == 'U' && mols[j]->seq=='U')  ){  distUU=mdist+distUU;numUU++;}else{
if((mols[i]->seq == 'U' && mols[j]->seq=='C') || (mols[i]->seq == 'C' && mols[j]->seq=='U')  ){  distCU=mdist+distCU;numCU++;}else{
if((mols[i]->seq == 'U' && mols[j]->seq=='G') || (mols[i]->seq == 'G' && mols[j]->seq=='U')  ){  distGU=mdist+distGU;numGU++;}else{
if((mols[i]->seq == 'C' && mols[j]->seq=='C') || (mols[i]->seq == 'C' && mols[j]->seq=='C')  ){  distCC=mdist+distCC;numCC++;}else{
if((mols[i]->seq == 'C' && mols[j]->seq=='G') || (mols[i]->seq == 'C' && mols[j]->seq=='G')  ){  distCG=mdist+distCG;numCG++;}else{
if((mols[i]->seq == 'G' && mols[j]->seq=='G') || (mols[i]->seq == 'G' && mols[j]->seq=='G')  ){  distGG=mdist+distGG;numGG++;}else{
}
}
}
}
}
}
}
}
}
}

}

}
}
//Foreach line write the coordinates of each c1' atom of the polymere

SaveFile <<  distAA << "\t" << distAU << "\t" << distAC << "\t" << distAG << "\t" << distUU << "\t" << distCU << "\t" << distGU << "\t" << distCC << "\t" << distCG << "\t" << distGG << "\t" << numAA << "\t" << numAU << "\t" << numAC << "\t" << numAG << "\t" << numUU << "\t" << numCU << "\t" << numGU << "\t" << numCC << "\t" << numCG << "\t" << numGG<< "\t" << score <<"\t" << id <<"\n";
SaveFile.close();
}

void Polymere::writeStructInfo(string file, double id,double round,double rmsd){
ofstream SaveFile("3dat",ios::app);
cout.precision(8);
//Foreach line write the coordinates of each c1' atom of the polymere
SaveFile <<  id << "\t" << round << "\t" << rmsd << "\n";
SaveFile.close();
}

void Polymere::writeStructInfo2(string file){
FILE * sfile = fopen(file.c_str(),"a");
cout.precision(8);
//Foreach line write the coordinates of each c1' atom of the polymere
fprintf(sfile,"%f\t%i\n",score,randtot);
fflush( sfile );
fclose(sfile);
}

void Polymere::writeStruct(string file, double id1){
//ofstream SaveFile(file,ios::app);
FILE * sfile = fopen(file.c_str(),"a");
cout.width(8);
//Foreach line write the coordinates of each c1' atom of the polymere
for(int w=0; w < num_mols; w++){
double x1 = mols[w]->x;
double y1 = mols[w]->y;
double z1 = mols[w]->z;
double bx1 = mols[w]->bx;
double by1 = mols[w]->by;
double bz1 = mols[w]->bz;
double b2x1 = mols[w]->b2x;
double b2y1 = mols[w]->b2y;
double b2z1 = mols[w]->b2z;
int group1 = mols[w]->group;
char seq1= mols[w]->seq;
int id2 = mols[w]->id;
//SaveFile <<  x1 << "\t" <<  y1 << "\t" <<  z1 << "\t" <<  bx1 << "\t" <<  by1 << "\t" <<  bz1 << "\t" <<  b2x1 << "\t" <<  b2y1 << "\t" <<  b2z1 << "\t"<<  group1 << "\t" <<  seq1 << "\t" << id1 << "\n";
fprintf(sfile,"%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%i\t%c\t%i\n",x1,y1,z1,bx1,by1,bz1,b2x1,b2y1,b2z1,group1,seq1,id2);
fflush( sfile );
}
//SaveFile.close();
fclose(sfile);
}
void Polymere::readMulti(string file,int start, int num){
	//Given a file with 2D coordinates generate from secondary structure it loads positions into polymere
	using namespace std;
	double x,y,z,bx,by,bz,b2x,b2y,b2z;
int group;
char seq;
int id;
	string field;
	std::ifstream input;
	input.open(file.c_str());
	string line;
	int a=0;
int count=0; 
while(getline(input, line) && count < (start+num)) 
{
if(count >= start && count < (start + num)){
//Read in file with each line being x \t y \t z \t group
	istringstream iss(line,istringstream::in);	
	getline(iss,field,'\t');
	x=strtodouble(field);
		getline(iss,field,'\t');
	y=strtodouble(field);
		getline(iss,field,'\t');
	z=strtodouble(field);
	getline(iss,field,'\t');
	bx=strtodouble(field);
		getline(iss,field,'\t');
	by=strtodouble(field);
		getline(iss,field,'\t');
	bz=strtodouble(field);
	getline(iss,field,'\t');
	b2x=strtodouble(field);
		getline(iss,field,'\t');
	b2y=strtodouble(field);
		getline(iss,field,'\t');
	b2z=strtodouble(field);
	getline(iss,field,'\t');
	group=strtoint(field);
		getline(iss,field,'\t');
	seq=strtochar(field);
		getline(iss,field,'\t');
	id=strtoint(field);
//cout << id << ":" << x << "," << y << "," << z << "," << bx << "," << by << "," << bz << "," << b2x << "," << b2y << "," << b2z<< "," << seq << "," <<group <<endl;
mols[a]->x = x;
mols[a]->y = y;
mols[a]->z = z;
mols[a]->bx = bx;
mols[a]->by = by;
mols[a]->bz = bz;
mols[a]->b2x = b2x;
mols[a]->b2y = b2y;
mols[a]->b2z = b2z;
mols[a]->group=group;
mols[a]->seq=seq;
mols[a]->id=id;
a++;
}
count++;
}
//cout << a << ": Molecules read in. " << endl;
num_mols=a;
}

void Polymere::readPDB(string file){
	//Given a file with 2D coordinates generate from secondary structure it loads positions into polymere
	using namespace std;
//Read PDB into PDB class	
Pdb * pdb = new Pdb(file);
//Get the model, chain, and residues
int a=0;
for(int resnums=0;resnums < pdb->models[0]->chains[0]->residues.size(); resnums++){
//Loop through all the resiudes 
string rname=pdb->models[0]->chains[0]->residues[resnums]->res_name;
printf("Resname: '%s'",rname.c_str());
if(rname.compare("  A")==0 || rname.compare("  G")==0){
vector<float> c1 = pdb->models[0]->chains[0]->residues[resnums]->getAtom(" C1'"); 
vector<float> c2 = pdb->models[0]->chains[0]->residues[resnums]->getAtom(" N9 "); 
vector<float> c3 = pdb->models[0]->chains[0]->residues[resnums]->getAtom(" C4 "); 
//Add coordinates 
mols[a]->x=c1[0];
mols[a]->y=c1[1];
mols[a]->z=c1[2];

mols[a]->bx=c2[0];
mols[a]->by=c2[1];
mols[a]->bz=c2[2];

mols[a]->b2x=c3[0];
mols[a]->b2y=c3[1];
mols[a]->b2z=c3[2];
if(rname.compare("  A")==0){
mols[a]->seq='A';
int bpos = a*81;
for(int atn = 0; atn < 26; atn++){
string at_name = pdb->models[0]->chains[0]->residues[resnums]->atomnames2[atn];
c1 = pdb->models[0]->chains[0]->residues[resnums]->getAtom(at_name);
full_atoms[bpos+ 0]=c1[0];
full_atoms[bpos+ 1]=c1[1];
full_atoms[bpos+ 2]=c1[2];
bpos = bpos+3;
}
}else{
if(rname.compare("  G")==0){
int bpos = a*81;
mols[a]->seq='G';
for(int atn = 0; atn < 27; atn++){
string at_name = pdb->models[0]->chains[0]->residues[resnums]->atomnames1[atn];
c1 = pdb->models[0]->chains[0]->residues[resnums]->getAtom(at_name);
full_atoms[bpos+ 0]=c1[0];
full_atoms[bpos+ 1]=c1[1];
full_atoms[bpos+ 2]=c1[2];
bpos = bpos+3;
}


}

}


}else{

if(rname.compare("  U")==0 || rname.compare("  C")==0){
vector<float> c1 = pdb->models[0]->chains[0]->residues[resnums]->getAtom(" C1'"); 
vector<float> c2 = pdb->models[0]->chains[0]->residues[resnums]->getAtom(" N1 "); 
vector<float> c3 = pdb->models[0]->chains[0]->residues[resnums]->getAtom(" C2 "); 

mols[a]->x=c1[0];
mols[a]->y=c1[1];
mols[a]->z=c1[2];

mols[a]->bx=c2[0];
mols[a]->by=c2[1];
mols[a]->bz=c2[2];

mols[a]->b2x=c3[0];
mols[a]->b2y=c3[1];
mols[a]->b2z=c3[2];
//Add atoms to full_atoms
if(rname.compare("  U")==0){
mols[a]->seq='U';
int bpos=a*81;
for(int atn = 0; atn < 23; atn++){
string at_name = pdb->models[0]->chains[0]->residues[resnums]->atomnames3[atn];
c1 = pdb->models[0]->chains[0]->residues[resnums]->getAtom(at_name);
full_atoms[bpos+ 0]=c1[0];
full_atoms[bpos+ 1]=c1[1];
full_atoms[bpos+ 2]=c1[2];
bpos = bpos+3;
}

}else{

if(rname.compare("  C")==0){
mols[a]->seq='C';
int bpos=a*81;
for(int atn = 0; atn < 24; atn++){
string at_name = pdb->models[0]->chains[0]->residues[resnums]->atomnames4[atn];
c1 = pdb->models[0]->chains[0]->residues[resnums]->getAtom(at_name);
full_atoms[bpos+ 0]=c1[0];
full_atoms[bpos+ 1]=c1[1];
full_atoms[bpos+ 2]=c1[2];
bpos = bpos+3;
}


}

}


}else{
//Failed to load residue
printf("Failed to load PDB: residue names dont' match(A C G U)");
//abort();

}
}
//Get the C1', N9 and C4

//Get the C1', N1 and C2

//Foreach residue in the file extract the C1' atoms coordinates and the base cooridnates and add to mols vector
a++;
}
//Loop through all residues and load all coordinates into the full_atoms vector 
num_mols=a;
//Delete PDB
for(int m=0; m<pdb->models.size(); m++){
for(int c=0; c<pdb->models[m]->chains.size(); c++){
for(int r=0; r<pdb->models[m]->chains[c]->residues.size(); r++){
for(int a=0; a<pdb->models[m]->chains[c]->residues[r]->atoms.size(); a++){
delete pdb->models[m]->chains[c]->residues[r]->atoms[a];
}
delete pdb->models[m]->chains[c]->residues[r];
}
delete pdb->models[m]->chains[c];
}
delete pdb->models[m];
} 
////////////////
delete pdb;
}


void Polymere::readStructNew(string file){
	//Given a file with 2D coordinates generate from secondary structure it loads positions into polymere
	using namespace std;
	double x,y,z,bx,by,bz,b2x,b2y,b2z;
int group;
char seq;
int id;
	string field;
	std::ifstream input;
	input.open(file.c_str());
	string line;
	int a=0;
while(getline(input, line)) 
{
//Read in file with each line being x \t y \t z \t group
	istringstream iss(line,istringstream::in);	
	getline(iss,field,'\t');
	x=strtodouble(field);
		getline(iss,field,'\t');
	y=strtodouble(field);
		getline(iss,field,'\t');
	z=strtodouble(field);
	getline(iss,field,'\t');
	bx=strtodouble(field);
		getline(iss,field,'\t');
	by=strtodouble(field);
		getline(iss,field,'\t');
	bz=strtodouble(field);
	getline(iss,field,'\t');
	b2x=strtodouble(field);
		getline(iss,field,'\t');
	b2y=strtodouble(field);
		getline(iss,field,'\t');
	b2z=strtodouble(field);
	getline(iss,field,'\t');
	group=strtoint(field);
		getline(iss,field,'\t');
	seq=strtochar(field);
		getline(iss,field,'\t');
	id=strtoint(field);
//cout << id << ":" << x << "," << y << "," << z << "," << bx << "," << by << "," << bz << "," << b2x << "," << b2y << "," << b2z<< "," << seq << "," <<group <<endl;
mols[a]->x = x;
mols[a]->y = y;
mols[a]->z = z;
mols[a]->bx = bx;
mols[a]->by = by;
mols[a]->bz = bz;
mols[a]->b2x = b2x;
mols[a]->b2y = b2y;
mols[a]->b2z = b2z;
mols[a]->group=group;
mols[a]->seq=seq;
mols[a]->id=id;
a++;
}
num_mols=a;
}

vector<double> Polymere::readDist(string file){
	//Given a file with 2D coordinates generate from secondary structure it loads positions into polymere
	using namespace std;
	double type,rho,theta,phi1,phi2; 
vector<double> params;
	string field;
	std::ifstream input;
	input.open(file.c_str());
	string line;
while(getline(input, line)) 
{
//Read in file with each line being x \t y \t z \t group
	istringstream iss(line,istringstream::in);	
	getline(iss,field,'\t');
	type=strtodouble(field);
		getline(iss,field,'\t');
	rho=strtodouble(field);
		getline(iss,field,'\t');
	theta=strtodouble(field);
	getline(iss,field,'\t');
	phi1=strtodouble(field);
		getline(iss,field,'\t');
	phi2=strtodouble(field);
params.push_back(type);
params.push_back(rho);
params.push_back(theta);
params.push_back(phi1);
params.push_back(phi2);
}
return params;
}
vector<Base *> Polymere::readStruct(string file){
	//Given a file with 2D coordinates generate from secondary structure it loads positions into polymere
	using namespace std;
	double x,y,z,bx,by,bz,b2x,b2y,b2z;
int group;
char seq;
int id;
vector<Base *> ols;
	string field;
	std::ifstream input;
	input.open(file.c_str());
	string line;
	int a=0;
while(getline(input, line)) 
{
//Read in file with each line being x \t y \t z \t group
	istringstream iss(line,istringstream::in);	
	getline(iss,field,'\t');
	x=strtodouble(field);
		getline(iss,field,'\t');
	y=strtodouble(field);
		getline(iss,field,'\t');
	z=strtodouble(field);
	getline(iss,field,'\t');
	bx=strtodouble(field);
		getline(iss,field,'\t');
	by=strtodouble(field);
		getline(iss,field,'\t');
	bz=strtodouble(field);
	getline(iss,field,'\t');
	b2x=strtodouble(field);
		getline(iss,field,'\t');
	b2y=strtodouble(field);
		getline(iss,field,'\t');
	b2z=strtodouble(field);
	getline(iss,field,'\t');
	group=strtoint(field);
		getline(iss,field,'\t');
	seq=strtochar(field);
		getline(iss,field,'\t');
	id=strtoint(field);
//cout << id << ":" << x << "," << y << "," << z << "," << bx << "," << by << "," << bz << "," << b2x << "," << b2y << "," << b2z<< "," << seq << "," <<group <<endl;
ols.push_back(new Base());
ols[a]->x = x;
ols[a]->y = y;
ols[a]->z = z;
ols[a]->bx = bx;
ols[a]->by = by;
ols[a]->bz = bz;
ols[a]->b2x = b2x;
ols[a]->b2y = b2y;
ols[a]->b2z = b2z;
ols[a]->group=group;
ols[a]->seq=seq;
ols[a]->id=id;
a++;
}
num_mols=a;
return ols;
}

Matrix * Polymere::coordSysF(int pos1){
//Given a position in a polymere it identifies the coordinate system relative to pos-1
int p1nn=pos1-3;
int p1n=pos1-2;
int p10=pos1-1;
int p11=pos1;
int p12=pos1+1;
int p13=pos1+2;
int p14=pos1+3;

double xp2=mols[p12]->x;
double yp2=mols[p12]->y;
double zp2=mols[p12]->z;
//cout << "X2 " << xp2 << " : " << yp2 << " : " << zp2 << endl; 
double x0=mols[p10]->x;
double y0=mols[p10]->y;
double z0=mols[p10]->z;
//cout << "X0 " << x0 << " : " << y0 << " : " << z0 << endl; 
double xp1=mols[p11]->x;
double yp1=mols[p11]->y;
double zp1=mols[p11]->z;
//cout << "X1 " << xp1 << " : " << yp1 << " : " << zp1 << endl; 

double a1=xp2-xp1;
double a2=yp2-yp1;
double a3=zp2-zp1;

//cout << "A: " << a1 << " : " << a2 << " : " << a3 << endl;
double b1=xp1-x0;
double b2=yp1-y0;
double b3=zp1-z0;

//cout << "B: " << b1 << " : " << b2 << " : " << b3 << endl;

//Coordinates of positive x-axis are b1...
//Z-axis is cross product
double z1=b2*a3-b3*a2;
double z2=b3*a1-b1*a3;
double z3=b1*a2-b2*a1;
//cout << "Z " << z1 << ":" << z2 << ":" << z3 << endl;
//y-axis is cross product zxb
double y1=z2*b3-z3*b2;
double y2=z3*b1-z1*b3;
double y3=z1*b2-z2*b1;
//Create normal vectors from y,z,b
double lb=sqrt(pow(b1,2)+pow(b2,2)+pow(b3,2));
double ly=sqrt(pow(y1,2)+pow(y2,2)+pow(y3,2));
double lz=sqrt(pow(z1,2)+pow(z2,2)+pow(z3,2));
y1=y1/ly;
y2=y2/ly;
y3=y3/ly;

z1=z1/lz;
z2=z2/lz;
z3=z3/lz;

b1=b1/lb;
b2=b2/lb;
b3=b3/lb;
//Create transformation matrix1
Matrix * m1 = new Matrix();
m1->a11=b1;
m1->a21=b2;
m1->a31=b3;

m1->a12=y1;
m1->a22=y2;
m1->a32=y3;

m1->a13=z1;
m1->a23=z2;
m1->a33=z3;
return m1;
}

Matrix * Polymere::coordSysG(int pos2){
int p21=pos2-1;
int p20=pos2+0;
int p2n=pos2+1;

double x21=mols[p21]->x;
double y21=mols[p21]->y;
double z21=mols[p21]->z;

double x20=mols[p20]->x;
double y20=mols[p20]->y;
double z20=mols[p20]->z;

double x2n=mols[p2n]->x;
double y2n=mols[p2n]->y;
double z2n=mols[p2n]->z;

double a1a=x2n-x20;
double a2a=y2n-y20;
double a3a=z2n-z20;

double b1a=x21-x20;
double b2a=y21-y20;
double b3a=z21-z20;

//Coordinates of positive x-axis are b1...
//Z-axis is cross product
double z1a=b2a*a3a-b3a*a2a;
double z2a=b3a*a1a-b1a*a3a;
double z3a=b1a*a2a-b2a*a1a;
//y-axis is cross product zxb
double y1a=z2a*b3a-z3a*b2a;
double y2a=z3a*b1a-z1a*b3a;
double y3a=z1a*b2a-z2a*b1a;
//Create normal vectors from y,z,b
double lba=sqrt(pow(b1a,2)+pow(b2a,2)+pow(b3a,2));
double lya=sqrt(pow(y1a,2)+pow(y2a,2)+pow(y3a,2));
double lza=sqrt(pow(z1a,2)+pow(z2a,2)+pow(z3a,2));
y1a=y1a/lya;
y2a=y2a/lya;
y3a=y3a/lya;

z1a=z1a/lza;
z2a=z2a/lza;
z3a=z3a/lza;

b1a=b1a/lba;
b2a=b2a/lba;
b3a=b3a/lba;
//Create transformation matrix1
Matrix * m2 = new Matrix();
m2->a11=b1a;
m2->a21=b2a;
m2->a31=b3a;

m2->a12=y1a;
m2->a22=y2a;
m2->a32=y3a;

m2->a13=z1a;
m2->a23=z2a;
m2->a33=z3a;
return m2;
}
Matrix * Polymere::coordSysE(int pos2){
int p21=pos2+1;
int p20=pos2+0;
int p2n=pos2-1;

double x21=mols[p21]->x;
double y21=mols[p21]->y;
double z21=mols[p21]->z;

double x20=mols[p20]->x;
double y20=mols[p20]->y;
double z20=mols[p20]->z;

double x2n=mols[p2n]->x;
double y2n=mols[p2n]->y;
double z2n=mols[p2n]->z;

double a1a=x2n-x20;
double a2a=y2n-y20;
double a3a=z2n-z20;

double b1a=x21-x20;
double b2a=y21-y20;
double b3a=z21-z20;

//Coordinates of positive x-axis are b1...
//Z-axis is cross product
double z1a=b2a*a3a-b3a*a2a;
double z2a=b3a*a1a-b1a*a3a;
double z3a=b1a*a2a-b2a*a1a;
//y-axis is cross product zxb
double y1a=z2a*b3a-z3a*b2a;
double y2a=z3a*b1a-z1a*b3a;
double y3a=z1a*b2a-z2a*b1a;
//Create normal vectors from y,z,b
double lba=sqrt(pow(b1a,2)+pow(b2a,2)+pow(b3a,2));
double lya=sqrt(pow(y1a,2)+pow(y2a,2)+pow(y3a,2));
double lza=sqrt(pow(z1a,2)+pow(z2a,2)+pow(z3a,2));
y1a=y1a/lya;
y2a=y2a/lya;
y3a=y3a/lya;

z1a=z1a/lza;
z2a=z2a/lza;
z3a=z3a/lza;

b1a=b1a/lba;
b2a=b2a/lba;
b3a=b3a/lba;
//Create transformation matrix1
Matrix * m2 = new Matrix();
m2->a11=b1a;
m2->a21=b2a;
m2->a31=b3a;

m2->a12=y1a;
m2->a22=y2a;
m2->a32=y3a;

m2->a13=z1a;
m2->a23=z2a;
m2->a33=z3a;
return m2;
}

int Polymere::getResidueStart(int pos){
//Each residue has 14 atoms in the base and 12 atoms in the backbone = 78 doubles
int posout = pos*81;
return posout;
}


void Polymere::coordSysAT(int pos,std::vector<Points> * pts){
//Adds backbone atoms to the RNA molecule
//coordSysNT must be called before this function in order to place the 
//base atoms properly
int p1=pos;
int p2=pos+1;
int p3=pos+2;

//Get C1' from the virtual bond model
double ax=mols[p1]->x;
double ay=mols[p1]->y;
double az=mols[p1]->z;

double bx=mols[p2]->x;
double by=mols[p2]->y;
double bz=mols[p2]->z;

double cx=mols[p3]->x;
double cy=mols[p3]->y;
double cz=mols[p3]->z;

//Create a forth point using the cross product of a-b X c-b + b
double A1=ax-bx;
double A2=ay-by;
double A3=az-bz;

double B1=cx-bx;
double B2=cy-by;
double B3=cz-bz;

double Z1=B2*A3-B3*A2;
double Z2=B3*A1-B1*A3;
double Z3=B1*A2-B2*A1;

//Forth point is cross product of A and B
double dx=Z1+bx;
double dy=Z2+by;
double dz=Z3+bz;

//Determine transformation Matrix 

Matrix * Q = new Matrix();
Q->a11=bx-ax;
Q->a21=cx-ax;
Q->a31=dx-ax;
Q->a12=by-ay;
Q->a22=cy-ay;
Q->a32=dy-ay;
Q->a13=bz-az;
Q->a23=cz-az;
Q->a33=dz-az;

//C1* coordinates from the full atom points
double ax2=pts->at(11).x;
double ay2=pts->at(11).y;
double az2=pts->at(11).z;

double bx2=pts->at(23).x;
double by2=pts->at(23).y;
double bz2=pts->at(23).z;

double cx2=pts->at(35).x;
double cy2=pts->at(35).y;
double cz2=pts->at(35).z;
//create the  point d
//Create a forth point using the cross product of a-b X c-b + b
 A1=ax2-bx2;
 A2=ay2-by2;
 A3=az2-bz2;

 B1=cx2-bx2;
 B2=cy2-by2;
 B3=cz2-bz2;

 Z1=B2*A3-B3*A2;
 Z2=B3*A1-B1*A3;
 Z3=B1*A2-B2*A1;

//Forth point is cross product of A and B
double dx2=Z1+bx2;
double dy2=Z2+by2;
double dz2=Z3+bz2;

 
Matrix * P = new Matrix();
P->a11=bx2-ax2;
P->a21=cx2-ax2;
P->a31=dx2-ax2;
P->a12=by2-ay2;
P->a22=cy2-ay2;
P->a32=dy2-ay2;
P->a13=bz2-az2;
P->a23=cz2-az2;
P->a33=dz2-az2;
/*

       (bx2-ax2) (by2-ay2) (bz2-az2) 
                                     
   Q = (cx2-ax2) (cy2-ay2) (cz2-az2) 
                                     
       (dx2-ax2) (dy2-ay2) (dz2-az2) 
*/
/*
cout << "Matrix P" << endl;
cout << P->a11 << ":" << P->a12<<":" << P->a13 << endl;
cout << P->a21 << ":" << P->a22<<":" << P->a23 << endl;
cout << P->a31 << ":" << P->a32<<":" << P->a33 << endl;

cout << "Inverse P" << endl;
*/
Matrix Ip = P->inverse();
/*
cout << Ip.a11 << ":" << Ip.a12<<":" << Ip.a13 << endl;
cout << Ip.a21 << ":" << Ip.a22<<":" << Ip.a23 << endl;
cout << Ip.a31 << ":" << Ip.a32<<":" << Ip.a33 << endl;

cout << "Matrix Q" << endl;
cout << Q->a11 << ":" << Q->a12<<":" << Q->a13 << endl;
cout << Q->a21 << ":" << Q->a22<<":" << Q->a23 << endl;
cout << Q->a31 << ":" << Q->a32<<":" << Q->a33 << endl;
*/
Matrix * T = new Matrix();

T->a11=Ip.a11*Q->a11+Ip.a12*Q->a21+Ip.a13*Q->a31;
T->a21=Ip.a21*Q->a11+Ip.a22*Q->a21+Ip.a23*Q->a31;
T->a31=Ip.a31*Q->a11+Ip.a32*Q->a21+Ip.a33*Q->a31;

T->a12=Ip.a11*Q->a12+Ip.a12*Q->a22+Ip.a13*Q->a32;
T->a22=Ip.a21*Q->a12+Ip.a22*Q->a22+Ip.a23*Q->a32;
T->a32=Ip.a31*Q->a12+Ip.a32*Q->a22+Ip.a33*Q->a32;

T->a13=Ip.a11*Q->a13+Ip.a12*Q->a23+Ip.a13*Q->a33;
T->a23=Ip.a21*Q->a13+Ip.a22*Q->a23+Ip.a23*Q->a33;
T->a33=Ip.a31*Q->a13+Ip.a32*Q->a23+Ip.a33*Q->a33;
/*
cout << "Matrix T" << endl;
cout << T->a11 << ":" << T->a12<<":" << T->a13 << endl;
cout << T->a21 << ":" << T->a22<<":" << T->a23 << endl;
cout << T->a31 << ":" << T->a32<<":" << T->a33 << endl;
*/
//Determine the translation vector
double ax3= ax-(ax2*T->a11+ay2*T->a21+az2*T->a31);
double ay3= ay-(ax2*T->a12+ay2*T->a22+az2*T->a32);
double az3= az-(ax2*T->a13+ay2*T->a23+az2*T->a33);

/*
The geometric condition on A,B,C,D not lying in a common plane amounts
to the matrix algebra condition that both P and Q are invertible, so:
         M  = P   * Q 
*/
//Convert each point to new coordinate system
//pts[1-3]=NT1-c1'
//pts[4-6]=NT2-c1'
//pts[7-9]=NT3-c1'
//pts[10-12]=NT1b-c
//pts[13-15]=NT2b-c
//pts[16-18]=NT3b-c
//pts[19-21]=NT1b-n
//pts[22-24]=NT2b-n
//pts[25-27]=NT3b-n

//Update points in base
//loop through 12 points at each position
//printf("\nMade transformation matrix, full_atoms size \%d\n",full_atoms.size());
for(int pos_num=0; pos_num < 3; pos_num++){
//Get the starting position of the backbone atoms
int bpos=getResidueStart((pos+pos_num));
for(int ps_num=0; ps_num < 12; ps_num++){
//printf("Loading atoms,%d\n",ps_num);
int bpos1=bpos+3*ps_num;
double x=pts->at(ps_num+(12*pos_num)).x;
double y=pts->at(ps_num+(12*pos_num)).y;
double z=pts->at(ps_num+(12*pos_num)).z;
double w1=x*T->a11+y*T->a21+z*T->a31+ax3;
double w2=x*T->a12+y*T->a22+z*T->a32+ay3;
double w3=x*T->a13+y*T->a23+z*T->a33+az3;
//Add transformed points to full atom polymere
full_atoms[bpos1]=w1;
full_atoms[bpos1+1]=w2;
full_atoms[bpos1+2]=w3;
}
}
//printf("Created backbone atoms\n");
//Update base atoms

//Loop through three positions
for(int pos_num=0; pos_num < 3; pos_num++){
vector <Points> p;
Points c1,b1,b2;
int type1;
c1.x=mols[pos+pos_num]->x;
c1.y=mols[pos+pos_num]->y;
c1.z=mols[pos+pos_num]->z;
b1.x=mols[pos+pos_num]->bx;
b1.y=mols[pos+pos_num]->by;
b1.z=mols[pos+pos_num]->bz;
b2.x=mols[pos+pos_num]->b2x;
b2.y=mols[pos+pos_num]->b2y;
b2.z=mols[pos+pos_num]->b2z;
p.push_back(c1);
p.push_back(b1);
p.push_back(b2);

if(mols[pos+pos_num]->seq=='G'){type1=1;}else{
if(mols[pos+pos_num]->seq=='A'){type1=2;}else{
if(mols[pos+pos_num]->seq=='U'){type1=3;}else{
if(mols[pos+pos_num]->seq=='C'){type1=4;}else{
}
}
}
}
//Create full atom base
FABase * base1 = new FABase(type1);  
//Align to coordinate system
base1->alignBase(p);
//Add to full_atom vector
int bpos=getResidueStart((pos+pos_num));
//Add 36 to bpos for the backbone atoms
bpos=bpos+36;
//Guanine: Total Number of atoms
if(type1==1){
full_atoms[bpos+ 0]=base1->N9.x;
full_atoms[bpos+ 1]=base1->N9.y;
full_atoms[bpos+ 2]=base1->N9.z;

full_atoms[bpos+ 3]=base1->C4.x;
full_atoms[bpos+ 4]=base1->C4.y;
full_atoms[bpos+ 5]=base1->C4.z;

full_atoms[bpos+ 6]=base1->N3.x;
full_atoms[bpos+ 7]=base1->N3.y;
full_atoms[bpos+ 8]=base1->N3.z;

full_atoms[bpos+ 9]=base1->C2.x;
full_atoms[bpos+10]=base1->C2.y;
full_atoms[bpos+11]=base1->C2.z;

full_atoms[bpos+12]=base1->N2.x;
full_atoms[bpos+13]=base1->N2.y;
full_atoms[bpos+14]=base1->N2.z;

full_atoms[bpos+15]=base1->N1.x;
full_atoms[bpos+16]=base1->N1.y;
full_atoms[bpos+17]=base1->N1.z;

full_atoms[bpos+18]=base1->C6.x;
full_atoms[bpos+19]=base1->C6.y;
full_atoms[bpos+20]=base1->C6.z;

full_atoms[bpos+21]=base1->O6.x;
full_atoms[bpos+22]=base1->O6.y;
full_atoms[bpos+23]=base1->O6.z;

full_atoms[bpos+24]=base1->C5.x;
full_atoms[bpos+25]=base1->C5.y;
full_atoms[bpos+26]=base1->C5.z;

full_atoms[bpos+27]=base1->N7.x;
full_atoms[bpos+28]=base1->N7.y;
full_atoms[bpos+29]=base1->N7.z;

full_atoms[bpos+30]=base1->C8.x;
full_atoms[bpos+31]=base1->C8.y;
full_atoms[bpos+32]=base1->C8.z;

full_atoms[bpos+33]=base1->H1.x;
full_atoms[bpos+34]=base1->H1.y;
full_atoms[bpos+35]=base1->H1.z;

full_atoms[bpos+36]=base1->H8.x;
full_atoms[bpos+37]=base1->H8.y;
full_atoms[bpos+38]=base1->H8.z;

full_atoms[bpos+39]=base1->H21.x;
full_atoms[bpos+40]=base1->H21.y;
full_atoms[bpos+41]=base1->H21.z;

full_atoms[bpos+42]=base1->H22.x;
full_atoms[bpos+43]=base1->H22.y;
full_atoms[bpos+44]=base1->H22.z;
}
//Adenine
if(type1==2){
full_atoms[bpos+ 0]=base1->N9.x;
full_atoms[bpos+ 1]=base1->N9.y;
full_atoms[bpos+ 2]=base1->N9.z;

full_atoms[bpos+ 3]=base1->C4.x;
full_atoms[bpos+ 4]=base1->C4.y;
full_atoms[bpos+ 5]=base1->C4.z;

full_atoms[bpos+ 6]=base1->N3.x;
full_atoms[bpos+ 7]=base1->N3.y;
full_atoms[bpos+ 8]=base1->N3.z;

full_atoms[bpos+ 9]=base1->C2.x;
full_atoms[bpos+10]=base1->C2.y;
full_atoms[bpos+11]=base1->C2.z;

full_atoms[bpos+12]=base1->N1.x;
full_atoms[bpos+13]=base1->N1.y;
full_atoms[bpos+14]=base1->N1.z;

full_atoms[bpos+15]=base1->C6.x;
full_atoms[bpos+16]=base1->C6.y;
full_atoms[bpos+17]=base1->C6.z;

full_atoms[bpos+18]=base1->N6.x;
full_atoms[bpos+19]=base1->N6.y;
full_atoms[bpos+20]=base1->N6.z;

full_atoms[bpos+21]=base1->C5.x;
full_atoms[bpos+22]=base1->C5.y;
full_atoms[bpos+23]=base1->C5.z;

full_atoms[bpos+24]=base1->N7.x;
full_atoms[bpos+25]=base1->N7.y;
full_atoms[bpos+26]=base1->N7.z;

full_atoms[bpos+27]=base1->C8.x;
full_atoms[bpos+28]=base1->C8.y;
full_atoms[bpos+29]=base1->C8.z;

full_atoms[bpos+30]=base1->H2.x;
full_atoms[bpos+31]=base1->H2.y;
full_atoms[bpos+32]=base1->H2.z;

full_atoms[bpos+33]=base1->H8.x;
full_atoms[bpos+34]=base1->H8.y;
full_atoms[bpos+35]=base1->H8.z;

full_atoms[bpos+36]=base1->H61.x;
full_atoms[bpos+37]=base1->H61.y;
full_atoms[bpos+38]=base1->H61.z;

full_atoms[bpos+39]=base1->H62.x;
full_atoms[bpos+40]=base1->H62.y;
full_atoms[bpos+41]=base1->H62.z;
}
//Uracil
if(type1==3){
full_atoms[bpos+ 0]=base1->N1.x;
full_atoms[bpos+ 1]=base1->N1.y;
full_atoms[bpos+ 2]=base1->N1.z;

full_atoms[bpos+ 3]=base1->C2.x;
full_atoms[bpos+ 4]=base1->C2.y;
full_atoms[bpos+ 5]=base1->C2.z;

full_atoms[bpos+ 6]=base1->C6.x;
full_atoms[bpos+ 7]=base1->C6.y;
full_atoms[bpos+ 8]=base1->C6.z;

full_atoms[bpos+ 9]=base1->O2.x;
full_atoms[bpos+10]=base1->O2.y;
full_atoms[bpos+11]=base1->O2.z;

full_atoms[bpos+12]=base1->N3.x;
full_atoms[bpos+13]=base1->N3.y;
full_atoms[bpos+14]=base1->N3.z;

full_atoms[bpos+15]=base1->C4.x;
full_atoms[bpos+16]=base1->C4.y;
full_atoms[bpos+17]=base1->C4.z;

full_atoms[bpos+18]=base1->O4.x;
full_atoms[bpos+19]=base1->O4.y;
full_atoms[bpos+20]=base1->O4.z;

full_atoms[bpos+21]=base1->C5.x;
full_atoms[bpos+22]=base1->C5.y;
full_atoms[bpos+23]=base1->C5.z;

full_atoms[bpos+24]=base1->H3.x;
full_atoms[bpos+25]=base1->H3.y;
full_atoms[bpos+26]=base1->H3.z;

full_atoms[bpos+27]=base1->H5.x;
full_atoms[bpos+28]=base1->H5.y;
full_atoms[bpos+29]=base1->H5.z;

full_atoms[bpos+30]=base1->H6.x;
full_atoms[bpos+31]=base1->H6.y;
full_atoms[bpos+32]=base1->H6.z;
}

if(type1==4){
//Cystosine
full_atoms[bpos+ 0]=base1->N1.x;
full_atoms[bpos+ 1]=base1->N1.y;
full_atoms[bpos+ 2]=base1->N1.z;

full_atoms[bpos+ 3]=base1->C2.x;
full_atoms[bpos+ 4]=base1->C2.y;
full_atoms[bpos+ 5]=base1->C2.z;

full_atoms[bpos+ 6]=base1->C6.x;
full_atoms[bpos+ 7]=base1->C6.y;
full_atoms[bpos+ 8]=base1->C6.z;

full_atoms[bpos+ 9]=base1->O2.x;
full_atoms[bpos+10]=base1->O2.y;
full_atoms[bpos+11]=base1->O2.z;

full_atoms[bpos+12]=base1->N3.x;
full_atoms[bpos+13]=base1->N3.y;
full_atoms[bpos+14]=base1->N3.z;

full_atoms[bpos+15]=base1->C4.x;
full_atoms[bpos+16]=base1->C4.y;
full_atoms[bpos+17]=base1->C4.z;

full_atoms[bpos+18]=base1->N4.x;
full_atoms[bpos+19]=base1->N4.y;
full_atoms[bpos+20]=base1->N4.z;

full_atoms[bpos+21]=base1->C5.x;
full_atoms[bpos+22]=base1->C5.y;
full_atoms[bpos+23]=base1->C5.z;

full_atoms[bpos+24]=base1->H5.x;
full_atoms[bpos+25]=base1->H5.y;
full_atoms[bpos+26]=base1->H5.z;

full_atoms[bpos+27]=base1->H6.x;
full_atoms[bpos+28]=base1->H6.y;
full_atoms[bpos+29]=base1->H6.z;

full_atoms[bpos+30]=base1->H41.x;
full_atoms[bpos+31]=base1->H41.y;
full_atoms[bpos+32]=base1->H41.z;

full_atoms[bpos+33]=base1->H42.x;
full_atoms[bpos+34]=base1->H42.y;
full_atoms[bpos+35]=base1->H42.z;
}
delete base1;
}

delete Q;
delete P;
delete T;
}


void Polymere::coordSysNt(int pos,std::vector<Points> * pts){
//Updates polymere with atoms describing base configuration
//pts[1-3]=NT1-c1'
//pts[4-6]=NT2-c1'
//pts[7-9]=NT3-c1'
//pts[10-12]=NT1b-c
//pts[13-15]=NT2b-c
//pts[16-18]=NT3b-c
//pts[19-21]=NT1b-n
//nts[22-24]=NT2b-n
//pts[25-27]=NT3b-n
//Takes in the NT base coordinates from the fragment database and 
//transforms them into the polymere coordinate system  
int p1=pos;
int p2=pos+1;
int p3=pos+2;

double ax=mols[p1]->x;
double ay=mols[p1]->y;
double az=mols[p1]->z;

double bx=mols[p2]->x;
double by=mols[p2]->y;
double bz=mols[p2]->z;

double cx=mols[p3]->x;
double cy=mols[p3]->y;
double cz=mols[p3]->z;

//Create a forth point using the cross product of a-b X c-b + b
double A1=ax-bx;
double A2=ay-by;
double A3=az-bz;

double B1=cx-bx;
double B2=cy-by;
double B3=cz-bz;

double Z1=B2*A3-B3*A2;
double Z2=B3*A1-B1*A3;
double Z3=B1*A2-B2*A1;

//Forth point is cross product of A and B
double dx=Z1+bx;
double dy=Z2+by;
double dz=Z3+bz;

//Determine transformation Matrix 

Matrix * Q = new Matrix();
Q->a11=bx-ax;
Q->a21=cx-ax;
Q->a31=dx-ax;
Q->a12=by-ay;
Q->a22=cy-ay;
Q->a32=dy-ay;
Q->a13=bz-az;
Q->a23=cz-az;
Q->a33=dz-az;
/*
   (bx-ax1) (by1-ay1) (bz1-az1) 
                                
 P=  (cx1-ax1) (cy1-ay1) (cz1-az1)
                                
   (dx1-ax1) (dy1-ay1) (dz1-az1)
*/

double ax2=pts->at(0).x;
double ay2=pts->at(0).y;
double az2=pts->at(0).z;

double bx2=pts->at(3).x;
double by2=pts->at(3).y;
double bz2=pts->at(3).z;

double cx2=pts->at(6).x;
double cy2=pts->at(6).y;
double cz2=pts->at(6).z;

//create the  point d
//Create a forth point using the cross product of a-b X c-b + b
 A1=ax2-bx2;
 A2=ay2-by2;
 A3=az2-bz2;

 B1=cx2-bx2;
 B2=cy2-by2;
 B3=cz2-bz2;

 Z1=B2*A3-B3*A2;
 Z2=B3*A1-B1*A3;
 Z3=B1*A2-B2*A1;

//Forth point is cross product of A and B
double dx2=Z1+bx2;
double dy2=Z2+by2;
double dz2=Z3+bz2;

 
Matrix * P = new Matrix();
P->a11=bx2-ax2;
P->a21=cx2-ax2;
P->a31=dx2-ax2;
P->a12=by2-ay2;
P->a22=cy2-ay2;
P->a32=dy2-ay2;
P->a13=bz2-az2;
P->a23=cz2-az2;
P->a33=dz2-az2;
/*

       (bx2-ax2) (by2-ay2) (bz2-az2) 
                                     
   Q = (cx2-ax2) (cy2-ay2) (cz2-az2) 
                                     
       (dx2-ax2) (dy2-ay2) (dz2-az2) 
*/
/*
cout << "Matrix P" << endl;
cout << P->a11 << ":" << P->a12<<":" << P->a13 << endl;
cout << P->a21 << ":" << P->a22<<":" << P->a23 << endl;
cout << P->a31 << ":" << P->a32<<":" << P->a33 << endl;

cout << "Inverse P" << endl;
*/
Matrix Ip = P->inverse();
/*
cout << Ip.a11 << ":" << Ip.a12<<":" << Ip.a13 << endl;
cout << Ip.a21 << ":" << Ip.a22<<":" << Ip.a23 << endl;
cout << Ip.a31 << ":" << Ip.a32<<":" << Ip.a33 << endl;

cout << "Matrix Q" << endl;
cout << Q->a11 << ":" << Q->a12<<":" << Q->a13 << endl;
cout << Q->a21 << ":" << Q->a22<<":" << Q->a23 << endl;
cout << Q->a31 << ":" << Q->a32<<":" << Q->a33 << endl;
*/
Matrix * T = new Matrix();

T->a11=Ip.a11*Q->a11+Ip.a12*Q->a21+Ip.a13*Q->a31;
T->a21=Ip.a21*Q->a11+Ip.a22*Q->a21+Ip.a23*Q->a31;
T->a31=Ip.a31*Q->a11+Ip.a32*Q->a21+Ip.a33*Q->a31;

T->a12=Ip.a11*Q->a12+Ip.a12*Q->a22+Ip.a13*Q->a32;
T->a22=Ip.a21*Q->a12+Ip.a22*Q->a22+Ip.a23*Q->a32;
T->a32=Ip.a31*Q->a12+Ip.a32*Q->a22+Ip.a33*Q->a32;

T->a13=Ip.a11*Q->a13+Ip.a12*Q->a23+Ip.a13*Q->a33;
T->a23=Ip.a21*Q->a13+Ip.a22*Q->a23+Ip.a23*Q->a33;
T->a33=Ip.a31*Q->a13+Ip.a32*Q->a23+Ip.a33*Q->a33;
/*
cout << "Matrix T" << endl;
cout << T->a11 << ":" << T->a12<<":" << T->a13 << endl;
cout << T->a21 << ":" << T->a22<<":" << T->a23 << endl;
cout << T->a31 << ":" << T->a32<<":" << T->a33 << endl;
*/
//Determine the translation vector
double ax3= ax-(ax2*T->a11+ay2*T->a21+az2*T->a31);
double ay3= ay-(ax2*T->a12+ay2*T->a22+az2*T->a32);
double az3= az-(ax2*T->a13+ay2*T->a23+az2*T->a33);

/*
The geometric condition on A,B,C,D not lying in a common plane amounts
to the matrix algebra condition that both P and Q are invertible, so:
         M  = P   * Q 
*/
//Convert each point to new coordinate system
//pts[1-3]=NT1-c1'
//pts[4-6]=NT2-c1'
//pts[7-9]=NT3-c1'
//pts[10-12]=NT1b-c
//pts[13-15]=NT2b-c
//pts[16-18]=NT3b-c
//pts[19-21]=NT1b-n
//pts[22-24]=NT2b-n
//pts[25-27]=NT3b-n

//Update points in base
double x=pts->at(0).x;
double y=pts->at(0).y;
double z=pts->at(0).z;
double w1=x*T->a11+y*T->a21+z*T->a31+ax3;
double w2=x*T->a12+y*T->a22+z*T->a32+ay3;
double w3=x*T->a13+y*T->a23+z*T->a33+az3;
mols[pos]->x=w1;
mols[pos]->y=w2;
mols[pos]->z=w3;

x=pts->at(1).x;
y=pts->at(1).y;
z=pts->at(1).z;
w1=x*T->a11+y*T->a21+z*T->a31+ax3;
w2=x*T->a12+y*T->a22+z*T->a32+ay3;
w3=x*T->a13+y*T->a23+z*T->a33+az3;
mols[pos]->bx=w1;
mols[pos]->by=w2;
mols[pos]->bz=w3;


x=pts->at(2).x;
y=pts->at(2).y;
z=pts->at(2).z;
w1=x*T->a11+y*T->a21+z*T->a31+ax3;
w2=x*T->a12+y*T->a22+z*T->a32+ay3;
w3=x*T->a13+y*T->a23+z*T->a33+az3;
mols[pos]->b2x=w1;
mols[pos]->b2y=w2;
mols[pos]->b2z=w3;

//Update points in base
pos=pos+1;
x=pts->at(3).x;
y=pts->at(3).y;
z=pts->at(3).z;
w1=x*T->a11+y*T->a21+z*T->a31+ax3;
w2=x*T->a12+y*T->a22+z*T->a32+ay3;
w3=x*T->a13+y*T->a23+z*T->a33+az3;
mols[pos]->x=w1;
mols[pos]->y=w2;
mols[pos]->z=w3;

x=pts->at(4).x;
y=pts->at(4).y;
z=pts->at(4).z;
w1=x*T->a11+y*T->a21+z*T->a31+ax3;
w2=x*T->a12+y*T->a22+z*T->a32+ay3;
w3=x*T->a13+y*T->a23+z*T->a33+az3;
mols[pos]->bx=w1;
mols[pos]->by=w2;
mols[pos]->bz=w3;


x=pts->at(5).x;
y=pts->at(5).y;
z=pts->at(5).z;
w1=x*T->a11+y*T->a21+z*T->a31+ax3;
w2=x*T->a12+y*T->a22+z*T->a32+ay3;
w3=x*T->a13+y*T->a23+z*T->a33+az3;
mols[pos]->b2x=w1;
mols[pos]->b2y=w2;
mols[pos]->b2z=w3;

//Update points in base
pos=pos+1;
x=pts->at(6).x;
y=pts->at(6).y;
z=pts->at(6).z;
w1=x*T->a11+y*T->a21+z*T->a31+ax3;
w2=x*T->a12+y*T->a22+z*T->a32+ay3;
w3=x*T->a13+y*T->a23+z*T->a33+az3;
mols[pos]->x=w1;
mols[pos]->y=w2;
mols[pos]->z=w3;

x=pts->at(7).x;
y=pts->at(7).y;
z=pts->at(7).z;
w1=x*T->a11+y*T->a21+z*T->a31+ax3;
w2=x*T->a12+y*T->a22+z*T->a32+ay3;
w3=x*T->a13+y*T->a23+z*T->a33+az3;
mols[pos]->bx=w1;
mols[pos]->by=w2;
mols[pos]->bz=w3;


x=pts->at(8).x;
y=pts->at(8).y;
z=pts->at(8).z;
w1=x*T->a11+y*T->a21+z*T->a31+ax3;
w2=x*T->a12+y*T->a22+z*T->a32+ay3;
w3=x*T->a13+y*T->a23+z*T->a33+az3;
mols[pos]->b2x=w1;
mols[pos]->b2y=w2;
mols[pos]->b2z=w3;

delete Q;
delete P;
delete T;
}

Matrix * Polymere::coordSysR(int pos2){
int p2nn=pos2+5;
int p2n=pos2+4;
int p20=pos2+3;
int p21=pos2+2;
int p22=pos2+1;
int p23=pos2+0;
int p24=pos2-1;

double xna=mols[p22]->x;
double yna=mols[p22]->y;
double zna=mols[p22]->z;

double xnna=mols[p21]->x;
double ynna=mols[p21]->y;
double znna=mols[p21]->z;

double x0a=mols[p20]->x;
double y0a=mols[p20]->y;
double z0a=mols[p20]->z;

double x1a=mols[p21]->x;
double yp1a=mols[p21]->y;
double zp1a=mols[p21]->z;

double a1a=xna-xnna;
double a2a=yna-ynna;
double a3a=zna-znna;

double b1a=x1a-x0a;
double b2a=yp1a-y0a;
double b3a=zp1a-z0a;

//Coordinates of positive x-axis are b1...
//Z-axis is cross product
double z1a=b2a*a3a-b3a*a2a;
double z2a=b3a*a1a-b1a*a3a;
double z3a=b1a*a2a-b2a*a1a;
//y-axis is cross product zxb
double y1a=z2a*b3a-z3a*b2a;
double y2a=z3a*b1a-z1a*b3a;
double y3a=z1a*b2a-z2a*b1a;
//Create normal vectors from y,z,b
double lba=sqrt(pow(b1a,2)+pow(b2a,2)+pow(b3a,2));
double lya=sqrt(pow(y1a,2)+pow(y2a,2)+pow(y3a,2));
double lza=sqrt(pow(z1a,2)+pow(z2a,2)+pow(z3a,2));
y1a=y1a/lya;
y2a=y2a/lya;
y3a=y3a/lya;

z1a=z1a/lza;
z2a=z2a/lza;
z3a=z3a/lza;

b1a=b1a/lba;
b2a=b2a/lba;
b3a=b3a/lba;
//Create transformation matrix1
Matrix * m2 = new Matrix();
m2->a11=b1a;
m2->a21=b2a;
m2->a31=b3a;

m2->a12=y1a;
m2->a22=y2a;
m2->a32=y3a;

m2->a13=z1a;
m2->a23=z2a;
m2->a33=z3a;
return m2;
}

void Polymere::chooseAllSingle(int pos1,Fragclass * frags,std::vector<int> * fidlist){
// Given a position if generates a set of moves that alters the given fragments
//Get frag1 type
int p10=pos1-1;
int p11=pos1;
int p12=pos1+1;
int p13=pos1+2;
int p14=pos1+3;
char s11=mols[p11]->seq;
char s12=mols[p12]->seq;
char s13=mols[p13]->seq;
int t1,t2,t3;
if(s11=='A' || s11=='G'){t1=1;}else{t1=0;}	
if(s12=='A' || s12=='G'){t2=1;}else{t2=0;}	
if(s13=='A' || s13=='G'){t3=1;}else{t3=0;}	
int type1=t1+2*t2+4*t3;
//get the fragment class 
if(type1==0){for(int j=0; j<frags->YYY1.size(); j++){if(frags->YYY1[j]->id!=-1){fidlist->push_back(frags->YYY1[j]->id);}}}
if(type1==1){for(int j=0; j<frags->RYY1.size(); j++){if(frags->RYY1[j]->id!=-1){fidlist->push_back(frags->RYY1[j]->id);}}}
if(type1==2){for(int j=0; j<frags->YRY1.size(); j++){if(frags->YRY1[j]->id!=-1){fidlist->push_back(frags->YRY1[j]->id);}}}
if(type1==3){for(int j=0; j<frags->RRY1.size(); j++){if(frags->RRY1[j]->id!=-1){fidlist->push_back(frags->RRY1[j]->id);}}}
if(type1==4){for(int j=0; j<frags->YYR1.size(); j++){if(frags->YYR1[j]->id!=-1){fidlist->push_back(frags->YYR1[j]->id);}}}
if(type1==5){for(int j=0; j<frags->RYR1.size(); j++){if(frags->RYR1[j]->id!=-1){fidlist->push_back(frags->RYR1[j]->id);}}}
if(type1==6){for(int j=0; j<frags->YRR1.size(); j++){if(frags->YRR1[j]->id!=-1){fidlist->push_back(frags->YRR1[j]->id);}}}
if(type1==7){for(int j=0; j<frags->RRR1.size(); j++){if(frags->RRR1[j]->id!=-1){fidlist->push_back(frags->RRR1[j]->id);}}}
printf("Fidsize: %di\n",fidlist->size());
}

int Polymere::chooseSingle(int pos1,Fragclass * frags){
// Given a position if generates a set of moves that alters the given fragments
//Get frag1 type
int p10=pos1-1;
int p11=pos1;
int p12=pos1+1;
int p13=pos1+2;
int p14=pos1+3;
char s11=mols[p11]->seq;
char s12=mols[p12]->seq;
char s13=mols[p13]->seq;
int t1,t2,t3;
if(s11=='A' || s11=='G'){t1=1;}else{t1=0;}	
if(s12=='A' || s12=='G'){t2=1;}else{t2=0;}	
if(s13=='A' || s13=='G'){t3=1;}else{t3=0;}	
int type1=t1+2*t2+4*t3;
int fid=frags->matchedSingle(type1);
return fid;
}

vector<Matches> Polymere::fixed3to5(int pos1, int pos2,double Lmin, double Lmax, Fragclass * frags){
int t1,t2,t3;
int p20=pos2+3;
int p21=pos2+2;
int p22=pos2+1;
int p23=pos2+0;
int p24=pos2-1;
char s21=mols[p23]->seq;
char s22=mols[p22]->seq;
char s23=mols[p21]->seq;
if(s21=='A' || s21=='G'){t1=1;}else{t1=0;}	
if(s22=='A' || s22=='G'){t2=1;}else{t2=0;}	
if(s23=='A' || s23=='G'){t3=1;}else{t3=0;}	
int type2=t1+2*t2+4*t3;
double x0a=mols[p20]->x;
double y0a=mols[p20]->y;
double z0a=mols[p20]->z;
Matrix * m2 = coordSysR(pos2);
//Create Vclass to translate
Vclass * o2=new Vclass();
o2->x=-1*x0a;
o2->y=-1*y0a;
o2->z=-1*z0a;
//Maybe return a vector of matches instead of pointers 
vector<Matches> mat;
mat=frags->matchedFMoves(type2,Lmin,Lmax,m2,o2,mols[pos1]->x,mols[pos1]->y,mols[pos1]->z,1,pos1,pos2);
delete o2;
delete m2;
return mat;
}

vector<Matches> Polymere::fixed5to3(int pos1, int pos2,double Lmin,double Lmax, Fragclass * frags){
//Given two positions pos1 and pos2 
//It identifies 5' -> 3' fragments starting at position 1
//that put the C1' atom of pos+3 between Lmin and Lmax angstroms
//Get frag1 type
int p10=pos1-1;
int p11=pos1;
int p12=pos1+1;
int p13=pos1+2;
int p14=pos1+3;
char s11=mols[p11]->seq;
char s12=mols[p12]->seq;
char s13=mols[p13]->seq;
int t1,t2,t3;
if(s11=='A' || s11=='G'){t1=1;}else{t1=0;}	
if(s12=='A' || s12=='G'){t2=1;}else{t2=0;}	
if(s13=='A' || s13=='G'){t3=1;}else{t3=0;}	
int type1=t1+2*t2+4*t3;

double x0=mols[p10]->x;
double y0=mols[p10]->y;
double z0=mols[p10]->z;
Matrix * m1 = coordSysF(pos1);
//Create Vclass to translate
Vclass * o1=new Vclass();
o1->x=-1*x0;
o1->y=-1*y0;
o1->z=-1*z0;
//Create transformation matrix to map frag vectors to global coordinates
//Positive X-axis p0 to p1
//Maybe return a vector of matches instead of pointers 
vector<Matches> mat;
mat=frags->matchedFMoves(type1,Lmin,Lmax,m1,o1,mols[pos2]->x,mols[pos2]->y,mols[pos2]->z,2,pos1,pos2);
delete o1;
delete m1;
return mat;
}

vector<Matches> Polymere::fixedMoves(int pos1, int pos2,Fragclass * frags,double L){
// Given a position if generates a set of moves that alters the given fragments
//Get frag1 type
int p10=pos1-1;
int p11=pos1;
int p12=pos1+1;
int p13=pos1+2;
int p14=pos1+3;
char s11=mols[p11]->seq;
char s12=mols[p12]->seq;
char s13=mols[p13]->seq;
int t1,t2,t3;
if(s11=='A' || s11=='G'){t1=1;}else{t1=0;}	
if(s12=='A' || s12=='G'){t2=1;}else{t2=0;}	
if(s13=='A' || s13=='G'){t3=1;}else{t3=0;}	
int type1=t1+2*t2+4*t3;
int p20=pos2+3;
int p21=pos2+2;
int p22=pos2+1;
int p23=pos2+0;
int p24=pos2-1;
char s21=mols[p23]->seq;
char s22=mols[p22]->seq;
char s23=mols[p21]->seq;
if(s21=='A' || s21=='G'){t1=1;}else{t1=0;}	
if(s22=='A' || s22=='G'){t2=1;}else{t2=0;}	
if(s23=='A' || s23=='G'){t3=1;}else{t3=0;}	
int type2=t1+2*t2+4*t3;
//cout<< "Type1: " << type1 <<":" << pos1 << " : " << " Type 2: " << type2 <<":" << pos2 << endl;
//Create transformation matrix to map frag vectors to global coordinates
//Positive X-axis p0 to p1

double x0=mols[p10]->x;
double y0=mols[p10]->y;
double z0=mols[p10]->z;
Matrix * m1 = coordSysF(pos1);
//Create Vclass to translate
Vclass * o1=new Vclass();
o1->x=-1*x0;
o1->y=-1*y0;
o1->z=-1*z0;
//Create transformation matrix to map frag vectors to global coordinates
//Positive X-axis p0 to p1

double x0a=mols[p20]->x;
double y0a=mols[p20]->y;
double z0a=mols[p20]->z;
Matrix * m2 = coordSysR(pos2);

//cout << "created 2" << endl;
//cout << m2->a11 << ":" << m2->a12<<":" << m2->a13 << endl;
//cout << m2->a21 << ":" << m2->a22<<":" << m2->a23 << endl;
//cout << m2->a31 << ":" << m2->a32<<":" << m2->a33 << endl;

//cout << "created 1" << endl;
//cout << m1->a11 << ":" << m1->a12<<":" << m1->a13 << endl;
//cout << m1->a21 << ":" << m1->a22<<":" << m1->a23 << endl;
//cout << m1->a31 << ":" << m1->a32<<":" << m1->a33 << endl;

//Create Vclass to translate
Vclass * o2=new Vclass();
o2->x=-1*x0a;
o2->y=-1*y0a;
o2->z=-1*z0a;
//Calculate original angle between fragments relative to pos1 coordinates 

double odv=0;
double odr=0;
double odz=0;
double Lx=mols[p24]->x - mols[p14]->x;
double Ly=mols[p24]->y - mols[p14]->y;
double Lz=mols[p24]->z - mols[p14]->z;
//Transform into position 1 coordinates
Vclass * v = new Vclass();
v->x = Lx;
v->y = Ly;
v->z = Lz;
Matrix * m3 = coordSysE((pos1+2));
//cout << "Original Matrix" << endl;	
//cout << m3->a11 << " : " << m3->a12 << " : " << m3->a13 << endl;
//cout << m3->a21 << " : " << m3->a22 << " : " << m3->a23 << endl;
//cout << m3->a31 << " : " << m3->a32 << " : " << m3->a33 << endl;
//cout << "Original Coords" << endl;
//cout << mols[p14]->x << " : " << mols[p14]->y << " : " << mols[p14]->z << " : " << mols[p24]->x << " : " << mols[p24]->y << " : " << mols[p24]->z << endl;
Matrix * m4 = coordSysG((pos2));
v->transform(m3);
Lx=v->x;
Ly=v->y;
Lz=v->z;
//cout << "L" << endl;
//cout << Lx << " : " << Ly << " : " << Lz << endl; 
odv = m3->a11*m4->a11+m3->a21*m4->a21+m3->a31*m4->a31;
odr = m3->a12*m4->a12+m3->a22*m4->a22+m3->a32*m4->a32;
odz = m3->a13*m4->a13+m3->a23*m4->a23+m3->a33*m4->a33;

L=pow(L,2);
//Maybe return a vector of matches instead of pointers 
vector<Matches> mat;
mat=frags->fixedMoves(type1,type2,L,m1,m2,o1,o2,pos1,pos2);//Returns all pairs of fragments that hold the ends within a particular distance
delete o2;
delete o1;
delete v;
delete m3;
delete m1;
delete m2;
delete m4;
return mat;
}

void Polymere::printMoves(int pos1, int pos2,Fragclass * frags){
// Given a position if generates a set of moves that alters the given fragments
//Get frag1 type
int p10=pos1-1;
int p11=pos1;
int p12=pos1+1;
int p13=pos1+2;
int p14=pos1+3;
char s11=mols[p11]->seq;
char s12=mols[p12]->seq;
char s13=mols[p13]->seq;
int t1,t2,t3;
if(s11=='A' || s11=='G'){t1=1;}else{t1=0;}	
if(s12=='A' || s12=='G'){t2=1;}else{t2=0;}	
if(s13=='A' || s13=='G'){t3=1;}else{t3=0;}	
int type1=t1+2*t2+4*t3;
int p20=pos2+3;
int p21=pos2+2;
int p22=pos2+1;
int p23=pos2+0;
int p24=pos2-1;
char s21=mols[p23]->seq;
char s22=mols[p22]->seq;
char s23=mols[p21]->seq;
if(s21=='A' || s21=='G'){t1=1;}else{t1=0;}	
if(s22=='A' || s22=='G'){t2=1;}else{t2=0;}	
if(s23=='A' || s23=='G'){t3=1;}else{t3=0;}	
int type2=t1+2*t2+4*t3;
//cout<< "Type1: " << type1 <<":" << pos1 << " : " << " Type 2: " << type2 <<":" << pos2 << endl;
//Create transformation matrix to map frag vectors to global coordinates
//Positive X-axis p0 to p1

double x0=mols[p10]->x;
double y0=mols[p10]->y;
double z0=mols[p10]->z;
Matrix * m1 = coordSysF(pos1);
//Create Vclass to translate
Vclass * o1=new Vclass();
o1->x=-1*x0;
o1->y=-1*y0;
o1->z=-1*z0;
//Create transformation matrix to map frag vectors to global coordinates
//Positive X-axis p0 to p1

double x0a=mols[p20]->x;
double y0a=mols[p20]->y;
double z0a=mols[p20]->z;
Matrix * m2 = coordSysR(pos2);

//cout << "created 2" << endl;
//cout << m2->a11 << ":" << m2->a12<<":" << m2->a13 << endl;
//cout << m2->a21 << ":" << m2->a22<<":" << m2->a23 << endl;
//cout << m2->a31 << ":" << m2->a32<<":" << m2->a33 << endl;

//cout << "created 1" << endl;
//cout << m1->a11 << ":" << m1->a12<<":" << m1->a13 << endl;
//cout << m1->a21 << ":" << m1->a22<<":" << m1->a23 << endl;
//cout << m1->a31 << ":" << m1->a32<<":" << m1->a33 << endl;

//Create Vclass to translate
Vclass * o2=new Vclass();
o2->x=-1*x0a;
o2->y=-1*y0a;
o2->z=-1*z0a;
//Calculate original angle between fragments relative to pos1 coordinates 

double odv=0;
double odr=0;
double odz=0;
double Lx=mols[p24]->x - mols[p14]->x;
double Ly=mols[p24]->y - mols[p14]->y;
double Lz=mols[p24]->z - mols[p14]->z;
//Transform into position 1 coordinates
Vclass * v = new Vclass();
v->x = Lx;
v->y = Ly;
v->z = Lz;
Matrix * m3 = coordSysE((pos1+2));
//cout << "Original Matrix" << endl;	
//cout << m3->a11 << " : " << m3->a12 << " : " << m3->a13 << endl;
//cout << m3->a21 << " : " << m3->a22 << " : " << m3->a23 << endl;
//cout << m3->a31 << " : " << m3->a32 << " : " << m3->a33 << endl;
//cout << "Original Coords" << endl;
//cout << mols[p14]->x << " : " << mols[p14]->y << " : " << mols[p14]->z << " : " << mols[p24]->x << " : " << mols[p24]->y << " : " << mols[p24]->z << endl;
Matrix * m4 = coordSysG((pos2));
v->transform(m3);
Lx=v->x;
Ly=v->y;
Lz=v->z;
//cout << "L" << endl;
//cout << Lx << " : " << Ly << " : " << Lz << endl; 
odv = m3->a11*m4->a11+m3->a21*m4->a21+m3->a31*m4->a31;
odr = m3->a12*m4->a12+m3->a22*m4->a22+m3->a32*m4->a32;
odz = m3->a13*m4->a13+m3->a23*m4->a23+m3->a33*m4->a33;

double L=pow(Lx,2)+pow(Ly,2)+pow(Lz,2);
//Maybe return a vector of matches instead of pointers 
frags->printMoves(type1,type2,L,odr,odv,odz,m1,m2,o1,o2,Lx,Ly,Lz,pos1,pos2);
delete o2;
delete o1;
delete v;
delete m3;
delete m1;
delete m2;
delete m4;
}

vector<Matches> Polymere::chooseMovesPairs(int pos1, int pos2,Fragclass * frags){
// Given a position if generates a set of moves that alters the given fragments
//Get frag1 type
int p10=pos1-1;
int p11=pos1;
int p12=pos1+1;
int p13=pos1+2;
int p14=pos1+3;
char s11=mols[p11]->seq;
char s12=mols[p12]->seq;
char s13=mols[p13]->seq;
int t1,t2,t3;
if(s11=='A' || s11=='G'){t1=1;}else{t1=0;}	
if(s12=='A' || s12=='G'){t2=1;}else{t2=0;}	
if(s13=='A' || s13=='G'){t3=1;}else{t3=0;}	
int type1=t1+2*t2+4*t3;
int p20=pos2+3;
int p21=pos2+2;
int p22=pos2+1;
int p23=pos2+0;
int p24=pos2-1;
char s21=mols[p23]->seq;
char s22=mols[p22]->seq;
char s23=mols[p21]->seq;
if(s21=='A' || s21=='G'){t1=1;}else{t1=0;}	
if(s22=='A' || s22=='G'){t2=1;}else{t2=0;}	
if(s23=='A' || s23=='G'){t3=1;}else{t3=0;}	
int type2=t1+2*t2+4*t3;
//cout<< "Type1: " << type1 <<":" << pos1 << " : " << " Type 2: " << type2 <<":" << pos2 << endl;
//Create transformation matrix to map frag vectors to global coordinates
//Positive X-axis p0 to p1

double x0=mols[p10]->x;
double y0=mols[p10]->y;
double z0=mols[p10]->z;
Matrix * m1 = coordSysF(pos1);
//Create Vclass to translate
Vclass * o1=new Vclass();
o1->x=-1*x0;
o1->y=-1*y0;
o1->z=-1*z0;
//Create transformation matrix to map frag vectors to global coordinates
//Positive X-axis p0 to p1

double x0a=mols[p20]->x;
double y0a=mols[p20]->y;
double z0a=mols[p20]->z;
Matrix * m2 = coordSysR(pos2);

//cout << "created 2" << endl;
//cout << m2->a11 << ":" << m2->a12<<":" << m2->a13 << endl;
//cout << m2->a21 << ":" << m2->a22<<":" << m2->a23 << endl;
//cout << m2->a31 << ":" << m2->a32<<":" << m2->a33 << endl;

//cout << "created 1" << endl;
//cout << m1->a11 << ":" << m1->a12<<":" << m1->a13 << endl;
//cout << m1->a21 << ":" << m1->a22<<":" << m1->a23 << endl;
//cout << m1->a31 << ":" << m1->a32<<":" << m1->a33 << endl;

//Create Vclass to translate
Vclass * o2=new Vclass();
o2->x=-1*x0a;
o2->y=-1*y0a;
o2->z=-1*z0a;
//Calculate original angle between fragments relative to pos1 coordinates 

double odv=0;
double odr=0;
double odz=0;
double Lx=mols[p24]->x - mols[p14]->x;
double Ly=mols[p24]->y - mols[p14]->y;
double Lz=mols[p24]->z - mols[p14]->z;
//Transform into position 1 coordinates
Vclass * v = new Vclass();
v->x = Lx;
v->y = Ly;
v->z = Lz;
Matrix * m3 = coordSysE((pos1+2));
//cout << "Original Matrix" << endl;	
//cout << m3->a11 << " : " << m3->a12 << " : " << m3->a13 << endl;
//cout << m3->a21 << " : " << m3->a22 << " : " << m3->a23 << endl;
//cout << m3->a31 << " : " << m3->a32 << " : " << m3->a33 << endl;
//cout << "Original Coords" << endl;
//cout << mols[p14]->x << " : " << mols[p14]->y << " : " << mols[p14]->z << " : " << mols[p24]->x << " : " << mols[p24]->y << " : " << mols[p24]->z << endl;
Matrix * m4 = coordSysG((pos2));
v->transform(m3);
Lx=v->x;
Ly=v->y;
Lz=v->z;
//cout << "L" << endl;
//cout << Lx << " : " << Ly << " : " << Lz << endl; 
odv = m3->a11*m4->a11+m3->a21*m4->a21+m3->a31*m4->a31;
odr = m3->a12*m4->a12+m3->a22*m4->a22+m3->a32*m4->a32;
odz = m3->a13*m4->a13+m3->a23*m4->a23+m3->a33*m4->a33;

double L=pow(Lx,2)+pow(Ly,2)+pow(Lz,2);
//Maybe return a vector of matches instead of pointers 
vector<Matches> mat;
mat=frags->matchedMoves(type1,type2,-1,-1,-1,-1,m1,m2,o1,o2,Lx,Ly,Lz,-1,pos1,pos2);
delete o2;
delete o1;
delete v;
delete m3;
delete m1;
delete m2;
delete m4;
return mat;
}

vector<Matches> Polymere::chooseMoves(int pos1, int pos2,Fragclass * frags,int t){
// Given a position if generates a set of moves that alters the given fragments
//Get frag1 type
int p10=pos1-1;
int p11=pos1;
int p12=pos1+1;
int p13=pos1+2;
int p14=pos1+3;
char s11=mols[p11]->seq;
char s12=mols[p12]->seq;
char s13=mols[p13]->seq;
int t1,t2,t3;
if(s11=='A' || s11=='G'){t1=1;}else{t1=0;}	
if(s12=='A' || s12=='G'){t2=1;}else{t2=0;}	
if(s13=='A' || s13=='G'){t3=1;}else{t3=0;}	
int type1=t1+2*t2+4*t3;
int p20=pos2+3;
int p21=pos2+2;
int p22=pos2+1;
int p23=pos2+0;
int p24=pos2-1;
char s21=mols[p23]->seq;
char s22=mols[p22]->seq;
char s23=mols[p21]->seq;
if(s21=='A' || s21=='G'){t1=1;}else{t1=0;}	
if(s22=='A' || s22=='G'){t2=1;}else{t2=0;}	
if(s23=='A' || s23=='G'){t3=1;}else{t3=0;}	
int type2=t1+2*t2+4*t3;
//cout<< "Type1: " << type1 <<":" << pos1 << " : " << " Type 2: " << type2 <<":" << pos2 << endl;
//Create transformation matrix to map frag vectors to global coordinates
//Positive X-axis p0 to p1

double x0=mols[p10]->x;
double y0=mols[p10]->y;
double z0=mols[p10]->z;
Matrix * m1 = coordSysF(pos1);
//Create Vclass to translate
Vclass * o1=new Vclass();
o1->x=-1*x0;
o1->y=-1*y0;
o1->z=-1*z0;
//Create transformation matrix to map frag vectors to global coordinates
//Positive X-axis p0 to p1

double x0a=mols[p20]->x;
double y0a=mols[p20]->y;
double z0a=mols[p20]->z;
Matrix * m2 = coordSysR(pos2);

//cout << "created 2" << endl;
//cout << m2->a11 << ":" << m2->a12<<":" << m2->a13 << endl;
//cout << m2->a21 << ":" << m2->a22<<":" << m2->a23 << endl;
//cout << m2->a31 << ":" << m2->a32<<":" << m2->a33 << endl;

//cout << "created 1" << endl;
//cout << m1->a11 << ":" << m1->a12<<":" << m1->a13 << endl;
//cout << m1->a21 << ":" << m1->a22<<":" << m1->a23 << endl;
//cout << m1->a31 << ":" << m1->a32<<":" << m1->a33 << endl;

//Create Vclass to translate
Vclass * o2=new Vclass();
o2->x=-1*x0a;
o2->y=-1*y0a;
o2->z=-1*z0a;
//Calculate original angle between fragments relative to pos1 coordinates 

double odv=0;
double odr=0;
double odz=0;
double Lx=mols[p24]->x - mols[p14]->x;
double Ly=mols[p24]->y - mols[p14]->y;
double Lz=mols[p24]->z - mols[p14]->z;
//Transform into position 1 coordinates
Vclass * v = new Vclass();
v->x = Lx;
v->y = Ly;
v->z = Lz;
Matrix * m3 = coordSysE((pos1+2));
//cout << "Original Matrix" << endl;	
//cout << m3->a11 << " : " << m3->a12 << " : " << m3->a13 << endl;
//cout << m3->a21 << " : " << m3->a22 << " : " << m3->a23 << endl;
//cout << m3->a31 << " : " << m3->a32 << " : " << m3->a33 << endl;
//cout << "Original Coords" << endl;
//cout << mols[p14]->x << " : " << mols[p14]->y << " : " << mols[p14]->z << " : " << mols[p24]->x << " : " << mols[p24]->y << " : " << mols[p24]->z << endl;
Matrix * m4 = coordSysG((pos2));
v->transform(m3);
Lx=v->x;
Ly=v->y;
Lz=v->z;
//cout << "L" << endl;
//cout << Lx << " : " << Ly << " : " << Lz << endl; 
odv = m3->a11*m4->a11+m3->a21*m4->a21+m3->a31*m4->a31;
odr = m3->a12*m4->a12+m3->a22*m4->a22+m3->a32*m4->a32;
odz = m3->a13*m4->a13+m3->a23*m4->a23+m3->a33*m4->a33;

double L=pow(Lx,2)+pow(Ly,2)+pow(Lz,2);
//Maybe return a vector of matches instead of pointers 
vector<Matches> mat;
mat=frags->matchedMoves(type1,type2,L,odr,odv,odz,m1,m2,o1,o2,Lx,Ly,Lz,t,pos1,pos2);
delete o2;
delete o1;
delete v;
delete m3;
delete m1;
delete m2;
delete m4;
return mat;
}


vector<Fragatoms *> Polymere::readAtoms(string file){
//Reads in the base positions in fragment database	
using namespace std;
vector<Fragatoms *> frags;
string field;
std::ifstream input;
input.open(file.c_str());
string line;
char nt1,nt2,nt3;
int a=0;
int w=0;
int id;
double x,y,z;
while(getline(input, line)) 
{
w=0;
//Read in file with each line being x \t y \t z \t group
//NT1_C1,NT1_B1,NT1_B2....
frags.push_back(new Fragatoms());
	istringstream iss(line,istringstream::in);	
		getline(iss,field,'\t');
	id=strtoint(field);
frags[a]->id=id;
		getline(iss,field,'\t');
	x=strtodouble(field);
		getline(iss,field,'\t');
	y=strtodouble(field);
		getline(iss,field,'\t');
	z=strtodouble(field);
frags[a]->atoms[w].x=x;
frags[a]->atoms[w].y=y;
frags[a]->atoms[w++].z=z;
		getline(iss,field,'\t');
	x=strtodouble(field);
		getline(iss,field,'\t');
	y=strtodouble(field);
		getline(iss,field,'\t');
	z=strtodouble(field);
frags[a]->atoms[w].x=x;
frags[a]->atoms[w].y=y;
frags[a]->atoms[w++].z=z;
		getline(iss,field,'\t');
	x=strtodouble(field);
		getline(iss,field,'\t');
	y=strtodouble(field);
		getline(iss,field,'\t');
	z=strtodouble(field);
frags[a]->atoms[w].x=x;
frags[a]->atoms[w].y=y;
frags[a]->atoms[w++].z=z;
		getline(iss,field,'\t');
	nt1=strtochar(field);
		getline(iss,field,'\t');
	x=strtodouble(field);
		getline(iss,field,'\t');
	y=strtodouble(field);
		getline(iss,field,'\t');
	z=strtodouble(field);
frags[a]->atoms[w].x=x;
frags[a]->atoms[w].y=y;
frags[a]->atoms[w++].z=z;
		getline(iss,field,'\t');
	x=strtodouble(field);
		getline(iss,field,'\t');
	y=strtodouble(field);
		getline(iss,field,'\t');
	z=strtodouble(field);
frags[a]->atoms[w].x=x;
frags[a]->atoms[w].y=y;
frags[a]->atoms[w++].z=z;
		getline(iss,field,'\t');
	x=strtodouble(field);
		getline(iss,field,'\t');
	y=strtodouble(field);
		getline(iss,field,'\t');
	z=strtodouble(field);
frags[a]->atoms[w].x=x;
frags[a]->atoms[w].y=y;
frags[a]->atoms[w++].z=z;
		getline(iss,field,'\t');
	nt2=strtochar(field);
		getline(iss,field,'\t');
	x=strtodouble(field);
		getline(iss,field,'\t');
	y=strtodouble(field);
		getline(iss,field,'\t');
	z=strtodouble(field);
frags[a]->atoms[w].x=x;
frags[a]->atoms[w].y=y;
frags[a]->atoms[w++].z=z;
		getline(iss,field,'\t');
	x=strtodouble(field);
		getline(iss,field,'\t');
	y=strtodouble(field);
		getline(iss,field,'\t');
	z=strtodouble(field);
frags[a]->atoms[w].x=x;
frags[a]->atoms[w].y=y;
frags[a]->atoms[w++].z=z;
		getline(iss,field,'\t');
	x=strtodouble(field);
		getline(iss,field,'\t');
	y=strtodouble(field);
		getline(iss,field,'\t');
	z=strtodouble(field);
frags[a]->atoms[w].x=x;
frags[a]->atoms[w].y=y;
frags[a]->atoms[w++].z=z;
		getline(iss,field,'\t');
	nt3=strtochar(field);
a++;
}
cout << "number of frags atoms: " << a << ":" << w <<endl;
return frags;
}
vector<Fragment *> Polymere::readFrags2(string file){
	//Given a file with 2D coordinates generate from secondary structure it loads positions into polymere
	using namespace std;
double len_1;
double len_2;
double len_3;
double len_4;
double tor_1;
double tor_2;
double bond_0;
double bond_1;
double bond_2;
int frag_type;
vector<Fragment *> frags;
string field;
	std::ifstream input;
	input.open(file.c_str());
	string line;
	int a=0;
while(getline(input, line)) 
{
//Read in file with each line being x \t y \t z \t group
	istringstream iss(line,istringstream::in);	
		getline(iss,field,'\t');
	len_1=strtodouble(field);
		getline(iss,field,'\t');
	len_2=strtodouble(field);
		getline(iss,field,'\t');
	len_3=strtodouble(field);
		getline(iss,field,'\t');
	len_4=strtodouble(field);
		getline(iss,field,'\t');
	tor_1=strtodouble(field);
		getline(iss,field,'\t');
	tor_2=strtodouble(field);
		getline(iss,field,'\t');
	bond_0=strtodouble(field);
		getline(iss,field,'\t');
	bond_1=strtodouble(field);
		getline(iss,field,'\t');
	bond_2=strtodouble(field);
		getline(iss,field,'\t');
	frag_type=strtoint(field);
frags.push_back(new Fragment());
frags[a]->len_1=  len_1;
frags[a]->len_2=  len_2;
frags[a]->len_3=  len_3;
frags[a]->len_4=  len_4;
frags[a]->tor_1=  tor_1;
frags[a]->tor_2= tor_2;
frags[a]->bond_0=bond_0;
frags[a]->bond_1=bond_1;
frags[a]->bond_2=bond_2;
frags[a]->frag_type=frag_type;
//cout << len_1 << ":" << len_2 << ":" << tor_1 << ":" << tor_2 << ":" << bond << endl;
a++;
}
//cout << "number of frags: " << a << endl;
return frags;
}

float Polymere::doublet(std::vector<Doublet> dlt){
//returns the sum of average doublet rmsd's for each position
float ave=0;
int tot=0;
for(int i =0; i < num_mols; i++){
float ave1=0;
int tot1=0;
for(int j =0; j < num_mols; j++){
if(j!=i){
double dist=moldist(i,j);
if(dist < 14){
double cop=coplanar(i,j);
if(cop > 2){
float d=rms_doublet(dlt,i,j);
if(d >= 0){
ave1=ave1+d;
tot1=tot1++;
}
}
}
}
}
if(tot1 > 0){
ave = ave+ave1/tot1;
}
tot = tot+tot1;
}
ave = ave/tot *100;

return ave;
}

//Make a scoring function that takes in different weights as its arguments
/*
float Polymere::final_score(std::vector<Doublet> dlt23,float hydro_weight,float dbl_weight){


} 
*/

float Polymere::rms_doubletA(std::vector<Doublet> dlt23, int i,int j){
//Returns the average of doublet fits
//AlignBase
vector <Points> p;
vector <Points> p2;
Points c1,b1,b2;
int type1,type2;
//Returns a score for hydrogen bonding
c1.x=mols[i]->x;
c1.y=mols[i]->y;
c1.z=mols[i]->z;
b1.x=mols[i]->bx;
b1.y=mols[i]->by;
b1.z=mols[i]->bz;
b2.x=mols[i]->b2x;
b2.y=mols[i]->b2y;
b2.z=mols[i]->b2z;
p.push_back(c1);
p.push_back(b1);
p.push_back(b2);

c1.x=mols[j]->x;
c1.y=mols[j]->y;
c1.z=mols[j]->z;
b1.x=mols[j]->bx;
b1.y=mols[j]->by;
b1.z=mols[j]->bz;
b2.x=mols[j]->b2x;
b2.y=mols[j]->b2y;
b2.z=mols[j]->b2z;
p2.push_back(c1);
p2.push_back(b1);
p2.push_back(b2);

if(mols[i]->seq=='G'){type1=1;}else{
if(mols[i]->seq=='A'){type1=2;}else{
if(mols[i]->seq=='U'){type1=3;}else{
if(mols[i]->seq=='C'){type1=4;}else{
}
}
}
}
if(mols[j]->seq=='G'){type2=1;}else{
if(mols[j]->seq=='A'){type2=2;}else{
if(mols[j]->seq=='U'){type2=3;}else{
if(mols[j]->seq=='C'){type2=4;}else{
}
}
}
}
FABase * base1 = new FABase(type1);  
FABase * base2 = new FABase(type2);  

if(type1==1 || type1==2){type1=1;}else{type1=0;}
if(type2==1 || type2==2){type2=1;}else{type2=0;}

//Move to coordinate system
/*
cout << "Base Coordinates: " << endl;
cout << base1->C1.x << ":" << base1->C1.y << ":"  << base1->C1.z << endl; 
cout << base1->C2.x << ":" << base1->C2.y << ":"  << base1->C2.z << endl; 
cout << base1->C4.x << ":" << base1->C4.y << ":"  << base1->C4.z << endl; 
cout << base1->C6.x << ":" << base1->C6.y << ":"  << base1->C6.z << endl; 
cout << "Base Coordinates: " << endl;
cout << base2->C1.x << ":" << base2->C1.y << ":"  << base2->C1.z << endl; 
cout << base2->C2.x << ":" << base2->C2.y << ":"  << base2->C2.z << endl; 
cout << base2->C4.x << ":" << base2->C4.y << ":"  << base2->C4.z << endl; 
cout << base2->C6.x << ":" << base2->C6.y << ":"  << base2->C6.z << endl; 
*/
base1->alignBase(p);
base2->alignBase(p2);
float c1x1 =(float) (base1->C1.x); 
float c2x1 =(float) (base1->C2.x); 
float c4x1 =(float) (base1->C4.x); 
float c6x1 =(float) (base1->C6.x);
float c1x2 =(float) (base2->C1.x); 
float c2x2 =(float) (base2->C2.x); 
float c4x2 =(float) (base2->C4.x);
float c6x2 =(float) (base2->C6.x); 
float c1y1 =(float) (base1->C1.y); 
float c2y1 =(float) (base1->C2.y); 
float c4y1 =(float) (base1->C4.y); 
float c6y1 =(float) (base1->C6.y); 
float c1y2 =(float) (base2->C1.y); 
float c2y2 =(float) (base2->C2.y); 
float c4y2 =(float) (base2->C4.y); 
float c6y2 =(float) (base2->C6.y); 
float c1z1 =(float) (base1->C1.z);
float c2z1 =(float) (base1->C2.z);
float c4z1 =(float) (base1->C4.z);
float c6z1 =(float) (base1->C6.z);
float c1z2 =(float) (base2->C1.z);
float c2z2 =(float) (base2->C2.z);
float c4z2 =(float) (base2->C4.z);
float c6z2 =(float) (base2->C6.z);
//If no pairwise distance < 4 return -1

double d1 =sqrt(pow((c1x1-c1x2),2)+pow((c1y1-c1y2),2)+pow((c1z1-c1z2),2));
double d2 =sqrt(pow((c1x1-c2x2),2)+pow((c1y1-c2y2),2)+pow((c1z1-c2z2),2));
double d3 =sqrt(pow((c1x1-c4x2),2)+pow((c1y1-c4y2),2)+pow((c1z1-c4z2),2));
double d4 =sqrt(pow((c1x1-c6x2),2)+pow((c1y1-c6y2),2)+pow((c1z1-c6z2),2));
double d5 =sqrt(pow((c2x1-c1x2),2)+pow((c2y1-c1y2),2)+pow((c2z1-c1z2),2));
double d6 =sqrt(pow((c2x1-c2x2),2)+pow((c2y1-c2y2),2)+pow((c2z1-c2z2),2));
double d7 =sqrt(pow((c2x1-c4x2),2)+pow((c2y1-c4y2),2)+pow((c2z1-c4z2),2));
double d8 =sqrt(pow((c2x1-c6x2),2)+pow((c2y1-c6y2),2)+pow((c2z1-c6z2),2));
double d9 =sqrt(pow((c4x1-c2x2),2)+pow((c4y1-c2y2),2)+pow((c4z1-c2z2),2));
double d10=sqrt(pow((c4x1-c4x2),2)+pow((c4y1-c4y2),2)+pow((c4z1-c4z2),2));
double d11=sqrt(pow((c4x1-c6x2),2)+pow((c4y1-c6y2),2)+pow((c4z1-c6z2),2));
double d12=sqrt(pow((c6x1-c1x2),2)+pow((c6y1-c1y2),2)+pow((c6z1-c1z2),2));
double d13=sqrt(pow((c6x1-c2x2),2)+pow((c6y1-c2y2),2)+pow((c6z1-c2z2),2));
double d14=sqrt(pow((c6x1-c4x2),2)+pow((c6y1-c4y2),2)+pow((c6z1-c4z2),2));
double d15=sqrt(pow((c6x1-c6x2),2)+pow((c6y1-c6y2),2)+pow((c6z1-c6z2),2));
double d16 =sqrt(pow((c4x1-c1x2),2)+pow((c4y1-c1y2),2)+pow((c4z1-c1z2),2));
int tot=0;
float sumd=0;
float mind=-1;
if(d1  > 15.5 && d2  > 15.5 && d3  > 15.5 && d4  > 15.5 && d5  > 15.5 && d6  > 15.5 && d7  > 15.5 && d8  > 15.5 && d16  > 15.5 && d9  > 15.5 && d10 > 15.5 && d11 > 15.5 && d12 > 15.5 && d13 > 15.5 && d14 > 15.5 && d15 > 15.5){
mind=-1;
}else{
//Extract coordinates from each base
for(int d_ind=0; d_ind < dlt23.size(); d_ind++){
if(dlt23[d_ind].type1==type1 && dlt23[d_ind].type2==type2){
float d=9999999;
float *v1,*v2;
float *mtx;
float x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16;
v1=NULL;
v2=NULL;
mtx=NULL;
mtx=(float*) malloc(16*sizeof(float));
  v1 =(float*) malloc(24 * sizeof(double));
  v2 =(float*) malloc(24 * sizeof(double));
if(v1==NULL || v2==NULL || mtx==NULL){
d=-1;
cout << "Allocation Failed" << endl;
}else{
x5 =c1x1;
y5 =c1y1;
z5 =c1z1;
v1[0]=x5;
v1[1]=y5;
v1[2]=z5;
x6 =c4x1;
y6 =c4y1;
z6 =c4z1;
v1[3]=x6;
v1[4]=y6;
v1[5]=z6;
x7 =c2x1;
y7 =c2y1;
z7 =c2z1;
v1[6]=x7;
v1[7]=y7;
v1[8]=z7;
x8 =c6x1;
y8 =c6y1;
z8 =c6z1;
v1[9]=x8;
v1[10]=y8;
v1[11]=z8;
x1 =c1x2;
y1 =c1y2;
z1 =c1z2;
v1[12]=x1;
v1[13]=y1;
v1[14]=z1;
x2 =c4x2;
y2 =c4y2;
z2 =c4z2;
v1[15]=x2;
v1[16]=y2;
v1[17]=z2;
x3 =c2x2;
y3 =c2y2;
z3 =c2z2;
v1[18]=x3;
v1[19]=y3;
v1[20]=z3;
x4 =c6x2;
y4 =c6y2;
z4 =c6z2;
v1[21]=x4;
v1[22]=y4;
v1[23]=z4;

///Doublet information 
x9 =dlt23[d_ind].C1_1x;
y9 =dlt23[d_ind].C1_1y;
z9 =dlt23[d_ind].C1_1z;
v2[0]=x9;
v2[1]=y9;
v2[2]=z9;
x10 =dlt23[d_ind].C4_1x;
y10 =dlt23[d_ind].C4_1y;
z10 =dlt23[d_ind].C4_1z;
v2[3]=x10;
v2[4]=y10;
v2[5]=z10;
x11 =dlt23[d_ind].C2_1x;
y11 =dlt23[d_ind].C2_1y;
z11 =dlt23[d_ind].C2_1z;
v2[6]=x11;
v2[7]=y11;
v2[8]=z11;
x12 =dlt23[d_ind].C6_1x;
y12 =dlt23[d_ind].C6_1y;
z12 =dlt23[d_ind].C6_1z;
 v2[9]=x12;
v2[10]=y12;
v2[11]=z12;
x13 =dlt23[d_ind].C1_2x;
y13 =dlt23[d_ind].C1_2y;
z13 =dlt23[d_ind].C1_2z;
v2[12]=x13;
v2[13]=y13;
v2[14]=z13;
x14 =dlt23[d_ind].C4_2x;
y14 =dlt23[d_ind].C4_2y;
z14 =dlt23[d_ind].C4_2z;
v2[15]=x14;
v2[16]=y14;
v2[17]=z14;
x15 =dlt23[d_ind].C2_2x;
y15 =dlt23[d_ind].C2_2y;
z15 =dlt23[d_ind].C2_2z;
v2[18]=x15;
v2[19]=y15;
v2[20]=z15;
x16 =dlt23[d_ind].C6_2x;
y16 =dlt23[d_ind].C6_2y;
z16 =dlt23[d_ind].C6_2z;
v2[21]=x16;
v2[22]=y16;
v2[23]=z16;
d=rmsd(v1,v2,8,mtx);
sumd = sumd+d;
tot++;
//cout<< "RMSD Doublet: " << d << endl;
  	free(v1);
  	free(mtx);
  	free(v2);
}
}
}
mind = sumd/tot;
}
return mind;
}

float Polymere::rms_doublet(std::vector<Doublet> dlt23, int i,int j){
//Returns the average of the best fits for each possible doublet
//AlignBase
vector <Points> p;
vector <Points> p2;
Points c1,b1,b2;
int type1,type2;
//Returns a score for hydrogen bonding
c1.x=mols[i]->x;
c1.y=mols[i]->y;
c1.z=mols[i]->z;
b1.x=mols[i]->bx;
b1.y=mols[i]->by;
b1.z=mols[i]->bz;
b2.x=mols[i]->b2x;
b2.y=mols[i]->b2y;
b2.z=mols[i]->b2z;
p.push_back(c1);
p.push_back(b1);
p.push_back(b2);

c1.x=mols[j]->x;
c1.y=mols[j]->y;
c1.z=mols[j]->z;
b1.x=mols[j]->bx;
b1.y=mols[j]->by;
b1.z=mols[j]->bz;
b2.x=mols[j]->b2x;
b2.y=mols[j]->b2y;
b2.z=mols[j]->b2z;
p2.push_back(c1);
p2.push_back(b1);
p2.push_back(b2);

if(mols[i]->seq=='G'){type1=1;}else{
if(mols[i]->seq=='A'){type1=2;}else{
if(mols[i]->seq=='U'){type1=3;}else{
if(mols[i]->seq=='C'){type1=4;}else{
}
}
}
}
if(mols[j]->seq=='G'){type2=1;}else{
if(mols[j]->seq=='A'){type2=2;}else{
if(mols[j]->seq=='U'){type2=3;}else{
if(mols[j]->seq=='C'){type2=4;}else{
}
}
}
}
FABase * base1 = new FABase(type1);  
FABase * base2 = new FABase(type2);  

if(type1==1 || type1==2){type1=1;}else{type1=0;}
if(type2==1 || type2==2){type2=1;}else{type2=0;}

//Move to coordinate system
/*
cout << "Base Coordinates: " << endl;
cout << base1->C1.x << ":" << base1->C1.y << ":"  << base1->C1.z << endl; 
cout << base1->C2.x << ":" << base1->C2.y << ":"  << base1->C2.z << endl; 
cout << base1->C4.x << ":" << base1->C4.y << ":"  << base1->C4.z << endl; 
cout << base1->C6.x << ":" << base1->C6.y << ":"  << base1->C6.z << endl; 
cout << "Base Coordinates: " << endl;
cout << base2->C1.x << ":" << base2->C1.y << ":"  << base2->C1.z << endl; 
cout << base2->C2.x << ":" << base2->C2.y << ":"  << base2->C2.z << endl; 
cout << base2->C4.x << ":" << base2->C4.y << ":"  << base2->C4.z << endl; 
cout << base2->C6.x << ":" << base2->C6.y << ":"  << base2->C6.z << endl; 
*/
base1->alignBase(p);
base2->alignBase(p2);
float c1x1 =(float) (base1->C1.x); 
float c2x1 =(float) (base1->C2.x); 
float c4x1 =(float) (base1->C4.x); 
float c6x1 =(float) (base1->C6.x);
float c1x2 =(float) (base2->C1.x); 
float c2x2 =(float) (base2->C2.x); 
float c4x2 =(float) (base2->C4.x);
float c6x2 =(float) (base2->C6.x); 
float c1y1 =(float) (base1->C1.y); 
float c2y1 =(float) (base1->C2.y); 
float c4y1 =(float) (base1->C4.y); 
float c6y1 =(float) (base1->C6.y); 
float c1y2 =(float) (base2->C1.y); 
float c2y2 =(float) (base2->C2.y); 
float c4y2 =(float) (base2->C4.y); 
float c6y2 =(float) (base2->C6.y); 
float c1z1 =(float) (base1->C1.z);
float c2z1 =(float) (base1->C2.z);
float c4z1 =(float) (base1->C4.z);
float c6z1 =(float) (base1->C6.z);
float c1z2 =(float) (base2->C1.z);
float c2z2 =(float) (base2->C2.z);
float c4z2 =(float) (base2->C4.z);
float c6z2 =(float) (base2->C6.z);
//If no pairwise distance < 4 return -1

double d1 =sqrt(pow((c1x1-c1x2),2)+pow((c1y1-c1y2),2)+pow((c1z1-c1z2),2));
double d2 =sqrt(pow((c1x1-c2x2),2)+pow((c1y1-c2y2),2)+pow((c1z1-c2z2),2));
double d3 =sqrt(pow((c1x1-c4x2),2)+pow((c1y1-c4y2),2)+pow((c1z1-c4z2),2));
double d4 =sqrt(pow((c1x1-c6x2),2)+pow((c1y1-c6y2),2)+pow((c1z1-c6z2),2));
double d5 =sqrt(pow((c2x1-c1x2),2)+pow((c2y1-c1y2),2)+pow((c2z1-c1z2),2));
double d6 =sqrt(pow((c2x1-c2x2),2)+pow((c2y1-c2y2),2)+pow((c2z1-c2z2),2));
double d7 =sqrt(pow((c2x1-c4x2),2)+pow((c2y1-c4y2),2)+pow((c2z1-c4z2),2));
double d8 =sqrt(pow((c2x1-c6x2),2)+pow((c2y1-c6y2),2)+pow((c2z1-c6z2),2));
double d9 =sqrt(pow((c4x1-c2x2),2)+pow((c4y1-c2y2),2)+pow((c4z1-c2z2),2));
double d10=sqrt(pow((c4x1-c4x2),2)+pow((c4y1-c4y2),2)+pow((c4z1-c4z2),2));
double d11=sqrt(pow((c4x1-c6x2),2)+pow((c4y1-c6y2),2)+pow((c4z1-c6z2),2));
double d12=sqrt(pow((c6x1-c1x2),2)+pow((c6y1-c1y2),2)+pow((c6z1-c1z2),2));
double d13=sqrt(pow((c6x1-c2x2),2)+pow((c6y1-c2y2),2)+pow((c6z1-c2z2),2));
double d14=sqrt(pow((c6x1-c4x2),2)+pow((c6y1-c4y2),2)+pow((c6z1-c4z2),2));
double d15=sqrt(pow((c6x1-c6x2),2)+pow((c6y1-c6y2),2)+pow((c6z1-c6z2),2));
double d16 =sqrt(pow((c4x1-c1x2),2)+pow((c4y1-c1y2),2)+pow((c4z1-c1z2),2));
float mind=99999;

if(d1  > 15.5 && d2  > 15.5 && d3  > 15.5 && d4  > 15.5 && d5  > 15.5 && d6  > 15.5 && d7  > 15.5 && d8  > 15.5 && d16  > 15.5 && d9  > 15.5 && d10 > 15.5 && d11 > 15.5 && d12 > 15.5 && d13 > 15.5 && d14 > 15.5 && d15 > 15.5){
mind=-1;

}else{

 
//Extract coordinates from each base
for(int d_ind=0; d_ind < dlt23.size(); d_ind++){
if(dlt23[d_ind].type1==type1 && dlt23[d_ind].type2==type2){
float d=9999999;
float *v1,*v2;
float *mtx;
float x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16;
v1=NULL;
v2=NULL;
mtx=NULL;
mtx=(float*) malloc(16*sizeof(float));
  v1 =(float*) malloc(24 * sizeof(double));
  v2 =(float*) malloc(24 * sizeof(double));
if(v1==NULL || v2==NULL || mtx==NULL){
d=-1;
cout << "Allocation Failed" << endl;
}else{
x5 =c1x1;
y5 =c1y1;
z5 =c1z1;
v1[0]=x5;
v1[1]=y5;
v1[2]=z5;
x6 =c4x1;
y6 =c4y1;
z6 =c4z1;
v1[3]=x6;
v1[4]=y6;
v1[5]=z6;
x7 =c2x1;
y7 =c2y1;
z7 =c2z1;
v1[6]=x7;
v1[7]=y7;
v1[8]=z7;
x8 =c6x1;
y8 =c6y1;
z8 =c6z1;
v1[9]=x8;
v1[10]=y8;
v1[11]=z8;
x1 =c1x2;
y1 =c1y2;
z1 =c1z2;
v1[12]=x1;
v1[13]=y1;
v1[14]=z1;
x2 =c4x2;
y2 =c4y2;
z2 =c4z2;
v1[15]=x2;
v1[16]=y2;
v1[17]=z2;
x3 =c2x2;
y3 =c2y2;
z3 =c2z2;
v1[18]=x3;
v1[19]=y3;
v1[20]=z3;
x4 =c6x2;
y4 =c6y2;
z4 =c6z2;
v1[21]=x4;
v1[22]=y4;
v1[23]=z4;

///Doublet information 
x9 =dlt23[d_ind].C1_1x;
y9 =dlt23[d_ind].C1_1y;
z9 =dlt23[d_ind].C1_1z;
v2[0]=x9;
v2[1]=y9;
v2[2]=z9;
x10 =dlt23[d_ind].C4_1x;
y10 =dlt23[d_ind].C4_1y;
z10 =dlt23[d_ind].C4_1z;
v2[3]=x10;
v2[4]=y10;
v2[5]=z10;
x11 =dlt23[d_ind].C2_1x;
y11 =dlt23[d_ind].C2_1y;
z11 =dlt23[d_ind].C2_1z;
v2[6]=x11;
v2[7]=y11;
v2[8]=z11;
x12 =dlt23[d_ind].C6_1x;
y12 =dlt23[d_ind].C6_1y;
z12 =dlt23[d_ind].C6_1z;
 v2[9]=x12;
v2[10]=y12;
v2[11]=z12;
x13 =dlt23[d_ind].C1_2x;
y13 =dlt23[d_ind].C1_2y;
z13 =dlt23[d_ind].C1_2z;
v2[12]=x13;
v2[13]=y13;
v2[14]=z13;
x14 =dlt23[d_ind].C4_2x;
y14 =dlt23[d_ind].C4_2y;
z14 =dlt23[d_ind].C4_2z;
v2[15]=x14;
v2[16]=y14;
v2[17]=z14;
x15 =dlt23[d_ind].C2_2x;
y15 =dlt23[d_ind].C2_2y;
z15 =dlt23[d_ind].C2_2z;
v2[18]=x15;
v2[19]=y15;
v2[20]=z15;
x16 =dlt23[d_ind].C6_2x;
y16 =dlt23[d_ind].C6_2y;
z16 =dlt23[d_ind].C6_2z;
v2[21]=x16;
v2[22]=y16;
v2[23]=z16;
d=rmsd(v1,v2,8,mtx);
//cout<< "RMSD Doublet: " << d << endl;
  	free(v1);
  	free(mtx);
  	free(v2);
}
if(d < mind){mind = d;}

}
}
}
return mind;
}

void Polymere::structParams(std::vector<Doublet> dlt23, int i,int j,int struct_id){
using namespace std;
float doublet=rms_doublet(dlt23,i,j);
if(doublet == -1){
//There are no heavy atoms within 5.3 Angstroms
cout << "Too far apart.." << endl;
}else{
float hyd = hydroP(i,j);
float minAA1=minAA(i,j);
float minDD1=minDD(i,j);
float doubletA=rms_doubletA(dlt23,i,j);
vector<double> params;
params.reserve(10);
for(int paidn=0; paidn<10; paidn++){params.push_back(0);}
getParams(i,j,&params);
int type1, type2;
if(mols[i]->seq=='G'){type1=1;}else{
if(mols[i]->seq=='A'){type1=2;}else{
if(mols[i]->seq=='U'){type1=3;}else{
if(mols[i]->seq=='C'){type1=4;}else{
}
}
}
}
if(mols[j]->seq=='G'){type2=1;}else{
if(mols[j]->seq=='A'){type2=2;}else{
if(mols[j]->seq=='U'){type2=3;}else{
if(mols[j]->seq=='C'){type2=4;}else{
}
}
}
}
/*
params->at(0)=type;
params->at(1)=rho;
params->at(2)=theta1;
params->at(3)=phi1;
params->at(4)=phi2;
params->at(5)=theta2;
params->at(6)=dx_i;
params->at(7)=dy_i;
params->at(8)=dz_i;
*/
//cout << struct_id << "," << i << "," << j << "," <<  type1 << "," << type2 << "," << hyd << "," << minAA1 << "," << minDD1 << "," << doublet << "," << doubletA << "," << params[6] << "," << params[7] << "," << params[8] << "," << params[1] << "," << params[2] << "," << params[3] << "," << params[5] << "," << params[4] << "," <<  params[9] << endl;
}

}

double Polymere::rmsf(Polymere * p2){
double d;
float *v1,*v2;
float *mtx;
v1=NULL;
v2=NULL;
mtx=NULL;
mtx=(float*) malloc(16*sizeof(float));
  v1 =(float*) malloc(num_mols * 9 * sizeof(float));
  v2 =(float*) malloc(num_mols * 9 * sizeof(float));
if(v1==NULL || v2==NULL || mtx==NULL){
d=-1;
cout << "Allocation Failed" << endl;
}else{
int i=0;
for(int j=0; j<num_mols; j++){
if(p2->mols[j]->group==1 || p2->mols[j]->group==2){
float x1 =(float) mols[i]->x;
float y1 = (float)mols[i]->y;
float z1 = (float)mols[i]->z;
v1[i*9+0]=x1;
v1[i*9+1]=y1;
v1[i*9+2]=z1;

float bx1 =(float) mols[i]->bx;
float by1 =(float) mols[i]->by;
float bz1 =(float) mols[i]->bz;
v1[i*9+3]=bx1;
v1[i*9+4]=by1;
v1[i*9+5]=bz1;

float b2x1 =(float) mols[i]->b2x;
float b2y1 =(float) mols[i]->b2y;
float b2z1 =(float) mols[i]->b2z;
v1[i*9+6]=b2x1;
v1[i*9+7]=b2y1;
v1[i*9+8]=b2z1;

float x2 =(float) p2->mols[i]->x;
float y2 =(float) p2->mols[i]->y;
float z2 =(float) p2->mols[i]->z;
v2[i*9+0]=x2;
v2[i*9+1]=y2;
v2[i*9+2]=z2;

float bx2 =(float) p2->mols[i]->bx;
float by2 =(float) p2->mols[i]->by;
float bz2 =(float) p2->mols[i]->bz;
v2[i*9+3]=bx2;
v2[i*9+4]=by2;
v2[i*9+5]=bz2;

float b2x2 =(float) p2->mols[i]->b2x;
float b2y2 =(float) p2->mols[i]->b2y;
float b2z2 =(float) p2->mols[i]->b2z;
v2[i*9+6]=b2x2;
v2[i*9+7]=b2y2;
v2[i*9+8]=b2z2;
i++;
}
}
//cout << "Groups Numbers: " << i << endl;
d=rmsd(v1,v2,i,mtx);
  free(v1);
  free(v2);
free(mtx);
}
if(d > 60){
//for(int j=0; j<num_mols; j++){cout << v1[(3*j)] << "," << v1[(3*j+1)] << "," << v1[(3*j+2)] << endl;}
//for(int j=0; j<num_mols; j++){cout << v2[(3*j)] << "," << v2[(3*j+1)] << "," << v2[(3*j+2)] << endl;}
}
return d;

}
double Polymere::moldist(int i,int j){
double mdist;
mdist =pow((mols[i]->x - mols[j]->x),2);
mdist =mdist + pow((mols[i]->y - mols[j]->y),2);
mdist =mdist + pow((mols[i]->z - mols[j]->z),2);
mdist=sqrt(mdist);
return mdist;
}

double Polymere::hdist(int i,int j){
double dist=0;
vector <Points> p;
vector <Points> p2;
Points c1,b1,b2;
int type1,type2;
//Returns a score for hydrogen bonding
c1.x=mols[i]->x;
c1.y=mols[i]->y;
c1.z=mols[i]->z;
b1.x=mols[i]->bx;
b1.y=mols[i]->by;
b1.z=mols[i]->bz;
b2.x=mols[i]->b2x;
b2.y=mols[i]->b2y;
b2.z=mols[i]->b2z;
p.push_back(c1);
p.push_back(b1);
p.push_back(b2);

c1.x=mols[j]->x;
c1.y=mols[j]->y;
c1.z=mols[j]->z;
b1.x=mols[j]->bx;
b1.y=mols[j]->by;
b1.z=mols[j]->bz;
b2.x=mols[j]->b2x;
b2.y=mols[j]->b2y;
b2.z=mols[j]->b2z;
p2.push_back(c1);
p2.push_back(b1);
p2.push_back(b2);

if(mols[i]->seq=='G'){type1=1;}else{
if(mols[i]->seq=='A'){type1=2;}else{
if(mols[i]->seq=='U'){type1=3;}else{
if(mols[i]->seq=='C'){type1=4;}else{
}
}
}
}
if(mols[j]->seq=='G'){type2=1;}else{
if(mols[j]->seq=='A'){type2=2;}else{
if(mols[j]->seq=='U'){type2=3;}else{
if(mols[j]->seq=='C'){type2=4;}else{
}
}
}
}
FABase * base1 = new FABase(type1);  
FABase * base2 = new FABase(type2);  
//Move to coordinate system
base1->alignBase(p);
base2->alignBase(p2);
//Calculate minimal pairwise distance between all atoms
return dist;
}
float Polymere::clusterDist(std::vector<Polymere *> centers){
float mindist=999;
for(int p=0; p < centers.size(); p++){


float dist = rms(centers[p]);
if(dist < mindist){
mindist=dist;
}
}
return mindist;
}

double Polymere::checkhbonds(){
//Checks if secondary structure is preserved
double out=1;
double mdist=0;
int j;
for(int i=0; i<num_mols; i++){
if(mols[i]->id >= 0){
if(mols[i]->id > i){
j=mols[i]->id;
mdist =pow((mols[i]->x - mols[j]->x),2);
mdist =mdist + pow((mols[i]->y - mols[j]->y),2);
mdist =mdist + pow((mols[i]->z - mols[j]->z),2);
//cout << i << ":" << j << ":" << mdist << endl;
if(mdist < 84){
out=0;
}
if(mdist > 144){
out=0;
}
float hbonds=hydroP(i,j);
//float hbonds=0;
if(hbonds < 1){
out=0;
}
}
}
}
return out;
}
double Polymere::checkbps(){
//Checks if secondary structure is preserved
//Returns 0 if bad, 1 if good
double out=0;
double mind=99999;
double xind=0;
double mdist;
//Check Base-Pair Distances
int j;
//cout << "Checking Base-Pairs" << endl;
for(int i=0; i<num_mols; i++){
//cout << "pos " << i <<" bp: " << mols[i]->id << endl;
if(mols[i]->id >= 0){
if(mols[i]->id > i){
j=mols[i]->id;
mdist =pow((mols[i]->x - mols[j]->x),2);
mdist =mdist + pow((mols[i]->y - mols[j]->y),2);
mdist =mdist + pow((mols[i]->z - mols[j]->z),2);
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
//mind = minimum distance between base-pairs 
//xind = maximum distance between base-pairs
if(mind < 81 || xind > 144){
out=0;
}else{
out=1;
//check all pairwise distances in movement region for c1 mindist
}
return out;
}
double Polymere::rmsAllg(const Polymere * p2,int movegroup){
double d;
float *v1,*v2;
float *mtx;
v1=NULL;
v2=NULL;
mtx=NULL;
mtx=(float*) malloc(16*sizeof(float));
  v1 =(float*) malloc(num_mols * 3 * sizeof(float));
  v2 =(float*) malloc(num_mols * 3 * sizeof(float));
if(v1==NULL || v2==NULL || mtx==NULL){
d=-1;
cout << "Allocation Failed" << endl;
}else{
int i=0;
for(int j=0; j<num_mols; j++){
if(mols[j]->group==movegroup){
float x1 =(float) mols[j]->x;
float y1 = (float)mols[j]->y;
float z1 = (float)mols[j]->z;
v1[i*3+0]=x1;
v1[i*3+1]=y1;
v1[i*3+2]=z1;

float x2 =(float) p2->mols[j]->x;
float y2 =(float) p2->mols[j]->y;
float z2 =(float) p2->mols[j]->z;
v2[i*3+0]=x2;
v2[i*3+1]=y2;
v2[i*3+2]=z2;
i++;
}
}
//cout << "Groups Numbers: " << i << endl;
d=(double)rmsd(v1,v2,i,mtx);
  free(v1);
  free(v2);
free(mtx);
v1=NULL;
v2=NULL;
mtx=NULL;
}
return d;
}

double Polymere::rmsAllc(const Polymere * p2){
double d;
float *v1,*v2;
float *mtx;
v1=NULL;
v2=NULL;
mtx=NULL;
mtx=(float*) malloc(16*sizeof(float));
  v1 =(float*) malloc(num_mols * 3 * sizeof(float));
  v2 =(float*) malloc(num_mols * 3 * sizeof(float));
if(v1==NULL || v2==NULL || mtx==NULL){
d=-1;
cout << "Allocation Failed" << endl;
}else{
int i=0;
for(int j=0; j<num_mols; j++){
float x1 =(float) mols[i]->x;
float y1 = (float)mols[i]->y;
float z1 = (float)mols[i]->z;
v1[i*3+0]=x1;
v1[i*3+1]=y1;
v1[i*3+2]=z1;

float x2 =(float) p2->mols[i]->x;
float y2 =(float) p2->mols[i]->y;
float z2 =(float) p2->mols[i]->z;
v2[i*3+0]=x2;
v2[i*3+1]=y2;
v2[i*3+2]=z2;
i++;
}
//cout << "Groups Numbers: " << i << endl;
d=(double)rmsd(v1,v2,i,mtx);
  free(v1);
  free(v2);
free(mtx);
v1=NULL;
v2=NULL;
mtx=NULL;
}
return d;
}
double Polymere::rmsAll(Polymere * p2){
double d;
float *v1,*v2;
float *mtx;
v1=NULL;
v2=NULL;
mtx=NULL;
mtx=(float*) malloc(16*sizeof(float));
  v1 =(float*) malloc(num_mols * 3 * sizeof(float));
  v2 =(float*) malloc(num_mols * 3 * sizeof(float));
if(v1==NULL || v2==NULL || mtx==NULL){
d=-1;
cout << "Allocation Failed" << endl;
}else{
int i=0;
for(int j=0; j<num_mols; j++){
float x1 =(float) mols[i]->x;
float y1 = (float)mols[i]->y;
float z1 = (float)mols[i]->z;
v1[i*3+0]=x1;
v1[i*3+1]=y1;
v1[i*3+2]=z1;

float x2 =(float) p2->mols[i]->x;
float y2 =(float) p2->mols[i]->y;
float z2 =(float) p2->mols[i]->z;
v2[i*3+0]=x2;
v2[i*3+1]=y2;
v2[i*3+2]=z2;
i++;
}
//cout << "Groups Numbers: " << i << endl;
d=(double)rmsd(v1,v2,i,mtx);
  free(v1);
  free(v2);
free(mtx);
v1=NULL;
v2=NULL;
mtx=NULL;
}
return d;
}



float Polymere::rms(Polymere * p2){
float d;
float *v1,*v2;
float *mtx;
v1=NULL;
v2=NULL;
mtx=NULL;
mtx=(float*) malloc(16*sizeof(float));
  v1 =(float*) malloc(num_mols * 3 * sizeof(float));
  v2 =(float*) malloc(num_mols * 3 * sizeof(float));
if(v1==NULL || v2==NULL || mtx==NULL){
d=-1;
cout << "Allocation Failed" << endl;
}else{
int i=0;
for(int j=0; j<num_mols; j++){
if(p2->mols[j]->group==1 || p2->mols[j]->group==2){
float x1 =(float) mols[i]->x;
float y1 = (float)mols[i]->y;
float z1 = (float)mols[i]->z;
v1[i*3+0]=x1;
v1[i*3+1]=y1;
v1[i*3+2]=z1;

float x2 =(float) p2->mols[i]->x;
float y2 =(float) p2->mols[i]->y;
float z2 =(float) p2->mols[i]->z;
v2[i*3+0]=x2;
v2[i*3+1]=y2;
v2[i*3+2]=z2;
i++;
}
}
//cout << "Groups Numbers: " << i << endl;
d=rmsd(v1,v2,i,mtx);
  free(v1);
  free(v2);
free(mtx);
}
if(d > 60){
//for(int j=0; j<num_mols; j++){cout << v1[(3*j)] << "," << v1[(3*j+1)] << "," << v1[(3*j+2)] << endl;}
//for(int j=0; j<num_mols; j++){cout << v2[(3*j)] << "," << v2[(3*j+1)] << "," << v2[(3*j+2)] << endl;}
}
return d;
}

void Polymere::localRMSD(vector<float> * rmsds,const Polymere * master,std::vector<int> * parent, std::vector<int> * child){
//Foreach residue in the structure look at all residues that are connected via base surface area graph
std::vector<float> fa1(full_atoms.begin(), full_atoms.end());
std::vector<float> fa2(master->full_atoms.begin(), master->full_atoms.end());
for(int nm=0; nm<num_mols; nm++){
//Foreach parent find the children and add to the RMSD vector
std::vector<int> resC;
resC.push_back(nm);
for(int pm=0; pm<parent->size(); pm++){
	if(nm==parent->at(pm)){resC.push_back(child->at(pm));}
}
/*
printf("Local: %d\n",nm);
for(int r=0; r<resC.size(); r++){
printf("%d ",resC[r]);
}
printf("\n");
*/
int msize=0;
for(int cm=0; cm<resC.size(); cm++){
int fsize=0;
if(master->mols[resC[cm]]->seq=='A'){fsize=22;}
if(master->mols[resC[cm]]->seq=='G'){fsize=23;}
if(master->mols[resC[cm]]->seq=='C'){fsize=21;}
if(master->mols[resC[cm]]->seq=='U'){fsize=20;}
if(resC[cm]==0){fsize=fsize-3;}
msize=msize+fsize;
}

//Initialize all vectors 
float d;
float *v1,*v2;
float *mtx;
v1=NULL;
v2=NULL;
mtx=NULL;
mtx=(float*) malloc(16*sizeof(float));
  v1 =(float*) malloc(msize * 3 * sizeof(float));
  v2 =(float*) malloc(msize * 3 * sizeof(float));
if(v1==NULL || v2==NULL || mtx==NULL){
d=-1;
cout << "Allocation Failed" << endl;
}else{
int i=0;
for(int cm=0; cm<resC.size(); cm++){
int fsize=0;
if(master->mols[resC[cm]]->seq=='A'){fsize=22;}
if(master->mols[resC[cm]]->seq=='G'){fsize=23;}
if(master->mols[resC[cm]]->seq=='C'){fsize=21;}
if(master->mols[resC[cm]]->seq=='U'){fsize=20;}

if(resC[cm]==0){fsize=fsize-3;}
//printf("Residue: %d\n",fsize);
int start =  getResidueStart(resC[cm]);
if(resC[cm]==0){start=9;}
int end = start + fsize;
int fi=start;
for(int fa=start; fa<end; fa++){
float x1 =fa1[fi+0];
float y1 =fa1[fi+1];
float z1 =fa1[fi+2];
//printf("%8.3f, %8.3f, %8.3f : ",x1,y1,z1);
v1[i*3+0]=x1;
v1[i*3+1]=y1;
v1[i*3+2]=z1;
float x2 =fa2[fi+0];
float y2 =fa2[fi+1];
float z2 =fa2[fi+2];
//printf("%8.3f, %8.3f, %8.3f\n",x2,y2,z2);
v2[i*3+0]=x2;
v2[i*3+1]=y2;
v2[i*3+2]=z2;
i++;
fi=fi+3;
}
}
//cout << "Groups Numbers: " << i << endl;
d=rmsd(v1,v2,i,mtx);
}
free(v1);
free(v2);
free(mtx);
rmsds->push_back(d);
//printf("rmsdl: %8.3f\n",d);
}
}

void Polymere::printVoro(){
Residue * Aresidue = new Residue;
// Set up the number of blocks that the container is divided into
const int n_x=20,n_y=20,n_z=20;
// Set the number of particles that are going to be randomly introduced
int particle_nums=0;
int i=0;
        double x,y,z,rad;
        // Create a container with the geometry given above, and make it
        // non-periodic in each of the three coordinates. Allocate space for
        // eight particles within each computational block

        // Randomly add particles into the container
double min_x=999;double max_x=-999;
double min_y=999;double max_y=-999;
double min_z=999;double max_z=-999;
vector<double> xps;
vector<double> yps;
vector<double> zps;

std::vector<int> res_ids;
std::vector<int> res_types;
std::vector<string> at_type;
std::vector<int> rtype;
for(int resn=0; resn<num_mols; resn++){
int bpos = resn*81;
if(mols[resn]->seq=='G'){
rtype.push_back(1);
for(int rn=0; rn < 27; rn++){
res_ids.push_back(resn);
if(rn < 12){res_types.push_back(-1);}else{res_types.push_back(1);}
string aname(Aresidue->atomnames1[rn]);
at_type.push_back(aname);
//create an atom x y z atom_name chain  res_name
int p1 = bpos+rn*3;
int p2 = bpos+rn*3+1;
int p3 = bpos+rn*3+2;
particle_nums++;
xps.push_back(full_atoms[p1]);
yps.push_back(full_atoms[p2]);
zps.push_back(full_atoms[p3]);

if(min_x > full_atoms[p1]){min_x=full_atoms[p1];}
if(min_y > full_atoms[p2]){min_y=full_atoms[p2];}
if(min_z > full_atoms[p3]){min_z=full_atoms[p3];}

if(max_x < full_atoms[p1]){max_x=full_atoms[p1];}
if(max_y < full_atoms[p2]){max_y=full_atoms[p2];}
if(max_z < full_atoms[p3]){max_z=full_atoms[p3];}

}
}else{
if(mols[resn]->seq=='A'){
rtype.push_back(2);
for(int rn=0; rn < 26; rn++){
res_ids.push_back(resn);
if(rn < 12){res_types.push_back(-2);}else{res_types.push_back(2);}
string aname(Aresidue->atomnames2[rn]);
at_type.push_back(aname);
//create an atom x y z atom_name chain  res_name
int p1 = bpos+rn*3;
int p2 = bpos+rn*3+1;
int p3 = bpos+rn*3+2;

particle_nums++;
xps.push_back(full_atoms[p1]);
yps.push_back(full_atoms[p2]);
zps.push_back(full_atoms[p3]);

if(min_x > full_atoms[p1]){min_x=full_atoms[p1];}
if(min_y > full_atoms[p2]){min_y=full_atoms[p2];}
if(min_z > full_atoms[p3]){min_z=full_atoms[p3];}

if(max_x < full_atoms[p1]){max_x=full_atoms[p1];}
if(max_y < full_atoms[p2]){max_y=full_atoms[p2];}
if(max_z < full_atoms[p3]){max_z=full_atoms[p3];}
}


}else{
if(mols[resn]->seq=='U'){
rtype.push_back(3);
for(int rn=0; rn < 23; rn++){
res_ids.push_back(resn);
if(rn < 12){res_types.push_back(-3);}else{res_types.push_back(3);}
string aname(Aresidue->atomnames3[rn]);
at_type.push_back(aname);
//create an atom x y z atom_name chain  res_name
int p1 = bpos+rn*3;
int p2 = bpos+rn*3+1;
int p3 = bpos+rn*3+2;

particle_nums++;
xps.push_back(full_atoms[p1]);
yps.push_back(full_atoms[p2]);
zps.push_back(full_atoms[p3]);

if(min_x > full_atoms[p1]){min_x=full_atoms[p1];}
if(min_y > full_atoms[p2]){min_y=full_atoms[p2];}
if(min_z > full_atoms[p3]){min_z=full_atoms[p3];}

if(max_x < full_atoms[p1]){max_x=full_atoms[p1];}
if(max_y < full_atoms[p2]){max_y=full_atoms[p2];}
if(max_z < full_atoms[p3]){max_z=full_atoms[p3];}
}

}else{
if(mols[resn]->seq=='C'){
rtype.push_back(4);
for(int rn=0; rn < 24; rn++){
res_ids.push_back(resn);
if(rn < 12){res_types.push_back(-4);}else{res_types.push_back(4);}
string aname(Aresidue->atomnames4[rn]);
at_type.push_back(aname);
//create an atom x y z atom_name chain  res_name
int p1 = bpos+rn*3;
int p2 = bpos+rn*3+1;
int p3 = bpos+rn*3+2;

particle_nums++;
xps.push_back(full_atoms[p1]);
yps.push_back(full_atoms[p2]);
zps.push_back(full_atoms[p3]);


if(min_x > full_atoms[p1]){min_x=full_atoms[p1];}
if(min_y > full_atoms[p2]){min_y=full_atoms[p2];}
if(min_z > full_atoms[p3]){min_z=full_atoms[p3];}

if(max_x < full_atoms[p1]){max_x=full_atoms[p1];}
if(max_y < full_atoms[p2]){max_y=full_atoms[p2];}
if(max_z < full_atoms[p3]){max_z=full_atoms[p3];}
}
}else{
}
}
}
}

}


min_x = min_x -9;
min_y = min_y -9;
min_z = min_z -9;

max_x = max_x +9;
max_y = max_y +9;
max_z = max_z +9;

//Voro container needs to be a cube
double minv = min_x;
if(minv > min_y){minv=min_y;}
if(minv > min_z){minv=min_z;}

double maxv = max_x;
if(maxv < max_y){maxv=max_y;}
if(maxv < max_z){maxv=max_z;}


const double x_min=minv;
const double x_max=maxv;
const double y_min=minv;
const double y_max=maxv;
const double z_min=minv;
const double z_max=maxv;

//printf("Ranges %8.3f - %8.3f : %8.3f - %8.3f : %8.3f - %8.3f\n",x_min,x_max,y_min,y_max,z_min,z_max);
//printf("Atom Information\n");
//for(int r=0; r<at_type.size(); r++){printf("%s:%d:%d\n",at_type[r].c_str(),res_types[r],res_ids[r]);}
container_poly con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,false,false,false,8);
//Residues are split by backbone and base atoms
for(int a = 0; a<at_type.size(); a++){
double x = xps[a];
double y = yps[a];
double z = zps[a];
double r=0.5;
if (at_type[a].find("H",0) != string::npos){
r=.31;
}else{
if (at_type[a].find("C",0) != string::npos){
r=.73;
}else{
if (at_type[a].find("O",0) != string::npos){
r=.66;
}else{
if (at_type[a].find("P",0) != string::npos){
r=1.07;
}else{
if (at_type[a].find("N",0) != string::npos){
r=0.71;
}else{
printf("Can't identify atom type: %s\n",at_type[a].c_str());
}
}
}
}
}
//printf("Adding: %8.3f %8.3f %8.3f %8.3f\n",x,y,z,r);
con.put(a,x,y,z,r);
}
con.draw_cells_rgl("voro.out",&res_ids);
}         

void Polymere::localgraph(vector<int> * parent, vector<int> * child){
Residue * Aresidue = new Residue;
// Set up the number of blocks that the container is divided into
const int n_x=20,n_y=20,n_z=20;
// Set the number of particles that are going to be randomly introduced
int particle_nums=0;
int i=0;
        double x,y,z,rad;
        // Create a container with the geometry given above, and make it
        // non-periodic in each of the three coordinates. Allocate space for
        // eight particles within each computational block

        // Randomly add particles into the container
double min_x=999;double max_x=-999;
double min_y=999;double max_y=-999;
double min_z=999;double max_z=-999;
vector<double> xps;
vector<double> yps;
vector<double> zps;

std::vector<int> res_ids;
std::vector<int> res_types;
std::vector<string> at_type;
std::vector<int> rtype;
for(int resn=0; resn<num_mols; resn++){
int bpos = resn*81;
if(mols[resn]->seq=='G'){
rtype.push_back(1);
for(int rn=0; rn < 27; rn++){
res_ids.push_back(resn);
if(rn < 12){res_types.push_back(-1);}else{res_types.push_back(1);}
string aname(Aresidue->atomnames1[rn]);
at_type.push_back(aname);
//create an atom x y z atom_name chain  res_name
int p1 = bpos+rn*3;
int p2 = bpos+rn*3+1;
int p3 = bpos+rn*3+2;
particle_nums++;
xps.push_back(full_atoms[p1]);
yps.push_back(full_atoms[p2]);
zps.push_back(full_atoms[p3]);

if(min_x > full_atoms[p1]){min_x=full_atoms[p1];}
if(min_y > full_atoms[p2]){min_y=full_atoms[p2];}
if(min_z > full_atoms[p3]){min_z=full_atoms[p3];}

if(max_x < full_atoms[p1]){max_x=full_atoms[p1];}
if(max_y < full_atoms[p2]){max_y=full_atoms[p2];}
if(max_z < full_atoms[p3]){max_z=full_atoms[p3];}

}
}else{
if(mols[resn]->seq=='A'){
rtype.push_back(2);
for(int rn=0; rn < 26; rn++){
res_ids.push_back(resn);
if(rn < 12){res_types.push_back(-2);}else{res_types.push_back(2);}
string aname(Aresidue->atomnames2[rn]);
at_type.push_back(aname);
//create an atom x y z atom_name chain  res_name
int p1 = bpos+rn*3;
int p2 = bpos+rn*3+1;
int p3 = bpos+rn*3+2;

particle_nums++;
xps.push_back(full_atoms[p1]);
yps.push_back(full_atoms[p2]);
zps.push_back(full_atoms[p3]);

if(min_x > full_atoms[p1]){min_x=full_atoms[p1];}
if(min_y > full_atoms[p2]){min_y=full_atoms[p2];}
if(min_z > full_atoms[p3]){min_z=full_atoms[p3];}

if(max_x < full_atoms[p1]){max_x=full_atoms[p1];}
if(max_y < full_atoms[p2]){max_y=full_atoms[p2];}
if(max_z < full_atoms[p3]){max_z=full_atoms[p3];}
}


}else{
if(mols[resn]->seq=='U'){
rtype.push_back(3);
for(int rn=0; rn < 23; rn++){
res_ids.push_back(resn);
if(rn < 12){res_types.push_back(-3);}else{res_types.push_back(3);}
string aname(Aresidue->atomnames3[rn]);
at_type.push_back(aname);
//create an atom x y z atom_name chain  res_name
int p1 = bpos+rn*3;
int p2 = bpos+rn*3+1;
int p3 = bpos+rn*3+2;

particle_nums++;
xps.push_back(full_atoms[p1]);
yps.push_back(full_atoms[p2]);
zps.push_back(full_atoms[p3]);

if(min_x > full_atoms[p1]){min_x=full_atoms[p1];}
if(min_y > full_atoms[p2]){min_y=full_atoms[p2];}
if(min_z > full_atoms[p3]){min_z=full_atoms[p3];}

if(max_x < full_atoms[p1]){max_x=full_atoms[p1];}
if(max_y < full_atoms[p2]){max_y=full_atoms[p2];}
if(max_z < full_atoms[p3]){max_z=full_atoms[p3];}
}

}else{
if(mols[resn]->seq=='C'){
rtype.push_back(4);
for(int rn=0; rn < 24; rn++){
res_ids.push_back(resn);
if(rn < 12){res_types.push_back(-4);}else{res_types.push_back(4);}
string aname(Aresidue->atomnames4[rn]);
at_type.push_back(aname);
//create an atom x y z atom_name chain  res_name
int p1 = bpos+rn*3;
int p2 = bpos+rn*3+1;
int p3 = bpos+rn*3+2;

particle_nums++;
xps.push_back(full_atoms[p1]);
yps.push_back(full_atoms[p2]);
zps.push_back(full_atoms[p3]);


if(min_x > full_atoms[p1]){min_x=full_atoms[p1];}
if(min_y > full_atoms[p2]){min_y=full_atoms[p2];}
if(min_z > full_atoms[p3]){min_z=full_atoms[p3];}

if(max_x < full_atoms[p1]){max_x=full_atoms[p1];}
if(max_y < full_atoms[p2]){max_y=full_atoms[p2];}
if(max_z < full_atoms[p3]){max_z=full_atoms[p3];}
}
}else{
}
}
}
}

}


min_x = min_x -9;
min_y = min_y -9;
min_z = min_z -9;

max_x = max_x +9;
max_y = max_y +9;
max_z = max_z +9;

//Voro container needs to be a cube
double minv = min_x;
if(minv > min_y){minv=min_y;}
if(minv > min_z){minv=min_z;}

double maxv = max_x;
if(maxv < max_y){maxv=max_y;}
if(maxv < max_z){maxv=max_z;}


const double x_min=minv;
const double x_max=maxv;
const double y_min=minv;
const double y_max=maxv;
const double z_min=minv;
const double z_max=maxv;

//printf("Ranges %8.3f - %8.3f : %8.3f - %8.3f : %8.3f - %8.3f\n",x_min,x_max,y_min,y_max,z_min,z_max);
//printf("Atom Information\n");
//for(int r=0; r<at_type.size(); r++){printf("%s:%d:%d\n",at_type[r].c_str(),res_types[r],res_ids[r]);}
container_poly con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,false,false,false,8);
//Residues are split by backbone and base atoms
for(int a = 0; a<at_type.size(); a++){
double x = xps[a];
double y = yps[a];
double z = zps[a];
double r=0.5;
if (at_type[a].find("H",0) != string::npos){
r=.31;
}else{
if (at_type[a].find("C",0) != string::npos){
r=.73;
}else{
if (at_type[a].find("O",0) != string::npos){
r=.66;
}else{
if (at_type[a].find("P",0) != string::npos){
r=1.07;
}else{
if (at_type[a].find("N",0) != string::npos){
r=0.71;
}else{
printf("Can't identify atom type: %s\n",at_type[a].c_str());
}
}
}
}
}
//printf("Adding: %8.3f %8.3f %8.3f %8.3f\n",x,y,z,r);
con.put(a,x,y,z,r);
}

vector<fpoint> volumes;
vector<fpoint> surface;
vector<fpoint> surfA;
vector<fpoint>surfG;
vector<fpoint>surfC;
vector<fpoint>surfU;
vector<int> orderA;
vector<int> orderG;
vector<int> orderC;
vector<int> orderU;
vector<int> unpacked;
vector<fpoint> bvolumes;
vector<fpoint> bsurface;
vector<fpoint> bsurfA;
vector<fpoint> bsurfG;
vector<fpoint> bsurfC;
vector<fpoint> bsurfU;
vector<int> borderA;
vector<int> borderG;
vector<int> borderC;
vector<int> borderU;
vector<int> bunpacked;
vector<int> order;

vector<int> atom0,atom1,atom2,atom3,atom4,atom5,atom6,atom7,atom8,atom9,atom10,atom11,atom12,atom13,atom14,atom15,atom16,atom17,atom18,atom19,atom20,atom21,atom22,atom23,atom24,atom25,atom26;
vector<fpoint> vatom0,vatom1,vatom2,vatom3,vatom4,vatom5,vatom6,vatom7,vatom8,vatom9,vatom10,vatom11,vatom12,vatom13,vatom14,vatom15,vatom16,vatom17,vatom18,vatom19,vatom20,vatom21,vatom22,vatom23,vatom24,vatom25,vatom26;

con.get_res_volume(
	&res_types,
	&res_ids,
	&at_type,
	&volumes,
	&unpacked,
	&surface,
	&surfA,
	&surfG,
	&surfC,
	&surfU,
	&orderA,
	&orderG,
	&orderC,
	&orderU,
	
	&bvolumes,
	&bunpacked,
	&bsurface,
	&bsurfA,
	&bsurfG,
	&bsurfC,
	&bsurfU,
	&borderA,
	&borderG,
	&borderC,
	&borderU,
	
	&vatom0,
	&vatom1,
	&vatom2,
	&vatom3,
	&vatom4,
	&vatom5,
	&vatom6,
	&vatom7,
	&vatom8,
	&vatom9,
	&vatom10,
	&vatom11,
	&vatom12,
	&vatom13,
	&vatom14,
	&vatom15,
	&vatom16,
	&vatom17,
	&vatom18,
	&vatom19,
	&vatom20,
	&vatom21,
	&vatom22,
	&vatom23,
	&vatom24,
	&vatom25,
	&vatom26,

	&atom0,
	&atom1,
	&atom2,
	&atom3,
	&atom4,
	&atom5,
	&atom6,
	&atom7,
	&atom8,
	&atom9,
	&atom10,
	&atom11,
	&atom12,
	&atom13,
	&atom14,
	&atom15,
	&atom16,
	&atom17,
	&atom18,
	&atom19,
	&atom20,
	&atom21,
	&atom22,
	&atom23,
	&atom24,
	&atom25,
	&atom26,
	
	parent,
	child,
	&order
);
}         

void Polymere::lowres(Score * score,std::vector<float> * resscores){
//Guassian distributions were determined from small_rna volume distributions
//Scores for graph orders based on negative log of frequency
std::vector<float> gorder;
gorder.push_back(5.712800);
gorder.push_back(4.273583);
gorder.push_back(3.585604);
gorder.push_back(2.778712);
gorder.push_back(2.276711);
gorder.push_back(1.422936);
gorder.push_back(1.350809);
gorder.push_back(1.874283);
gorder.push_back(2.540195);
gorder.push_back(3.213504);
gorder.push_back(4.504488);
gorder.push_back(5.552457);
gorder.push_back(6.283345);
gorder.push_back(7.056534);
gorder.push_back(9.541441);



for(int i=0; i<score->s1.size(); i++){
//No hydrogens
float u=0;
float sig=1;
if(score->s5[i]==1){u=152.6237;sig=32.46298;}
if(score->s5[i]==2){u=141.3321;sig=36.57051;}
if(score->s5[i]==3){u=110.2345;sig=27.54622;}
if(score->s5[i]==4){u=112.3131;sig=27.03934;}
float x = score->s6[i];
float gaussian=-1*log(1/(sqrt(2*PI*pow(sig,2)))*exp(-1*pow((x-u),2)/(2*pow(sig,2))));
float goscore=0;
if(score->s0[i]<15){
goscore=gorder[score->s0[i]];
}else{
goscore=10;
}
resscores->push_back((gaussian+goscore-5*score->s82[i]));
}
}


int Polymere::atomClashes(){
//returns the total number of heavy atom clashes
//pairs of heavy atoms within 1.483 Angstroms of each other
int uatoms[23]={1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0};
int gatoms[27]={1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0};
int aatoms[26]={1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0};
int catoms[24]={1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0};
int bpos1=0;
int bpos2=0;
int m1tot=0;
int m2tot=0;
int atom1=0;
int atom2=0;
float pdist=0;
int clashes=0;
float x1,x2,y1,y2,z1,z2,dist;
for(int m=0; m<(num_mols-1); m++){
bpos1=m*81;
m1tot=0;
if(mols[m]->seq == 'G'){m1tot=27;}
if(mols[m]->seq == 'A'){m1tot=26;}
if(mols[m]->seq == 'U'){m1tot=23;}
if(mols[m]->seq == 'C'){m1tot=24;}

for(int m2=(m+1); m2<(num_mols); m2++){
bpos2=m2*81;
m2tot=0;
pdist = (pow((mols[m]->x+mols[m2]->x),2) + pow((mols[m]->y+mols[m2]->y),2) + pow((mols[m]->z+mols[m2]->z),2));
if(pdist < 169){  //If bases C1' atoms are within 13 Angstroms check for clashes
if(mols[m2]->seq == 'G'){m2tot=27;}
if(mols[m2]->seq == 'A'){m2tot=26;}
if(mols[m2]->seq == 'U'){m2tot=23;}
if(mols[m2]->seq == 'C'){m2tot=24;}

for(int i=0; i<m1tot; i++){
atom1=0;
if(m1tot==27){atom1=gatoms[i];}
if(m1tot==23){atom1=uatoms[i];}
if(m1tot==26){atom1=aatoms[i];}
if(m1tot==24){atom1=catoms[i];}
 x1=full_atoms[(bpos1+i*3)];
 y1=full_atoms[(bpos1+i*3+1)];
 z1=full_atoms[(bpos1+i*3+2)];
for(int j=0; j<m2tot; j++){
atom2=0;
if(m2tot==27){atom2=gatoms[j];}
if(m2tot==23){atom2=uatoms[j];}
if(m2tot==26){atom2=aatoms[j];}
if(m2tot==24){atom2=catoms[j];}
x2=full_atoms[(bpos2+j*3)];
y2=full_atoms[(bpos2+j*3+1)];
z2=full_atoms[(bpos2+j*3+2)];
if((m==0 && i > 2) || m > 0){
if(atom1==1 && atom2==1){//Only look at heavy atom clashes
dist=pow((x1-x2),2)+pow((y1-y2),2)+pow((z1-z2),2);
//This was the minimal distance^2 observed between any two heavy atoms in a structure
if(dist < 4){
//printf("Clash %d-%d, %8.3f: Atoms: %d %d\n", m,m2,dist,i,j);
clashes++;
}
}
else{
//Hydrogen/Heavy Atom 1.09Ang
if((atom1+atom2)==1){
if(dist < 1.21){
clashes++;
}
}else{
if(atom1==0 && atom2==0){//Hydrogen-Hydrogen distance = .74Ang
if(dist < .64){
clashes++;
}
}
}
}





}
}
}

} 
}

}

return clashes;
}
/*
void Polymere::pairwiseDist(){
char uatoms[23][1]={'P','O','O','O','C','C','O','C','O','C','O','C','N','C','C','O','N','C','O','C','H','H','H'};
char gatoms[27][1]={'P','O','O','O','C','C','O','C','O','C','O','C','N','C','N','C','N','N','C','O','C','N','C','H','H','H','H'};
char aatoms[26][1]={'P','O','O','O','C','C','O','C','O','C','O','C','N','C','N','C','N','C','N','C','N','C','H','H','H','H'};
char catoms[24][1]={'P','O','O','O','C','C','O','C','O','C','O','C','N','C','C','O','N','C','N','C','H','H','H','H'};

//Counts the number of backbone P-O bonds that are 
// < 1.4 and > 1.8 Angstroms in distance
 
int counts=0;
for(int m=0; m<(num_mols-1); m++){
int bpos2=(m+1)*81;
int bpos1=m*81+3*8;
float x1=full_atoms[bpos1];
float y1=full_atoms[bpos1+1];
float z1=full_atoms[bpos1+2];

float x2=full_atoms[bpos2];
float y2=full_atoms[bpos2+1];
float z2=full_atoms[bpos2+2];
float dist=sqrt(pow((x1-x2),2)+pow((y1-y2),2)+pow((z1-z2),2));
if(dist < 1.4 || dist > 1.8){
counts++;
}


//Below is for creating the pair-wise distance scores

printf("BPO: %8.3f\n",dist);
if(dist < 30){
//All other pairwise distances
for(int m2=(m+1); m2<(num_mols); m2++){
int bpos1=m*81;
int bpos2=m2*81;
int m1tot=0;
int m2tot=0;

if(mols[m]->seq == 'G'){m1tot=27;}
if(mols[m]->seq == 'A'){m1tot=26;}
if(mols[m]->seq == 'U'){m1tot=23;}
if(mols[m]->seq == 'C'){m1tot=24;}

if(mols[m2]->seq == 'G'){m2tot=27;}
if(mols[m2]->seq == 'A'){m2tot=26;}
if(mols[m2]->seq == 'U'){m2tot=23;}
if(mols[m2]->seq == 'C'){m2tot=24;}

for(int i=0; i<m1tot; i++){
char atom1;
if(m1tot==27){atom1=gatoms[i][1];}
if(m1tot==23){atom1=uatoms[i][1];}
if(m1tot==26){atom1=aatoms[i][1];}
if(m1tot==24){atom1=catoms[i][1];}
float x1=full_atoms[(bpos1+i*3)];
float y1=full_atoms[(bpos1+i*3+1)];
float z1=full_atoms[(bpos1+i*3+2)];
for(int j=0; j<m2tot; j++){
char atom2;
if(m1tot==27){atom2=gatoms[j][1];}
if(m1tot==23){atom2=uatoms[j][1];}
if(m1tot==26){atom2=aatoms[j][1];}
if(m1tot==24){atom2=catoms[j][1];}
float x2=full_atoms[(bpos2+j*3)];
float y2=full_atoms[(bpos2+j*3+1)];
float z2=full_atoms[(bpos2+j*3+2)];
float dist=sqrt(pow((x1-x2),2)+pow((y1-y2),2)+pow((z1-z2),2));
if(dist < 2.5 & dist > 0){
printf("PW:%c:%c:%8.3f\n",atom1,atom2,dist); 
}
}
}

} 
}

}
printf("Total Bad Backbone Bonds: %d\n",counts);
return 0;



}
*/
///////////////////////////////////////////////////////////////
//   Creating a set of checks and scores that work on a 
//   particular group of the polymer.
///////////////////////////////////////////////////////////////
int Polymere::groupCheckClashes(int group_id){
//returns the total number of heavy atom clashes
//pairs of heavy atoms within 1.483 Angstroms of each other
int uatoms[23]={1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0};
int gatoms[27]={1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0};
int aatoms[26]={1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0};
int catoms[24]={1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0};

int bpos1=0;
int bpos2=0;
int m1tot=0;
int m2tot=0;
int atom1=0;
int atom2=0;
float pdist=0;
int clashes=0;
float x1,x2,y1,y2,z1,z2,dist;
for(int m=0; m<this->mols.size(); m++){
if(mols[m]->group==group_id){
bpos1=m*81;
m1tot=27;
if(mols[m]->seq == 'G'){m1tot=27;}
if(mols[m]->seq == 'A'){m1tot=26;}
if(mols[m]->seq == 'U'){m1tot=23;}
if(mols[m]->seq == 'C'){m1tot=24;}

for(int m2=(m+1); m2<mols.size(); m2++){
if(mols[m2]->group==group_id){
bpos2=m2*81;
m2tot=27;
pdist = (pow((mols[m]->x-mols[m2]->x),2) + pow((mols[m]->y-mols[m2]->y),2) + pow((mols[m]->z-mols[m2]->z),2));
if(pdist < 169){  //If bases C1' atoms are within 13 Angstroms check for clashes
if(mols[m2]->seq == 'G'){m2tot=27;}
if(mols[m2]->seq == 'A'){m2tot=26;}
if(mols[m2]->seq == 'U'){m2tot=23;}
if(mols[m2]->seq == 'C'){m2tot=24;}

for(int i=0; i<m1tot; i++){
atom1=1;
if(m1tot==27){atom1=gatoms[i];}
if(m1tot==23){atom1=uatoms[i];}
if(m1tot==26){atom1=aatoms[i];}
if(m1tot==24){atom1=catoms[i];}
 x1=full_atoms[(bpos1+i*3)];
 y1=full_atoms[(bpos1+i*3+1)];
 z1=full_atoms[(bpos1+i*3+2)];
for(int j=0; j<m2tot; j++){
atom2=1;
if(m2tot==27){atom2=gatoms[j];}
if(m2tot==23){atom2=uatoms[j];}
if(m2tot==26){atom2=aatoms[j];}
if(m2tot==24){atom2=catoms[j];}
x2=full_atoms[(bpos2+j*3)];
y2=full_atoms[(bpos2+j*3+1)];
z2=full_atoms[(bpos2+j*3+2)];
if((m==0 && i > 2) || m > 0){
dist=pow((x1-x2),2)+pow((y1-y2),2)+pow((z1-z2),2);
if(atom1==1 && atom2==1){//Only look at heavy atom clashes
//This was the minimal distance^2 observed between any two heavy atoms in a structure
if(dist < 2.5){//Heavy-Heavy = 1.48
//printf("Clash %d-%d, %8.3f: Atoms: %d %d\n", m,m2,dist,i,j);
clashes++;
}
}
else{
//Hydrogen/Heavy Atom 1.09Ang
if((atom1+atom2)==1){
if(dist < 1.21){
clashes++;
}
}else{
if(atom1==0 && atom2==0){//Hydrogen-Hydrogen distance = .74Ang
if(dist < .64){
clashes++;
}
}
}

}
}
}
}

} 
}
}


}

}

return clashes;
}
int Polymere::groupCheckBackBone(int group_id,float p_o_lower, float p_o_upper){
//Counts the number of backbone P-O bonds that are 
// < 1.4 and > 1.8 Angstroms in distance
int counts=0;
int bpos2,bpos1;
float x1,y1,z1,x2,y2,z2,dist;
for(int m=0; m<(num_mols-1); m++){
if(mols[m]->group==group_id && mols[(m+1)]->group==group_id){
bpos2=(m+1)*81;
bpos1=m*81+3*8;
x1=full_atoms[bpos1];
y1=full_atoms[bpos1+1];
z1=full_atoms[bpos1+2];

x2=full_atoms[bpos2];
y2=full_atoms[bpos2+1];
z2=full_atoms[bpos2+2];
dist=(pow((x1-x2),2)+pow((y1-y2),2)+pow((z1-z2),2));
if(dist < p_o_lower || dist > p_o_upper){ //Keep P-O3' bond length 1.4 - 1.952
//if(dist < 1.6 || dist > 4){ //Keep P-O3' bond length 1.4 - 2.0
//if(dist < 1 || dist > 16){
counts++;
//printf("BackBone Position: %d\n",m);
}
}
}
return counts;
}

int Polymere::checkBackBone(){
//Counts the number of backbone P-O bonds that are 
// < 1.4 and > 1.8 Angstroms in distance
int counts=0;
int bpos2,bpos1;
float x1,y1,z1,x2,y2,z2,dist;
for(int m=0; m<(num_mols-1); m++){
bpos2=(m+1)*81;
bpos1=m*81+3*8;
x1=full_atoms[bpos1];
y1=full_atoms[bpos1+1];
z1=full_atoms[bpos1+2];

x2=full_atoms[bpos2];
y2=full_atoms[bpos2+1];
z2=full_atoms[bpos2+2];
dist=(pow((x1-x2),2)+pow((y1-y2),2)+pow((z1-z2),2));
//if(dist < 1.96 || dist > 3.81){ //Keep P-O3' bond length 1.4 - 1.952
if(dist < 1.6 || dist > 4){ //Keep P-O3' bond length 1.4 - 2.0
//printf("  dist: %8.3f\n", dist);
counts++;
}
}
return counts;
}

void Polymere::group_resscore(Score * score,vector<double> * dparams,int group_id){
vector<float> hydros;
vector<float> doublets;
//Residue Level Scores
doubletVg(dparams,3.0,3.0,3.0,0.2,0.2,&doublets,group_id);
hydroVg(&hydros,group_id);
//Phosphate to O3` distances


Residue * Aresidue = new Residue;
// Set up the number of blocks that the container is divided into
const int n_x=20,n_y=20,n_z=20;
// Set the number of particles that are going to be randomly introduced
int particle_nums=0;
int i=0;
        double x,y,z,rad;
        // Create a container with the geometry given above, and make it
        // non-periodic in each of the three coordinates. Allocate space for
        // eight particles within each computational block

        // Randomly add particles into the container
double min_x=999;double max_x=-999;
double min_y=999;double max_y=-999;
double min_z=999;double max_z=-999;
vector<double> xps;
vector<double> yps;
vector<double> zps;

std::vector<int> res_ids;
std::vector<int> res_types;
std::vector<string> at_type;
std::vector<int> rtype;
for(int resn=0; resn<num_mols; resn++){
if(mols[resn]->group==group_id){
int bpos = resn*81;
if(mols[resn]->seq=='G'){
rtype.push_back(1);
for(int rn=0; rn < 23; rn++){
res_ids.push_back(resn);
if(rn < 12){res_types.push_back(-1);}else{res_types.push_back(1);}
string aname(Aresidue->atomnames1[rn]);
at_type.push_back(aname);
//create an atom x y z atom_name chain  res_name
int p1 = bpos+rn*3;
int p2 = bpos+rn*3+1;
int p3 = bpos+rn*3+2;
particle_nums++;
xps.push_back(full_atoms[p1]);
yps.push_back(full_atoms[p2]);
zps.push_back(full_atoms[p3]);

if(min_x > full_atoms[p1]){min_x=full_atoms[p1];}
if(min_y > full_atoms[p2]){min_y=full_atoms[p2];}
if(min_z > full_atoms[p3]){min_z=full_atoms[p3];}

if(max_x < full_atoms[p1]){max_x=full_atoms[p1];}
if(max_y < full_atoms[p2]){max_y=full_atoms[p2];}
if(max_z < full_atoms[p3]){max_z=full_atoms[p3];}

}
}else{
if(mols[resn]->seq=='A'){
rtype.push_back(2);
for(int rn=0; rn < 22; rn++){
res_ids.push_back(resn);
if(rn < 12){res_types.push_back(-2);}else{res_types.push_back(2);}
string aname(Aresidue->atomnames2[rn]);
at_type.push_back(aname);
//create an atom x y z atom_name chain  res_name
int p1 = bpos+rn*3;
int p2 = bpos+rn*3+1;
int p3 = bpos+rn*3+2;

particle_nums++;
xps.push_back(full_atoms[p1]);
yps.push_back(full_atoms[p2]);
zps.push_back(full_atoms[p3]);

if(min_x > full_atoms[p1]){min_x=full_atoms[p1];}
if(min_y > full_atoms[p2]){min_y=full_atoms[p2];}
if(min_z > full_atoms[p3]){min_z=full_atoms[p3];}

if(max_x < full_atoms[p1]){max_x=full_atoms[p1];}
if(max_y < full_atoms[p2]){max_y=full_atoms[p2];}
if(max_z < full_atoms[p3]){max_z=full_atoms[p3];}
}


}else{
if(mols[resn]->seq=='U'){
rtype.push_back(3);
for(int rn=0; rn < 20; rn++){
res_ids.push_back(resn);
if(rn < 12){res_types.push_back(-3);}else{res_types.push_back(3);}
string aname(Aresidue->atomnames3[rn]);
at_type.push_back(aname);
//create an atom x y z atom_name chain  res_name
int p1 = bpos+rn*3;
int p2 = bpos+rn*3+1;
int p3 = bpos+rn*3+2;

particle_nums++;
xps.push_back(full_atoms[p1]);
yps.push_back(full_atoms[p2]);
zps.push_back(full_atoms[p3]);

if(min_x > full_atoms[p1]){min_x=full_atoms[p1];}
if(min_y > full_atoms[p2]){min_y=full_atoms[p2];}
if(min_z > full_atoms[p3]){min_z=full_atoms[p3];}

if(max_x < full_atoms[p1]){max_x=full_atoms[p1];}
if(max_y < full_atoms[p2]){max_y=full_atoms[p2];}
if(max_z < full_atoms[p3]){max_z=full_atoms[p3];}
}

}else{
if(mols[resn]->seq=='C'){
rtype.push_back(4);
for(int rn=0; rn < 20; rn++){
res_ids.push_back(resn);
if(rn < 12){res_types.push_back(-4);}else{res_types.push_back(4);}
string aname(Aresidue->atomnames4[rn]);
at_type.push_back(aname);
//create an atom x y z atom_name chain  res_name
int p1 = bpos+rn*3;
int p2 = bpos+rn*3+1;
int p3 = bpos+rn*3+2;

particle_nums++;
xps.push_back(full_atoms[p1]);
yps.push_back(full_atoms[p2]);
zps.push_back(full_atoms[p3]);


if(min_x > full_atoms[p1]){min_x=full_atoms[p1];}
if(min_y > full_atoms[p2]){min_y=full_atoms[p2];}
if(min_z > full_atoms[p3]){min_z=full_atoms[p3];}

if(max_x < full_atoms[p1]){max_x=full_atoms[p1];}
if(max_y < full_atoms[p2]){max_y=full_atoms[p2];}
if(max_z < full_atoms[p3]){max_z=full_atoms[p3];}
}
}else{
}
}
}
}

}
}

min_x = min_x -9;
min_y = min_y -9;
min_z = min_z -9;

max_x = max_x +9;
max_y = max_y +9;
max_z = max_z +9;

//Voro container needs to be a cube
double minv = min_x;
if(minv > min_y){minv=min_y;}
if(minv > min_z){minv=min_z;}

double maxv = max_x;
if(maxv < max_y){maxv=max_y;}
if(maxv < max_z){maxv=max_z;}


const double x_min=minv;
const double x_max=maxv;
const double y_min=minv;
const double y_max=maxv;
const double z_min=minv;
const double z_max=maxv;

//printf("Ranges %8.3f - %8.3f : %8.3f - %8.3f : %8.3f - %8.3f\n",x_min,x_max,y_min,y_max,z_min,z_max);
//printf("Atom Information\n");
//for(int r=0; r<at_type.size(); r++){printf("%s:%d:%d\n",at_type[r].c_str(),res_types[r],res_ids[r]);}
container_poly con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,false,false,false,8);
//Residues are split by backbone and base atoms
for(int a = 0; a<at_type.size(); a++){
double x = xps[a];
double y = yps[a];
double z = zps[a];
double r=0.5;
if (at_type[a].find("H",0) != string::npos){
r=.31;
}else{
if (at_type[a].find("C",0) != string::npos){
r=.73;
}else{
if (at_type[a].find("O",0) != string::npos){
r=.66;
}else{
if (at_type[a].find("P",0) != string::npos){
r=1.07;
}else{
if (at_type[a].find("N",0) != string::npos){
r=0.71;
}else{
printf("Can't identify atom type: %s\n",at_type[a].c_str());
}
}
}
}
}
//printf("Adding: %8.3f %8.3f %8.3f %8.3f\n",x,y,z,r);
con.put(a,x,y,z,r);
}

vector<fpoint> volumes;
vector<fpoint> surface;
vector<fpoint> surfA;
vector<fpoint>surfG;
vector<fpoint>surfC;
vector<fpoint>surfU;
vector<int> orderA;
vector<int> orderG;
vector<int> orderC;
vector<int> orderU;
vector<int> unpacked;
vector<fpoint> bvolumes;
vector<fpoint> bsurface;
vector<fpoint> bsurfA;
vector<fpoint> bsurfG;
vector<fpoint> bsurfC;
vector<fpoint> bsurfU;
vector<int> borderA;
vector<int> borderG;
vector<int> borderC;
vector<int> borderU;
vector<int> bunpacked;
vector<int> parent;
vector<int> child;
vector<int> order;

vector<int> atom0,atom1,atom2,atom3,atom4,atom5,atom6,atom7,atom8,atom9,atom10,atom11,atom12,atom13,atom14,atom15,atom16,atom17,atom18,atom19,atom20,atom21,atom22,atom23,atom24,atom25,atom26;
vector<fpoint> vatom0,vatom1,vatom2,vatom3,vatom4,vatom5,vatom6,vatom7,vatom8,vatom9,vatom10,vatom11,vatom12,vatom13,vatom14,vatom15,vatom16,vatom17,vatom18,vatom19,vatom20,vatom21,vatom22,vatom23,vatom24,vatom25,vatom26;

con.get_res_volume(
	&res_types,
	&res_ids,
	&at_type,
	&volumes,
	&unpacked,
	&surface,
	&surfA,
	&surfG,
	&surfC,
	&surfU,
	&orderA,
	&orderG,
	&orderC,
	&orderU,
	
	&bvolumes,
	&bunpacked,
	&bsurface,
	&bsurfA,
	&bsurfG,
	&bsurfC,
	&bsurfU,
	&borderA,
	&borderG,
	&borderC,
	&borderU,
	
	&vatom0,
	&vatom1,
	&vatom2,
	&vatom3,
	&vatom4,
	&vatom5,
	&vatom6,
	&vatom7,
	&vatom8,
	&vatom9,
	&vatom10,
	&vatom11,
	&vatom12,
	&vatom13,
	&vatom14,
	&vatom15,
	&vatom16,
	&vatom17,
	&vatom18,
	&vatom19,
	&vatom20,
	&vatom21,
	&vatom22,
	&vatom23,
	&vatom24,
	&vatom25,
	&vatom26,

	&atom0,
	&atom1,
	&atom2,
	&atom3,
	&atom4,
	&atom5,
	&atom6,
	&atom7,
	&atom8,
	&atom9,
	&atom10,
	&atom11,
	&atom12,
	&atom13,
	&atom14,
	&atom15,
	&atom16,
	&atom17,
	&atom18,
	&atom19,
	&atom20,
	&atom21,
	&atom22,
	&atom23,
	&atom24,
	&atom25,
	&atom26,
	
	&parent,
	&child,
	&order
);
//printf("New Score Function: %d\n",volumes.size());
int vind=0;

vector<int> gorder;
vector<int> gtype1;
vector<int> gtype2;
vector<int> gtype3;
vector<int> gtype4;

for(int r = 0; r < num_mols; r++){
if(mols[r]->group==group_id){ 
	int tmporder=0;
	int g1=0;
	int g2=0;
	int g3=0;
	int g4=0;
	for(int g = 0; g< parent.size(); g++){
	if(parent[g]==r){
		tmporder++;
	if(rtype[child[g]]==1){g1++;}	
	if(rtype[child[g]]==2){g2++;}	
	if(rtype[child[g]]==3){g3++;}	
	if(rtype[child[g]]==4){g4++;}	
		}	
	}
gorder.push_back(tmporder);
gtype1.push_back(g1);
gtype2.push_back(g2);
gtype3.push_back(g3);
gtype4.push_back(g4);
}
}
//printf("Graph Information\n");



//for(int g=0; g< num_mols; g++){printf("%c %d %d %d %d %d\n",mols[g]->seq,gorder[g],gtype1[g],gtype2[g],gtype3[g],gtype4[g]);}
//printf("End Score Function\n");
//Write everything to the score matrix
vind=0;
for(int vgind=0; vgind<(num_mols); vgind++){
if(mols[vgind]->group==group_id){
score->s0[vind]=gorder[vind];
score->s1[vind]=gtype1[vind];
score->s2[vind]=gtype2[vind];
score->s3[vind]=gtype3[vind];
score->s4[vind]=gtype4[vind];
score->s5[vind]=rtype[vind];
score->s6[vind]=volumes[vgind];
score->s7[vind]=unpacked[vgind];
score->s8[vind]=surface[vgind];
score->s9[vind]=bvolumes[vgind];
score->s10[vind]=bunpacked[vgind];
score->s11[vind]=bsurface[vgind];
score->s12[vind]=surfG[vgind];
score->s13[vind]=surfA[vgind];
score->s14[vind]=surfU[vgind];
score->s15[vind]=surfC[vgind];
score->s16[vind]=bsurfG[vgind];
score->s17[vind]=bsurfA[vgind];
score->s18[vind]=bsurfU[vgind];
score->s19[vind]=bsurfC[vgind];
score->s20[vind]=orderG[vgind];
score->s21[vind]=orderA[vgind];
score->s22[vind]=orderU[vgind];
score->s23[vind]=orderC[vgind];
score->s24[vind]=borderG[vgind];
score->s25[vind]=borderA[vgind];
score->s26[vind]=borderU[vgind];
score->s27[vind]=borderC[vgind];
score->s28[vind]=vatom0[vgind];
score->s29[vind]=vatom1[vgind];
score->s30[vind]=vatom2[vgind];
score->s31[vind]=vatom3[vgind];
score->s32[vind]=vatom4[vgind];
score->s33[vind]=vatom5[vgind];
score->s34[vind]=vatom6[vgind];
score->s35[vind]=vatom7[vgind];
score->s36[vind]=vatom8[vgind];
score->s37[vind]=vatom9[vgind];
score->s38[vind]=vatom10[vgind];
score->s39[vind]=vatom11[vgind];
score->s40[vind]=vatom12[vgind];
score->s41[vind]=vatom13[vgind];
score->s42[vind]=vatom14[vgind];
score->s43[vind]=vatom15[vgind];
score->s44[vind]=vatom16[vgind];
score->s45[vind]=vatom17[vgind];
score->s46[vind]=vatom18[vgind];
score->s47[vind]=vatom19[vgind];
score->s48[vind]=vatom20[vgind];
score->s49[vind]=vatom21[vgind];
score->s50[vind]=vatom22[vgind];
score->s51[vind]=vatom23[vgind];
score->s52[vind]=vatom24[vgind];
score->s53[vind]=vatom25[vgind];
score->s54[vind]=vatom26[vgind];
score->s55[vind]=atom0[vgind];
score->s56[vind]=atom1[vgind];
score->s57[vind]=atom2[vgind];
score->s58[vind]=atom3[vgind];
score->s59[vind]=atom4[vgind];
score->s60[vind]=atom5[vgind];
score->s61[vind]=atom6[vgind];
score->s62[vind]=atom7[vgind];
score->s63[vind]=atom8[vgind];
score->s64[vind]=atom9[vgind];
score->s65[vind]=atom10[vgind];
score->s66[vind]=atom11[vgind];
score->s67[vind]=atom12[vgind];
score->s68[vind]=atom13[vgind];
score->s69[vind]=atom14[vgind];
score->s70[vind]=atom15[vgind];
score->s71[vind]=atom16[vgind];
score->s72[vind]=atom17[vgind];
score->s73[vind]=atom18[vgind];
score->s74[vind]=atom19[vgind];
score->s75[vind]=atom20[vgind];
score->s76[vind]=atom21[vgind];
score->s77[vind]=atom22[vgind];
score->s78[vind]=atom23[vgind];
score->s79[vind]=atom24[vgind];
score->s80[vind]=atom25[vgind];
score->s81[vind]=atom26[vgind];
score->s82[vind]=hydros[vind];
score->s83[vind]=doublets[vind];
vind++;
}      
}  
}         

void Polymere::resscore(Score * score,vector<double> * dparams){
vector<float> hydros;
vector<int> hpairs;
vector<float> doublets;
//Residue Level Scores
doubletV(dparams,3.0,3.0,3.0,0.2,0.2,&doublets);
hydroPairs(&hydros,&hpairs);
//Phosphate to O3` distances


Residue * Aresidue = new Residue;
// Set up the number of blocks that the container is divided into
const int n_x=20,n_y=20,n_z=20;
// Set the number of particles that are going to be randomly introduced
int particle_nums=0;
int i=0;
        double x,y,z,rad;
        // Create a container with the geometry given above, and make it
        // non-periodic in each of the three coordinates. Allocate space for
        // eight particles within each computational block

        // Randomly add particles into the container
double min_x=999;double max_x=-999;
double min_y=999;double max_y=-999;
double min_z=999;double max_z=-999;
vector<double> xps;
vector<double> yps;
vector<double> zps;

std::vector<int> res_ids;
std::vector<int> res_types;
std::vector<string> at_type;
std::vector<int> rtype;
for(int resn=0; resn<(num_mols); resn++){
int bpos = resn*81;
if(mols[resn]->seq=='G'){
rtype.push_back(1);
for(int rn=0; rn < 23; rn++){
res_ids.push_back(resn);
if(rn < 12){res_types.push_back(-1);}else{res_types.push_back(1);}
string aname(Aresidue->atomnames1[rn]);
at_type.push_back(aname);
//create an atom x y z atom_name chain  res_name
int p1 = bpos+rn*3;
int p2 = bpos+rn*3+1;
int p3 = bpos+rn*3+2;
particle_nums++;
xps.push_back(full_atoms[p1]);
yps.push_back(full_atoms[p2]);
zps.push_back(full_atoms[p3]);

if(min_x > full_atoms[p1]){min_x=full_atoms[p1];}
if(min_y > full_atoms[p2]){min_y=full_atoms[p2];}
if(min_z > full_atoms[p3]){min_z=full_atoms[p3];}

if(max_x < full_atoms[p1]){max_x=full_atoms[p1];}
if(max_y < full_atoms[p2]){max_y=full_atoms[p2];}
if(max_z < full_atoms[p3]){max_z=full_atoms[p3];}

}
}else{
if(mols[resn]->seq=='A'){
rtype.push_back(2);
for(int rn=0; rn < 22; rn++){
res_ids.push_back(resn);
if(rn < 12){res_types.push_back(-2);}else{res_types.push_back(2);}
string aname(Aresidue->atomnames2[rn]);
at_type.push_back(aname);
//create an atom x y z atom_name chain  res_name
int p1 = bpos+rn*3;
int p2 = bpos+rn*3+1;
int p3 = bpos+rn*3+2;

particle_nums++;
xps.push_back(full_atoms[p1]);
yps.push_back(full_atoms[p2]);
zps.push_back(full_atoms[p3]);

if(min_x > full_atoms[p1]){min_x=full_atoms[p1];}
if(min_y > full_atoms[p2]){min_y=full_atoms[p2];}
if(min_z > full_atoms[p3]){min_z=full_atoms[p3];}

if(max_x < full_atoms[p1]){max_x=full_atoms[p1];}
if(max_y < full_atoms[p2]){max_y=full_atoms[p2];}
if(max_z < full_atoms[p3]){max_z=full_atoms[p3];}
}


}else{
if(mols[resn]->seq=='U'){
rtype.push_back(3);
for(int rn=0; rn < 20; rn++){
res_ids.push_back(resn);
if(rn < 12){res_types.push_back(-3);}else{res_types.push_back(3);}
string aname(Aresidue->atomnames3[rn]);
at_type.push_back(aname);
//create an atom x y z atom_name chain  res_name
int p1 = bpos+rn*3;
int p2 = bpos+rn*3+1;
int p3 = bpos+rn*3+2;

particle_nums++;
xps.push_back(full_atoms[p1]);
yps.push_back(full_atoms[p2]);
zps.push_back(full_atoms[p3]);

if(min_x > full_atoms[p1]){min_x=full_atoms[p1];}
if(min_y > full_atoms[p2]){min_y=full_atoms[p2];}
if(min_z > full_atoms[p3]){min_z=full_atoms[p3];}

if(max_x < full_atoms[p1]){max_x=full_atoms[p1];}
if(max_y < full_atoms[p2]){max_y=full_atoms[p2];}
if(max_z < full_atoms[p3]){max_z=full_atoms[p3];}
}

}else{
if(mols[resn]->seq=='C'){
rtype.push_back(4);
for(int rn=0; rn < 20; rn++){
res_ids.push_back(resn);
if(rn < 12){res_types.push_back(-4);}else{res_types.push_back(4);}
string aname(Aresidue->atomnames4[rn]);
at_type.push_back(aname);
//create an atom x y z atom_name chain  res_name
int p1 = bpos+rn*3;
int p2 = bpos+rn*3+1;
int p3 = bpos+rn*3+2;

particle_nums++;
xps.push_back(full_atoms[p1]);
yps.push_back(full_atoms[p2]);
zps.push_back(full_atoms[p3]);


if(min_x > full_atoms[p1]){min_x=full_atoms[p1];}
if(min_y > full_atoms[p2]){min_y=full_atoms[p2];}
if(min_z > full_atoms[p3]){min_z=full_atoms[p3];}

if(max_x < full_atoms[p1]){max_x=full_atoms[p1];}
if(max_y < full_atoms[p2]){max_y=full_atoms[p2];}
if(max_z < full_atoms[p3]){max_z=full_atoms[p3];}
}
}else{
}
}
}
}

}


min_x = min_x -9;
min_y = min_y -9;
min_z = min_z -9;

max_x = max_x +9;
max_y = max_y +9;
max_z = max_z +9;

//Voro container needs to be a cube
double minv = min_x;
if(minv > min_y){minv=min_y;}
if(minv > min_z){minv=min_z;}

double maxv = max_x;
if(maxv < max_y){maxv=max_y;}
if(maxv < max_z){maxv=max_z;}


const double x_min=minv;
const double x_max=maxv;
const double y_min=minv;
const double y_max=maxv;
const double z_min=minv;
const double z_max=maxv;

//printf("Ranges %8.3f - %8.3f : %8.3f - %8.3f : %8.3f - %8.3f\n",x_min,x_max,y_min,y_max,z_min,z_max);
//printf("Atom Information\n");
//for(int r=0; r<at_type.size(); r++){printf("%s:%d:%d\n",at_type[r].c_str(),res_types[r],res_ids[r]);}
container_poly con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,false,false,false,8);
//Residues are split by backbone and base atoms
for(int a = 0; a<at_type.size(); a++){
double x = xps[a];
double y = yps[a];
double z = zps[a];
double r=0.5;
if (at_type[a].find("H",0) != string::npos){
r=.31;
}else{
if (at_type[a].find("C",0) != string::npos){
r=.73;
}else{
if (at_type[a].find("O",0) != string::npos){
r=.66;
}else{
if (at_type[a].find("P",0) != string::npos){
r=1.07;
}else{
if (at_type[a].find("N",0) != string::npos){
r=0.71;
}else{
printf("Can't identify atom type: %s\n",at_type[a].c_str());
}
}
}
}
}
//printf("Adding: %8.3f %8.3f %8.3f %8.3f\n",x,y,z,r);
con.put(a,x,y,z,r);
}

vector<fpoint> volumes;
vector<fpoint> surface;
vector<fpoint> surfA;
vector<fpoint>surfG;
vector<fpoint>surfC;
vector<fpoint>surfU;
vector<int> orderA;
vector<int> orderG;
vector<int> orderC;
vector<int> orderU;
vector<int> unpacked;
vector<fpoint> bvolumes;
vector<fpoint> bsurface;
vector<fpoint> bsurfA;
vector<fpoint> bsurfG;
vector<fpoint> bsurfC;
vector<fpoint> bsurfU;
vector<int> borderA;
vector<int> borderG;
vector<int> borderC;
vector<int> borderU;
vector<int> bunpacked;
vector<int> parent;
vector<int> child;
vector<int> order;

vector<int> atom0,atom1,atom2,atom3,atom4,atom5,atom6,atom7,atom8,atom9,atom10,atom11,atom12,atom13,atom14,atom15,atom16,atom17,atom18,atom19,atom20,atom21,atom22,atom23,atom24,atom25,atom26;
vector<fpoint> vatom0,vatom1,vatom2,vatom3,vatom4,vatom5,vatom6,vatom7,vatom8,vatom9,vatom10,vatom11,vatom12,vatom13,vatom14,vatom15,vatom16,vatom17,vatom18,vatom19,vatom20,vatom21,vatom22,vatom23,vatom24,vatom25,vatom26;

con.get_res_volume(
	&res_types,
	&res_ids,
	&at_type,
	&volumes,
	&unpacked,
	&surface,
	&surfA,
	&surfG,
	&surfC,
	&surfU,
	&orderA,
	&orderG,
	&orderC,
	&orderU,
	
	&bvolumes,
	&bunpacked,
	&bsurface,
	&bsurfA,
	&bsurfG,
	&bsurfC,
	&bsurfU,
	&borderA,
	&borderG,
	&borderC,
	&borderU,
	
	&vatom0,
	&vatom1,
	&vatom2,
	&vatom3,
	&vatom4,
	&vatom5,
	&vatom6,
	&vatom7,
	&vatom8,
	&vatom9,
	&vatom10,
	&vatom11,
	&vatom12,
	&vatom13,
	&vatom14,
	&vatom15,
	&vatom16,
	&vatom17,
	&vatom18,
	&vatom19,
	&vatom20,
	&vatom21,
	&vatom22,
	&vatom23,
	&vatom24,
	&vatom25,
	&vatom26,

	&atom0,
	&atom1,
	&atom2,
	&atom3,
	&atom4,
	&atom5,
	&atom6,
	&atom7,
	&atom8,
	&atom9,
	&atom10,
	&atom11,
	&atom12,
	&atom13,
	&atom14,
	&atom15,
	&atom16,
	&atom17,
	&atom18,
	&atom19,
	&atom20,
	&atom21,
	&atom22,
	&atom23,
	&atom24,
	&atom25,
	&atom26,
	
	&parent,
	&child,
	&order
);
//printf("New Score Function: %d\n",volumes.size());
int vind=0;
/*
for(int v=0; v<volumes.size(); v++){
printf("%d   %c  %8.3f %d %8.3f %8.3f %d %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %d %d %d %d %d %d %d %d %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.f %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %8.3f %8.3f\n",
vind,
mols[vind]->seq,
volumes[vind],
unpacked[vind],
surface[vind],
bvolumes[vind],
bunpacked[vind],
bsurface[vind],
surfG[vind],
surfA[vind], //10
surfU[vind],
surfC[vind], 
bsurfG[vind],
bsurfA[vind],
bsurfU[vind],
bsurfC[vind],
orderG[vind],
orderA[vind],
orderU[vind],
orderC[vind], //20
borderG[vind],
borderA[vind],
borderU[vind],
borderC[vind],
vatom0[vind],
vatom1[vind],
vatom2[vind],
vatom3[vind],
vatom4[vind],
vatom5[vind],//30
vatom6[vind],
vatom7[vind],
vatom8[vind],
vatom9[vind],
vatom10[vind],
vatom11[vind],
vatom12[vind],
vatom13[vind],
vatom14[vind],
vatom15[vind],//40
vatom16[vind],
vatom17[vind],
vatom18[vind],
vatom19[vind],
vatom20[vind],
vatom21[vind],
vatom22[vind],
vatom23[vind],
vatom24[vind],
vatom25[vind],//50
vatom26[vind],
atom0[vind],
atom1[vind],
atom2[vind],
atom3[vind],
atom4[vind],
atom5[vind],
atom6[vind],
atom7[vind],
atom8[vind],
atom9[vind],
atom10[vind],//60
atom11[vind],
atom12[vind],
atom13[vind],
atom14[vind],
atom15[vind],
atom16[vind],
atom17[vind],
atom18[vind],
atom19[vind],
atom20[vind],//70
atom21[vind],
atom22[vind],
atom23[vind],
atom24[vind],
atom25[vind],
atom26[vind],
hydros[vind],
doublets[vind]
);

vind++;
}
*/

vector<int> gorder;
vector<int> gtype1;
vector<int> gtype2;
vector<int> gtype3;
vector<int> gtype4;

for(int r = 0; r < num_mols; r++){ 
	int tmporder=0;
	int g1=0;
	int g2=0;
	int g3=0;
	int g4=0;
	for(int g = 0; g< parent.size(); g++){
	if(parent[g]==r){
		tmporder++;
	if(rtype[child[g]]==1){g1++;}	
	if(rtype[child[g]]==2){g2++;}	
	if(rtype[child[g]]==3){g3++;}	
	if(rtype[child[g]]==4){g4++;}	
		}	
	}
gorder.push_back(tmporder);
gtype1.push_back(g1);
gtype2.push_back(g2);
gtype3.push_back(g3);
gtype4.push_back(g4);
}

//printf("Graph Information\n");



//for(int g=0; g< num_mols; g++){printf("%c %d %d %d %d %d\n",mols[g]->seq,gorder[g],gtype1[g],gtype2[g],gtype3[g],gtype4[g]);}
//printf("End Score Function\n");
//Write everything to the score matrix
for(int vind=0; vind<num_mols; vind++){
score->s0[vind]=gorder[vind];
score->s1[vind]=gtype1[vind];
score->s2[vind]=gtype2[vind];
score->s3[vind]=gtype3[vind];
score->s4[vind]=gtype4[vind];
score->s5[vind]=rtype[vind];
score->s6[vind]=volumes[vind];
score->s7[vind]=unpacked[vind];
score->s8[vind]=surface[vind];
score->s9[vind]=bvolumes[vind];
score->s10[vind]=bunpacked[vind];
score->s11[vind]=bsurface[vind];
score->s12[vind]=surfG[vind];
score->s13[vind]=surfA[vind];
score->s14[vind]=surfU[vind];
score->s15[vind]=surfC[vind];
score->s16[vind]=bsurfG[vind];
score->s17[vind]=bsurfA[vind];
score->s18[vind]=bsurfU[vind];
score->s19[vind]=bsurfC[vind];
score->s20[vind]=orderG[vind];
score->s21[vind]=orderA[vind];
score->s22[vind]=orderU[vind];
score->s23[vind]=orderC[vind];
score->s24[vind]=borderG[vind];
score->s25[vind]=borderA[vind];
score->s26[vind]=borderU[vind];
score->s27[vind]=borderC[vind];
score->s28[vind]=vatom0[vind];
score->s29[vind]=vatom1[vind];
score->s30[vind]=vatom2[vind];
score->s31[vind]=vatom3[vind];
score->s32[vind]=vatom4[vind];
score->s33[vind]=vatom5[vind];
score->s34[vind]=vatom6[vind];
score->s35[vind]=vatom7[vind];
score->s36[vind]=vatom8[vind];
score->s37[vind]=vatom9[vind];
score->s38[vind]=vatom10[vind];
score->s39[vind]=vatom11[vind];
score->s40[vind]=vatom12[vind];
score->s41[vind]=vatom13[vind];
score->s42[vind]=vatom14[vind];
score->s43[vind]=vatom15[vind];
score->s44[vind]=vatom16[vind];
score->s45[vind]=vatom17[vind];
score->s46[vind]=vatom18[vind];
score->s47[vind]=vatom19[vind];
score->s48[vind]=vatom20[vind];
score->s49[vind]=vatom21[vind];
score->s50[vind]=vatom22[vind];
score->s51[vind]=vatom23[vind];
score->s52[vind]=vatom24[vind];
score->s53[vind]=vatom25[vind];
score->s54[vind]=vatom26[vind];
score->s55[vind]=atom0[vind];
score->s56[vind]=atom1[vind];
score->s57[vind]=atom2[vind];
score->s58[vind]=atom3[vind];
score->s59[vind]=atom4[vind];
score->s60[vind]=atom5[vind];
score->s61[vind]=atom6[vind];
score->s62[vind]=atom7[vind];
score->s63[vind]=atom8[vind];
score->s64[vind]=atom9[vind];
score->s65[vind]=atom10[vind];
score->s66[vind]=atom11[vind];
score->s67[vind]=atom12[vind];
score->s68[vind]=atom13[vind];
score->s69[vind]=atom14[vind];
score->s70[vind]=atom15[vind];
score->s71[vind]=atom16[vind];
score->s72[vind]=atom17[vind];
score->s73[vind]=atom18[vind];
score->s74[vind]=atom19[vind];
score->s75[vind]=atom20[vind];
score->s76[vind]=atom21[vind];
score->s77[vind]=atom22[vind];
score->s78[vind]=atom23[vind];
score->s79[vind]=atom24[vind];
score->s80[vind]=atom25[vind];
score->s81[vind]=atom26[vind];
score->s82[vind]=hydros[vind];
score->s83[vind]=doublets[vind];
score->s84[vind]=hpairs[vind];
}        
}         
         
void Polymere::scoreit(int nmax, int movegroup, const std::vector< Fragment *>  frags,std::vector< Fragatoms *> atoms,std::vector<Fragall *> batoms,Fragclass * dat,Polymere * pol,std::vector<double> * dparams){
//Uses a biased Metropolis Monte-Carlo Method to solve the structure
//The proportion of the score contributed by a particular fragment is measured
//Movements are focused on regions performing poorly
//Hydrogen Bonding 
//Doublet Library 
using namespace std;
Polymere  * pol2 = new Polymere();
Polymere  * outpol = new Polymere();
//interations = Maximum number of iterations to perform
//max_rand    = Maximum amount of randomization 
//frags       = Fragment database
//dat         = Indexed fragments
//group       = Movement group to randomize
//pol         = This is the master polymere that everything is compared to
std::vector<Matches> mat; 
mat.reserve(7000);
double id=0;
double pid=0;
//uniqpol is a vector of structures that results from a single closed moves search for a given structure
std::vector<Polymere *> lastround;
std::vector<Polymere *> last5;
std::vector<Polymere *> uniqpol;
uniqpol.reserve(1000);

std::vector< std::vector<int> > groups;
int maxg=0;
//Read in the movement group from the ct file 
for(int g=0; g<pol->num_mols; g++){
int gc = pol->mols[g]->group;
//cout << "group: " << gc << " bp1: " << pol->mols[g]->id << endl;
if(gc>maxg){
maxg=gc;
}
}
maxg++;
//cout << "Max Groups: " << maxg << endl;
for(int i=1; i < maxg; i++){
std::vector<int> group;
for(int g=0; g<pol->num_mols; g++){
int gc = pol->mols[g]->group;
if(gc==i){
group.push_back(g);
}
}
groups.push_back(group);
}

int sig = groups[movegroup].size();

//Set counts to 0 for each position
for(int g=0; g<groups[movegroup].size(); g++){
mols[groups[movegroup][g]]->group=0;
mols[groups[movegroup][g]]->prob=0;
}

//cout << "Starting randomization" << endl;
int tot;

//Vectors containing the indices of valid movement groups
std::vector<int> groupOrder1;
std::vector<int> groupOrder2;
std::vector<int> groupOrder0;

for(int i=1; i < sig-2; i++){
groupOrder1.push_back(i);
groupOrder2.push_back(i);
}
int nt=(sig-3)*(sig-3);
for(int i=0; i <nt; i++){groupOrder0.push_back(i);} 
//Set initial Score and RMSD threshold
float score_threshold=9999;
int rand_threshold=0;
//Foreach starting structure score them based on some criteria for each move.  Need to loop through 
//a certain number of cycles or until the simulation converges
//Add loop to perform multiple generations of movements
//printf("Nmax: %d\n",nmax);
for(int iteration=0; iteration < nmax; iteration++){
//cout << "Iteration: " << iteration ;
////////////////////////////////////////////////////////////
// Move current polymere, score, add improvements to uniq //
////////////////////////////////////////////////////////////
//Generate Random order for movements
int n = groupOrder1.size();
int k;
int temp;
while(n>1){
k=(rand() % n);
--n;
temp=groupOrder1[n];
groupOrder1[n]=groupOrder1[k];
groupOrder1[k]=temp;
}

n=groupOrder2.size();
while(n>1){
k=(rand() % n);
--n;
temp=groupOrder2[n];
groupOrder2[n]=groupOrder2[k];
groupOrder2[k]=temp;
}
n=groupOrder1.size()*groupOrder2.size();
while(n>1){
k=(rand() % n);
--n;
temp=groupOrder0[n];
groupOrder0[n]=groupOrder0[k];
groupOrder0[k]=temp;
}
//cout << "Group1: ";
for(int i=0; i< groupOrder1.size(); i++){
//cout << groupOrder1[i] << " ";
}
//cout << endl;
int maxout=(sig-3)*(sig-3);
int breakout = 0;
int newtot=0;

//Reset the structure list
//cout << "Last Size: " << lastround.size() << endl;
for(int u=0; u<lastround.size(); u++){
for(int m=0; m<lastround[u]->mols.size(); m++){
delete lastround[u]->mols[m];
}
lastround[u]->mols.clear();
}
lastround.clear();

for(int t=0; t < maxout; t++){
int i=groupOrder0[t] % (sig-3);
int j= (groupOrder0[t]-i)/(sig-3);
//Clear Uniqpols From the previous Run
for(int u=0; u<uniqpol.size(); u++){
for(int m=0; m<uniqpol[u]->mols.size(); m++){
delete uniqpol[u]->mols[m];
}
uniqpol[u]->mols.clear();
}
uniqpol.clear();
//Check to make sure the groups are in focus
//Decide if the group should be moved
float prob=rand()/(float)(RAND_MAX);
float totalprob = mols[groups[movegroup][groupOrder1[i]]]->prob + mols[ groups[movegroup][groupOrder2[j]] ]->prob;
totalprob = totalprob/10;
//totalprob=totalprob*(float)groups[movegroup].size()/2;
if((groupOrder1[i]+3) < groupOrder2[j] && totalprob < prob){
//Loop through possible position combinations
//Only use positions that aren't set to -1 
//If position is less than -1 that is the cluster number to use
//As soon as move is identified make the move
//Find center of each group movement and save to branches
//cout << "Group 1: " << groups[movegroup][groupOrder1[i]] << " Group 2: " << groups[movegroup][groupOrder2[j]] << endl;
mat=chooseMoves(groups[movegroup][groupOrder1[i]],groups[movegroup][groupOrder2[j]],dat,-1);
tot=mat.size();
//cout << "tot: " << tot << " ";
//Identify the center of cluster of moves  
double id=0;
double pid=0;
for(int w=0; w < tot; w++){
//cout << "Creating new structure " << endl;
//Make a copy of the molecule array
//Use Group to track randomization
int nmol=pol->num_mols;
pol2->copy(this);
if(batoms[mat[w].frag1]->check()==1 && batoms[mat[w].frag2]->check()==1){
pol2->changeFrag2(mat[w].pos1,mat[w].frag1,frags);
pol2->changeFrag2(mat[w].pos2,mat[w].frag2,frags);
//Check base-pair differences
int c1dist1=pol2->c1dist();
if(c1dist1!=0){id=pol2->checkhbonds();}else{id=0;
//cout << "Failed C1dist: " << c1dist1 << endl;
}
//Make a copy of pol2 and save
if(id!=0){
//Generate a score 
//Update full atom coordinates
pol2->coordSysNt(mat[w].pos1,&atoms[mat[w].frag1]->atoms);
pol2->coordSysNt(mat[w].pos2,&atoms[mat[w].frag2]->atoms);
pol2->updateFull();
pol2->coordSysAT(mat[w].pos1,&batoms[mat[w].frag1]->atoms);
pol2->coordSysAT((mat[w].pos2),&batoms[mat[w].frag2]->atoms);
//Need to make it so that the probabilities are adjusted 
std::vector<int> group_inds=groups[movegroup];
//float s=pol2->hydroG(&group_inds);

//GLOBAL SCORES 
Score * pscore = new Score(pol2->num_mols);
pol2->resscore(pscore,dparams);
//float stacked = pol2->gaussDbl(dparams,3,3,3,0.2,0.2,&group_inds);
//float hydro = pol2->hydroG(&group_inds);
//pol2->score = stacked-hydro;
float g1=pscore->gscore(1);//Graph Order
float g2=pscore->gscore(2);//Base Packing
float g3=pscore->gscore(3);//Backbone Packing
float g4=pscore->gscore(4);//Hydrogen Bond
float g5=pscore->gscore(5);//Surface Area
float g6=pscore->gscore(6);//Low resolution score
float ogscore1=(g1+g2+g3+g5)/pol->num_mols;
float ogscore2=g6/pol->num_mols;
float ogscore3=g4/pol->num_mols;

float overall_score = (10000+ogscore1) + (-1*(1/(1+exp(-1*(ogscore1-20)))-1))*(ogscore2-1000) +  ((-1*(1/(1+exp(-1*(ogscore2-15)))-1))*(-1*(1/(1+exp(-1*(ogscore1-20)))-1)))*(-10*ogscore3-100);
printf("Overall Score: %8.3f\n",overall_score);

if(overall_score < -10000){
pscore->write_globals("bad_globals.txt");
overall_score=9999; 
}
delete pscore;
pol2->score=overall_score;
//GLOBAL SCORES 


//accept the move using metropolis criteria
int accept=0;
newtot++;
//To prevent us from constantly using the same mols remove lastround that are same as current
float curr_rms=rmsAll(pol2);
if(curr_rms>0.05){
lastround.push_back(new Polymere());
lastround.back()->copy(pol2);
}else{
accept=-1;
}
if(accept==0){
if(pol2->score < score){
accept=1;
}else{
float prob1=exp(score-(pol2->score));
printf("\nCurrent Score: %8.3f New Score: %8.3f Prob: %8.3f\n",score,pol2->score,prob1);
float rand1 = rand()/((float)RAND_MAX);
if(rand1 < prob){
accept=1;
}else{
}
}
}
if(accept==1){
float minrmsdist=999;
for(int i=0; i<last5.size(); i++){
float dist = pol2->rmsAll(last5[i]);
if(dist < minrmsdist){
minrmsdist=dist;
}
}
if(minrmsdist < 0.1){
accept=0;
//Remove last entry from lastround
for(int m=0; m<lastround.back()->mols.size(); m++){
delete lastround.back()->mols[m];
}
lastround.pop_back();
}
}
//cout << "Matches: " << mat[w].pos1 << ":" << mat[w].pos2 << endl;
if(accept==1){
uniqpol.push_back(new Polymere());
uniqpol.back()->copy(pol2);
}
}
}else{
//printf("Bad Frags: %d %d\n",mat[w].frag1,mat[w].frag2);
}
}
//cout << "Uniqpol Size: " << uniqpol.size();
if(uniqpol.size()>0){
////////////////////////////////////////////////////////////
//Identify the center of structures in uniqpol
////////////////////////////////////////////////////////////
//Find structure with minimal average rmsd from all others
int center_ind;
float minrmsd=99999;
if(uniqpol.size() > 2){
for(int i=0; i<uniqpol.size(); i++){
float avermsd=0;
for(int j=0; j<uniqpol.size(); j++){
float rmsd_diff = uniqpol[i]->rmsAll(uniqpol[j]);
avermsd=avermsd+rmsd_diff;
}
avermsd=avermsd/uniqpol.size();

if(avermsd<minrmsd){
center_ind=i;
minrmsd=avermsd;
}

}
}else{
//If only two of them exist take the first one
center_ind=0;
}
minrmsd=rmsAll(pol);

//Set mols to uniqpol[center_ind]
if(uniqpol[center_ind]->score < score_threshold){
//Save structure to outPol
score_threshold=uniqpol[center_ind]->score;
outpol->copy(uniqpol[center_ind]);
}
last5.push_back(new Polymere);
copy(uniqpol[center_ind]);
last5.back()->copy(uniqpol[center_ind]);
//Break out of loop
cout << " Moved: " << center_ind << " RMSD: " << minrmsd << " Randomization: " << uniqpol[center_ind]->randtot << " Score: " << uniqpol[center_ind]->score;

Score * pscore = new Score(uniqpol[center_ind]->num_mols);
uniqpol[center_ind]->resscore(pscore,dparams);
//float stacked = pol2->gaussDbl(dparams,3,3,3,0.2,0.2,&group_inds);
//float hydro = pol2->hydroG(&group_inds);
//pol2->score = stacked-hydro;
float g1=pscore->gscore(1);//Graph Order
float g2=pscore->gscore(2);//Base Packing
float g3=pscore->gscore(3);//Backbone Packing
float g4=pscore->gscore(4);//Hydrogen Bond
float g5=pscore->gscore(5);//Hydrogen Bond
float overall_score=(g1+g2-g3+g5)/uniqpol[center_ind]->num_mols;
cout << " Current: "<< g1 << " " << g2 << " " << g3 << " " << g4 << " " << g5 << " " << overall_score << " P50: ";
//float stacked = pol2->gaussDbl(dparams,3,3,3,0.2,0.2,&group_inds);
//float hydro = pol2->hydroG(&group_inds);
//pol2->score = stacked-hydro;
pol->resscore(pscore,dparams);
g1=pscore->gscore(1);//Graph Order
g2=pscore->gscore(2);//Base Packing
g3=pscore->gscore(3);//Backbone Packing
g4=pscore->gscore(4);//Hydrogen Bond
g5=pscore->gscore(5);//Hydrogen Bond
overall_score=(g1+g2-g3+g5)/uniqpol[center_ind]->num_mols;
cout << g1 << " " << g2 << " " << g3 << " " << g4 << " " << g5 << " " <<  overall_score << endl;
delete pscore;
breakout=1;
break;
}else{
}


}

}
if(breakout==0){
cout << endl;
//Select member from lastround list with best score if lastround.size() > 0
if(lastround.size() > 1){
//Set current polymere to lastround[0]
std::sort(lastround.begin(), lastround.end(), PairComparePol ());
//cout << "Total Structures: " << lastround.size() << " : " << lastround[0]->score << " : " << lastround.back()->score <<  endl;
last5.push_back(new Polymere);
copy(lastround.back());
last5.back()->copy(lastround.back());
}else{
//Set last element in last5 to have randtot = 1
//Search for First element in list from the end that isn't randtot = 1
if(last5.size()>0){
last5.back()->randtot = 1;
int i=last5.size()-1;
int jump=0;
while(i>0){
if(last5[i]->randtot == 0){
jump=1;
cout << "Jumped Back to " << i << endl;
copy(last5[i]);
i=0;
}

i--;
}
if(jump==0){
break; 
cout << "Structure can not be moved with current constraints" << endl;
}
}else{
cout << "Structure can not be moved with given constraints" << endl;
break;
}
}
  


}
}

//Changes Values of current polymere to outpol
copy(outpol);
cout << "Minimal Structure Score: " << score << endl; 

for(int m=0; m<pol2->mols.size(); m++){
delete pol2->mols[m];
}
pol2->mols.clear();
//Remove outpol
for(int m=0; m<outpol->mols.size(); m++){
delete outpol->mols[m];
}
outpol->mols.clear();
for(int u=0; u<lastround.size(); u++){
for(int m=0; m<lastround[u]->mols.size(); m++){
delete lastround[u]->mols[m];
}
lastround[u]->mols.clear();
}
lastround.clear();
for(int u=0; u<last5.size(); u++){
for(int m=0; m<last5[u]->mols.size(); m++){
delete last5[u]->mols[m];
}
last5[u]->mols.clear();
}
last5.clear();
for(int u=0; u<uniqpol.size(); u++){
for(int m=0; m<uniqpol[u]->mols.size(); m++){
delete uniqpol[u]->mols[m];
}
uniqpol[u]->mols.clear();
}
uniqpol.clear();
delete outpol;
delete pol2;
//Clear uniqpol;
}


//Strands are represent by the id variable in the base class
//need to make sure backbone check and frag replacements only occur
//within a strand
/*
void Polymere::buildHelix(){



}

void Polymere::connectStrands(int strand_id1, int strand_id2,std::vector< std::pair<int,int> > * cons){




}

void Polymere::buildLoop(int loop_id,const std::vector< Fragment *>  frags,Fragclass * dat,std::vector< Fragatoms *> atoms,std::vector<Fragall *> batoms,std::vector< std::pair<int,int> > * cons,std::vector<int> * fixed){
//Builds an RNA tertiary structure using the secondary structure constraints
//This ignores all atoms not in the loop
//Requirements for a loop
//	Minimum stand size is 5 for each section of the loop
//	Extra sequence can be added after the base pair connecting strands to meet this requirement
//	Minimum loop size ignoring the extra sequence beyond the bases connecting strands is 12 

//Loop_id should be mapped to the group variable in the base
//Base->id is used to distinguish between strands in the loop
//fixed identifies loop_id's that should have fixed conformations
//           With fixed conformations you should calculate the distance from one chain connection to the other.  The total base length is x/5.5 for calculating the mid point of the loop
//           When moving the fixed conformations only moves 	
//Check that the constraints are sufficient to connect all elements in the loop
//Constraints use polymer molecular index to describe the constraints

//1) Each strand should have two connections


//Using loop_id and base_id create a new polymer that only contains the loop strands




//Create a single stand by connecting all but on of the bases in the constraints

//Identify the position that is equal distance from each of the bases in the unpaired constraint(this becomes the midpoint) base-pairs in the loop sequence count as two sequential bases

//Identify a single fragment(mid-2,mid-1,mid,mid+1 that brings the distance of mid+3 and mid-3 within the helical realm 9.5-11.5 C1' distance with no clashes
//If the midpoint is a half base mid + 3.5 and mid-3.5 

//Using paired movement groups symmetric around the midpoint, Identify closed moves that optimize packing (mid+2,mid-2),(mid+3,mid-3)... of those bases in the fragments. Depending the direction of the strands going away from the midpoint use the correct positions (pos1 = 5'-3' frag, pos2 = 3'->5' frag)
//If the midpoint is a half base mid + 3.5 and mid-3.5 


}
*/
/*
void Polymere::buildMultiChain(std::vector<Polymere * > polList,string seq,const std::vector< Fragment *>  frags,Fragclass * dat,std::vector< Fragatoms *> atoms,std::vector<Fragall *> batoms,std::vector<double> * dparams){
//Connects a list of polymers at specific postions by creating a loop 
//that brings the two strands together
//polList = vector of polymers that are connected in the order given end to end


//Initialize current polymer to the first polymere in the list
int cur_pos=0;
this->copyOver(polList[0],cur_pos,0,polList[0]->num_mols);


for(int pl=1; pl < pairPos.size(); pl++){
//Add 6 nt's plus pol2 to end of pol1
std::string seq;
seq.append(polList[pl-1]->mols[polList[pl-1]->num_mols]);
seq.append("GAUU");
seq.append(polList[pl]->mols[0]);
Polymere * add6 = new Polymere(seq,frags,dat,atoms,batoms,dparams);
this->copyOver(add6,cur_pos,0,6);
this->copyOver(polList[pl],cur_pos,0,polList[0]->num_mols);
//Link polymeres togther by adding a pseudo loop and moving frags until no clashes exist
//find loop that brings pol1 end into a base pair with start of pol1





} 




}

int Polymere::buildMultiBuldge(int s1, int s2,std::vector<int> * small_side,const std::vector< Fragment *>  frags,Fragclass * dat,std::vector< Fragatoms *> atoms,std::vector<Fragall *> batoms,int s_side,std::string outfile){
//Builds a helical section of the polymere given the start and end of two strands 5'-3' 
//and a set of constraints for each movement group
//genPairedStacks - 	uses stacked fragments only
//			identifies pairs of fragments that don't clash
//			and abide by particular contraints
// RULES:
// 1) Any assymetric buldge can only have a ratio less than 2:1
// 2) Small side is a vector describing the pairings
//Generate constraints map
std::vector< std::vector< std::pair<int,int> > > cons;
//Check that lenght small_side.size + small_size==1 = e2-s2
//Small Side
//	 1 = 1:2 ratio
//	 0 = 1:1 ratio
//	-1 = 2:1 ratio
//	 2 = 2:3 ratio (small side has a base pair in it)
//      -2 = 3:2 ratio 



int lcheck=0;
int success;
int m=0;
if(s_side==1){//Small side is position 1
for(int j=0; j<(small_side->size()-2); j++){
std::vector< std::pair<int,int> > bpcons;
int k=0;
int icr=1;	
int icr1=1;	
	if(small_side->at(j)==0){
		std::pair<int,int> scon((s1-j+0),(s2+m+0));	
		bpcons.push_back(scon);
		icr=1;
		icr1=1;
		}else{
		std::pair<int,int> scon1((s1-j+0),(s2+m+0));	
		bpcons.push_back(scon1);
		std::pair<int,int> scon2((s1-j+0),(s2+m+1));	
		bpcons.push_back(scon2);
		icr=2;
		icr1=2;
		}
	if(small_side->at(j+1)==0){
		std::pair<int,int> scon((s1-j-1),(s2+m+icr));	
		bpcons.push_back(scon);
		icr++;
		}else{
		std::pair<int,int> scon1((s1-j-1),(s2+m+icr));	
		bpcons.push_back(scon1);
		std::pair<int,int> scon2((s1-j-1),(s2+m+icr+1));	
		bpcons.push_back(scon2);
		icr=icr+2;
	}
	if(small_side->at(j+2)==0){
		std::pair<int,int> scon((s1-j-2),(s2+m+icr));	
		bpcons.push_back(scon);
		}else{
		std::pair<int,int> scon1((s1-j-2),(s2+m+icr));	
		bpcons.push_back(scon1);
		std::pair<int,int> scon2((s1-j-2),(s2+m+icr+1));	
		bpcons.push_back(scon2);
	}
m=m+icr1;
cons.push_back(bpcons);
}
}else{
if(s_side==2){
for(int j=0; j<(small_side->size()-2); j++){
std::vector< std::pair<int,int> > bpcons;
int k=0;
int icr=1;	
int icr1=1;	
	if(small_side->at(j)==0){
		std::pair<int,int> scon((s1-m+0),(s2+j+0));	
		bpcons.push_back(scon);
		icr=1;
		icr1=1;
		}else{
		std::pair<int,int> scon1((s1-m+0),(s2+j+0));	
		bpcons.push_back(scon1);
		std::pair<int,int> scon2((s1-m-1),(s2+j+0));	
		bpcons.push_back(scon2);
		icr=2;
		icr1=2;
		}
	if(small_side->at(j+1)==0){
		std::pair<int,int> scon((s1-m-icr),(s2+j+1));	
		bpcons.push_back(scon);
		icr++;
		}else{
		std::pair<int,int> scon1((s1-m-icr),(s2+j+1));	
		bpcons.push_back(scon1);
		std::pair<int,int> scon2((s1-m-icr-1),(s2+j+1));	
		bpcons.push_back(scon2);
		icr=icr+2;
	}
	if(small_side->at(j+2)==0){
		std::pair<int,int> scon((s1-m-icr),(s2+j+2));	
		bpcons.push_back(scon);
		}else{
		std::pair<int,int> scon1((s1-m-icr),(s2+j+2));	
		bpcons.push_back(scon1);
		std::pair<int,int> scon2((s1-m-icr-1),(s2+j+2));	
		bpcons.push_back(scon2);
	}
m=m+icr1;
cons.push_back(bpcons);
}


}else{
printf("Small side must be either 1 or 2 for buildBuldge\n");
return 0;
}

}
//Print out constraints for debugging purpose
for(int j=0; j<cons.size(); j++){
	printf("pos %d\n",j);
	for(int m=0; m<cons[j].size(); m++){
		printf("\t%d-%d\n",cons[j][m].first,cons[j][m].second);
	}
}

//Using cons build the helix.
//The first pair of positions in cons for each step is the starting 
//positions for the fragments
//5->3 starts at cons[i][0].second
//3->5 starts at cons[i][0].first - 2
std::vector<std::vector< std::pair<int,int> > > path_tree;
std::vector<int> path_pos;
int c=0;
while(path_tree.size() < cons.size() && c > -1){//Continue to extend the tree until all constraints are met
int pos1=cons[c][0].first-2;
int pos2=cons[c][0].second;
//Update the clash group to include those in this movement groups set of constraints
int clashGroup=5;
int unclashGroup=0;
std::vector< std::pair<int,int> > branch;
int goforward=this->genPairedCons(clashGroup,&cons[c],frags,dat,atoms,batoms,&branch);
if(goforward==1){
printf("Moving Forward: %d\n",branch.size());
path_pos.push_back((branch.size()-1)); //set current position to the last pair
path_tree.push_back(branch);
c++;//Move to the next set of constraints
}else{
printf("Moving backward:");
//Remove any clash constraints 
this->mols[pos1+0]->group=unclashGroup;
this->mols[pos1+1]->group=unclashGroup;
this->mols[pos1+2]->group=unclashGroup;
//Set group clashGroup to include 
this->mols[pos2+0]->group=unclashGroup;
this->mols[pos2+1]->group=unclashGroup;
this->mols[pos2+2]->group=unclashGroup;
//Remove clash constraints from the current group
//move back in path_tree
int cur_branch=c-1;
c=-1;
while(cur_branch>=0){
printf("Branch:%d  Position:%d\n",cur_branch,path_pos[cur_branch]);
if(path_pos[cur_branch]>=0){
path_pos[cur_branch]=path_pos[cur_branch]-1;
//move the positions of the polymere
int pos1=cons[cur_branch][0].first-2;
int pos2=cons[cur_branch][0].second;
this->mols[pos1+0]->group=clashGroup;
this->mols[pos1+1]->group=clashGroup;
this->mols[pos1+2]->group=clashGroup;
//Set group clashGroup to include 
this->mols[pos2+0]->group=clashGroup;
this->mols[pos2+1]->group=clashGroup;
this->mols[pos2+2]->group=clashGroup;
int fid1=path_tree[cur_branch][path_pos[cur_branch]].first;
int fid2=path_tree[cur_branch][path_pos[cur_branch]].second;
this->changeFrag2(pos1,fid1,frags);
this->changeFrag2(pos2,fid2,frags);
this->coordSysNt(pos1,&atoms[fid1]->atoms);
this->coordSysNt(pos2,&atoms[fid2]->atoms);
this->updateFull();
this->coordSysAT(pos1,&batoms[fid1]->atoms);
this->coordSysAT(pos2,&batoms[fid2]->atoms);
c=cur_branch+1;
cur_branch=-1;
}else{
//Remove any clash constraints 
int pos1=cons[cur_branch][0].first-2;
int pos2=cons[cur_branch][0].second;
this->mols[pos1+0]->group=unclashGroup;
this->mols[pos1+1]->group=unclashGroup;
this->mols[pos1+2]->group=unclashGroup;
//Set group clashGroup to include 
this->mols[pos2+0]->group=unclashGroup;
this->mols[pos2+1]->group=unclashGroup;
this->mols[pos2+2]->group=unclashGroup;
//Take a step back in the tree
cur_branch=cur_branch-1;
c=-1;
//Remove the branch
path_tree.pop_back();
path_pos.pop_back();
}
}


}
}
if(c==-1){printf("Failed to find perfect structure\n"); success=0;}else{
success=1;
Pdb * pdb = new Pdb(this,'j');
pdb->write(outfile);
pdb->clear();
delete pdb;
}
return success;
}

*/
int Polymere::generalCon(std::pair<int,int>  gencons,float * hydro1){
//Returns the constraint type that was found

//ctype = -1 Look at all constraint types
//ctype = 0-5 only consider the particular constraint
//Uses the doublet library to create a tertiary structure
//Finds fragments that score the best for a particular doublet library member 
//Returns 1 if all constraints are met and 0 if any of the constraints fail

//Type -2 - Didn't match any constraints
//Type -1 - C1 distance < 14A - used for filtering fragments
//Type 0 - Planar orientation but no defined hydrogen bonding
//W = Watson side(0)
//H = Hoogsteen side(1)
//S = Sugar side(2)

//CIS
//Type 1 - W:W 
//Type 2 - W:H 
//Type 3 - W:S 
//Type 4 - H:W 
//Type 5 - H:H 
//Type 6 - H:S 
//Type 7 - S:W 
//Type 8 - S:H 
//Type 9 - S:S 
//TRANS 10-18

//First position list 
//Type 19 - inward stacking (bases coordinate systems are in same direction) 
//Type 20 - outward stacking (bases coordinate systems are facing each other)
//Type 21 - upward stacking(bases coordinate systems are in opposite directions)

std::vector<double> params1;
std::vector<double> params2;
params2.push_back(0);//type;
params2.push_back(1);//rho;
params2.push_back(2);//theta1;
params2.push_back(3);//phi1;
params2.push_back(4);//phi2;
params2.push_back(5);//theta2;
params2.push_back(6);//dx_i;
params2.push_back(7);//dy_i;
params2.push_back(8);//dz_i;
params2.push_back(9);//typeRY;
params2.push_back(10);//zzang;
params1.push_back(0);//type;
params1.push_back(1);//rho;
params1.push_back(2);//theta1;
params1.push_back(3);//phi1;
params1.push_back(4);//phi2;
params1.push_back(5);//theta2;
params1.push_back(6);//dx_i;
params1.push_back(7);//dy_i;
params1.push_back(8);//dz_i;
params1.push_back(9);//typeRY;
params1.push_back(10);//zzang;
this->getParams(gencons.first,gencons.second,&params1);
this->getParams(gencons.second,gencons.first,&params2);

int cis=-1;
int inward=-1;
int conType=-2;
if(params1[1] < 14){
conType=-1;
	//Planar Orientation
	if(abs(params1[10]-PI)<0.35 || abs(params1[10])<0.35){
	conType=0;
//cis
//	1 = z-axis oriented in same direction 
//	0 = z-axis oriented in opposite directions 
		if(abs(params1[10]-PI)<0.35){cis=0;}else{cis=1;}
//inward
//	0 = z-axis are not pointed toward each other
//	1 = z-axis are pointed toward each other
		if(cis==0){
			if(params1[8] > 0 && params2[8] > 0){inward = 1;}
			if(params1[8] < 0 && params2[8] < 0){inward = 0;}
		}
//If z distance is within [2A,6A] and dx < 6A 
if(abs(params1[8])>1.5 && abs(params2[8])>1.5 && abs(params1[8])<6 && abs(params2[8])<6 && params1[6]<6){
//Inward stacking
if(cis==0 && inward==1){conType=19;
if(params1[9]==params2[9]){*hydro1=0;}
}
//Outward stacking
if(cis==0 && inward==0){conType=20;
if(params1[9]==params2[9]){*hydro1=0;}

}
//Upward stacking
if(cis==1){ conType=21;*hydro1=0;
if(params1[9]==params2[9]){*hydro1=0;}

}
}else{
	//check for hydrogen bonding
std::pair<int,int> faces;
faces.first=-1;
faces.second=-1;
float hydro=this->hydroType(gencons.first,gencons.second,&faces);
*hydro1 = hydro;
//printf("%8.3f: %d, %d\n",hydro,faces.first,faces.second);
if(faces.first==-1 || faces.second==-1){
conType=0;//planar with no hydroType
}else{

if(abs(params2[8])<1.5 && abs(params1[8])<1.5){
//Set minimal requirements for each type of interaction
if(faces.first==0 && faces.second==0){conType=1;}
if(faces.first==0 && faces.second==1){conType=2;}
if(faces.first==0 && faces.second==2){conType=3;}
if(faces.first==1 && faces.second==0){conType=4;}
if(faces.first==1 && faces.second==1){conType=5;}
if(faces.first==1 && faces.second==2){conType=6;}
if(faces.first==2 && faces.second==0){conType=7;}
if(faces.first==2 && faces.second==1){conType=8;}
if(faces.first==2 && faces.second==2){conType=9;}
if(cis==1){conType=conType+9;}
}else{
conType=0;
}
}
}
}else{
	if(abs(params1[10]-PI)<1.5 || abs(params2[10])<1.5){
	conType=0;
//cis
//	1 = z-axis oriented in same direction 
//	0 = z-axis oriented in opposite directions 
		if(abs(params1[10]-PI)<1.5){cis=0;}else{cis=1;}
//inward
//	0 = z-axis are not pointed toward each other
//	1 = z-axis are pointed toward each other
		if(cis==0){
			if(params1[8] > 0 && params2[8] > 0){inward = 1;}
			if(params1[8] < 0 && params2[8] < 0){inward = 0;}
		}
//If z distance is within [2A,6A] and dx < 6A 
if(abs(params1[8])>1.5 && abs(params2[8])>1.5 && abs(params1[8])<6 && abs(params2[8])<6 && params1[6]<6){
//Inward stacking
if(cis==0 && inward==1){conType=219;
if(params1[9]==params2[9]){*hydro1=0;}
}
//Outward stacking
if(cis==0 && inward==0){conType=220;
if(params1[9]==params2[9]){*hydro1=0;}

}
//Upward stacking
if(cis==1){ conType=221;*hydro1=0;
if(params1[9]==params2[9]){*hydro1=0;}

}
}else{
	//check for hydrogen bonding
std::pair<int,int> faces;
faces.first=-1;
faces.second=-1;
float hydro=this->hydroType(gencons.first,gencons.second,&faces);
*hydro1 = hydro;
//printf("%8.3f: %d, %d\n",hydro,faces.first,faces.second);
if(faces.first==-1 || faces.second==-1){
conType=0;//planar with no hydroType
}else{

if(abs(params2[8])<3 && abs(params1[8])<3){
//Set minimal requirements for each type of interaction
if(faces.first==0 && faces.second==0){conType=301;}
if(faces.first==0 && faces.second==1){conType=302;}
if(faces.first==0 && faces.second==2){conType=303;}
if(faces.first==1 && faces.second==0){conType=304;}
if(faces.first==1 && faces.second==1){conType=305;}
if(faces.first==1 && faces.second==2){conType=306;}
if(faces.first==2 && faces.second==0){conType=307;}
if(faces.first==2 && faces.second==1){conType=308;}
if(faces.first==2 && faces.second==2){conType=309;}
if(cis==1){conType=conType+9;}
}else{
conType=0;
}
}
}
}




}
}
return conType;
}

int Polymere::generalCons(std::vector<std::pair<std::pair<int,int> ,int> > gencons,float * tot_hydro){
//loops through all constraints to determine if they match
int allgood=1;
float hydro1;
for(int g = 0; g < gencons.size(); g++){
int conType=this->generalCon(std::pair<int,int>(gencons[g].first.first,gencons[g].first.second),&hydro1);
int conTypeDG=conType;
*tot_hydro=*tot_hydro+hydro1;
//If the constraint isn't met set allgood to 0
if(gencons[g].second==-1 && conType==-2){allgood=0;}
if(gencons[g].second==-9 && (conType < 1 || conType > 18)){allgood=0;}
if(gencons[g].second==-8 && (conType < 1)){allgood=0;}
if(gencons[g].second==-8 && (conType > 18 && conType < 301)){allgood=0;}
if(gencons[g].second==-8 && (conType > 318)){allgood=0;}
if(gencons[g].second>0 && gencons[g].second!=conType){allgood=0;}

}

return allgood;
}

int Polymere::consBuild1(int movepos,std::vector<std::pair<std::pair<int,int>,int> > gencons, const std::vector< Fragment *>  frags,Fragclass * dat,std::vector< Fragatoms *> atoms,std::vector<Fragall *> batoms,std::vector<std::vector<int> > loopO,std::vector<int> clashGroup,int max_iterations,int max_leaves,int group,std::vector<float> * con_score,std::vector<int> * constree,std::vector< std::vector<int> > * fragtree,int step,std::string ofile,int write_steps,float p_o_lower,float p_o_upper){
//Assembles even loops that bring together movepos-1 and movepos+6
//movepos 	- Position of the start of the loop
//frags   	- course-grain fragments 
//atoms   	- all atom fragments
//batoms  	- all atom fragments
//loopO   	- Database of looping fragments
//clashGroup 	- Group of nucleotides that are considered when looking at clashes
//group         - The number associated with the clash group
//con_score     - score of a particular loop conformation
//constree      - keeps track of the current number of attempts for a given step 
//fragtree      - fragments associated with a given step 
//step          - step 
//ofile         - prefix for the output file if each constraint step is to be written to a file
//write_step    - 0 = don't write, 1=write each step   

Polymere * pol2=new Polymere();
Polymere * pol1=new Polymere();
pol1->copy(this);

//set the clash group up
for(int m=0; m<pol1->mols.size(); m++){pol1->mols[m]->group=0;}
for(int m=0; m < pol1->mols.size() && m < clashGroup.size(); m++){
if(clashGroup[m]==0){pol1->mols[m]->group=0;}else{pol1->mols[m]->group=group;}
}

//Set the maximum iterations to the total size of the loopO database if it is <0 or larger than the database size
if(max_iterations < 0 || max_iterations > loopO.size()){max_iterations=loopO.size();}
int maxl=0;
int tot_it=0;
while(maxl<max_leaves && tot_it < max_iterations){
pol2->copy(pol1);
int rand1 = rand() % (int) loopO.size();
//printf("%d %d\n",rand1,tot_it);
pol2->changeFrag2(movepos-1,loopO[rand1][0],frags);
pol2->changeFrag2((movepos+1),loopO[rand1][1],frags);
pol2->changeFrag2((movepos+3),loopO[rand1][2],frags);
pol2->changeFrag2((movepos+5),loopO[rand1][3],frags);
pol2->coordSysNt((movepos-1),&atoms[loopO[rand1][0]]->atoms);
pol2->coordSysNt((movepos+1),&atoms[loopO[rand1][1]]->atoms);
pol2->coordSysNt((movepos+3),&atoms[loopO[rand1][2]]->atoms);
pol2->coordSysNt((movepos+5),&atoms[loopO[rand1][3]]->atoms);
pol2->updateFull();
pol2->coordSysAT((movepos-1),&batoms[loopO[rand1][0]]->atoms);
pol2->coordSysAT((movepos+1),&batoms[loopO[rand1][1]]->atoms);
pol2->coordSysAT((movepos+3),&batoms[loopO[rand1][2]]->atoms);
pol2->coordSysAT((movepos+5),&batoms[loopO[rand1][3]]->atoms);

//IF no clashes exist add to fragtree
if(pol2->groupCheckBackBone(group,p_o_lower,p_o_upper)<1 && pol2->groupCheckClashes(group)<1){
//printf("Conecting %d %d\n",movepos[at][0]-1,movepos[at][0]+7);
//Add to fragtree
float hydroL=0;
std::vector<std::pair<std::pair<int,int>,int> > tmpcons;
std::pair<std::pair<int,int>,int> icons(std::pair<int,int>((movepos-1),movepos+7),gencons[0].second);
tmpcons.push_back(icons);
for(int a1=1; a1<gencons.size(); a1++){tmpcons.push_back(gencons[a1]);}
int conout=pol2->generalCons(tmpcons,&hydroL);
if(conout==1){
con_score->push_back(hydroL);//All loops treated equal
constree->push_back(step);
fragtree->push_back(loopO[rand1]);
maxl++;
//If write_step == 1 print out PDB 
if(write_steps==1){
std::string ofile2=ofile;
ofile2+=boost::lexical_cast<std::string>(step);
ofile2.append(".pdb");
Pdb * npdb = new Pdb(pol2,'J');
npdb->write(ofile.c_str());
npdb->clear();
delete npdb;
}
}
}
tot_it++;
}
return maxl;
}

int Polymere::consBuild2(int movepos,std::vector<std::pair<std::pair<int,int>,int> > gencons, const std::vector< Fragment *>  frags,Fragclass * dat,std::vector< Fragatoms *> atoms,std::vector<Fragall *> batoms,std::vector<std::vector<int> > loopE,std::vector<int> clashGroup,int max_iterations,int max_leaves,int group,std::vector<float> * con_score,std::vector<int> * constree,std::vector< std::vector<int> > * fragtree,int step,std::string ofile,int write_steps,float p_o_lower,float p_o_upper){
//Assembles even loops that bring together movepos-1 and movepos+6
//movepos 	- Position of the start of the loop
//frags   	- course-grain fragments 
//atoms   	- all atom fragments
//batoms  	- all atom fragments
//loopO   	- Database of looping fragments
//clashGroup 	- Group of nucleotides that are considered when looking at clashes
//group         - The number associated with the clash group
//con_score     - score of a particular loop conformation
//constree      - keeps track of the current number of attempts for a given step 
//fragtree      - fragments associated with a given step 
//step          - step 
//ofile         - prefix for the output file if each constraint step is to be written to a file
//write_step    - 0 = don't write, 1=write each step   
Polymere * pol2=new Polymere();
Polymere * pol1=new Polymere();
pol1->copy(this);
//set the clash group up
for(int m=0; m<pol1->mols.size(); m++){pol1->mols[m]->group=0;}
for(int m=0; m<pol1->mols.size() && m < clashGroup.size(); m++){
if(clashGroup[m]==0){pol1->mols[m]->group=0;}else{pol1->mols[m]->group=group;}
}


for(int m=0; m<pol1->mols.size(); m++){printf("%d",pol1->mols[m]->group);}
printf("\n");
//Set the maximum iterations to the total size of the loopO database if it is <0 or larger than the database size
if(max_iterations < 0 || max_iterations > loopE.size()){max_iterations=loopE.size();}
int maxl=0;
int tot_it=0;
while(maxl<max_leaves && tot_it < max_iterations){
pol2->copy(pol1);
int rand1 = rand() % (int) loopE.size();
//printf("%d %d\n",rand1,tot_it);
pol2->changeFrag2(movepos-1,loopE[rand1][0],frags);
pol2->changeFrag2((movepos+1),loopE[rand1][1],frags);
pol2->changeFrag2((movepos+3),loopE[rand1][2],frags);
pol2->changeFrag2((movepos+5),loopE[rand1][3],frags);
pol2->coordSysNt((movepos-1),&atoms[loopE[rand1][0]]->atoms);
pol2->coordSysNt((movepos+1),&atoms[loopE[rand1][1]]->atoms);
pol2->coordSysNt((movepos+3),&atoms[loopE[rand1][2]]->atoms);
pol2->coordSysNt((movepos+5),&atoms[loopE[rand1][3]]->atoms);
pol2->updateFull();
pol2->coordSysAT((movepos-1),&batoms[loopE[rand1][0]]->atoms);
pol2->coordSysAT((movepos+1),&batoms[loopE[rand1][1]]->atoms);
pol2->coordSysAT((movepos+3),&batoms[loopE[rand1][2]]->atoms);
pol2->coordSysAT((movepos+5),&batoms[loopE[rand1][3]]->atoms);

//IF no clashes exist add to fragtree
if(pol2->groupCheckBackBone(group,p_o_lower,p_o_upper)<1 && pol2->groupCheckClashes(group)<1){
//Add to fragtree
float hydroL=0;
std::vector<std::pair<std::pair<int,int>,int> > tmpcons;
std::pair<std::pair<int,int>,int> icons(std::pair<int,int>((movepos-1),movepos+6),gencons[0].second);
tmpcons.push_back(icons);
for(int a1=1; a1<gencons.size(); a1++){tmpcons.push_back(gencons[a1]);}
int conout=pol2->generalCons(tmpcons,&hydroL);
if(conout==1){
con_score->push_back(hydroL);
constree->push_back(step);
fragtree->push_back(loopE[rand1]);
maxl++;

if(write_steps==1){
std::string ofile2=ofile;
ofile2+=boost::lexical_cast<std::string>(step);
ofile2.append(".pdb");
Pdb * npdb = new Pdb(pol2,'J');
npdb->write(ofile.c_str());
npdb->clear();
delete npdb;
}
}
}
tot_it++;
}
return maxl;
}

int Polymere::consBuild3(int movepos1,int movepos2,std::vector<std::pair<std::pair<int,int>,int> > gencons, const std::vector< Fragment *>  frags,Fragclass * dat,std::vector< Fragatoms *> atoms,std::vector<Fragall *> batoms,std::vector<int> clashGroup,int max_iterations,int max_leaves,int group,std::vector<float> * con_score,std::vector<int> * constree,std::vector< std::vector<int> > * fragtree,int step,std::string ofile,int write_steps,float p_o_lower, float p_o_upper){
//For use with two movement positions. Builds helices
Polymere * pol2=new Polymere();
Polymere * pol1=new Polymere();
pol1->copy(this);
//update Clash group
for(int m=0; m<pol1->mols.size(); m++){pol1->mols[m]->group=0;}
for(int m=0; m<pol1->mols.size() && m < clashGroup.size(); m++){
if(clashGroup[m]==0){pol1->mols[m]->group=0;}else{pol1->mols[m]->group=group;}
printf("%d",pol1->mols[m]->group);
}
printf("\n");
for(int m=0; m<pol1->mols.size() && m < clashGroup.size(); m++){
if(m==movepos1 || m==movepos2){
printf("|");}else{
printf("%c",pol1->mols[m]->seq);
}
}

printf("\nGetting Frags\n");
std::vector<int> pos1Fids;
std::vector<int> pos2Fids;
std::vector<Matches> mat; 
mat=pol1->chooseMovesPairs(movepos2,movepos1,dat);
printf("Total Pair Frags: %d\n",mat.size());
//pol1->chooseAllSingle(movepos1,dat,&pos1Fids);
//pol1->chooseAllSingle(movepos2,dat,&pos2Fids);
int totpairs=mat.size();
printf("Total Possible Fragment Pairs %d-%d: %d %d\n",movepos1,movepos2,mat.size());
printf("Initial clashes %d %d\n",pol1->groupCheckBackBone(group,p_o_lower,p_o_upper), pol1->groupCheckClashes(group));

//Filter for pairs that get the distance correct without clashes
int tot=0;
int total_out=0;
int max_its=max_iterations;
if(max_iterations > mat.size()){max_its=mat.size();}

std::vector<int> already_used;
while(tot < max_its  && total_out < max_leaves){
pol2->copy(pol1);

int rand1=rand() % (int) mat.size();
int go=1;
int id=rand1;
if(already_used.size()>0){ 
if(binary_search(already_used.begin(),already_used.end(),id)){
go=0;
}else{
already_used.push_back(id);
sort(already_used.begin(),already_used.end());
}
}else{
already_used.push_back(id);
}
//Pick random pair
//Change Frag and check for clashes and distance
if(go==1 && batoms[mat[rand1].frag2]->check()==1 && batoms[mat[rand1].frag1]->check()==1){   //If fragment is any good use it
pol2->changeFrag2(movepos1,mat[rand1].frag2,frags);
pol2->coordSysNt(movepos1,&atoms[mat[rand1].frag2]->atoms);
pol2->updateFull();
pol2->coordSysAT(movepos1,&batoms[mat[rand1].frag2]->atoms);

pol2->changeFrag2(movepos2,mat[rand1].frag1,frags);
pol2->coordSysNt(movepos2,&atoms[mat[rand1].frag1]->atoms);
pol2->updateFull();
pol2->coordSysAT(movepos2,&batoms[mat[rand1].frag1]->atoms);
//printf("clash: %d, backbone: %d\n",pol2->groupCheckClashes(group), pol2->groupCheckBackBone(group));
if(pol2->groupCheckBackBone(group,p_o_lower,p_o_upper) < 1 && pol2->groupCheckClashes(group)<1){
printf(".");
//If constraints are met then save to frag_set
int conout=1;
float hydro_score=0;
if(gencons[0].second!=0){conout=pol2->generalCons(gencons,&hydro_score);}
if(conout==1){
printf("\nfound one: %8.3f\n",hydro_score);
//add to fragset
con_score->push_back(hydro_score);//save moves score
std::vector<int> fragbranch;
fragbranch.push_back(mat[rand1].frag2);
fragbranch.push_back(mat[rand1].frag1);
fragtree->push_back(fragbranch);
constree->push_back(step);
total_out++;
//Write to temporaty output file 
if(write_steps==1){
std::string ofile2=ofile;
ofile2+=boost::lexical_cast<std::string>(step);
ofile2.append(".pdb");
Pdb * npdb = new Pdb(pol2,'J');
npdb->write(ofile2.c_str());
npdb->clear();
delete npdb;
}
}

}
}
tot++;
}
return total_out; 
}

int Polymere::consBuild4(std::vector<int> movepos,std::vector<std::pair<std::pair<int,int>,int> > gencons, const std::vector< Fragment *>  frags,Fragclass * dat,std::vector< Fragatoms *> atoms,std::vector<Fragall *> batoms,std::vector<int> clashGroup,int max_iterations,int max_leaves,int group,std::vector<float> * con_score,std::vector<int> * constree,std::vector< std::vector<int> > * fragtree,int step,std::string ofile,int write_steps,float min_res, int max_clust,float p_o_lower, float p_o_upper){
//All moves are made in 5' to 3' direction
//The clash group includes the original clash group plus all fragments 3' of current fragment
Polymere * pol2=new Polymere();
Polymere * pol1=new Polymere();
pol1->copy(this);
sort(movepos.begin(),movepos.end());
printf("Movement positions: ");
for(int j=0; j<movepos.size(); j++){printf(" %d",movepos[j]);}
printf("\nInitial clashes %d %d\n",pol1->groupCheckBackBone(group,p_o_lower,p_o_upper), pol1->groupCheckClashes(group));
//Filter for pairs that get the distance correct without clashes
std::vector<std::vector<int> > frag_set;
int mg=0;
int total_out=0;
while(mg < movepos.size()){
pol1->copy(this);
for(int g=0; g<pol1->mols.size(); g++){pol1->mols[g]->group=0;}
printf("Getting Frags: Pos %d\n",movepos[mg]);
//Update clashgroup
if(mg < movepos.size()-1){
for(int g=0; g<pol1->num_mols && g<clashGroup.size(); g++){
if(clashGroup[g]!=0 && g < movepos[mg+1]){pol1->mols[g]->group=group;}
} 
}else{
for(int g=0; g<pol1->num_mols && g<clashGroup.size(); g++){
if(clashGroup[g]==1){pol1->mols[g]->group=group;}
} 
}
for(int g=0; g<=mg; g++){
pol1->mols[movepos[g]+0]->group=group;
pol1->mols[movepos[g]+1]->group=group;
pol1->mols[movepos[g]+2]->group=group;
}
printf("Clash Group: ");
for(int m=0; m<pol1->mols.size(); m++){
printf("%d.",pol1->mols[m]->group);
}
printf("\n");

for(int m=0; m<pol1->mols.size(); m++){
if(m==movepos[mg]){
printf("|.");}else{
printf("%c.",pol1->mols[m]->seq);
}

}


///
if(mg > 0 && mg < (movepos.size()-1) && frag_set.size() > 0){
Dynaclust * dy = new Dynaclust(max_clust,min_res);
std::vector<int> pos1Fids;
std::vector<std::vector<int> > frag_set1;
pol1->chooseAllSingle(movepos[mg],dat,&pos1Fids);
int totpairs;
totpairs=pos1Fids.size()*frag_set.size();
int max_atmps=max_iterations;
if(max_iterations > totpairs){max_atmps=totpairs;}
printf("Total Possible Fragment Pairs Clustering:%d/%d/%8.3f\n",max_atmps,totpairs,min_res);
printf("Initial clashes %d %d\n",pol1->groupCheckBackBone(group,p_o_lower,p_o_upper), pol1->groupCheckClashes(group));
//Filter for pairs that get the distance correct without clashes
int tot=0;
int total_atms=0;
int iter=0;
int addclust=0;
std::vector<int> already_used;
while(tot < max_clust && total_atms < max_atmps){
total_atms++;
pol2->copy(pol1);
//Pick random pair
int rand1=rand() % (int)  pos1Fids.size();
int rand2=rand() % (int)  frag_set.size();
int rand_couple = rand1*1000000+rand2;
//make sure you didn't use either from pos1 or fragset
int go=1;
if(already_used.size()>0){ 
if(binary_search(already_used.begin(),already_used.end(),rand_couple)){
go=0;
}else{
already_used.push_back(rand_couple);
sort(already_used.begin(),already_used.end());
}
}else{
already_used.push_back(rand_couple);
}
//Change Frag and check for clashes and distance
if(go==1 && batoms[pos1Fids[rand1]]->check()==1){
for(int p=0; p<frag_set[rand2].size(); p++){
pol2->changeFrag2(movepos[p],frag_set[rand2][p],frags);
pol2->coordSysNt(movepos[p],&atoms[frag_set[rand2][p]]->atoms);
pol2->updateFull();
pol2->coordSysAT(movepos[p],&batoms[frag_set[rand2][p]]->atoms);
}

pol2->changeFrag2(movepos[mg],pos1Fids[rand1],frags);
pol2->coordSysNt(movepos[mg],&atoms[pos1Fids[rand1]]->atoms);
pol2->updateFull();
pol2->coordSysAT(movepos[mg],&batoms[pos1Fids[rand1]]->atoms);

if(pol2->groupCheckBackBone(group,p_o_lower,p_o_upper) < 1 && pol2->groupCheckClashes(group)<1){
//For all pairs of contraints that are currently in the clashGroup vector
//Begin constraint clustering
//Make sure new polymere is aligned to the first movement positon
pol2->alignStruct2(movepos[0]);
//Get cooridinates of the last 
float xc = pol2->mols[(movepos[mg]+3)]->x;
float yc = pol2->mols[(movepos[mg]+3)]->y;
float zc = pol2->mols[(movepos[mg]+3)]->z;
//printf("%d %d: %8.3f\t%8.3f\t%8.3f\n",pos1Fids[rand1],frag_set[rand2][0],xc,yc,zc);
addclust=dy->add(xc,yc,zc);
if(addclust==1){
//reset tot so more can be identified
tot++;
//Add the fragments to the frag_set
std::vector<int> frags_out;
for(int j=0; j<frag_set[rand2].size(); j++){frags_out.push_back(frag_set[rand2][j]);}
frags_out.push_back(pos1Fids[rand1]);
frag_set1.push_back(frags_out);
}
else{
if(addclust==-1){tot=max_clust;}
}
}
}
}//end while loop
frag_set.clear();
for(int j=0; j<frag_set1.size(); j++){frag_set.push_back(frag_set1[j]);}
frag_set1.clear();
dy->pos_nc.clear();
delete dy;
}else{
//////////////////////////////////////////////////////////////
//without clustering
//////////////////////////////////////////////////////////////
if(mg==0 && mg < (movepos.size()-1)){
frag_set.clear();
printf("Getting Frags\n");
std::vector<int> pos1Fids;
pol1->chooseAllSingle(movepos[mg],dat,&pos1Fids);
int totpairs=pos1Fids.size();
printf("Total Possible Fragment Pairs: %d\n",totpairs);
printf("Initial clashes %d %d\n",pol1->groupCheckBackBone(group,p_o_lower,p_o_upper), pol1->groupCheckClashes(group));
//Filter for pairs that get the distance correct without clashes
int tot=0;
std::vector<int> already_used;
while(tot < pos1Fids.size()){
pol2->copy(pol1);
int rand1=rand() % (int) pos1Fids.size();
int go=1;
if(already_used.size()>0){ 
if(binary_search(already_used.begin(),already_used.end(),rand1)){
go=0;
}else{
already_used.push_back(rand1);
sort(already_used.begin(),already_used.end());
}
}else{
already_used.push_back(rand1);
}
//Pick random pair
//Change Frag and check for clashes and distance
if(go==1 && batoms[pos1Fids[rand1]]->check()==1){   //If fragment is any good use it
pol2->changeFrag2(movepos[mg],pos1Fids[rand1],frags);
pol2->coordSysNt(movepos[mg],&atoms[pos1Fids[rand1]]->atoms);
pol2->updateFull();
pol2->coordSysAT(movepos[mg],&batoms[pos1Fids[rand1]]->atoms);

if(pol2->groupCheckBackBone(group,p_o_lower,p_o_upper) < 1 && pol2->groupCheckClashes(group)<1){
std::vector<int> frags_out;
frags_out.push_back(pos1Fids[rand1]);
frag_set.push_back(frags_out);
}
}
tot++;
}
}else{
//Final Constraints
std::vector<std::vector<int> > frag_set1;
std::vector<int> pos1Fids;
pol1->chooseAllSingle(movepos[mg],dat,&pos1Fids);
int totpairs;
if(mg==0){
totpairs=pos1Fids.size();
printf("Total Possible Fragment Subs Final: %d\n",totpairs);
printf("Initial clashes %d %d\n",pol1->groupCheckBackBone(group,p_o_lower,p_o_upper), pol1->groupCheckClashes(group));
//Filter for pairs that get the distance correct without clashes
int tot=0;
for(int m2=0; m2 < pos1Fids.size(); m2++){
pol2->copy(pol1);
int rand1=rand() % (int) pos1Fids.size();
//Change Frag and check for clashes and distance
if(batoms[pos1Fids[rand1]]->check()==1){
pol2->changeFrag2(movepos[mg],pos1Fids[rand1],frags);
pol2->coordSysNt(movepos[mg],&atoms[pos1Fids[rand1]]->atoms);
pol2->updateFull();
pol2->coordSysAT(movepos[mg],&batoms[pos1Fids[rand1]]->atoms);
if(pol2->groupCheckBackBone(group,p_o_lower,p_o_upper) < 1 && pol2->groupCheckClashes(group) < 1){
//For all pairs of contraints that are currently in the clashGroup vector
//Begin constraint clustering
//Make sure all hydroPairs are > 1
//Write Everything to a file 
//Find all that meet the contraints
int conout=1;
float hydro_score=0;
if(gencons[0].second!=0){conout=pol2->generalCons(gencons,&hydro_score);}

if(conout==1){
printf("\nfound one: %8.3f\n",hydro_score);
//add to fragset
con_score->push_back(hydro_score);//save moves score
std::vector<int> fragbranch;
fragbranch.push_back(pos1Fids[rand1]);
fragtree->push_back(fragbranch);
constree->push_back(step);
tot++;
total_out++;
if(write_steps==1){
std::string ofile2=ofile;
ofile2+=boost::lexical_cast<std::string>(step);
ofile2.append(".pdb");
Pdb * npdb = new Pdb(pol2,'J');
npdb->write(ofile2.c_str());
npdb->clear();
delete npdb;
}
}
}
}
if(tot>max_leaves){break;}
}
frag_set.clear();
}else{
totpairs=pos1Fids.size()*frag_set.size();
if(totpairs > 0){
printf("Total Possible Fragment Pairs Final: %d\n",totpairs);
printf("Initial clashes %d %d\n",pol1->groupCheckBackBone(group,p_o_lower,p_o_upper), pol1->groupCheckClashes(group));
//Filter for pairs that get the distance correct without clashes
int tot=0;
int total_atms=0;
std::vector<int> already_used;
int total_iteration=max_iterations;
if(max_iterations > totpairs){total_iteration=totpairs;}
while(tot < max_leaves  && total_atms<total_iteration){
total_atms++;
pol2->copy(pol1);
//Pick random pair
int rand1=rand() % (int) pos1Fids.size();
int rand2=rand() % (int) frag_set.size();
int rand_couple = rand1*1000000+rand2;
//make sure you didn't use either from pos1 or fragset
int go=1;
if(already_used.size()>0){ 
if(binary_search(already_used.begin(),already_used.end(),rand_couple)){
go=0;
}else{
already_used.push_back(rand_couple);
sort(already_used.begin(),already_used.end());
}
}else{
already_used.push_back(rand_couple);
}

//Change Frag and check for clashes and distance
if(go==1 && batoms[pos1Fids[rand1]]->check()==1){

for(int p=0; p<frag_set[rand2].size(); p++){
pol2->changeFrag2(movepos[p],frag_set[rand2][p],frags);
pol2->coordSysNt(movepos[p],&atoms[frag_set[rand2][p]]->atoms);
pol2->updateFull();
pol2->coordSysAT(movepos[p],&batoms[frag_set[rand2][p]]->atoms);
}
pol2->changeFrag2(movepos[mg],pos1Fids[rand1],frags);
pol2->coordSysNt(movepos[mg],&atoms[pos1Fids[rand1]]->atoms);
pol2->updateFull();
pol2->coordSysAT(movepos[mg],&batoms[pos1Fids[rand1]]->atoms);


if(pol2->groupCheckBackBone(group,p_o_lower,p_o_upper) < 1 && pol2->groupCheckClashes(group) < 1){
//For all pairs of contraints that are currently in the clashGroup vector
//Begin constraint clustering
//Make sure all hydroPairs are > 1
//Write Everything to a file 
//Find all that meet the contraints
int conout=1;
float hydro_score=0;
if(gencons[0].second!=0){conout=pol2->generalCons(gencons,&hydro_score);}

if(conout==1){
printf("\nfound one: %8.3f\n",hydro_score);
tot++;
std::vector<int> fragbranch;
for(int p=0; p<frag_set[rand2].size(); p++){
fragbranch.push_back(frag_set[rand2][p]);
}
fragbranch.push_back(pos1Fids[rand1]);
fragtree->push_back(fragbranch);
constree->push_back(step);
con_score->push_back(hydro_score);
total_out++;
if(write_steps==1){
std::string ofile2=ofile;
ofile2+=boost::lexical_cast<std::string>(step);
ofile2.append(".pdb");
Pdb * npdb = new Pdb(pol2,'J');
npdb->write(ofile2.c_str());
npdb->clear();
delete npdb;
}



}
}
}


}//end of while loop
}
frag_set.clear();
}

}
}
mg++;
}
return total_out; 
}

int Polymere::assemble(std::vector<std::vector<int> > movepos,std::vector<std::vector<std::pair<std::pair<int,int>,int> > > gencons,std::vector<std::vector<int> > clashgroups, const std::vector< Fragment *>  frags,Fragclass * dat,std::vector< Fragatoms *> atoms,std::vector<Fragall *> batoms,std::vector<std::vector<int> > loopO,std::vector<std::vector<int> > loopE, int max_iterations,int max_leaves, int group,std::string ofile,int write_steps,float min_res, int max_clust,int start_con,float p_o_lower,float p_o_upper){

std::vector<int> constree;
std::vector<std::vector<int> > fragtree;
std::vector<float> con_score;
Polymere * pol1=new Polymere();
pol1->copy(this);
//Set group variable on pol1
//Only consider clashes for postions that appear in previous 
//movement groups and the current fragment set
//Contains sets of fragments that have no clashes or backbone breaks
int at =start_con;
//If the starting constraint is greater than 0 move through the PDB unpdating the class group
int all_atmps=0;
while(at<movepos.size() && all_atmps < 1000){
//if the program tries to move backward stop it 
if(at < start_con){break;}
all_atmps++;
printf("Constraint Set: %d\n",at);
//Add new frags to clashGroup
//Types of movement groups
//Odd Loop - Include all loop NTs in clash group
//Even Loop - Include all loop NTs in cluster
//Clusters - Include only the members of the fragments updated as clash group
if(gencons[at][0].first.first==-2 && gencons[at][0].first.second==-2){
pol1->consBuild2(movepos[at][0],gencons[at],frags,dat,atoms,batoms,loopE,clashgroups[at],max_iterations,max_leaves,5,&con_score,&constree,&fragtree,at,ofile,write_steps,p_o_lower,p_o_upper);
}else{
if(gencons[at][0].first.first==-3 && gencons[at][0].first.second==-3){
//Add everything between and including postion and position+6
pol1->consBuild1(movepos[at][0],gencons[at],frags,dat,atoms,batoms,loopO,clashgroups[at],max_iterations,max_leaves,5,&con_score,&constree,&fragtree,at,ofile,write_steps,p_o_lower,p_o_upper);
}else{
//Add all fragments in the movement group to the clash group
if(movepos[at].size()==2){
pol1->consBuild3(movepos[at][0],movepos[at][1],gencons[at],frags,dat,atoms,batoms,clashgroups[at],max_iterations,max_leaves,5,&con_score,&constree,&fragtree,at,ofile,write_steps,p_o_lower,p_o_upper);
}else{
pol1->consBuild4(movepos[at],gencons[at],frags,dat,atoms,batoms,clashgroups[at],max_iterations,max_leaves,5,&con_score,&constree,&fragtree,at,ofile,write_steps,min_res,max_clust,p_o_lower,p_o_upper);
}
}
}

//print out current fragtree and contree
/*
for(int a=0; a<constree.size(); a++){
printf("%8.3f %d: ",con_score[a],constree[a]);
for(int t=0; t<movepos[constree[a]].size();t++){printf("%d ",movepos[constree[a]][t]);}
for(int t=0; t<fragtree[a].size();t++){printf("%d ",fragtree[a][t]);}
printf("\n");
}
*/
//Remove it from all the constraints
if(constree.size()>0){
int at1=constree.back();
printf("at1: %d,  at: %d\n",at1,at);
at=at1;
//Look at all constraints for the current step and find the one that has best score
int max_pos=(int)constree.size()-1;
float max_score =con_score[max_pos];

int cur_atmps=0;
for(int ia=(constree.size()-1); ia >=0; ia--){
if(constree[ia]!=at){break;}else{
cur_atmps++;
if(max_score < con_score[ia]){max_score = con_score[ia]; max_pos=ia;
}
}
}

int a = max_pos;
printf("%d Using:%d %d %d %8.3f %d: ",cur_atmps,con_score.size(), constree.size(),fragtree.size(),con_score[a],constree[a]);
for(int t=0; t<movepos[constree[a]].size();t++){printf("%d ",movepos[constree[a]][t]);}
for(int t=0; t<fragtree[a].size();t++){printf("%d ",fragtree[a][t]);}
printf("\n");
std::vector<int> update_frags=fragtree[max_pos];
if((gencons[at][0].first.first==-2 && gencons[at][0].first.second==-2) || (gencons[at][0].first.first==-3 && gencons[at][0].first.second==-3)){

for(int p=0; p<update_frags.size(); p++){pol1->changeFrag2((movepos[at][0]+(2*p)-1),update_frags[p],frags);}
for(int p=0; p<update_frags.size(); p++){pol1->coordSysNt((movepos[at][0]+(2*p)-1),&atoms[update_frags[p]]->atoms);}
pol1->updateFull();
for(int p=0; p<update_frags.size(); p++){pol1->coordSysAT((movepos[at][0]+(2*p)-1),&batoms[update_frags[p]]->atoms);}

}else{
for(int p=0; p<update_frags.size(); p++){
pol1->changeFrag2(movepos[at][p],update_frags[p],frags);
pol1->coordSysNt(movepos[at][p],&atoms[update_frags[p]]->atoms);
pol1->updateFull();
pol1->coordSysAT(movepos[at][p],&batoms[update_frags[p]]->atoms);
}
}
this->copy(pol1);
//Remove last item from each vector
constree.erase(constree.begin()+max_pos);
fragtree.erase(fragtree.begin()+max_pos);
con_score.erase(con_score.begin()+max_pos);
at++;
}else{
at=0;
}

}//End of while loop
this->copy(pol1);
int out=1;
if(all_atmps>=1000){
out=0;
}

return out;

}

int Polymere::buildBuldge(int s1, int s2,std::vector<int> * small_side,const std::vector< Fragment *>  frags,Fragclass * dat,std::vector< Fragatoms *> atoms,std::vector<Fragall *> batoms,int s_side,std::string outfile){
//Builds a helical section of the polymere given the start and end of two strands 5'-3' 
//and a set of constraints for each movement group
//genPairedStacks - 	uses stacked fragments only
//			identifies pairs of fragments that don't clash
//			and abide by particular contraints
// RULES:
// 1) Any assymetric buldge can only have a ratio less than 2:1
// 2) Small side is a vector describing the pairings
//Generate constraints map
std::vector< std::vector< std::pair<int,int> > > cons;
//Check that lenght small_side.size + small_size==1 = e2-s2
int lcheck=0;
int success;
int m=0;
if(s_side==1){//Small side is position 1
for(int j=0; j<(small_side->size()-2); j++){
std::vector< std::pair<int,int> > bpcons;
int k=0;
int icr=1;	
int icr1=1;	
	if(small_side->at(j)==0){
		std::pair<int,int> scon((s1-j+0),(s2+m+0));	
		bpcons.push_back(scon);
		icr=1;
		icr1=1;
		}else{
		std::pair<int,int> scon1((s1-j+0),(s2+m+0));	
		bpcons.push_back(scon1);
		std::pair<int,int> scon2((s1-j+0),(s2+m+1));	
		bpcons.push_back(scon2);
		icr=2;
		icr1=2;
		}
	if(small_side->at(j+1)==0){
		std::pair<int,int> scon((s1-j-1),(s2+m+icr));	
		bpcons.push_back(scon);
		icr++;
		}else{
		std::pair<int,int> scon1((s1-j-1),(s2+m+icr));	
		bpcons.push_back(scon1);
		std::pair<int,int> scon2((s1-j-1),(s2+m+icr+1));	
		bpcons.push_back(scon2);
		icr=icr+2;
	}
	if(small_side->at(j+2)==0){
		std::pair<int,int> scon((s1-j-2),(s2+m+icr));	
		bpcons.push_back(scon);
		}else{
		std::pair<int,int> scon1((s1-j-2),(s2+m+icr));	
		bpcons.push_back(scon1);
		std::pair<int,int> scon2((s1-j-2),(s2+m+icr+1));	
		bpcons.push_back(scon2);
	}
m=m+icr1;
cons.push_back(bpcons);
}
}else{
if(s_side==2){
for(int j=0; j<(small_side->size()-2); j++){
std::vector< std::pair<int,int> > bpcons;
int k=0;
int icr=1;	
int icr1=1;	
	if(small_side->at(j)==0){
		std::pair<int,int> scon((s1-m+0),(s2+j+0));	
		bpcons.push_back(scon);
		icr=1;
		icr1=1;
		}else{
		std::pair<int,int> scon1((s1-m+0),(s2+j+0));	
		bpcons.push_back(scon1);
		std::pair<int,int> scon2((s1-m-1),(s2+j+0));	
		bpcons.push_back(scon2);
		icr=2;
		icr1=2;
		}
	if(small_side->at(j+1)==0){
		std::pair<int,int> scon((s1-m-icr),(s2+j+1));	
		bpcons.push_back(scon);
		icr++;
		}else{
		std::pair<int,int> scon1((s1-m-icr),(s2+j+1));	
		bpcons.push_back(scon1);
		std::pair<int,int> scon2((s1-m-icr-1),(s2+j+1));	
		bpcons.push_back(scon2);
		icr=icr+2;
	}
	if(small_side->at(j+2)==0){
		std::pair<int,int> scon((s1-m-icr),(s2+j+2));	
		bpcons.push_back(scon);
		}else{
		std::pair<int,int> scon1((s1-m-icr),(s2+j+2));	
		bpcons.push_back(scon1);
		std::pair<int,int> scon2((s1-m-icr-1),(s2+j+2));	
		bpcons.push_back(scon2);
	}
m=m+icr1;
cons.push_back(bpcons);
}


}else{
printf("Small side must be either 1 or 2 for buildBuldge\n");
return 0;
}

}
//Print out constraints for debugging purpose
//for(int j=0; j<cons.size(); j++){
//	printf("pos %d\n",j);
//	for(int m=0; m<cons[j].size(); m++){
//		printf("\t%d-%d\n",cons[j][m].first,cons[j][m].second);
//	}
//}

//Using cons build the helix.
//The first pair of positions in cons for each step is the starting 
//positions for the fragments
//5->3 starts at cons[i][0].second
//3->5 starts at cons[i][0].first - 2
std::vector<std::vector< std::pair<int,int> > > path_tree;
std::vector<int> path_pos;
int c=0;
while(path_tree.size() < cons.size() && c > -1){//Continue to extend the tree until all constraints are met
int pos1=cons[c][0].first-2;
int pos2=cons[c][0].second;
//Update the clash group to include those in this movement groups set of constraints
int clashGroup=5;
int unclashGroup=0;
std::vector< std::pair<int,int> > branch;
int goforward=this->genPairedCons(clashGroup,&cons[c],frags,dat,atoms,batoms,&branch);
if(goforward==1){
printf("Moving Forward: %d\n",branch.size());
path_pos.push_back((branch.size()-1)); //set current position to the last pair
path_tree.push_back(branch);
c++;//Move to the next set of constraints
}else{
printf("Moving backward:");
//Remove any clash constraints 
this->mols[pos1+0]->group=unclashGroup;
this->mols[pos1+1]->group=unclashGroup;
this->mols[pos1+2]->group=unclashGroup;
//Set group clashGroup to include 
this->mols[pos2+0]->group=unclashGroup;
this->mols[pos2+1]->group=unclashGroup;
this->mols[pos2+2]->group=unclashGroup;
//Remove clash constraints from the current group
//move back in path_tree
int cur_branch=c-1;
c=-1;
while(cur_branch>=0){
printf("Branch:%d  Position:%d\n",cur_branch,path_pos[cur_branch]);
if(path_pos[cur_branch]>=0){
path_pos[cur_branch]=path_pos[cur_branch]-1;
//move the positions of the polymere
int pos1=cons[cur_branch][0].first-2;
int pos2=cons[cur_branch][0].second;
this->mols[pos1+0]->group=clashGroup;
this->mols[pos1+1]->group=clashGroup;
this->mols[pos1+2]->group=clashGroup;
//Set group clashGroup to include 
this->mols[pos2+0]->group=clashGroup;
this->mols[pos2+1]->group=clashGroup;
this->mols[pos2+2]->group=clashGroup;
int fid1=path_tree[cur_branch][path_pos[cur_branch]].first;
int fid2=path_tree[cur_branch][path_pos[cur_branch]].second;
this->changeFrag2(pos1,fid1,frags);
this->changeFrag2(pos2,fid2,frags);
this->coordSysNt(pos1,&atoms[fid1]->atoms);
this->coordSysNt(pos2,&atoms[fid2]->atoms);
this->updateFull();
this->coordSysAT(pos1,&batoms[fid1]->atoms);
this->coordSysAT(pos2,&batoms[fid2]->atoms);
c=cur_branch+1;
cur_branch=-1;
}else{
//Remove any clash constraints 
int pos1=cons[cur_branch][0].first-2;
int pos2=cons[cur_branch][0].second;
this->mols[pos1+0]->group=unclashGroup;
this->mols[pos1+1]->group=unclashGroup;
this->mols[pos1+2]->group=unclashGroup;
//Set group clashGroup to include 
this->mols[pos2+0]->group=unclashGroup;
this->mols[pos2+1]->group=unclashGroup;
this->mols[pos2+2]->group=unclashGroup;
//Take a step back in the tree
cur_branch=cur_branch-1;
c=-1;
//Remove the branch
path_tree.pop_back();
path_pos.pop_back();
}
}


}
}
if(c==-1){printf("Failed to find perfect structure\n"); success=0;}else{
success=1;
Pdb * pdb = new Pdb(this,'j');
pdb->write(outfile);
pdb->clear();
delete pdb;
}
return success;
}



int Polymere::genPairedCons(int clashGroup,std::vector<std::pair<int,int> > * cons, const std::vector< Fragment *>  frags,Fragclass * dat,std::vector< Fragatoms *> atoms,std::vector<Fragall *> batoms,std::vector<std::pair<int,int> > * branch)
{
//It if works it returns 1, if it fails it returns 0
// 1) Get all singles two fragment positions 5' to 3'
// 2) Identify all stacked configurations for position 1 that don't 
//    causes in the initial clash group
// 3) Identify all stacked configurations for position 2 that don't 
//    clash with itself plus the initial clashGroup
// 4) Loop through all possible combinations until pairs that are stacked
//    are identified. Return all pairs that are identified 
Polymere * pol1 = new Polymere();
Polymere * pol2 = new Polymere();
Polymere * outpol = new Polymere();
pol1->copy(this);
std::vector<int> pos1Fids,pos2Fids;
int pos1=cons->at(0).first-2; //pos1 is actually in the 5'-3' direction
int pos2=cons->at(0).second;
pol1->chooseAllSingle(pos1,dat,&pos1Fids);
pol1->chooseAllSingle(pos2,dat,&pos2Fids);
//Test singles from pos1 for clahses in clashGroup and stacking
//Determine the number of steps each fragment will have to move
int step1=0; //Step is the number of NT's in a given strain that are considered in the current movement group
int cmagb=cons->at(0).first;
for(int c=0; c<cons->size(); c++){ 
if(cons->at(c).first > cmagb-3){
step1++;
}
}

int step2=0;
cmagb=cons->at(0).second;
for(int c=0; c<cons->size(); c++){ 
if(cons->at(c).second < cmagb+3){
step2++;
}
}
step1=step1-1;
step2=step2-1;
printf("Step1: %d Step2: %d\n",step1,step2);

std::vector<int> stackprobs,stackFrags1,stackFrags2;
//Select fragments with the appropriate step-sizes
int accept;
float cutoff=0;
if(step1 != step2){
//anything goesi
cutoff=4;
}else{cutoff=.35;} 
printf("cutoff: %8.3f\n",cutoff);


for(int j1=0; j1<pos1Fids.size(); j1++){
int j = rand() % pos1Fids.size();
accept=0;
pol2->copy(this);
int pos1=cons->at(0).first-2;
//Set group clashGroup to include
for(int c=0; c<cons->size(); c++){ 
if(cons->at(c).first > cons->at(0).first-3){
pol2->mols[cons->at(c).first]->group=clashGroup;
}
}

if(batoms[pos1Fids[j]]->check()==1){
pol2->changeFrag2(pos1,pos1Fids[j],frags);
pol2->coordSysNt(pos1,&atoms[pos1Fids[j]]->atoms);
pol2->updateFull();
pol2->coordSysAT(pos1,&batoms[pos1Fids[j]]->atoms);
//printf("backbone:clashes %d:%d\n",pol2->groupCheckBackBone(clashGroup),pol2->groupCheckClashes(clashGroup));
if(pol2->groupCheckBackBone(clashGroup) < 1 && pol2->groupCheckClashes(clashGroup) < 1){
//Check for axis alignment and total z-distance 
//minimum angular distance = 10 degrees
float maxang=0;
for(int x=0; x<3; x++){
	for(int y=x+1; y<3; y++){
	float cang=pol2->zaligned((pos1+x),(pos1+y));
if(abs(cang)>maxang){maxang=abs(cang);}	
	}
}
if(maxang < cutoff){
//Appropriately Stacked
//Determine if its z-rise is appropriate for the step size required
/*
float dist_c1 = pol2->c1distf(pos1,pos1+1)+pol2->c1distf(pos1+1,pos1+2);
if(step1 > step2 && dist_c1 > 11){
accept=1;
}else{
if(step1 < step2 && dist_c1 < 11){
accept =1;
}else{
if(step1==step2){
accept=1;
}
}
} 
*/
accept=1;
if(accept==1){
stackFrags1.push_back(pos1Fids[j]);
}

}
}
}
}

for(int j2=0; j2<pos2Fids.size(); j2++){
int j = rand() % pos2Fids.size();
accept=0;
pol2->copy(this);
int pos2=cons->at(0).second;
//Set group clashGroup to include
for(int c=0; c<cons->size(); c++){
//Only use the first 3 in the clash group
if(cons->at(c).second < cons->at(0).second+3){
pol2->mols[cons->at(c).second]->group=clashGroup;}
}

if(batoms[pos2Fids[j]]->check()==1){
pol2->changeFrag2(pos2,pos2Fids[j],frags);
pol2->coordSysNt(pos2,&atoms[pos2Fids[j]]->atoms);
pol2->updateFull();
pol2->coordSysAT(pos2,&batoms[pos2Fids[j]]->atoms);
//printf("backbone:clashes %d:%d\n",pol2->groupCheckBackBone(clashGroup),pol2->groupCheckClashes(clashGroup));
if(pol2->groupCheckBackBone(clashGroup) < 1 && pol2->groupCheckClashes(clashGroup) < 1){
//Check for axis alignment and total z-distance 
//minimum angular distance = 10 degrees
float maxang=0;
for(int x=0; x<3; x++){
	for(int y=x+1; y<3; y++){
	float cang=pol2->zaligned((pos2+x),(pos2+y));
if(abs(cang)>maxang){maxang=abs(cang);}	
	}
}
if(maxang < cutoff){
//Appropriately Stacked
//Determine if its z-rise is appropriate for the step size required
/*
float dist_c1 = pol2->c1distf(pos2,pos2+1)+pol2->c1distf(pos2+1,pos2+2);
if(step1 > step2 && dist_c1 > 11){
accept=1;
}else{
if(step1 < step2 && dist_c1 < 11){
accept =1;
}else{
if(step1==step2){
accept=1;
}
}
}
*/ 
accept=1;
if(accept==1){
stackFrags2.push_back(pos2Fids[j]);
}

}


}
}
}


printf("Total Stacked Frags: %d - %d\n",stackFrags1.size(), stackFrags2.size());
printf("Testing all combinations for aligned pairs..\n");
int bclash=0;
int aclash=0;
int paired=0;
//Set up constraint pairs
std::vector< std::pair<int,int> > cons3; //If there are pairs 
std::vector< std::pair<int,int> > cons1;
std::vector< std::pair<int,int> > cons2;
int con_type=1;
for(int i=0; i<cons->size(); i++){
if(cons->at(i).second < (pos2+3)  && cons->at(i).first > pos1-1){
cons3.push_back(std::pair<int,int>(cons->at(i).first,cons->at(i).second));
}
printf("Reduced Cons: %d-%d\n",cons->at(i).first,cons->at(i).second);
}

if(cons3[0].first==cons3[1].first || cons3[0].second==cons3[1].second){
con_type=2;
cons1.push_back(cons3[0]);
cons1.push_back(cons3[2]);
cons2.push_back(cons3[1]);
cons2.push_back(cons3[2]);
printf("Using semi-cons: \n");
printf("c1: %d -- %d\n",cons1[0].first,cons1[0].second);
printf("c1: %d -- %d\n",cons1[1].first,cons1[1].second);
printf("c2: %d -- %d\n",cons2[0].first,cons2[0].second);
printf("c2: %d -- %d\n",cons2[1].first,cons2[1].second);
}else{
if(cons3[1].first==cons3[2].first || cons3[1].second==cons3[2].second){
con_type=2;
cons1.push_back(cons3[0]);
cons1.push_back(cons3[1]);

cons2.push_back(cons3[0]);
cons2.push_back(cons3[2]);
printf("Using semi-cons: \n");
printf("c1: %d -- %d\n",cons1[0].first,cons1[0].second);
printf("c1: %d -- %d\n",cons1[1].first,cons1[1].second);
printf("c2: %d -- %d\n",cons2[0].first,cons2[0].second);
printf("c2: %d -- %d\n",cons2[1].first,cons2[1].second);
}else{
//Just use cons3
printf("Using all-cons: \n");
printf("c3: %d -- %d\n",cons3[0].first,cons3[0].second);
printf("c3: %d -- %d\n",cons3[1].first,cons3[1].second);
printf("c3: %d -- %d\n",cons3[2].first,cons3[2].second);
}
}



int min_score = 9999;
outpol->copy(this);
printf("Number of molecules: %d\n",this->num_mols);
for(int w=0; w<stackFrags1.size(); w++){
for(int x=0; x<stackFrags2.size(); x++){
pol2->copy(this);
pol2->mols[pos1+0]->group=clashGroup;
pol2->mols[pos1+1]->group=clashGroup;
pol2->mols[pos1+2]->group=clashGroup;
//Set group clashGroup to include 
pol2->mols[pos2+0]->group=clashGroup;
pol2->mols[pos2+1]->group=clashGroup;
pol2->mols[pos2+2]->group=clashGroup;

pol2->changeFrag2(pos1,stackFrags1[w],frags);
pol2->changeFrag2(pos2,stackFrags2[x],frags);
pol2->coordSysNt(pos1,&atoms[stackFrags1[w]]->atoms);
pol2->coordSysNt(pos2,&atoms[stackFrags2[x]]->atoms);
pol2->updateFull();
pol2->coordSysAT(pos1,&batoms[stackFrags1[w]]->atoms);
pol2->coordSysAT(pos2,&batoms[stackFrags2[x]]->atoms);
//if(pol2->groupCheckBackBone(clashGroup) > 0){bclash++;}
//if(pol2->groupCheckClashes(clashGroup) > 0){aclash++;}
if(pol2->groupCheckBackBone(clashGroup) < 1 && pol2->groupCheckClashes(clashGroup) < 1 && ((con_type==1 && pol2->hardConstraints(&cons3)<1) || (con_type==2 && (pol2->hardConstraints(&cons1)<1 || pol2->hardConstraints(&cons2)<1)) ) ){
//Check for pairs that keep the frags within the constraints
outpol->copy(pol2);
printf("Perfect Pair\n");
min_score=0;
std::pair<int,int> fragpair(stackFrags1[w],stackFrags2[x]);
branch->push_back(fragpair);
if(branch->size() > 10){
break;
}
}else{
if(pol2->groupCheckBackBone(clashGroup) < 1 && pol2->groupCheckClashes(clashGroup) < 1){
//Check for pairs that keep the frags within the constraints
if(min_score > pol2->hardConstraints(&cons3)){
outpol->copy(pol2);
min_score=pol2->hardConstraints(&cons3);
printf("min_cons: %d\n",min_score);
}
}
}
}
if(branch->size()>10){break;}
}
int goforward;
if(branch->size() > 0){
this->copy(outpol);
goforward=1;
}else{
goforward=0;
}
pol1->clearMols();
pol2->clearMols();
outpol->clearMols();
//printf("B-Clashes: %d A-Clashes: %d\n",bclash,aclash);
delete pol1;
delete pol2;
delete outpol;
return goforward;
}

void Polymere::genPairedStacks(int pos1,  int pos2,int clashGroup, const std::vector< Fragment *>  frags,Fragclass * dat,std::vector< Fragatoms *> atoms,std::vector<Fragall *> batoms)
{
// 1) Get all singles two fragment positions 5' to 3'
// 2) Identify all stacked configurations for position 1 that don't 
//    causes in the initial clash group
// 3) Identify all stacked configurations for position 2 that don't 
//    clash with itself plus the initial clashGroup
// 4) Loop through all possible combinations until pairs that are stacked
//    are identified. Return all pairs that are identified 
Polymere * pol1 = new Polymere();
Polymere * pol2 = new Polymere();
Polymere * outpol = new Polymere();
pol1->copy(this);
std::vector<int> pos1Fids,pos2Fids;
pol1->chooseAllSingle(pos1,dat,&pos1Fids);
pol1->chooseAllSingle(pos2,dat,&pos2Fids);
//Test singles from pos1 for clahses in clashGroup and stacking
std::vector<int> stackprobs,stackFrags1,stackFrags2;
for(int j=0; j<pos1Fids.size(); j++){
pol2->copy(this);
//Set group clashGroup to include 
pol2->mols[pos1+0]->group=clashGroup;
pol2->mols[pos1+1]->group=clashGroup;
pol2->mols[pos1+2]->group=clashGroup;
pol2->mols[pos1+3]->group=clashGroup;

if(batoms[pos1Fids[j]]->check()==1){
pol2->changeFrag2(pos1,pos1Fids[j],frags);
pol2->coordSysNt(pos1,&atoms[pos1Fids[j]]->atoms);
pol2->updateFull();
pol2->coordSysAT(pos1,&batoms[pos1Fids[j]]->atoms);
//printf("backbone:clashes %d:%d\n",pol2->groupCheckBackBone(clashGroup),pol2->groupCheckClashes(clashGroup));
if(pol2->groupCheckBackBone(clashGroup) < 1 && pol2->groupCheckClashes(clashGroup) < 1){
stackprobs.clear();
float iscore=pol2->adjStacked(&stackprobs);
//Check to see if the current position becomes stacked
//printf("stackprobs1: %d %d %d %d\n",stackprobs[pos1-1],stackprobs[pos1],stackprobs[pos1+1],stackprobs[pos1+2]);
if((stackprobs[pos1]+stackprobs[pos1+1]+stackprobs[pos1+2])>0){
//save the fragment id
stackFrags1.push_back(pos1Fids[j]);
}
}
}


}

for(int j=0; j<pos2Fids.size(); j++){
pol2->copy(this);
//Set group clashGroup to include 
pol2->mols[pos2+0]->group=clashGroup;
pol2->mols[pos2+1]->group=clashGroup;
pol2->mols[pos2+2]->group=clashGroup;
pol2->mols[pos2+3]->group=clashGroup;

if(batoms[pos2Fids[j]]->check()==1){
pol2->changeFrag2(pos2,pos2Fids[j],frags);
pol2->coordSysNt(pos2,&atoms[pos2Fids[j]]->atoms);
pol2->updateFull();
pol2->coordSysAT(pos2,&batoms[pos2Fids[j]]->atoms);
if(pol2->groupCheckBackBone(clashGroup) < 1 && pol2->groupCheckClashes(clashGroup) < 1){
stackprobs.clear();
float iscore=pol2->adjStacked(&stackprobs);
//printf("stackprobs2: %d %d %d %d\n",stackprobs[pos2-1],stackprobs[pos2],stackprobs[pos2+1],stackprobs[pos2+2]);
//Check to see if the current position becomes stacked
if((stackprobs[pos2]+stackprobs[pos2+1]+stackprobs[pos2-1])>0){
//save the fragment id
stackFrags2.push_back(pos2Fids[j]);
}
}
}
}


printf("Total Stacked Frags: %d - %d\n",stackFrags1.size(), stackFrags2.size());
printf("Testing all combinations for aligned pairs..\n");
int paired=0;
//Set up constraint pairs
std::vector< std::pair<int,int> > cons;
cons.push_back(std::pair<int,int>(pos1,pos2+2));
cons.push_back(std::pair<int,int>(pos1+1,pos2+1));
cons.push_back(std::pair<int,int>(pos1+2,pos2+0));
int min_score = 9999;
outpol->copy(this);
for(int w=0; w<stackFrags1.size(); w++){
for(int x=0; x<stackFrags2.size(); x++){
pol2->copy(this);
pol2->mols[pos1+0]->group=clashGroup;
pol2->mols[pos1+1]->group=clashGroup;
pol2->mols[pos1+2]->group=clashGroup;
pol2->mols[pos1+3]->group=clashGroup;
//Set group clashGroup to include 
pol2->mols[pos2+0]->group=clashGroup;
pol2->mols[pos2+1]->group=clashGroup;
pol2->mols[pos2+2]->group=clashGroup;
pol2->mols[pos2+3]->group=clashGroup;

pol2->changeFrag2(pos1,stackFrags1[w],frags);
pol2->changeFrag2(pos2,stackFrags2[x],frags);
pol2->coordSysNt(pos1,&atoms[stackFrags1[w]]->atoms);
pol2->coordSysNt(pos2,&atoms[stackFrags2[x]]->atoms);
pol2->updateFull();
pol2->coordSysAT(pos1,&batoms[stackFrags1[w]]->atoms);
pol2->coordSysAT(pos2,&batoms[stackFrags2[x]]->atoms);
if(pol2->groupCheckBackBone(clashGroup) < 1 && pol2->groupCheckClashes(clashGroup) < 1 && pol2->hardConstraints(&cons)<1){
//Check for pairs that keep the frags within the constraints
outpol->copy(pol2);
printf("Perfect Pair\n");
Pdb * pdb = new Pdb(pol2,'j');
pdb->write("helix0.pdb");
pdb->clear();
delete pdb;
min_score=0;
break;
}else{
if(pol2->groupCheckBackBone(clashGroup) < 1 && pol2->groupCheckClashes(clashGroup) < 1){
//Check for pairs that keep the frags within the constraints
if(min_score > pol2->hardConstraints(&cons)){
outpol->copy(pol2);
min_score=pol2->hardConstraints(&cons);
printf("min_cons: %d\n",min_score);
}
}
}
}
if(min_score==0){break;}
}

this->copy(outpol);
pol1->clearMols();
pol2->clearMols();
outpol->clearMols();
delete pol1;
delete pol2;
delete outpol;
}

void Polymere::genStacked(const std::vector< Fragment *>  frags,Fragclass * dat,std::vector< Fragatoms *> atoms,std::vector<Fragall *> batoms)
{
//WARNING this removes all movement probability information contained in the polymere
//first = mid.point start-3, second = midpoint start+1
//For each movement group identify the best pair of fragments that
//by score,clashes, and backbone connections
//Create an initial group that contains all elements between the frist 
//two positions in the movement group

for(int g=2; g<(num_mols-3); g++){mols[g]->prob=0;}
printf("Setting up groups and probs: %d,%d\n",this->groupCheckBackBone(5),this->groupCheckClashes(5));
//Create a temporary polymere to hold the tested positions
Polymere * pol2=new Polymere();
Polymere * pol1=new Polymere();
Polymere * outpol=new Polymere();
//Use single moves between first pair of nt's to identify fragments
//that align the first pair nt's
pol1->copy(this);//reset the starting position
outpol->copy(this);
for(int g=0; g<pol1->mols.size(); g++){pol1->mols[g]->group=0;}
int best_score=99999; //Best Possible Score is -1*(num_mols - 1)
std::vector<int> stackprobs;
std::vector<int> movepos;
for(int s=0; s<100; s++){
//printf("Attempt: %d\n",s);
//Reset group information
pol1->copy(this);//reset the starting position
for(int g=0; g<pol1->mols.size(); g++){pol1->mols[g]->group=0;} // Only consider the new fragment that is being added for clashes
pol1->mols[3]->group=5;
pol1->mols[4]->group=5;
pol1->mols[5]->group=5;

//Pick random positions from start loop
//Try to minimize score for starting NT's
//Get a list of positions that still need to be stacked
pol2->copy(pol1);
//For each unstacked position replace it with a stacked fragment
for(int rpos=3; rpos<(num_mols-3); rpos=rpos+2){
//printf("Moving %d\n",rpos);
//Move it until it stacks or 100 moves
//reset pol2 to pol1
pol2->copy(pol1);
int madeit=0;
int iter=0;
while(madeit==0){
if(iter>1000){
//printf("failed: %d\n",rpos);
madeit=2;
break;
//Reset and try again
}
int fid=pol2->chooseSingle(rpos,dat);
if(batoms[fid]->check()==1){
pol2->changeFrag2(rpos,fid,frags);
pol2->coordSysNt(rpos,&atoms[fid]->atoms);
pol2->updateFull();
pol2->coordSysAT(rpos,&batoms[fid]->atoms);
if(pol2->groupCheckBackBone(5) < 1 && pol2->groupCheckClashes(5) < 1){
stackprobs.clear();
float iscore=-1*pol2->adjStacked(&stackprobs);
//Check to see if the current position becomes stacked

if( (rpos > 3 && (stackprobs[rpos-1]+stackprobs[rpos]+stackprobs[rpos+1]+stackprobs[rpos+2])==4) || rpos ==3 && (stackprobs[rpos]+stackprobs[rpos+1]+stackprobs[rpos+2])==3){
pol2->mols[rpos]->prob=1;
pol2->mols[rpos+3]->group=5;
pol2->mols[rpos+4]->group=5;
//printf("madit: %d\n",rpos);
madeit=1;
pol1->copy(pol2);
}
}
}
iter++;
}//Loop madeit==0
if(madeit==2){break;}
}//Loop rpos
//Update stacking probs
stackprobs.clear();
float tot;
tot=pol1->adjStacked(&stackprobs);
int unstacked=0;
for(int m=4; m<(pol1->num_mols-5); m++){
if(stackprobs[m-1]==1 && stackprobs[m]==1 && stackprobs[m+1]==1 && stackprobs[m+2]==1){pol1->mols[m]->prob=1;}else{
//printf("%d - unstacked\n",m);
unstacked++;}
}
if(unstacked < best_score){
outpol->copy(pol1);
//printf("Total Unstacked: %d\n",unstacked);
}
if(unstacked==0){
outpol->copy(pol1);
break;
}
}
this->copy(outpol);
//Clean-up
pol1->clearMols();
pol2->clearMols();
outpol->clearMols();
delete pol1; 
delete pol2; 
delete outpol; 
}


int Polymere::findLoop(int start_pos,int type,std::vector< std::vector<int> > * frag_list,std::vector<Fragment *> frags, Fragclass * dat,std::vector<Fragatoms* > atoms, std::vector<Fragall *> batoms){
//loop type = 0 => start_pos is matched to start_pos+7
//loop type = 1 => start_pos is matched to start_pos+8
//return 1 if loop is found or 0 if attempts failed
int success=0;
Polymere * pol1 = new Polymere();
Polymere * pol2 = new Polymere();
pol1->copy(this);
std::vector< std::pair<int,int> > cons;
if(type==1){
cons.push_back(std::pair<int,int>(start_pos,start_pos+8));
}else{
if(type==0){
cons.push_back(std::pair<int,int>(start_pos,start_pos+7));
}else{
printf("bad loop type: %d",type);
printf("should be 0 or 1\n");
return success;
}
}

//Set clash detection to start - 2 to start+8
for(int m=0; m<pol1->mols.size(); m++){
if(m < (start_pos) || m > start_pos+7){
pol1->mols[m]->group=0;
}else{
pol1->mols[m]->group=5;
}
}
//Foreach random member of frag_list update the polymere
for(int i =0; i < frag_list->size(); i++){
pol2->copy(pol1);
int rand1 = rand() % frag_list->size();
pol2->changeFrag2(start_pos+0,frag_list->at(rand1)[0],frags);
pol2->changeFrag2((start_pos+2),frag_list->at(rand1)[1],frags);
pol2->changeFrag2((start_pos+4),frag_list->at(rand1)[2],frags);
pol2->changeFrag2((start_pos+6),frag_list->at(rand1)[3],frags);
pol2->coordSysNt((start_pos+0),&atoms[frag_list->at(rand1)[0]]->atoms);
pol2->coordSysNt((start_pos+2),&atoms[frag_list->at(rand1)[1]]->atoms);
pol2->coordSysNt((start_pos+4),&atoms[frag_list->at(rand1)[2]]->atoms);
pol2->coordSysNt((start_pos+6),&atoms[frag_list->at(rand1)[3]]->atoms);
pol2->updateFull();
pol2->coordSysAT((start_pos+0),&batoms[frag_list->at(rand1)[0]]->atoms);
pol2->coordSysAT((start_pos+2),&batoms[frag_list->at(rand1)[1]]->atoms);
pol2->coordSysAT((start_pos+4),&batoms[frag_list->at(rand1)[2]]->atoms);
pol2->coordSysAT((start_pos+6),&batoms[frag_list->at(rand1)[3]]->atoms);
if(pol2->groupCheckBackBone(5) < 1 && pol2->groupCheckClashes(5) < 1){
if(pol2->hardConstraints(&cons) < 1){
this->copy(pol2);
success=1;
break;
}
}
}
return success;
}

void Polymere::buildLoopO(int start_pos,int end_pos,int nmax,int iter_max,const std::vector< Fragment *>  frags,Fragclass * dat,std::vector< Fragatoms *> atoms,std::vector<Fragall *> batoms,string gfile,string loopseq)
{
FILE * sfile = fopen(gfile.c_str(),"a");
//Start pos is the first position in the loop
//	This is the position of the first base-pair in the helix + 1

//Workflow:
// 1) Find all pairs of fragments(pos1,pos1+2) that:
//	don't clash
// 	get cons within appropriate distance
// 2) From the possible pairs filter for those that:
//	preserve angle and stacking
// 3) If results from 2 is the empty set 
//	use no clashing pairs (pos1,pos1+2),pos1+4)
//	repeat filters
// 	If results are empty set return -1
//	Else return list for all fragment 
//Create a temporary polymere to hold the tested positions

Polymere * pol2=new Polymere();
Polymere * pol1=new Polymere();
Polymere * outpol=new Polymere();
std::vector< std::pair<int,int> > cons;
cons.push_back(std::pair<int,int>((start_pos),(start_pos+8)));
pol1->copy(this);
//Set group variable on pol1 
for(int m=0; m < mols.size(); m++){
if(m >= (start_pos-1) && m <=start_pos+4){pol1->mols[m]->group=5;}else{pol1->mols[m]->group=0;} 
}
printf("Getting Frags\n");
std::vector<int> pos1Fids,pos2Fids,pos3Fids;
pol1->chooseAllSingle(start_pos,dat,&pos1Fids);
pol1->chooseAllSingle((start_pos+2),dat,&pos2Fids);
int totpairs=pos1Fids.size()*pos2Fids.size();
printf("Total Possible Fragment Pairs: %d\n",totpairs);
printf("Initial clashes %d %d\n",pol1->groupCheckBackBone(5), pol1->groupCheckClashes(5));
//Filter for pairs that get the distance correct without clashes
int tot=0;
int total_out=0;
std::vector< std::pair<int,int> > frag_set1;
while(tot < totpairs && total_out < 5000){
pol2->copy(pol1);
//Pick random pair
int rand1=rand() % pos1Fids.size();
int rand2=rand() % pos2Fids.size();
//Change Frag and check for clashes and distance
if(batoms[pos1Fids[rand1]]->check()==1 && batoms[pos2Fids[rand2]]->check()==1){
pol2->changeFrag2(start_pos,pos1Fids[rand1],frags);
pol2->changeFrag2((start_pos+2),pos2Fids[rand2],frags);
pol2->coordSysNt(start_pos,&atoms[pos1Fids[rand1]]->atoms);
pol2->coordSysNt((start_pos+2),&atoms[pos2Fids[rand2]]->atoms);
pol2->updateFull();
pol2->coordSysAT(start_pos,&batoms[pos1Fids[rand1]]->atoms);
pol2->coordSysAT((start_pos+2),&batoms[pos2Fids[rand2]]->atoms);
if(pol2->groupCheckBackBone(5) < 1 && pol2->groupCheckClashes(5) < 1){
//Check distance
float c1d=pol2->c1distf((start_pos),(start_pos+4));
if(c1d < 22 && c1d > 6.5){
if(total_out % 100 == 0){printf("%d\n",total_out);}
frag_set1.push_back(std::pair<int,int>(pos1Fids[rand1],pos2Fids[rand2]));
//Check all constraints
total_out++;
//printf("\ndist: %8.3f\n",c1d);
}
}
}
tot++;
}
printf("Total Pairs of Moves: %d\n",frag_set1.size());
std::vector< std::vector<int> > frag_list;

for(int m=0; m < mols.size(); m++){
if(m >= (start_pos-1) && m <=start_pos+6){pol1->mols[m]->group=5;}else{pol1->mols[m]->group=0;} 
}
//For each pair of frags that didn't clash or break the backbone
//use a third position and look for frags that work
pol1->chooseAllSingle(start_pos+4,dat,&pos3Fids);
tot=0;
total_out=0;
totpairs=pos3Fids.size()*frag_set1.size();
while(tot < totpairs && total_out < 5000){
pol2->copy(pol1);
//Pick random pair
int rand1=rand() % frag_set1.size();
int rand2=rand() % pos3Fids.size();
//Change Frag and check for clashes and distance
if(batoms[pos3Fids[rand2]]->check()==1){
pol2->changeFrag2(start_pos,frag_set1[rand1].first,frags);
pol2->changeFrag2((start_pos+2),frag_set1[rand1].second,frags);
pol2->changeFrag2((start_pos+4),pos3Fids[rand2],frags);
pol2->coordSysNt((start_pos+0),&atoms[frag_set1[rand1].first]->atoms);
pol2->coordSysNt((start_pos+2),&atoms[frag_set1[rand1].second]->atoms);
pol2->coordSysNt((start_pos+4),&atoms[pos3Fids[rand2]]->atoms);
pol2->updateFull();
pol2->coordSysAT(start_pos,&batoms[frag_set1[rand1].first]->atoms);
pol2->coordSysAT((start_pos+2),&batoms[frag_set1[rand1].second]->atoms);
pol2->coordSysAT((start_pos+4),&batoms[pos3Fids[rand2]]->atoms);
if(pol2->groupCheckBackBone(5) < 1 && pol2->groupCheckClashes(5) < 1){
float c1d=pol2->c1distf((start_pos),(start_pos+6));
if(c1d < 20 && c1d > 6.5){
if(total_out % 100 == 0){printf("%d\n",total_out);}
std::vector<int> frags_o;
frags_o.push_back(frag_set1[rand1].first);
frags_o.push_back(frag_set1[rand1].second);
frags_o.push_back(pos3Fids[rand2]);
frag_list.push_back(frags_o);
total_out++;
}
}
}
tot++;
}

printf("Total Pairs of Moves: %d\n",frag_list.size());
for(int m=0; m < mols.size(); m++){
if(m >= (start_pos-1) && m <=start_pos+8){pol1->mols[m]->group=5;}else{pol1->mols[m]->group=0;} 
}
//For each pair of frags that didn't clash or break the backbone
//use a third position and look for frags that work
pos3Fids.clear();
pol1->chooseAllSingle(start_pos+6,dat,&pos3Fids);
tot=0;
total_out=0;
totpairs=pos3Fids.size()*frag_list.size();
while(tot < totpairs && total_out < 5000){
//if(tot%1000==0){printf("tot_both: %d, %d\n",tot,total_out);}
pol2->copy(pol1);
//Pick random pair
int rand1=rand() % frag_list.size();
int rand2=rand() % pos3Fids.size();
//Change Frag and check for clashes and distance
if(batoms[pos3Fids[rand2]]->check()==1){
pol2->changeFrag2(start_pos+0,frag_list[rand1][0],frags);
pol2->changeFrag2((start_pos+2),frag_list[rand1][1],frags);
pol2->changeFrag2((start_pos+4),frag_list[rand1][2],frags);
pol2->changeFrag2((start_pos+6),pos3Fids[rand2],frags);
pol2->coordSysNt((start_pos+0),&atoms[frag_list[rand1][0]]->atoms);
pol2->coordSysNt((start_pos+2),&atoms[frag_list[rand1][1]]->atoms);
pol2->coordSysNt((start_pos+4),&atoms[frag_list[rand1][2]]->atoms);
pol2->coordSysNt((start_pos+6),&atoms[pos3Fids[rand2]]->atoms);
pol2->updateFull();
pol2->coordSysAT(start_pos+0,&batoms[frag_list[rand1][0]]->atoms);
pol2->coordSysAT((start_pos+2),&batoms[frag_list[rand1][1]]->atoms);
pol2->coordSysAT((start_pos+4),&batoms[frag_list[rand1][2]]->atoms);
pol2->coordSysAT((start_pos+6),&batoms[pos3Fids[rand2]]->atoms);
if(pol2->groupCheckBackBone(5) < 1 && pol2->groupCheckClashes(5) < 1){
float c1d=pol2->c1distf((start_pos),(start_pos+8));
if(c1d < 14 && c1d > 5){
if(pol2->hardConstraints(&cons) < 1){
if(total_out % 100 == 0){printf("%d\n",total_out);}
int seq_size=num_mols-6;
Polymere *tmp=new Polymere();
tmp->copy(pol2);
tmp->mols.erase(tmp->mols.begin(),tmp->mols.begin()+9);
tmp->mols.erase(tmp->mols.begin()+(seq_size-11),tmp->mols.end());
tmp->full_atoms.erase(tmp->full_atoms.begin(),tmp->full_atoms.begin()+(81*9));
tmp->full_atoms.erase(tmp->full_atoms.begin()+(81*(seq_size-11)),tmp->full_atoms.end());
//Loop Database
fprintf(sfile,"%s %d %d %d %d\n",loopseq.c_str(),frag_list[rand1][0],frag_list[rand1][1],frag_list[rand1][2],pos3Fids[rand2]);
fflush( sfile );
Pdb * pdb =  new Pdb(tmp,'J');
string lpout=gfile;
lpout.append(".pdb");
pdb->write(lpout);
pdb->clear();
delete pdb;
for(int m=0; m<tmp->mols.size(); m++){
delete tmp->mols[m];
}
tmp->full_atoms.clear();
tmp->mols.clear();
delete tmp;
total_out++;
}
}
}
}
tot++;
}


fflush( sfile );
fclose(sfile);
pol1->clearMols();
pol2->clearMols();
outpol->clearMols();
delete pol1;
delete pol2;
delete outpol;
}

void Polymere::buildLoop(int start_pos,int end_pos,int nmax,int iter_max,const std::vector< Fragment *>  frags,Fragclass * dat,std::vector< Fragatoms *> atoms,std::vector<Fragall *> batoms,string gfile,string loopseq)
{
FILE * sfile = fopen(gfile.c_str(),"a");
//Start pos is the first position in the loop
//	This is the position of the first base-pair in the helix + 1

//Workflow:
// 1) Find all pairs of fragments(pos1,pos1+2) that:
//	don't clash
// 	get cons within appropriate distance
// 2) From the possible pairs filter for those that:
//	preserve angle and stacking
// 3) If results from 2 is the empty set 
//	use no clashing pairs (pos1,pos1+2),pos1+4)
//	repeat filters
// 	If results are empty set return -1
//	Else return list for all fragment 
//Create a temporary polymere to hold the tested positions

Polymere * pol2=new Polymere();
Polymere * pol1=new Polymere();
Polymere * outpol=new Polymere();
std::vector< std::pair<int,int> > cons;
cons.push_back(std::pair<int,int>((start_pos),(start_pos+7)));
pol1->copy(this);
//Set group variable on pol1 
for(int m=0; m < mols.size(); m++){
if(m >= (start_pos-1) && m <=start_pos+4){pol1->mols[m]->group=5;}else{pol1->mols[m]->group=0;} 
}
printf("Getting Frags\n");
std::vector<int> pos1Fids,pos2Fids,pos3Fids;
pol1->chooseAllSingle(start_pos,dat,&pos1Fids);
pol1->chooseAllSingle((start_pos+2),dat,&pos2Fids);
int totpairs=pos1Fids.size()*pos2Fids.size();
printf("Total Possible Fragment Pairs: %d\n",totpairs);
printf("Initial clashes %d %d\n",pol1->groupCheckBackBone(5), pol1->groupCheckClashes(5));
//Filter for pairs that get the distance correct without clashes
int tot=0;
int total_out=0;
std::vector< std::pair<int,int> > frag_set1;
while(tot < totpairs && total_out < 5000){
pol2->copy(pol1);
//Pick random pair
int rand1=rand() % pos1Fids.size();
int rand2=rand() % pos2Fids.size();
//Change Frag and check for clashes and distance
if(batoms[pos1Fids[rand1]]->check()==1 && batoms[pos2Fids[rand2]]->check()==1){
pol2->changeFrag2(start_pos,pos1Fids[rand1],frags);
pol2->changeFrag2((start_pos+2),pos2Fids[rand2],frags);
pol2->coordSysNt(start_pos,&atoms[pos1Fids[rand1]]->atoms);
pol2->coordSysNt((start_pos+2),&atoms[pos2Fids[rand2]]->atoms);
pol2->updateFull();
pol2->coordSysAT(start_pos,&batoms[pos1Fids[rand1]]->atoms);
pol2->coordSysAT((start_pos+2),&batoms[pos2Fids[rand2]]->atoms);
if(pol2->groupCheckBackBone(5) < 1 && pol2->groupCheckClashes(5) < 1){
//Check distance
float c1d=pol2->c1distf((start_pos),(start_pos+4));
if(c1d < 22 && c1d > 6.5){
if(total_out % 100 == 0){printf("%d\n",total_out);}
frag_set1.push_back(std::pair<int,int>(pos1Fids[rand1],pos2Fids[rand2]));
//Check all constraints
total_out++;
//printf("\ndist: %8.3f\n",c1d);
}
}
}
tot++;
}
printf("Total Pairs of Moves: %d\n",frag_set1.size());
std::vector< std::vector<int> > frag_list;

for(int m=0; m < mols.size(); m++){
if(m >= (start_pos-1) && m <=start_pos+6){pol1->mols[m]->group=5;}else{pol1->mols[m]->group=0;} 
}
//For each pair of frags that didn't clash or break the backbone
//use a third position and look for frags that work
pol1->chooseAllSingle(start_pos+4,dat,&pos3Fids);
tot=0;
total_out=0;
totpairs=pos3Fids.size()*frag_set1.size();
while(tot < totpairs && total_out < 5000){
pol2->copy(pol1);
//Pick random pair
int rand1=rand() % frag_set1.size();
int rand2=rand() % pos3Fids.size();
//Change Frag and check for clashes and distance
if(batoms[pos3Fids[rand2]]->check()==1){
pol2->changeFrag2(start_pos,frag_set1[rand1].first,frags);
pol2->changeFrag2((start_pos+2),frag_set1[rand1].second,frags);
pol2->changeFrag2((start_pos+4),pos3Fids[rand2],frags);
pol2->coordSysNt((start_pos+0),&atoms[frag_set1[rand1].first]->atoms);
pol2->coordSysNt((start_pos+2),&atoms[frag_set1[rand1].second]->atoms);
pol2->coordSysNt((start_pos+4),&atoms[pos3Fids[rand2]]->atoms);
pol2->updateFull();
pol2->coordSysAT(start_pos,&batoms[frag_set1[rand1].first]->atoms);
pol2->coordSysAT((start_pos+2),&batoms[frag_set1[rand1].second]->atoms);
pol2->coordSysAT((start_pos+4),&batoms[pos3Fids[rand2]]->atoms);
if(pol2->groupCheckBackBone(5) < 1 && pol2->groupCheckClashes(5) < 1){
float c1d=pol2->c1distf((start_pos),(start_pos+6));
if(c1d < 20 && c1d > 6.5){
if(total_out % 100 == 0){printf("%d\n",total_out);}
std::vector<int> frags_o;
frags_o.push_back(frag_set1[rand1].first);
frags_o.push_back(frag_set1[rand1].second);
frags_o.push_back(pos3Fids[rand2]);
frag_list.push_back(frags_o);
total_out++;
}
}
}
tot++;
}

printf("Total Pairs of Moves: %d\n",frag_list.size());
for(int m=0; m < mols.size(); m++){
if(m >= (start_pos-1) && m <=start_pos+7){pol1->mols[m]->group=5;}else{pol1->mols[m]->group=0;} 
}
//For each pair of frags that didn't clash or break the backbone
//use a third position and look for frags that work
pos3Fids.clear();
pol1->chooseAllSingle(start_pos+6,dat,&pos3Fids);
tot=0;
total_out=0;
totpairs=pos3Fids.size()*frag_list.size();
while(tot < totpairs && total_out < 5000){
//if(tot%1000==0){printf("tot_both: %d, %d\n",tot,total_out);}
pol2->copy(pol1);
//Pick random pair
int rand1=rand() % frag_list.size();
int rand2=rand() % pos3Fids.size();
//Change Frag and check for clashes and distance
if(batoms[pos3Fids[rand2]]->check()==1){
pol2->changeFrag2(start_pos+0,frag_list[rand1][0],frags);
pol2->changeFrag2((start_pos+2),frag_list[rand1][1],frags);
pol2->changeFrag2((start_pos+4),frag_list[rand1][2],frags);
pol2->changeFrag2((start_pos+6),pos3Fids[rand2],frags);
pol2->coordSysNt((start_pos+0),&atoms[frag_list[rand1][0]]->atoms);
pol2->coordSysNt((start_pos+2),&atoms[frag_list[rand1][1]]->atoms);
pol2->coordSysNt((start_pos+4),&atoms[frag_list[rand1][2]]->atoms);
pol2->coordSysNt((start_pos+6),&atoms[pos3Fids[rand2]]->atoms);
pol2->updateFull();
pol2->coordSysAT(start_pos+0,&batoms[frag_list[rand1][0]]->atoms);
pol2->coordSysAT((start_pos+2),&batoms[frag_list[rand1][1]]->atoms);
pol2->coordSysAT((start_pos+4),&batoms[frag_list[rand1][2]]->atoms);
pol2->coordSysAT((start_pos+6),&batoms[pos3Fids[rand2]]->atoms);
if(pol2->groupCheckBackBone(5) < 1 && pol2->groupCheckClashes(5) < 1){
float c1d=pol2->c1distf((start_pos),(start_pos+7));
if(c1d < 14 && c1d > 5){
if(pol2->hardConstraints(&cons) < 1){
if(total_out % 100 == 0){printf("%d\n",total_out);}
int seq_size=num_mols-6;
Polymere *tmp=new Polymere();
tmp->copy(pol2);
tmp->mols.erase(tmp->mols.begin(),tmp->mols.begin()+9);
tmp->mols.erase(tmp->mols.begin()+(seq_size-12),tmp->mols.end());
tmp->full_atoms.erase(tmp->full_atoms.begin(),tmp->full_atoms.begin()+(81*9));
tmp->full_atoms.erase(tmp->full_atoms.begin()+(81*(seq_size-12)),tmp->full_atoms.end());
//Loop Database
fprintf(sfile,"%s %d %d %d %d\n",loopseq.c_str(),frag_list[rand1][0],frag_list[rand1][1],frag_list[rand1][2],pos3Fids[rand2]);
fflush( sfile );
string lpout=gfile;
lpout.append(".pdb");
Pdb * pdb =  new Pdb(tmp,'J');
pdb->write(lpout);
pdb->clear();
delete pdb;
for(int m=0; m<tmp->mols.size(); m++){
delete tmp->mols[m];
}
tmp->full_atoms.clear();
tmp->mols.clear();
delete tmp;
total_out++;
}
}
}
}
tot++;
}


fflush( sfile );
fclose(sfile);
pol1->clearMols();
pol2->clearMols();
outpol->clearMols();
delete pol1;
delete pol2;
delete outpol;
}


void Polymere::midPointBuild(std::vector<std::pair<int,int> > movegroups,std::vector< std::pair<float,float> > cons,const std::vector< Fragment *>  frags,Fragclass * dat,std::vector< Fragatoms *> atoms,std::vector<Fragall *> batoms,std::vector<double> * dparams)
{
//WARNING this removes all group information contained in the polymere
//first = mid.point start-3, second = midpoint start+1
//For each movement group identify the best pair of fragments that
//by score,clashes, and backbone connections

//Create an initial group that contains all elements between the frist 
//two positions in the movement group
for(int g=0; g<mols.size(); g++){
mols[g]->group=0;
} 
int current_score;
double L = 8;
//Create a temporary polymere to hold the tested positions
Polymere * pol2=new Polymere();
Polymere * pol1=new Polymere();
Polymere * outpol=new Polymere();

//Use single moves between first pair of nt's to identify fragments
//that align the first pair nt's

int pos1=0;
int pos2=0;
std::vector<int> startloop;
pos1=movegroups[0].first+1;
if(pos1 > (movegroups[0].second)){
printf("Starting positions are to close\n"); 
abort();
}
while(pos1 < (movegroups[0].second+1)){
startloop.push_back(pos1);
pos1++;
}
//Set group for clash checking
for(int g=0; g<num_mols; g++){mols[g]->group=0;}
for(int g=movegroups[0].first-2; g< movegroups[0].second+5; g++){mols[g]->group=5;}
std::vector< std::pair<int,int> > constart;
constart.push_back(std::pair<int,int>((movegroups[0].first),(movegroups[0].second+2)));
constart.push_back(std::pair<int,int>((movegroups[0].first-1),(movegroups[0].second+3)));
score=99999;
outpol->copy(this);
pol1->copy(this);
float best_score=9999;
float start_temp=2;
for(int s=0; s<100; s++){
printf("\nSimulation: %d Current Best: %8.3f: Current Temp: %8.3f\n",s,pol1->score,start_temp);
pol1->copy(this);//reset the starting position
for(int n=0; n<20000; n++){
//Pick random positions from start loop
//Try to minimize score for starting NT's
pol2->copy(pol1);
//pick a random position from startloop
int rpos=rand() % startloop.size();
int fid=chooseSingle(startloop[rpos],dat);
if(batoms[fid]->check()==1){
pol2->changeFrag2(startloop[rpos],fid,frags);
pol2->coordSysNt(startloop[rpos],&atoms[fid]->atoms);
pol2->updateFull();
pol2->coordSysAT(startloop[rpos],&batoms[fid]->atoms);
if(pol2->groupCheckBackBone(5) < 1 && pol2->groupCheckClashes(5) < 1){
//float iscore=-1*(pol2->constraints(&constart,5000,start_temp));
float iscore=(pol2->hardConstraints(&constart));
if(iscore <= pol1->score){
printf("score: %8.3f\n",iscore);
pol1->score = iscore;
pol2->score=iscore;
pol1->copy(pol2);
}
}
if(pol1->score < best_score){
outpol->copy(pol1);
best_score = pol1->score;
}
}
}
int ocons=outpol->hardConstraints(&constart);
if(ocons==0){
break;
}
printf("Outpol contraints: %d\n\n",ocons);
}
//Use best of quick 10
this->copy(outpol);
float prevscore=0;
std::vector<Matches> mat; 
std::vector< std::pair<int,int> > conshard;
int tot_groups=0;
conshard.push_back(std::pair<int,int>((movegroups[0].first+3),(movegroups[0].second-1)));
conshard.push_back(std::pair<int,int>((movegroups[0].first+2),(movegroups[0].second+0)));
conshard.push_back(std::pair<int,int>((movegroups[0].first+1),(movegroups[0].second+1)));
conshard.push_back(std::pair<int,int>((movegroups[0].first+0),(movegroups[0].second+2)));
for(int pmove=0; pmove < movegroups.size(); pmove++){
printf("Adding Constraints: \n");
conshard.push_back(std::pair<int,int>((movegroups[pmove].first-1),(movegroups[pmove].second+3)));
printf("%d-%d\n",(movegroups[pmove].first+0),(movegroups[pmove].second+2));
outpol->copy(this);
tot_groups=0;
pos1=movegroups[pmove].first;
pos2=movegroups[pmove].second;
//Update the movement groups
for(int g=(pos1+1); g<(pos2+2); g++){
mols[g]->group=5;
tot_groups++;
}
//The first fragment is moving in the 5'-> 3' direction, the 
//second fragment is moveing in the 3'->5' direction
mat=this->fixedMoves(pos2,pos1,dat,L);
current_score=this->hardConstraints(&conshard);
printf("%d Total Paired Moves:%d Current Score: %d\n",pmove,mat.size(),current_score);
//Limit the total sampled to 10,000
float prob = (float)100000/mat.size();
for(int m=0; m<mat.size(); m++){
if(rand()/(float)RAND_MAX < prob){
if(batoms[mat[m].frag2]->check()==1){
pol2->copy(this);
pol2->changeFrag2(pos1,mat[m].frag2,frags);
pol2->changeFrag2(pos2,mat[m].frag1,frags);
pol2->coordSysNt(pos1,&atoms[mat[m].frag2]->atoms);
pol2->coordSysNt(pos2,&atoms[mat[m].frag1]->atoms);
pol2->updateFull();
pol2->coordSysAT(pos1,&batoms[mat[m].frag2]->atoms);
pol2->coordSysAT(pos2,&batoms[mat[m].frag1]->atoms);
//Write to PDB file
if(pol2->groupCheckBackBone(5) < 1 && pol2->groupCheckClashes(5) < 1){
int iscore=pol2->hardConstraints(&conshard);
if(iscore < current_score){
pol2->score=(float)iscore;
//Copy pol2 into this
printf("Overall Score: %d\n",iscore);
outpol->copy(pol2);
current_score = iscore;
if(iscore==0){break;
this->copy(pol2);
}
}
//Check to see if th changeFrag1 is the same as changeFrag2
}
}
}
}
this->copy(outpol);
}

}






void Polymere::simulate(int rmax, int nmax, int movegroup, const std::vector< Fragment *>  frags,Fragclass * dat,const Polymere * pol,std::vector< Fragatoms *> atoms,std::vector<Fragall *> batoms,std::vector< std::pair<int,int> > * cons){
using namespace std;
Polymere  * pol2 = new Polymere();     // Polymer used in the simulation. Everything is copied into and out of this class
std::vector< std::vector<int> > groups;// Movement groups of the polymere
int maxg=0;
for(int g=0; g<num_mols; g++){         // Determine the total number of movement groups. The user must define groups without any
int gc = pol->mols[g]->group;          // gaps. For example, 0 1 2 3 is ok, whereas 0 2 3 isn't
if(gc>maxg){
maxg=gc;
}
}
maxg++;
cout << "Max Groups: " << maxg << endl; 

for(int i=1; i < maxg; i++){         // Set-up the groups vector.  Stores vectors of positions for each group
std::vector<int> group;
for(int g=0; g<pol->num_mols; g++){
int gc = pol->mols[g]->group;
if(gc==i){
group.push_back(g);
}
}
groups.push_back(group);
}
int sig = groups[movegroup].size();   // Total number of positions in the movement group that will be used
for(int g=0; g<groups[movegroup].size(); g++){mols[groups[movegroup][g]]->prob=0;}  //Sets the (1-movement probability) for the each molecule to zero
pol2->copy(this);
this->score = (float)(pol2->atomClashes()+pol2->checkBackBone()); // Set the initial score
printf("Initial Clashes: %8.3f\n",score);
//rmax is from the -r option.  Sets the maximum number of randomizations to do before begining dynamic constraints
float cur_score = 0;
float temp = 70;
float ratio=0;
for(int iteration=0; iteration < rmax; iteration++){
pol2->copy(this);//Copy current conformation into pol2
int krand=(rand() % sig);//Pick a random position
float prob=rand()/(float)(RAND_MAX);
//float totalprob = pol2->mols[groups[movegroup][krand]]->prob*(groups[movegroup].size()/2)/((float)(iteration)); //probability that the move should be excepted.  This allows for equal sampling across the polymer
float totalprob=0;
if(totalprob < prob){
int fid=chooseSingle(groups[movegroup][krand],dat);
if(batoms[fid]->check()==1){
pol2->changeFrag2(groups[movegroup][krand],fid,frags);
pol2->coordSysNt(groups[movegroup][krand],&atoms[fid]->atoms);
pol2->updateFull();
pol2->coordSysAT(groups[movegroup][krand],&batoms[fid]->atoms);

//Minimize the number of clashes and bad backbone atoms
cur_score = (float)(pol2->atomClashes() + pol2->checkBackBone());
if(cur_score <= score){
prob=0; 
ratio=1;
}else{
ratio = exp(-1*(cur_score - score)/temp);
prob = rand()/(float)(RAND_MAX);
}
//printf("%8.3f - %8.3f: %8.3f, %8.3f, %8.3f,%d,%d\n",prob,ratio,cur_score,score,temp,pol2->checkBackBone(),pol2->atomClashes());
if(prob < ratio){
printf("    score: %8.3f\n",cur_score);
//Adjust the probability for the krand position
pol2->mols[groups[movegroup][krand]]->prob=(pol2->mols[groups[movegroup][krand]]->prob + 1);
pol2->score = cur_score;
this->copy(pol2);
if(cur_score == 0){
Pdb * pdb =  new Pdb(pol2,'J');
pdb->write("sim1.pdb");
pdb->clear();
delete pdb;
}
}
}
}
temp = temp*.96; // Slowing reduces the temperature 
}

printf("Clash Constraints %8.3f\n",this->score);
//Identifies conformations with maximal closed moves entropy that also abide by secondary structure constraints
//Vectors containing the indices of valid movement groups
std::vector<int> groupOrder1;
std::vector<int> groupOrder2;
std::vector<int> groupOrder0;

for(int i=1; i < sig-2; i++){
groupOrder1.push_back(i);
groupOrder2.push_back(i);
}
int nt=(sig-3)*(sig-3);
for(int i=0; i <nt; i++){groupOrder0.push_back(i);} 
int tot=0;
//Reading scoring function
int bout;
float ctemp=10;
float vcstr;
vcstr = this->constraints(cons,0,ctemp);
this->score=vcstr;
for(int iteration=0; iteration < nmax; iteration++){
bout=0;
if(score > 8.5){
printf("Constraint Level: %8.3f - %8.3f\n", ctemp,this->score);
ctemp=this->updateConstraints(cons);
printf("New Constraint Level: %8.3f\n",ctemp);
score=5;
}
////////////////////////////////////////////////////////////
// Move current polymere, score, add improvements to uniq //
////////////////////////////////////////////////////////////
pol2->copy(this);//Copy current conformation into pol2
//Create a list of random pairs
//Instead of continuing to pick random pairs
//We create the list of all possible pairs in a random order
//If the search exhausts all possiblities we stop the simulation 
std::vector<Matches> mat; 
int n = groupOrder1.size();
int k;
int temp;
while(n>1){
k=(rand() % n);
--n;
temp=groupOrder1[n];
groupOrder1[n]=groupOrder1[k];
groupOrder1[k]=temp;
}

n=groupOrder2.size();
while(n>1){
k=(rand() % n);
--n;
temp=groupOrder2[n];
groupOrder2[n]=groupOrder2[k];
groupOrder2[k]=temp;
}
n=groupOrder1.size()*groupOrder2.size();
while(n>1){
k=(rand() % n);
--n;
temp=groupOrder0[n];
groupOrder0[n]=groupOrder0[k];
groupOrder0[k]=temp;
}
int maxout=(sig-3)*(sig-3);
//Randomly go through pairs of group members
for(int tout=0; tout < maxout; tout++){
bout=2;
int i=groupOrder0[tout] % (sig-3);
int j= (groupOrder0[tout]-i)/(sig-3);
float rprob=(rand()/(float)RAND_MAX);
//printf("Probs: %8.3f - %8.3f,%8.3f\n",rprob,mols[groupOrder2[j]]->prob,mols[groupOrder1[i]]->prob);
//Find center of each group movement and save to branches
if(groups[movegroup][groupOrder1[i]] > groups[movegroup][groupOrder2[j]] +3){
mat=chooseMoves(groups[movegroup][groupOrder1[i]],groups[movegroup][groupOrder2[j]],dat,-1);
tot=mat.size();
for(int w=0; w < tot; w++){
if(batoms[mat[w].frag1]->check()==1 && batoms[mat[w].frag2]->check()==1){
pol2->copy(this);
pol2->changeFrag2(mat[w].pos1,mat[w].frag1,frags);
pol2->changeFrag2(mat[w].pos2,mat[w].frag2,frags);
pol2->coordSysNt(mat[w].pos1,&atoms[mat[w].frag1]->atoms);
pol2->coordSysNt(mat[w].pos2,&atoms[mat[w].frag2]->atoms);
pol2->updateFull();
pol2->coordSysAT(mat[w].pos1,&batoms[mat[w].frag1]->atoms);
pol2->coordSysAT(mat[w].pos2,&batoms[mat[w].frag2]->atoms);
vcstr=pol2->constraints(cons,0,ctemp); //Returns 1=Breaks constraints, 0=Satisfies Constraints
pol2->score=vcstr;
if(vcstr > 0){
if(vcstr > score){
this->copy(pol2);
printf("Constraint: %8.3f\t\tScore: %8.3f\n",ctemp,score);
Pdb * pdb =  new Pdb(pol2,'J');
pdb->write("sim1.pdb");
pdb->clear();
delete pdb;
bout=1;
break;
}else{
float prob1=exp(-1*(pol2->score)/0.3);
float prob2=exp(-1*score/0.3);
float ratio = prob2/prob1;
float rand1 = rand()/((float)RAND_MAX);
if(rand1 < ratio){
printf("Accepted Constraint: %8.3f\t\tScore: %8.3f \t\t",ctemp,score);
printf("Ratio: %8.3f:%8.3f\n",pol2->score,score);
this->copy(pol2);
Pdb * pdb =  new Pdb(pol2,'J');
pdb->write("sim1.pdb");
pdb->clear();
delete pdb;
bout=1;
break;
}
}
}

}
}
if(bout==1){break;}
}
}
if(bout==2){break;}
}
/* DELETE ALL POINTERS */
pol2->clearMols();
delete pol2;
//Clear uniqpol;
}




float Polymere::updateConstraints(std::vector< std::pair<int,int> > * constraints){
std::vector<double> params;
params.push_back(0);//type;
params.push_back(1);//rho;
params.push_back(2);//theta1;
params.push_back(3);//phi1;
params.push_back(4);//phi2;
params.push_back(5);//theta2;
params.push_back(6);//dx_i;
params.push_back(7);//dy_i;
params.push_back(8);//dz_i;
params.push_back(9);//typeRY;
int illegals=0;
float newcons,ang1,ang2,dd1,gang1,gang2,graa;
newcons=0;
for(int j=0; j<constraints->size(); j++){
	int pos1 = constraints->at(j).first;		
	int pos2 = constraints->at(j).second;		
        this->getParams(pos1, pos2,&params);
	float ang1=(float)atanh(cos(params[4]));
	float ang2=(float)atanh(cos(params[3]));
	float dd1=(float)log(pow(params[1],3));
//Get the score of the 3D gaussian for the distance and angle parameters
	float gang1=-1*pow((ang1-0.68),2);
	float gang2=-1*pow((ang2-0.68),2);
//	float gang1=0;
//	float gang2=0;
	float grho =-1*pow((dd1-6.5),2);
	float graa = (gang1+gang2+grho)/log(0.5);
		if(graa > newcons ){newcons = graa;} 
}
return newcons;
}

int Polymere::hardConstraintsL(std::vector< std::pair<int,int> > * cons){
//Determines if all pairs of NT's in the list meet the set of constraints
std::vector<double> params;
params.push_back(0);//type;
params.push_back(1);//rho;
params.push_back(2);//theta1;
params.push_back(3);//phi1;
params.push_back(4);//phi2;
params.push_back(5);//theta2;
params.push_back(6);//dx_i;
params.push_back(7);//dy_i;
params.push_back(8);//dz_i;
params.push_back(9);//typeRY;

int illegals=0;
for(int j=0; j<cons->size(); j++){
	int pos1 = cons->at(j).first;		
	int pos2 = cons->at(j).second;		
        this->getParams(pos1, pos2,&params);
	if(params[1] > 16 || params[4] < 4){illegals++;}
	if(params[4] > -0.5 || params[4] < -1.5){illegals++;}
	if(params[3] > -0.5 || params[3] < -1.5){illegals++;}
	if(params[2] > 2.0 || params[2] < 1.2){illegals++;}
	if(params[5] > 2.0 || params[5] < 1.2){illegals++;}
	if(params[8] > 3.0 || params[8] < -3.0){illegals++;}
	}	
return illegals;
}

int Polymere::hardConstraints(std::vector< std::pair<int,int> > * cons){
//Determines if all pairs of NT's in the list meet the set of constraints
std::vector<double> params;
params.push_back(0);//type;
params.push_back(1);//rho;
params.push_back(2);//theta1;
params.push_back(3);//phi1;
params.push_back(4);//phi2;
params.push_back(5);//theta2;
params.push_back(6);//dx_i;
params.push_back(7);//dy_i;
params.push_back(8);//dz_i;
params.push_back(9);//typeRY;


int illegals=0;
for(int j=0; j<cons->size(); j++){
	int pos1 = cons->at(j).first;		
	int pos2 = cons->at(j).second;		
        this->getParams(pos1, pos2,&params);
	if(params[1] < 11){
	//coplanar1
	if(params[1] > 8){
	//coplanar
	if(params[4] > -0.5 || params[4] < -1.5){illegals++;}
	if(params[3] > -0.5 || params[3] < -1.5){illegals++;}
	if(params[2] > 2.0 || params[2] < 1.2){illegals++;}
	if(params[5] > 2.0 || params[5] < 1.2){illegals++;}
	if(params[8] > 1.5 || params[8] < -1.5){illegals++;}
	}else{
	if(params[1] > 5){
	if(params[4] > 1 || params[4] < -1.5){illegals++;}
	if(params[3] > -0.5 || params[3] < -1.5){illegals++;}
	if(params[2] > 2.0 || params[2] < 1.2){illegals++;}
	if(params[5] > 2.0 || params[5] < 1.2){illegals++;}
	if(params[8] > 1.5 || params[8] < -1.5){illegals++;}
	}else{
	//stacked only stack if they are purine/purine combos or py/py
	if(params[9]==1 || params[9]==2){
	if(params[8] < 1.5 && params[8] > -1.5){illegals++;}
	if(abs(params[2]-params[5]) > 0.15){illegals++;}
}else{
illegals++;
}
	}	
	}
}else{
illegals=illegals+2;
}
if(illegals>0){
//break;
}
}

return illegals;
}

float Polymere::constraints(std::vector< std::pair<int, int> > * cons,int clash_level,float temp){
//std::pair stores the base-pairings
//clash_level - number of allowable steric clashes
//Temp is the current temperature used in the annealing process [100-0]
// LD > x > 7
//where LD = Seq_Dist
//initial distance is set to 1/2 sequence distance	
std::vector<double> params;
params.push_back(0);//type;
params.push_back(1);//rho;
params.push_back(2);//theta1;
params.push_back(3);//phi1;
params.push_back(4);//phi2;
params.push_back(5);//theta2;
params.push_back(6);//dx_i;
params.push_back(7);//dy_i;
params.push_back(8);//dz_i;
params.push_back(9);//typeRY;
float illegals=0;
float ang1,ang2,dd1,gang1,gang2,graa;
//illegals=this->checkBackBone();
if(illegals <= clash_level){
//illegals=this->atomClashes();
if(illegals <= clash_level){
illegals=11;
for(int j=0; j<cons->size(); j++){
	int pos1 = cons->at(j).first;		
	int pos2 = cons->at(j).second;		
        this->getParams(pos1, pos2,&params);
	float ang1=(float)atanh(cos(params[4]));
	float ang2=(float)atanh(cos(params[3]));
	float dd1=(float)log(pow(params[1],3));
//	printf("Rho: %8.3f\n",dd1);
//Get the score of the 3D gaussian for the distance and angle parameters
	float gang1=-1*pow((ang1-0.68),2)/(temp);
	float gang2=-1*pow((ang2-0.68),2)/(temp);
	float grho =-1*pow((dd1-6.5),2)/(temp);
	float graa = 10*exp(gang1+gang2+grho);	
	if(graa < illegals){illegals = graa;}
}
return illegals; 
}else{
//printf("Failed Atom Clashes %d\n",illegals);
illegals=-1;
}
}else{
//printf("Failed B-Bone Clashes %d\n",illegals);
illegals=-1;
}
if(illegals>=0){
//printf(" Dist: %8.3f\n",illegals);
} 
return illegals;
}


void Polymere::dynamix(int nmax, int movegroup, const std::vector< Fragment *>  frags,Fragclass * dat,const Polymere * pol,std::vector< Fragatoms *> atoms,std::vector<Fragall *> batoms,int nt1, int nt2){

using namespace std;
Polymere  * pol2 = new Polymere();     // Polymer used in the simulation. Everything is copied into and out of this class
std::vector< std::vector<int> > groups;// Movement groups of the polymere
int maxg=0;


for(int g=0; g<num_mols; g++){         // Determine the total number of movement groups. The user must define groups without any
int gc = pol->mols[g]->group;          // gaps. For example, 0 1 2 3 is ok, whereas 0 2 3 isn't
if(gc>maxg){
maxg=gc;
}
}
maxg++;
cout << "Max Groups: " << maxg << endl;

for(int i=1; i < maxg; i++){
std::vector<int> group;
for(int g=0; g<pol->num_mols; g++){
int gc = pol->mols[g]->group;
if(gc==i){
group.push_back(g);
}
}
groups.push_back(group);
}

int sig = groups[movegroup].size();   // Total number of positions in the movement group that will be used

for(int g=0; g<groups[movegroup].size(); g++){mols[groups[movegroup][g]]->group=0;}

score = 900000;
//Maximum attempts before giving up
int last_curr=nmax;
int escore=0;
int escore1=0;
for(int iteration=0; iteration < nmax; iteration++){
if(iteration > (last_curr+500)){if(escore < 0){break;}else{if(iteration > (last_curr+1000)){score=90000;}}}
//Copy into pol2
pol2->num_mols=num_mols;
pol2->randtot=0;
for(int t=0; t<pol->num_mols;  t++){
pol2->mols[t]-> z     =mols[t]->z;
pol2->mols[t]-> x     =mols[t]->x;
pol2->mols[t]-> y     =mols[t]->y; 
pol2->mols[t]-> seq   =mols[t]->seq;
//Adjust the probability for the movement group
pol2->mols[t]-> group  = mols[t]->group;
pol2->mols[t]-> bx    =mols[t]->bx;
pol2->mols[t]-> by    =mols[t]->by;
pol2->mols[t]-> bz    =mols[t]->bz;
pol2->mols[t]-> b2x   =mols[t]->b2x;
pol2->mols[t]-> b2y   =mols[t]->b2y;
pol2->mols[t]-> b2z   =mols[t]->b2z;
pol2->mols[t]-> id    =mols[t]->id;
if(pol2->mols[t]->group < 0){pol2->randtot=pol2->randtot++;}
}
for(int t=0; t<full_atoms.size(); t++){pol2->full_atoms[t]=full_atoms[t];}
pol2->score =score;
//Pick a random position
int krand=(rand() % sig);
//Check to make sure the groups are in focus
//Decide if the group should be moved
float prob=rand()/(float)(RAND_MAX);
float totalprob = (float)(-1*(pol2->mols[groups[movegroup][krand]]->group)*groups[movegroup].size())/((float)(iteration));
//printf("Iteration: %d:%8.3f:%8.3f  Group: %d K: %d\n",iteration,totalprob,prob,groups[movegroup][krand],krand);

if(totalprob < prob){
int fid=chooseSingle(groups[movegroup][krand],dat);
pol2->changeFrag2(groups[movegroup][krand],fid,frags);
pol2->coordSysNt(groups[movegroup][krand],&atoms[fid]->atoms);
pol2->updateFull();
pol2->coordSysAT(groups[movegroup][krand],&batoms[fid]->atoms);
int accept = 0;
//Check base-pair differences
int c1dist1=pol2->c1dist();
if(c1dist1!=0){
if(nt1!=nt2){
//Score is the distance between the two base pairs
double dist = pol2->moldist(nt1,nt2);
double score1=0;
if(dist < 14 & dist > 7){
//Check for acceptor donor classhes
double stmp=pol2->hydroG(&groups[movegroup]);
if(stmp > 0){
score1=(pol2->hydroP(nt1,nt2)+pol2->hydroP(nt2,nt1))+(pol2->hydroP((nt1+1),(nt2-1))+pol2->hydroP((nt2-1),(nt1+1)));
double h1=pol2->hydroP(nt1,nt2)+pol2->hydroP(nt2,nt1);
double h2=pol2->hydroP((nt1+1),(nt2-1))+pol2->hydroP((nt2-1),(nt1+1));
if(h1 > 2 & h2 > 2){
escore1=pol2->entropy(200,movegroup,frags,dat,pol,atoms,batoms);
score1=-50*escore1-stmp-10*score1;
escore1=-1*escore1;
if(escore1 < -10){printf("escore: %d\n",escore1);}
}else{
score1=-1*stmp-10*score1;
}

}else{
score1=-1*stmp;
}

}else{
score1 = pow((dist-10.5),2);
}
pol2->score = score1;
//Metropolis method to bring two bases close together
if(pol2->score <= score){
accept=1;
}else{
float prob1=exp(-1*(pol2->score));
float prob2=exp(-1*score);
float ratio = prob1/prob2;
float rand1 = rand()/((float)RAND_MAX);
if(rand1 < ratio){
accept=1;
}
}
}else{
double stmp=pol2->hydroG(&groups[movegroup]);
pol2->score = -3*stmp;
escore=0;
last_curr=iteration;
if(pol2->score <= score){
accept=1;
}else{
double prob1=exp(-1*(pol2->score));
double prob2=exp(-1*score);
double ratio = prob1/prob2;
double rand1 = rand()/((double)RAND_MAX);
if(rand1 < ratio){
accept=1;
}
}

}
}else{
//cout << "Failed C1dist: " << c1dist1 << endl;
accept=0;
}

if(accept==1){
escore=escore1;
last_curr=iteration;
randtot=0;
for(int t=0; t<num_mols;  t++){
mols[t]-> z     =pol2->mols[t]->z;
mols[t]-> x     =pol2->mols[t]->x;
mols[t]-> y     =pol2->mols[t]->y; 
mols[t]-> seq   =pol2->mols[t]->seq;
//Adjust the probability for the movement group
if(t==groups[movegroup][krand]){
mols[t]-> group   = pol2->mols[t]->group - 1;
}else{mols[t]-> group   =pol2->mols[t]->group;}
mols[t]-> bx    =pol2->mols[t]->bx;
mols[t]-> by    =pol2->mols[t]->by;
mols[t]-> bz    =pol2->mols[t]->bz;
mols[t]-> b2x   =pol2->mols[t]->b2x;
mols[t]-> b2y   =pol2->mols[t]->b2y;
mols[t]-> b2z   =pol2->mols[t]->b2z;
mols[t]-> id    =pol2->mols[t]->id;
if(mols[t]->group < 0){randtot++;}
}
for(int t=0; t<full_atoms.size(); t++){full_atoms[t]=pol2->full_atoms[t];}
score =pol2->score;
double rmsd1 = pol2->rmsAllc(pol);
printf("Current Score %d: %8.3f RMSD: %8.3f\n",iteration,pol2->score,rmsd1);
}

}

}


/* DELETE ALL POINTERS */

for(int m=0; m<pol2->mols.size(); m++){
delete pol2->mols[m];
}
pol2->full_atoms.clear();
pol2->mols.clear();
delete pol2;
//Clear uniqpol;
}

int Polymere::entropy(int nmax,int movegroup, const std::vector< Fragment *>  frags,Fragclass * dat,const Polymere * pol,std::vector< Fragatoms *> atoms,std::vector<Fragall *> batoms){
using namespace std;
//Returns the total number of possible closed moves for a given group

std::vector<Matches> mat; 
mat.reserve(7000);
//uniqpol is a vector of structures that results from a single closed moves search for a given structure

std::vector< std::vector<int> > groups;
int maxg=0;
//Read in the movement group from the ct file 
for(int g=0; g<num_mols; g++){
int gc = pol->mols[g]->group;
//cout << "group: " << gc << " bp1: " << pol->mols[g]->id << endl;
if(gc>maxg){
maxg=gc;
}
}
maxg++;
//cout << "Max Groups: " << maxg << endl;
for(int i=1; i < maxg; i++){
std::vector<int> group;
for(int g=0; g<num_mols; g++){
int gc = pol->mols[g]->group;
if(gc==i){
group.push_back(g);
}
}
groups.push_back(group);
}

int sig = groups[movegroup].size();
//cout << "Group Size: " << movegroup << ":" << sig << endl;
//cout << "Total Movement Groups: " << groups.size() << endl;

//Vectors containing the indices of valid movement groups
std::vector<int> groupOrder1;
std::vector<int> groupOrder2;
std::vector<int> groupOrder0;

for(int i=1; i < sig-2; i++){
groupOrder1.push_back(i);
groupOrder2.push_back(i);
}
int nt=(sig-3)*(sig-3);
for(int i=0; i <nt; i++){groupOrder0.push_back(i);} 
//Set initial Score and RMSD threshold
//cout << "Starting randomization: " << nmax << ":" << sig <<  endl;
//Foreach starting structure score them based on some criteria for each move.  Need to loop through 
//a certain number of cycles or until the simulation converges
//Add loop to perform multiple generations of movements
////////////////////////////////////////////////////////////
// Move current polymere, score, add improvements to uniq //
////////////////////////////////////////////////////////////
//Generate Random order for movements
int n = groupOrder1.size();
int k;
int temp;
while(n>1){
k=(rand() % n);
--n;
temp=groupOrder1[n];
groupOrder1[n]=groupOrder1[k];
groupOrder1[k]=temp;
}

n=groupOrder2.size();
while(n>1){
k=(rand() % n);
--n;
temp=groupOrder2[n];
groupOrder2[n]=groupOrder2[k];
groupOrder2[k]=temp;
}
n=groupOrder1.size()*groupOrder2.size();
while(n>1){
k=(rand() % n);
--n;
temp=groupOrder0[n];
groupOrder0[n]=groupOrder0[k];
groupOrder0[k]=temp;
}
//cout << endl;

int maxout=(sig-3)*(sig-3);
if(maxout > nmax){maxout=nmax;}

int entropy=0;
Polymere * pol2= new Polymere();
for(int tout=0; tout < maxout; tout++){
int tot;
int i=groupOrder0[tout] % (sig-3);
int j= (groupOrder0[tout]-i)/(sig-3);
//Loop through possible position combinations
//Only use positions that aren't set to -1 
//If position is less than -1 that is the cluster number to use
//As soon as move is identified make the move
//Find center of each group movement and save to branches
mat=chooseMoves(groups[movegroup][groupOrder1[i]],groups[movegroup][groupOrder2[j]],dat,-1);
tot=mat.size();
//printf("Identified Moves: %d\n",tot);
//Identify the center of cluster of moves  
for(int w=0; w<tot; w++){
pol2->num_mols=num_mols;
pol2->randtot=0;
for(int t=0; t< num_mols; t++){
pol2->mols[t]-> z     =mols[t]->z;
pol2->mols[t]-> x     =mols[t]->x;
pol2->mols[t]-> y     =mols[t]->y; 
pol2->mols[t]-> seq   =mols[t]->seq;
//Adjust the probability for the movement group
pol2->mols[t]->group=0;
pol2->mols[t]-> bx    =mols[t]->bx;
pol2->mols[t]-> by    =mols[t]->by;
pol2->mols[t]-> bz    =mols[t]->bz;
pol2->mols[t]-> b2x   =mols[t]->b2x;
pol2->mols[t]-> b2y   =mols[t]->b2y;
pol2->mols[t]-> b2z   =mols[t]->b2z;
pol2->mols[t]-> id    =mols[t]->id;
}
pol2->randtot=0;
for(int t=0; t<full_atoms.size(); t++){pol2->full_atoms[t]=full_atoms[t];}
pol2->changeFrag2(mat[w].pos1,mat[w].frag1,frags);
pol2->changeFrag2(mat[w].pos2,mat[w].frag2,frags);
pol2->coordSysNt(mat[w].pos1,&atoms[mat[w].frag1]->atoms);
pol2->coordSysNt(mat[w].pos2,&atoms[mat[w].frag2]->atoms);
pol2->updateFull();
pol2->coordSysAT(mat[w].pos1,&batoms[mat[w].frag1]->atoms);
pol2->coordSysAT(mat[w].pos2,&batoms[mat[w].frag2]->atoms);

/*
printf("c1 batoms: %8.3f %8.3f %8.3f\n",batoms[mat[w].frag2]->atoms[11].x,batoms[mat[w].frag2]->atoms[23].x,batoms[mat[w].frag2]->atoms[35].x); 
printf("c1 atoms: %8.3f %8.3f %8.3f\n",atoms[mat[w].frag2]->atoms[0].x,atoms[mat[w].frag2]->atoms[3].x,atoms[mat[w].frag2]->atoms[6].x ); 
*/




//Check base-pair differences
int c1dist1=pol2->c1dist();
double id=0;
if(c1dist1!=0){id=pol2->checkhbonds();}else{id=0;}
if(id>0){entropy++; break;}
}



}
return entropy;
}

int Polymere::modelit(int cycles,int capture, int movegroup, const std::vector< Fragment *>  frags,Fragclass * dat,const Polymere * pol,std::vector< Fragatoms *> atoms,std::vector<Fragall *> batoms,Parser * model,vector<double> * dparams,string fileout,vector<int> * parent, vector<int> * child,string out_id){
//Returns star_center if it is below the rmsd_res already, otherwise it returns -1 and the current polymere is the new center that should be added to the file
using namespace std;
int bout=1;
Polymere  * pol2 = new Polymere();
pol2->copy(this);
Polymere  * outpol = new Polymere();
double minscore = 99999999;
score = 99999999;
for(int g=0; g<pol2->num_mols; g++){pol2->mols[g]->prob=1;mols[g]->prob=1;}
//Need to capture the last "capture" polymeres
int curr_size=0;
std::vector<Polymere  *> captured;

std::vector<Matches> mat; 
mat.reserve(7000);
double id=0;

//uniqpol is a vector of structures that results from a single closed moves search for a given structure
std::vector< std::vector<int> > groups;
int maxg=0;
//Read in the movement group from the ct file 
for(int g=0; g<num_mols; g++){
int gc = pol->mols[g]->group;
//cout << "group: " << gc << " bp1: " << pol->mols[g]->id << endl;
if(gc>maxg){
maxg=gc;
}
}
maxg++;

for(int i=1; i < maxg; i++){
std::vector<int> group;
for(int g=0; g<num_mols; g++){
int gc = pol->mols[g]->group;
if(gc==i){
group.push_back(g);
}
}
groups.push_back(group);
}

int sig = groups[movegroup].size();
//Set counts to 0 for each position
for(int g=0; g<groups[movegroup].size(); g++){
mols[groups[movegroup][g]]->group=0;
}

int tot;

//Initializing vectors contianing movement group information
//Vectors containing the indices of valid movement groups
std::vector<int> groupOrder1;
std::vector<int> groupOrder2;
std::vector<int> groupOrder0;

for(int i=1; i < sig-2; i++){
groupOrder1.push_back(i);
groupOrder2.push_back(i);
}
int nt=(sig-3)*(sig-3);
for(int i=0; i <nt; i++){groupOrder0.push_back(i);} 


//Reading scoring function

for(int iteration=0; iteration < cycles; iteration++){
////////////////////////////////////////////////////////////
// Move current polymere, score, add improvements to uniq //
////////////////////////////////////////////////////////////
int n = groupOrder1.size();
int k;
int temp;
while(n>1){
k=(rand() % n);
--n;
temp=groupOrder1[n];
groupOrder1[n]=groupOrder1[k];
groupOrder1[k]=temp;
}

n=groupOrder2.size();
while(n>1){
k=(rand() % n);
--n;
temp=groupOrder2[n];
groupOrder2[n]=groupOrder2[k];
groupOrder2[k]=temp;
}
n=groupOrder1.size()*groupOrder2.size();
while(n>1){
k=(rand() % n);
--n;
temp=groupOrder0[n];
groupOrder0[n]=groupOrder0[k];
groupOrder0[k]=temp;
}

int maxout=(sig-3)*(sig-3);
//Randomly go through pairs of group members
bout=1;
for(int tout=0; tout < maxout; tout++){
int i=groupOrder0[tout] % (sig-3);
int j= (groupOrder0[tout]-i)/(sig-3);
float rprob=(rand()/(float)RAND_MAX);
//printf("Probs: %8.3f - %8.3f,%8.3f\n",rprob,mols[groupOrder2[j]]->prob,mols[groupOrder1[i]]->prob);
if((groupOrder1[i]+3) < groupOrder2[j] & 2*rprob < (mols[groupOrder2[j]]->prob+mols[groupOrder1[i]]->prob)){
//Find center of each group movement and save to branches
mat=chooseMoves(groups[movegroup][groupOrder1[i]],groups[movegroup][groupOrder2[j]],dat,-1);
tot=mat.size();
//printf("Identified Moves: %d\n",tot);
//Identify the center of cluster of moves  
for(int w=0; w < tot; w++){
id=0;
//cout << "Creating new structure " << endl;
//Make a copy of the molecule array
//Use Group to track randomization
pol2->num_mols=num_mols;
pol2->randtot=0;
for(int t=0; t< num_mols; t++){
pol2->mols[t]-> z     =mols[t]->z;
pol2->mols[t]-> x     =mols[t]->x;
pol2->mols[t]-> y     =mols[t]->y; 
pol2->mols[t]-> seq   =mols[t]->seq;
//Adjust the probability for the movement group
if(t==groups[movegroup][groupOrder1[i]] || t==groups[movegroup][groupOrder2[j]]){
pol2->mols[t]-> group   = mols[t]->group - 1;
}else{pol2->mols[t]-> group   =mols[t]->group;}
pol2->mols[t]-> bx    =mols[t]->bx;
pol2->mols[t]-> by    =mols[t]->by;
pol2->mols[t]-> bz    =mols[t]->bz;
pol2->mols[t]-> b2x   =mols[t]->b2x;
pol2->mols[t]-> b2y   =mols[t]->b2y;
pol2->mols[t]-> b2z   =mols[t]->b2z;
pol2->mols[t]-> id    =mols[t]->id;
pol2->mols[t]-> prob    =mols[t]->prob;
if(pol2->mols[t]->group < 0){pol2->randtot=pol2->randtot++;}
}
for(int t=0; t<full_atoms.size(); t++){pol2->full_atoms[t]=full_atoms[t];}

pol2->changeFrag2(mat[w].pos1,mat[w].frag1,frags);
pol2->changeFrag2(mat[w].pos2,mat[w].frag2,frags);
pol2->coordSysNt(mat[w].pos1,&atoms[mat[w].frag1]->atoms);
pol2->coordSysNt(mat[w].pos2,&atoms[mat[w].frag2]->atoms);
pol2->updateFull();
pol2->coordSysAT(mat[w].pos1,&batoms[mat[w].frag1]->atoms);
pol2->coordSysAT(mat[w].pos2,&batoms[mat[w].frag2]->atoms);

int c1dist1=pol2->c1dist();
//Check to make sure it is different than the last 5 moves
double ormsd=999999;
double ormsd1=rmsAll(pol2);
if(ormsd1<ormsd){ormsd=ormsd1;}
//printf("rmsd: %8.3f\n",ormsd1);
if(ormsd > 0.1){
if(c1dist1!=0){
id=pol2->checkhbonds();
}else{id=0;}
//Determine minimal RMSD from any 
if(id>0){
//ADD SCORING SYSTEM 
Score * scoreclass = new Score(pol2->num_mols);
pol2->resscore(scoreclass,dparams);
    // create a parser object
std::vector<double> results;
model->score(scoreclass,&results);
std::vector<float> labmeeting;
pol2->localRMSD(&labmeeting,pol,parent,child);
double gscore = 0;
float mgscore=0; 
for(int g=0; g<labmeeting.size(); g++){gscore = gscore + (double)labmeeting[g]; if(labmeeting[g] > mgscore){mgscore=labmeeting[g];}}
pol2->score = 10*pol2->rmsAllc(pol);
delete scoreclass;
for(int g=0; g< pol2->num_mols; g++){
pol2->mols[g]->prob=labmeeting[g]/mgscore;
//printf("prob %d: %8.3f\n",g,pol2->mols[g]->prob);
}
//Set each group to the residue level score.  Try to normalize the score across the structure
if(pol2->score <= score){
id=1;
}else{
id=0;
double ratio=exp(-1*(pol2->score))/exp(-1*score);
double rand1 = rand()/((double)RAND_MAX);
printf("%8.3f : %8.3f\n",ratio,rand1);
if(rand1 < ratio){
id=1;
}
}
}
}else{
id=0;
}
//Make a copy of pol2 and save
if(id!=0){
printf("minscore: %8.3f \n",pol2->score);
if(pol2->score < minscore){
printf("MINSCORE: %8.3f \n",pol2->score);
minscore = pol2->score;
//Copy to outpol
//printf("Min Score: %8.3f\n", minscore);
outpol->copy(pol2);
}
//Add to last five
if(curr_size < capture){
captured.push_back(new Polymere());
curr_size++;
}else{
//Delete first polymere and push back new polymere
captured[0]->clearMols();
captured.erase(captured.begin());
captured.push_back(new Polymere());
//After you have captured all that are required quit simulating
}
captured.back()->copy(pol2);
//Update current conformation 
copy(pol2);
bout=0;
break;
}
}//Loop through all matches
}//If groups are 3 apart
if(bout==0){break;}
}//Loop through random group pairs
if(bout==1){
break;
}

}//Iteration
copy(outpol);


//Global Scores
string gfile=fileout;
gfile.append("/gscore/gscore_");
gfile.append(out_id);
gfile.append(".txt");
FILE * sfile = fopen(gfile.c_str(),"a");
for(int c=0; c<captured.size(); c++){
double grmsd=captured[c]->rmsAllc(pol);
fprintf(sfile,"%s %d %8.3f %8.3f\n",out_id.c_str(),c,captured[c]->score,grmsd);
}
fflush( sfile );
fclose(sfile);
//Local Scores
//Write out the score for each captured structure
for(int c=0; c<captured.size(); c++){
Score * scoreclass = new Score(captured[c]->num_mols);
captured[c]->resscore(scoreclass,dparams);
vector<float> lrmsds;
//printf("Initial RMSD: %8.3f\n",r0);
captured[c]->localRMSD(&lrmsds,pol,parent,child);
for(int s =0; s<scoreclass->s0.size(); s++){
scoreclass->local_rmsd[s]=lrmsds[s];
}
scoreclass->sim_id=c;
string ofile=fileout;
ofile.append("/lscore/lscore_");
ofile.append(out_id);
ofile.append(".txt");
scoreclass->write(ofile);
delete scoreclass;
}

//Write captured structures to output file
for(int c=0; c<captured.size(); c++){
Pdb * pdb = new Pdb(captured[c],'j');
string cfile=fileout;
cfile.append("/captured/pdb_");
cfile.append(out_id);
cfile.append(".pdb");
pdb->write(cfile);
pdb->clear();
delete pdb;
}

//Write minimal score file
Pdb * pdb = new Pdb(outpol,'j');
string mfile=fileout;
mfile.append("/minimal/min_");
mfile.append(out_id);
mfile.append(".pdb");
pdb->write(mfile);
pdb->clear();
delete pdb;

for(int m=0; m<pol2->mols.size(); m++){
delete pol2->mols[m];
}
pol2->full_atoms.clear();
pol2->mols.clear();
delete pol2;

for(int m=0; m<outpol->mols.size(); m++){
delete outpol->mols[m];
}
outpol->full_atoms.clear();
outpol->mols.clear();
delete outpol;
return 0;
}

float Polymere::monteCons(int nmax, int movegroup, const std::vector< Fragment *>  frags,Fragclass * dat,const Polymere * pol,std::vector< Fragatoms *> atoms,std::vector<Fragall *> batoms,std::vector<double> * dparams,int stype,std::vector< std::vector <std::pair<int,int> > > allcons,int sim_id,std::string out_prefix){
//Returns star_center if it is below the rmsd_res already, otherwise it returns -1 and the current polymere is the new center that should be added to the file
using namespace std;
FILE * sfile = fopen(out_prefix.c_str(),"a");
int bout=1;
int time=0;

Polymere  * pol2 = new Polymere();
pol2->copy(this);
Polymere  * outpol = new Polymere();
double minscore = 99999999;
score = 99999999;
outpol->copy(this);
vector<Polymere  *> last5;
std::vector<Matches> mat; 
mat.reserve(7000);
double id=0;
//uniqpol is a vector of structures that results from a single closed moves search for a given structure
std::vector< std::vector<int> > groups;
int maxg=0;
//Read in the movement group from the ct file 
for(int g=0; g<num_mols; g++){
int gc = pol2->mols[g]->group;
//cout << "group: " << gc << " bp1: " << pol->mols[g]->id << endl;
if(gc>maxg){
maxg=gc;
}
}
maxg++;
//cout << "Max Groups: " << maxg << endl;
for(int i=0; i < maxg; i++){
std::vector<int> group;
for(int g=0; g<num_mols; g++){
int gc = pol2->mols[g]->group;
if(gc==i){
group.push_back(g);
}
}
groups.push_back(group);
}

int sig = groups[movegroup].size();
printf("move group size: %d\n",sig);

int tot;
//Vectors containing the indices of valid movement groups
std::vector<int> groupOrder1;
std::vector<int> groupOrder2;
std::vector<int> groupOrder0;

for(int i=1; i < sig-2; i++){
groupOrder1.push_back(i);
groupOrder2.push_back(i);
}
int nt=(sig-3)*(sig-3);
for(int i=0; i <nt; i++){groupOrder0.push_back(i);} 

for(int iteration=0; iteration < nmax; iteration++){
//cout << "Iteration: " << iteration << endl;
////////////////////////////////////////////////////////////
// Move current polymere, score, add improvements to uniq //
////////////////////////////////////////////////////////////
//Generate Random order for movements
int n = groupOrder1.size();
int k;
int temp;
while(n>1){
k=(rand() % n);
--n;
temp=groupOrder1[n];
groupOrder1[n]=groupOrder1[k];
groupOrder1[k]=temp;
}

n=groupOrder2.size();
while(n>1){
k=(rand() % n);
--n;
temp=groupOrder2[n];
groupOrder2[n]=groupOrder2[k];
groupOrder2[k]=temp;
}
n=groupOrder1.size()*groupOrder2.size();
while(n>1){
k=(rand() % n);
--n;
temp=groupOrder0[n];
groupOrder0[n]=groupOrder0[k];
groupOrder0[k]=temp;
}
//cout << "Group1: ";
for(int i=0; i< groupOrder1.size(); i++){
//cout << groupOrder1[i] << " ";
}
//cout << endl;
int maxout=(sig-3)*(sig-3);

//Randomly go through pairs of group members
bout=1;
for(int tout=0; tout < maxout; tout++){
int i=groupOrder0[tout] % (sig-3);
int j= (groupOrder0[tout]-i)/(sig-3);


if((groupOrder1[i]+3) < groupOrder2[j]){
//Find center of each group movement and save to branches
mat=chooseMoves(groups[movegroup][groupOrder1[i]],groups[movegroup][groupOrder2[j]],dat,-1);
tot=mat.size();
//Max out total at 1000
//printf("Total Moves: %d\n",tot);
if(tot > 500){
//printf("Total Moves: %d\n",tot);
tot=500;
}
printf("Identified Moves: %d\n",tot);
//Identify the center of cluster of moves  
//printf("Total: %d\n",tot);
for(int w=0; w < tot; w++){
float g0,g1,g2,g3,g4,g5,g6;
id=0;
if(batoms[mat[w].frag1]->check()==1 && batoms[mat[w].frag2]->check()==1){
//cout << "Creating new structure " << endl;
//Make a copy of the molecule array
//Use Group to track randomization
pol2->copy(this);
pol2->changeFrag2(mat[w].pos1,mat[w].frag1,frags);
pol2->changeFrag2(mat[w].pos2,mat[w].frag2,frags);
pol2->coordSysNt(mat[w].pos1,&atoms[mat[w].frag1]->atoms);
pol2->coordSysNt(mat[w].pos2,&atoms[mat[w].frag2]->atoms);
pol2->updateFull();
pol2->coordSysAT(mat[w].pos1,&batoms[mat[w].frag1]->atoms);
pol2->coordSysAT(mat[w].pos2,&batoms[mat[w].frag2]->atoms);
int c1dist1=1;
//Check to make sure it is different than the last 5 moves
double ormsd=999999;
//for(int j=0; j<curr_size; j++){
double ormsd1=this->rmsAllg(pol2,movegroup);
if(ormsd1<ormsd){ormsd=ormsd1;}
//printf("  rmsd: %8.3f\n",ormsd1); 
//}

if(ormsd > 0.001){
if(c1dist1!=0){
int aclash=pol2->groupCheckClashes(movegroup);
int bclash=pol2->groupCheckBackBone(movegroup);
//printf("Clashes: %d %d\n",aclash,bclash);
if(aclash > 500 || bclash > 0){id=0;}else{id=1;}
}else{id=0;}
//Determine minimal RMSD from any 
if(id>0){
//Add scoring function
float overall_score=0;
if(stype<4){
Score * pscore = new Score(sig);
pol2->group_resscore(pscore,dparams,movegroup);
//float stacked = pol2->gaussDbl(dparams,3,3,3,0.2,0.2,&group_inds);
//float hydro = pol2->hydroG(&group_inds);
//pol2->score = stacked-hydro;
g0=pscore->gscore(0);//Doublet
g1=pscore->gscore(1);//Graph Order
g2=pscore->gscore(2);//Base Packing
g3=pscore->gscore(3);//Backbone Packing
g4=pscore->gscore(4);//Hydrogen Bond
g5=pscore->gscore(5);//Surface Area
g6=pscore->gscore(6);//Low resolution score
if(g6 < -10000){
pscore->write_globals("bad_globals.txt");
overall_score=9999; 
}
if(stype==0){
//Set score to total hydrogen
overall_score=overall_score - 10*(g4/sig);
}
if(stype==1){
//doublet 
overall_score=overall_score - 10*(g0/sig); 
}
if(stype==2){
//packing
float ogscore1=(g1+g2+g3+g5)/sig;
float ogscore2=g6/sig;
overall_score = (ogscore1) + ((1/(1+exp(-1*(ogscore1-20)))))*(12) + (-1*(1/(1+exp(-1*(ogscore1-20)))-1))*(ogscore2);
}
if(stype==3){
//combination score
float ogscore1=(g1+g2+g3+g5)/sig;
float ogscore2=g6/sig;
float ogscore3=g4/sig;
overall_score = (ogscore1) + ((1/(1+exp(-1*(ogscore1-20)))))*(12) + (-1*(1/(1+exp(-1*(ogscore1-20)))-1))*(ogscore2) +  ((-1*(1/(1+exp(-1*(ogscore2-15)))-1))*(-1*(1/(1+exp(-1*(ogscore1-20)))-1)))*(-10*ogscore3);
}
delete pscore;
}else{
if(stype==4){
//Set score to total hydrogen
overall_score=overall_score + (float)(pol2->checkBackBone() + pol2->atomClashes());
}
}
//printf("Overall Score: %8.3f\n",overall_score);

pol2->score=overall_score;
//printf("     score: %8.3f - %8.3f\n",pol2->score,score);
if(pol2->score <= score || iteration < 20 ){
//printf("     accept score: %8.3f - %8.3f\n",pol2->score,score);
id=1;
}else{
id=0;
double ratio=exp(-1*(pol2->score))/exp(-1*score);
double rand1 = rand()/((double)RAND_MAX);
if(rand1 < ratio){
id=1;
}
}
}
}else{
id=0;
}
//Make a copy of pol2 and save
if(id!=0){
//printf("minscore: %8.3f \n",pol2->score);
std::vector< std::pair<int,int> > mcons= pol2->confMatch(allcons);
//printf("mcons[0]: %d\n",mcons[0].second);
if(mcons[0].second < 1){
if(pol2->score < minscore){
minscore = pol2->score;
outpol->copy(pol2);
}
double nat_rmsd=pol2->rmsAllg(pol,movegroup);
//Write to simpath file
//sim_id,time,conf_id1,conf_id2,conf_id3,conf_id4,frag1,frag2,pos1,pos2,rmsd_native,g0,g1..g6,overall_score
for(int mc=mcons.size(); mc<5; mc++){
mcons.push_back(std::pair<int,int>(-1,-1));
}
fprintf(sfile,"%d %d %d %d %d %d %d %d %d %d %d %8.3f  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f\n",
sim_id,
time++,
mcons[0].first, 
mcons[1].first, 
mcons[2].first, 
mcons[3].first, 
mcons[4].first, 
mat[w].frag1,
mat[w].frag2,
mat[w].pos1,
mat[w].pos2,
nat_rmsd,
g0,g1,g2,g3,g4,g5,g6,pol2->score);
fflush(sfile);
//Add to last five
this->copy(pol2);
//Pdb * npdb=new Pdb(pol2,'J');
//npdb->write("paths.ent");
//npdb->clear();
//delete npdb;
bout=0;
//if(score < (-5*rmsd_res*1.5)){bout=5;}
//Stop the simulation if something happens: set bout=5
break;
}
}
}//Loop through all matches
}

}//If groups are 3 apart
if(bout==0 || bout==5){break;}
}//Loop through random group pairs
if(bout==1 || bout==5){
//No more moves are possible 
break;
}

//if(minscore < (-5*(rmsd_res*2))){break;} 
}//Iteration

this->copy(outpol);
fclose(sfile);
printf("Final Structure: %8.3f Replaced: %d\n",outpol->score, outpol->randtot);

for(int m=0; m<pol2->mols.size(); m++){
delete pol2->mols[m];
}
pol2->full_atoms.clear();
pol2->mols.clear();
delete pol2;

for(int m=0; m<outpol->mols.size(); m++){
delete outpol->mols[m];
}
outpol->full_atoms.clear();
outpol->mols.clear();
delete outpol;
minscore=minscore/5;
return minscore;
} 


float Polymere::monte(int nmax, int movegroup, const std::vector< Fragment *>  frags,Fragclass * dat,const Polymere * pol,std::vector< Fragatoms *> atoms,std::vector<Fragall *> batoms,std::vector<double> * dparams,int stype){
//Returns star_center if it is below the rmsd_res already, otherwise it returns -1 and the current polymere is the new center that should be added to the file
using namespace std;
int bout=1;
Polymere  * pol2 = new Polymere();
pol2->copy(this);
Polymere  * outpol = new Polymere();
double minscore = 99999999;
score = 99999999;
outpol->copy(this);
vector<Polymere  *> last5;
std::vector<Matches> mat; 
mat.reserve(7000);
double id=0;
//uniqpol is a vector of structures that results from a single closed moves search for a given structure
std::vector< std::vector<int> > groups;
int maxg=0;
//Read in the movement group from the ct file 
for(int g=0; g<num_mols; g++){
int gc = pol->mols[g]->group;
//cout << "group: " << gc << " bp1: " << pol->mols[g]->id << endl;
if(gc>maxg){
maxg=gc;
}
}
maxg++;
//cout << "Max Groups: " << maxg << endl;
for(int i=1; i < maxg; i++){
std::vector<int> group;
for(int g=0; g<num_mols; g++){
int gc = pol->mols[g]->group;
if(gc==i){
group.push_back(g);
}
}
groups.push_back(group);
}

int sig = groups[movegroup].size();

//Set counts to 0 for each position
for(int g=0; g<groups[movegroup].size(); g++){
mols[groups[movegroup][g]]->group=0;
}

int tot;
//Vectors containing the indices of valid movement groups
std::vector<int> groupOrder1;
std::vector<int> groupOrder2;
std::vector<int> groupOrder0;

for(int i=1; i < sig-2; i++){
groupOrder1.push_back(i);
groupOrder2.push_back(i);
}
int nt=(sig-3)*(sig-3);
for(int i=0; i <nt; i++){groupOrder0.push_back(i);} 

for(int iteration=0; iteration < nmax; iteration++){
//cout << "Iteration: " << iteration << endl;
////////////////////////////////////////////////////////////
// Move current polymere, score, add improvements to uniq //
////////////////////////////////////////////////////////////
//Generate Random order for movements
int n = groupOrder1.size();
int k;
int temp;
while(n>1){
k=(rand() % n);
--n;
temp=groupOrder1[n];
groupOrder1[n]=groupOrder1[k];
groupOrder1[k]=temp;
}

n=groupOrder2.size();
while(n>1){
k=(rand() % n);
--n;
temp=groupOrder2[n];
groupOrder2[n]=groupOrder2[k];
groupOrder2[k]=temp;
}
n=groupOrder1.size()*groupOrder2.size();
while(n>1){
k=(rand() % n);
--n;
temp=groupOrder0[n];
groupOrder0[n]=groupOrder0[k];
groupOrder0[k]=temp;
}
//cout << "Group1: ";
for(int i=0; i< groupOrder1.size(); i++){
//cout << groupOrder1[i] << " ";
}
//cout << endl;
int maxout=(sig-3)*(sig-3);

//Randomly go through pairs of group members
bout=1;
for(int tout=0; tout < maxout; tout++){
int i=groupOrder0[tout] % (sig-3);
int j= (groupOrder0[tout]-i)/(sig-3);


if((groupOrder1[i]+3) < groupOrder2[j]){
//Find center of each group movement and save to branches
mat=chooseMoves(groups[movegroup][groupOrder1[i]],groups[movegroup][groupOrder2[j]],dat,-1);
tot=mat.size();
//Max out total at 1000
if(tot > 500){
printf("Total Moves: %d\n",tot);
tot=500;
}
//printf("Identified Moves: %d\n",tot);
//Identify the center of cluster of moves  
//printf("Total: %d\n",tot);
for(int w=0; w < tot; w++){
id=0;
if(batoms[mat[w].frag1]->check()==1 && batoms[mat[w].frag2]->check()==1){
//cout << "Creating new structure " << endl;
//Make a copy of the molecule array
//Use Group to track randomization
pol2->num_mols=num_mols;
pol2->randtot=0;
for(int t=0; t< num_mols; t++){

pol2->mols[t]-> z     =mols[t]->z;
pol2->mols[t]-> x     =mols[t]->x;
pol2->mols[t]-> y     =mols[t]->y; 
pol2->mols[t]-> seq   =mols[t]->seq;
//Adjust the probability for the movement group
if(t==groups[movegroup][groupOrder1[i]] || t==groups[movegroup][groupOrder2[j]]){
pol2->mols[t]-> group   = mols[t]->group - 1;
}else{pol2->mols[t]-> group   =mols[t]->group;}
pol2->mols[t]-> bx    =mols[t]->bx;
pol2->mols[t]-> by    =mols[t]->by;
pol2->mols[t]-> bz    =mols[t]->bz;
pol2->mols[t]-> b2x   =mols[t]->b2x;
pol2->mols[t]-> b2y   =mols[t]->b2y;
pol2->mols[t]-> b2z   =mols[t]->b2z;
pol2->mols[t]-> id    =mols[t]->id;
if(pol2->mols[t]->group < 0){pol2->randtot=pol2->randtot++;}
}
for(int t=0; t<full_atoms.size(); t++){pol2->full_atoms[t]=full_atoms[t];}
pol2->changeFrag2(mat[w].pos1,mat[w].frag1,frags);
pol2->changeFrag2(mat[w].pos2,mat[w].frag2,frags);
pol2->coordSysNt(mat[w].pos1,&atoms[mat[w].frag1]->atoms);
pol2->coordSysNt(mat[w].pos2,&atoms[mat[w].frag2]->atoms);
pol2->updateFull();
pol2->coordSysAT(mat[w].pos1,&batoms[mat[w].frag1]->atoms);
pol2->coordSysAT(mat[w].pos2,&batoms[mat[w].frag2]->atoms);
int c1dist1=pol2->c1dist();
//Check to make sure it is different than the last 5 moves
double ormsd=999999;
//for(int j=0; j<curr_size; j++){
double ormsd1=rmsAllc(pol2);
if(ormsd1<ormsd){ormsd=ormsd1;}
//}

if(ormsd > 0.001){
if(c1dist1!=0){
id=pol2->checkhbonds();
}else{id=0;}
//Determine minimal RMSD from any 
if(id>0){
//Add scoring function
float overall_score=0;
if(stype<4){
Score * pscore = new Score(pol2->num_mols);
pol2->resscore(pscore,dparams);
//float stacked = pol2->gaussDbl(dparams,3,3,3,0.2,0.2,&group_inds);
//float hydro = pol2->hydroG(&group_inds);
//pol2->score = stacked-hydro;
float g0=pscore->gscore(0);//Doublet
float g1=pscore->gscore(1);//Graph Order
float g2=pscore->gscore(2);//Base Packing
float g3=pscore->gscore(3);//Backbone Packing
float g4=pscore->gscore(4);//Hydrogen Bond
float g5=pscore->gscore(5);//Surface Area
float g6=pscore->gscore(6);//Low resolution score
if(g6 < -10000){
pscore->write_globals("bad_globals.txt");
overall_score=9999; 
}
if(stype==0){
//Set score to total hydrogen
overall_score=overall_score - 10*(g4/pol->num_mols);
}
if(stype==1){
//doublet 
overall_score=overall_score - 10*(g0/pol->num_mols); 
}
if(stype==2){
//packing
float ogscore1=(g1+g2+g3+g5)/pol->num_mols;
float ogscore2=g6/pol->num_mols;
overall_score = (ogscore1) + ((1/(1+exp(-1*(ogscore1-20)))))*(12) + (-1*(1/(1+exp(-1*(ogscore1-20)))-1))*(ogscore2);
}
if(stype==3){
//combination score
float ogscore1=(g1+g2+g3+g5)/pol->num_mols;
float ogscore2=g6/pol->num_mols;
float ogscore3=g4/pol->num_mols;
overall_score = (ogscore1) + ((1/(1+exp(-1*(ogscore1-20)))))*(12) + (-1*(1/(1+exp(-1*(ogscore1-20)))-1))*(ogscore2) +  ((-1*(1/(1+exp(-1*(ogscore2-15)))-1))*(-1*(1/(1+exp(-1*(ogscore1-20)))-1)))*(-10*ogscore3);
}
delete pscore;
}else{
if(stype==4){
//Set score to total hydrogen
overall_score=overall_score + (float)(pol2->checkBackBone() + pol2->atomClashes());
}
}
//printf("Overall Score: %8.3f\n",overall_score);

pol2->score=overall_score;

if(pol2->score <= score || iteration < 5 ){
id=1;
}else{
id=0;
double ratio=exp(-1*(pol2->score))/exp(-1*score);
double rand1 = rand()/((double)RAND_MAX);
if(rand1 < ratio){
id=1;
}
}
}
}else{
id=0;
}
//Make a copy of pol2 and save
if(id!=0){
//printf("minscore: %8.3f \n",pol2->score);
if(pol2->score < minscore){
minscore = pol2->score;
//Copy to outpol
//printf("Min Score: %8.3f\n", minscore);
for(int w=0; w<pol2->num_mols; w++){
outpol->mols[w]-> z     =pol2->mols[w]->z;
outpol->mols[w]-> x     =pol2->mols[w]->x;
outpol->mols[w]-> y     =pol2->mols[w]->y; 
outpol->mols[w]-> group =pol2->mols[w]->group;
outpol->mols[w]-> seq   =pol2->mols[w]->seq;
outpol->mols[w]-> bx     =pol2->mols[w]-> bx  ;
outpol->mols[w]-> by     =pol2->mols[w]-> by  ;
outpol->mols[w]-> bz    =pol2->mols[w]-> bz  ;
outpol->mols[w]-> b2x =pol2->mols[w]-> b2x ;
outpol->mols[w]-> b2y   =pol2->mols[w]-> b2y ;
outpol->mols[w]-> b2z   =pol2->mols[w]-> b2z ;
outpol->mols[w]-> id   =pol2->mols[w]-> id  ;
}
outpol->num_mols=pol2->num_mols;
outpol->randtot=pol2->randtot;
outpol->score = pol2->score;
for(int t=0; t<pol2->full_atoms.size(); t++){outpol->full_atoms[t]=pol2->full_atoms[t];}
}

//Add to last five


for(int w=0; w<pol2->num_mols; w++){
mols[w]-> z     =pol2->mols[w]->z;
mols[w]-> x     =pol2->mols[w]->x;
mols[w]-> y     =pol2->mols[w]->y; 
mols[w]-> group =pol2->mols[w]->group;
mols[w]-> seq   =pol2->mols[w]->seq;
mols[w]-> bx     =pol2->mols[w]-> bx  ;
mols[w]-> by     =pol2->mols[w]-> by  ;
mols[w]-> bz    =pol2->mols[w]-> bz  ;
mols[w]-> b2x =pol2->mols[w]-> b2x ;
mols[w]-> b2y   =pol2->mols[w]-> b2y ;
mols[w]-> b2z   =pol2->mols[w]-> b2z ;
mols[w]-> id   =pol2->mols[w]-> id  ;
}
num_mols=pol2->num_mols;
randtot=pol2->randtot;
score = pol2->score;
for(int t=0; t<pol2->full_atoms.size(); t++){full_atoms[t]=pol2->full_atoms[t];}
bout=0;
float minrmsd = pol2->rmsAllc(pol);
cout << " RMSD: " << minrmsd  << " Score: " << pol2->score << endl;

//if(score < (-5*rmsd_res*1.5)){bout=5;}
//Stop the simulation if something happens: set bout=5
break;
}

}//Loop through all matches
}

}//If groups are 3 apart
if(bout==0 || bout==5){break;}
}//Loop through random group pairs
if(bout==1 || bout==5){
//No more moves are possible 
break;
}

//if(minscore < (-5*(rmsd_res*2))){break;} 
}//Iteration

for(int t=0; t< num_mols; t++){
mols[t]-> z     =outpol->mols[t]->z;
mols[t]-> x     =outpol->mols[t]->x;
mols[t]-> y     =outpol->mols[t]->y; 
mols[t]-> seq   =outpol->mols[t]->seq;
mols[t]-> group   =outpol->mols[t]->group;
mols[t]-> bx    =outpol->mols[t]->bx;
mols[t]-> by    =outpol->mols[t]->by;
mols[t]-> bz    =outpol->mols[t]->bz;
mols[t]-> b2x   =outpol->mols[t]->b2x;
mols[t]-> b2y   =outpol->mols[t]->b2y;
mols[t]-> b2z   =outpol->mols[t]->b2z;
mols[t]-> id    =outpol->mols[t]->id;
}
for(int t=0; t<outpol->full_atoms.size(); t++){full_atoms[t]=outpol->full_atoms[t];}
randtot = outpol->randtot;
score = outpol->score;
printf("Final Structure: %8.3f Replaced: %d\n",outpol->score, outpol->randtot);

for(int m=0; m<pol2->mols.size(); m++){
delete pol2->mols[m];
}
pol2->full_atoms.clear();
pol2->mols.clear();
delete pol2;

for(int m=0; m<outpol->mols.size(); m++){
delete outpol->mols[m];
}
outpol->full_atoms.clear();
outpol->mols.clear();
delete outpol;
minscore=minscore/5;
return minscore;
} 

float Polymere::drimix(double rmsd_res,int start_center,int nmax, int movegroup, const std::vector< Fragment *>  frags,Fragclass * dat,const Polymere * pol,std::vector<Polymere *> * centers,std::vector< Fragatoms *> atoms,std::vector<Fragall *> batoms){
//Returns star_center if it is below the rmsd_res already, otherwise it returns -1 and the current polymere is the new center that should be added to the file
using namespace std;
int bout=1;
Polymere  * pol2 = new Polymere();
pol2->copy(centers->at(start_center));
int InitialClashTotal=pol2->checkBackBone();// + pol2->atomClashes();
Polymere  * outpol = new Polymere();
double minscore = 99999999;
score = 99999999;
vector<Polymere  *> last5;
int curr_size=1;
std::vector<Matches> mat; 
mat.reserve(7000);
double id=0;
//Read file and create a vector of centers
//printf("Total Centers Loaded: %d\n",centers->size());
//Doesn't read in the polymere that is currently being used 
//If it is not in the file set start_center < 0
//Calculate Current RMSD from all the centers
//If min rmsd is less than rmsd_res  
double minrmsd=99999;
for(int m=0; m<centers->size(); m++){
	if(m!=start_center){
		double nrmsd=rmsAllc(centers->at(m));	
		if(nrmsd < minrmsd){minrmsd=nrmsd;}
	}
}
if(minrmsd < 1){
//printf("Below Resolution: %8.3f\n",minrmsd);
return 0;
} 


//uniqpol is a vector of structures that results from a single closed moves search for a given structure
std::vector< std::vector<int> > groups;
int maxg=0;
//Read in the movement group from the ct file 
for(int g=0; g<num_mols; g++){
int gc = pol->mols[g]->group;
//cout << "group: " << gc << " bp1: " << pol->mols[g]->id << endl;
if(gc>maxg){
maxg=gc;
}
}
maxg++;
//cout << "Max Groups: " << maxg << endl;
for(int i=1; i < maxg; i++){
std::vector<int> group;
for(int g=0; g<num_mols; g++){
int gc = pol->mols[g]->group;
if(gc==i){
group.push_back(g);
}
}
groups.push_back(group);
}

int sig = groups[movegroup].size();

//Set counts to 0 for each position
for(int g=0; g<groups[movegroup].size(); g++){
mols[groups[movegroup][g]]->group=0;
}

int tot;
//Vectors containing the indices of valid movement groups
std::vector<int> groupOrder1;
std::vector<int> groupOrder2;
std::vector<int> groupOrder0;

for(int i=1; i < sig-2; i++){
groupOrder1.push_back(i);
groupOrder2.push_back(i);
}
int nt=(sig-3)*(sig-3);
for(int i=0; i <nt; i++){groupOrder0.push_back(i);} 
for(int iteration=0; iteration < nmax; iteration++){
//cout << "Iteration: " << iteration << endl;
////////////////////////////////////////////////////////////
// Move current polymere, score, add improvements to uniq //
////////////////////////////////////////////////////////////
//Generate Random order for movements
int n = groupOrder1.size();
int k;
int temp;
while(n>1){
k=(rand() % n);
--n;
temp=groupOrder1[n];
groupOrder1[n]=groupOrder1[k];
groupOrder1[k]=temp;
}

n=groupOrder2.size();
while(n>1){
k=(rand() % n);
--n;
temp=groupOrder2[n];
groupOrder2[n]=groupOrder2[k];
groupOrder2[k]=temp;
}
n=groupOrder1.size()*groupOrder2.size();
while(n>1){
k=(rand() % n);
--n;
temp=groupOrder0[n];
groupOrder0[n]=groupOrder0[k];
groupOrder0[k]=temp;
}
//cout << "Group1: ";
for(int i=0; i< groupOrder1.size(); i++){
//cout << groupOrder1[i] << " ";
}
//cout << endl;
int maxout=(sig-3)*(sig-3);

//Randomly go through pairs of group members
bout=1;
for(int tout=0; tout < maxout; tout++){
int i=groupOrder0[tout] % (sig-3);
int j= (groupOrder0[tout]-i)/(sig-3);


if((groupOrder1[i]+3) < groupOrder2[j]){
//Find center of each group movement and save to branches
mat=chooseMoves(groups[movegroup][groupOrder1[i]],groups[movegroup][groupOrder2[j]],dat,-1);
tot=mat.size();
if(tot > 500){
printf("Identified Moves: %d\n",tot);
tot=500;
}
//Identify the center of cluster of moves  
for(int w=0; w < tot; w++){
id=0;
if(batoms[mat[w].frag1]->check()==1 && batoms[mat[w].frag2]->check()==1){
//cout << "Creating new structure " << endl;
//Make a copy of the molecule array
//Use Group to track randomization
pol2->num_mols=num_mols;
pol2->randtot=0;
for(int t=0; t< num_mols; t++){

pol2->mols[t]-> z     =mols[t]->z;
pol2->mols[t]-> x     =mols[t]->x;
pol2->mols[t]-> y     =mols[t]->y; 
pol2->mols[t]-> seq   =mols[t]->seq;
//Adjust the probability for the movement group
if(t==groups[movegroup][groupOrder1[i]] || t==groups[movegroup][groupOrder2[j]]){
pol2->mols[t]-> group   = mols[t]->group - 1;
}else{pol2->mols[t]-> group   =mols[t]->group;}
pol2->mols[t]-> bx    =mols[t]->bx;
pol2->mols[t]-> by    =mols[t]->by;
pol2->mols[t]-> bz    =mols[t]->bz;
pol2->mols[t]-> b2x   =mols[t]->b2x;
pol2->mols[t]-> b2y   =mols[t]->b2y;
pol2->mols[t]-> b2z   =mols[t]->b2z;
pol2->mols[t]-> id    =mols[t]->id;
if(pol2->mols[t]->group < 0){pol2->randtot=pol2->randtot++;}
}
for(int t=0; t<full_atoms.size(); t++){pol2->full_atoms[t]=full_atoms[t];}


pol2->changeFrag2(mat[w].pos1,mat[w].frag1,frags);
pol2->changeFrag2(mat[w].pos2,mat[w].frag2,frags);
pol2->coordSysNt(mat[w].pos1,&atoms[mat[w].frag1]->atoms);
pol2->coordSysNt(mat[w].pos2,&atoms[mat[w].frag2]->atoms);
pol2->updateFull();
pol2->coordSysAT(mat[w].pos1,&batoms[mat[w].frag1]->atoms);
pol2->coordSysAT(mat[w].pos2,&batoms[mat[w].frag2]->atoms);

int c1dist1=pol2->c1dist();//check c1' clashes
int currentClashes=99;
if(c1dist1!=0){
currentClashes=pol2->checkBackBone(); // pol2->atomClashes();//check backbone bond distances
}
//Check to make sure it is different than the last 5 moves
double ormsd=999999;
//for(int j=0; j<curr_size; j++){
double ormsd1=rmsAllc(pol2);
if(ormsd1<ormsd){ormsd=ormsd1;}
//}

if(ormsd > 0.01){
if(c1dist1!=0 && currentClashes <= InitialClashTotal){
InitialClashTotal=currentClashes;


//printf("Backbones: %d Atom Clashes: %d\n",backbone,atomclash);
id=pol2->checkhbonds();
}else{id=0;}
//Determine minimal RMSD from any 
if(id>0){
double minrmsd12=9999999;
for(int c=0; c<centers->size(); c++){
double nrmsd12=pol2->rmsAllc(centers->at(c));
if(nrmsd12 < minrmsd12){minrmsd12=nrmsd12;}
}
pol2->score = -5*minrmsd12;
//printf("Score: %8.3f: %d\n",pol2->score,pol2->randtot);
if(pol2->score <= score){
id=1;
}else{
id=0;
double ratio=exp(-1*(pol2->score))/exp(-1*score);
float rand1 = rand()/((double)RAND_MAX);
if(rand1 < ratio){
id=1;
}
}
}
}else{
id=0;
}
//Make a copy of pol2 and save
if(id!=0){
printf("minscore: %8.3f: %d \n",pol2->score,currentClashes);
if(pol2->score < minscore){
minscore = pol2->score;
//Copy to outpol
//printf("Min Score: %8.3f\n", minscore);
for(int w=0; w<pol2->num_mols; w++){
outpol->mols[w]-> z     =pol2->mols[w]->z;
outpol->mols[w]-> x     =pol2->mols[w]->x;
outpol->mols[w]-> y     =pol2->mols[w]->y; 
outpol->mols[w]-> group =pol2->mols[w]->group;
outpol->mols[w]-> seq   =pol2->mols[w]->seq;
outpol->mols[w]-> bx     =pol2->mols[w]-> bx  ;
outpol->mols[w]-> by     =pol2->mols[w]-> by  ;
outpol->mols[w]-> bz    =pol2->mols[w]-> bz  ;
outpol->mols[w]-> b2x =pol2->mols[w]-> b2x ;
outpol->mols[w]-> b2y   =pol2->mols[w]-> b2y ;
outpol->mols[w]-> b2z   =pol2->mols[w]-> b2z ;
outpol->mols[w]-> id   =pol2->mols[w]-> id  ;
}
outpol->num_mols=pol2->num_mols;
outpol->randtot=pol2->randtot;
outpol->score = pol2->score;
for(int t=0; t<pol2->full_atoms.size(); t++){outpol->full_atoms[t]=pol2->full_atoms[t];}
}

//Add to last five



for(int w=0; w<pol2->num_mols; w++){
mols[w]-> z     =pol2->mols[w]->z;
mols[w]-> x     =pol2->mols[w]->x;
mols[w]-> y     =pol2->mols[w]->y; 
mols[w]-> group =pol2->mols[w]->group;
mols[w]-> seq   =pol2->mols[w]->seq;
mols[w]-> bx     =pol2->mols[w]-> bx  ;
mols[w]-> by     =pol2->mols[w]-> by  ;
mols[w]-> bz    =pol2->mols[w]-> bz  ;
mols[w]-> b2x =pol2->mols[w]-> b2x ;
mols[w]-> b2y   =pol2->mols[w]-> b2y ;
mols[w]-> b2z   =pol2->mols[w]-> b2z ;
mols[w]-> id   =pol2->mols[w]-> id  ;
}
num_mols=pol2->num_mols;
randtot=pol2->randtot;
score = pol2->score;
for(int t=0; t<pol2->full_atoms.size(); t++){full_atoms[t]=pol2->full_atoms[t];}
bout=0;
if(score < (-5*rmsd_res*1.5)){bout=5;}
break;
}

}//Loop through all matches
}

}//If groups are 3 apart
if(bout==0 || bout==5){break;}
}//Loop through random group pairs
if(bout==1 || bout==5){
//No more moves are possible 
break;
}

//if(minscore < (-5*(rmsd_res*2))){break;} 
}//Iteration

for(int t=0; t< num_mols; t++){
mols[t]-> z     =outpol->mols[t]->z;
mols[t]-> x     =outpol->mols[t]->x;
mols[t]-> y     =outpol->mols[t]->y; 
mols[t]-> seq   =outpol->mols[t]->seq;
mols[t]-> group   =outpol->mols[t]->group;
mols[t]-> bx    =outpol->mols[t]->bx;
mols[t]-> by    =outpol->mols[t]->by;
mols[t]-> bz    =outpol->mols[t]->bz;
mols[t]-> b2x   =outpol->mols[t]->b2x;
mols[t]-> b2y   =outpol->mols[t]->b2y;
mols[t]-> b2z   =outpol->mols[t]->b2z;
mols[t]-> id    =outpol->mols[t]->id;
}
for(int t=0; t<outpol->full_atoms.size(); t++){full_atoms[t]=outpol->full_atoms[t];}
outpol->randtot=(this->checkBackBone()+this->atomClashes());
outpol->score  =minscore;
randtot = outpol->randtot;
score = outpol->score;
printf("Final Structure: %8.3f Clashes: %d",outpol->score, outpol->randtot);

for(int m=0; m<pol2->mols.size(); m++){
delete pol2->mols[m];
}
pol2->full_atoms.clear();
pol2->mols.clear();
delete pol2;

for(int m=0; m<outpol->mols.size(); m++){
delete outpol->mols[m];
}
outpol->full_atoms.clear();
outpol->mols.clear();
delete outpol;
minscore=minscore/5;
float output_score;
if(minscore < -1){
output_score=minscore;
}else{
output_score = 0;
}
return output_score; 
} 

float Polymere::trimix(double rmsd_res,int start_center,int nmax, int movegroup, const std::vector< Fragment *>  frags,Fragclass * dat,const Polymere * pol,std::vector<Polymere *> * centers,std::vector< Fragatoms *> atoms,std::vector<Fragall *> batoms){
//Returns star_center if it is below the rmsd_res already, otherwise it returns -1 and the current polymere is the new center that should be added to the file
using namespace std;
int bout=1;
Polymere  * pol2 = new Polymere();
pol2->copy(centers->at(start_center));
int prevbone=pol2->checkBackBone();
Polymere  * outpol = new Polymere();
double minscore = 99999999;
score = 99999999;
vector<Polymere  *> last5;
int curr_size=1;
std::vector<Matches> mat; 
mat.reserve(7000);
double id=0;
//Read file and create a vector of centers
//printf("Total Centers Loaded: %d\n",centers->size());
//Doesn't read in the polymere that is currently being used 
//If it is not in the file set start_center < 0
//Calculate Current RMSD from all the centers
//If min rmsd is less than rmsd_res  
double minrmsd=99999;
for(int m=0; m<centers->size(); m++){
	if(m!=start_center){
		double nrmsd=rmsAllc(centers->at(m));	
		if(nrmsd < minrmsd){minrmsd=nrmsd;}
	}
}
if(minrmsd < 1){
//printf("Below Resolution: %8.3f\n",minrmsd);
return 0;
} 


//uniqpol is a vector of structures that results from a single closed moves search for a given structure
std::vector< std::vector<int> > groups;
int maxg=0;
//Read in the movement group from the ct file 
for(int g=0; g<num_mols; g++){
int gc = pol->mols[g]->group;
//cout << "group: " << gc << " bp1: " << pol->mols[g]->id << endl;
if(gc>maxg){
maxg=gc;
}
}
maxg++;
//cout << "Max Groups: " << maxg << endl;
for(int i=1; i < maxg; i++){
std::vector<int> group;
for(int g=0; g<num_mols; g++){
int gc = pol->mols[g]->group;
if(gc==i){
group.push_back(g);
}
}
groups.push_back(group);
}

int sig = groups[movegroup].size();

//Set counts to 0 for each position
for(int g=0; g<groups[movegroup].size(); g++){
mols[groups[movegroup][g]]->group=0;
}

int tot;
//Vectors containing the indices of valid movement groups
std::vector<int> groupOrder1;
std::vector<int> groupOrder2;
std::vector<int> groupOrder0;

for(int i=1; i < sig-2; i++){
groupOrder1.push_back(i);
groupOrder2.push_back(i);
}
int nt=(sig-3)*(sig-3);
for(int i=0; i <nt; i++){groupOrder0.push_back(i);} 

for(int iteration=0; iteration < nmax; iteration++){
//cout << "Iteration: " << iteration << endl;
////////////////////////////////////////////////////////////
// Move current polymere, score, add improvements to uniq //
////////////////////////////////////////////////////////////
//Generate Random order for movements
int n = groupOrder1.size();
int k;
int temp;
while(n>1){
k=(rand() % n);
--n;
temp=groupOrder1[n];
groupOrder1[n]=groupOrder1[k];
groupOrder1[k]=temp;
}

n=groupOrder2.size();
while(n>1){
k=(rand() % n);
--n;
temp=groupOrder2[n];
groupOrder2[n]=groupOrder2[k];
groupOrder2[k]=temp;
}
n=groupOrder1.size()*groupOrder2.size();
while(n>1){
k=(rand() % n);
--n;
temp=groupOrder0[n];
groupOrder0[n]=groupOrder0[k];
groupOrder0[k]=temp;
}
//cout << "Group1: ";
for(int i=0; i< groupOrder1.size(); i++){
//cout << groupOrder1[i] << " ";
}
//cout << endl;
int maxout=(sig-3)*(sig-3);

//Randomly go through pairs of group members
bout=1;
for(int tout=0; tout < maxout; tout++){
int i=groupOrder0[tout] % (sig-3);
int j= (groupOrder0[tout]-i)/(sig-3);


if((groupOrder1[i]+3) < groupOrder2[j]){
//Find center of each group movement and save to branches
mat=chooseMoves(groups[movegroup][groupOrder1[i]],groups[movegroup][groupOrder2[j]],dat,-1);
tot=mat.size();
if(tot > 500){
printf("Identified Moves: %d\n",tot);
tot=500;
}
//Identify the center of cluster of moves  
for(int w=0; w < tot; w++){
id=0;
if(batoms[mat[w].frag1]->check()==1 && batoms[mat[w].frag2]->check()==1){
//cout << "Creating new structure " << endl;
//Make a copy of the molecule array
//Use Group to track randomization
pol2->num_mols=num_mols;
pol2->randtot=0;
for(int t=0; t< num_mols; t++){

pol2->mols[t]-> z     =mols[t]->z;
pol2->mols[t]-> x     =mols[t]->x;
pol2->mols[t]-> y     =mols[t]->y; 
pol2->mols[t]-> seq   =mols[t]->seq;
//Adjust the probability for the movement group
if(t==groups[movegroup][groupOrder1[i]] || t==groups[movegroup][groupOrder2[j]]){
pol2->mols[t]-> group   = mols[t]->group - 1;
}else{pol2->mols[t]-> group   =mols[t]->group;}
pol2->mols[t]-> bx    =mols[t]->bx;
pol2->mols[t]-> by    =mols[t]->by;
pol2->mols[t]-> bz    =mols[t]->bz;
pol2->mols[t]-> b2x   =mols[t]->b2x;
pol2->mols[t]-> b2y   =mols[t]->b2y;
pol2->mols[t]-> b2z   =mols[t]->b2z;
pol2->mols[t]-> id    =mols[t]->id;
if(pol2->mols[t]->group < 0){pol2->randtot=pol2->randtot++;}
}
for(int t=0; t<full_atoms.size(); t++){pol2->full_atoms[t]=full_atoms[t];}


pol2->changeFrag2(mat[w].pos1,mat[w].frag1,frags);
pol2->changeFrag2(mat[w].pos2,mat[w].frag2,frags);
pol2->coordSysNt(mat[w].pos1,&atoms[mat[w].frag1]->atoms);
pol2->coordSysNt(mat[w].pos2,&atoms[mat[w].frag2]->atoms);
pol2->updateFull();
pol2->coordSysAT(mat[w].pos1,&batoms[mat[w].frag1]->atoms);
pol2->coordSysAT(mat[w].pos2,&batoms[mat[w].frag2]->atoms);

int c1dist1=pol2->c1dist();//check c1' clashes
int backbone=99;
if(c1dist1!=0){
backbone=pol2->checkBackBone();//check backbone bond distances
}
//Check to make sure it is different than the last 5 moves
double ormsd=999999;
//for(int j=0; j<curr_size; j++){
double ormsd1=rmsAllc(pol2);
if(ormsd1<ormsd){ormsd=ormsd1;}
//}

if(ormsd > 0.01){
if(c1dist1!=0 && backbone <= prevbone){
int atomclash=pol2->atomClashes();
//printf("Backbones: %d Atom Clashes: %d\n",backbone,atomclash);
id=pol2->checkhbonds();
}else{id=0;}
//Determine minimal RMSD from any 
if(id>0){
double minrmsd12=9999999;
for(int c=0; c<centers->size(); c++){
double nrmsd12=pol2->rmsAllc(centers->at(c));
if(nrmsd12 < minrmsd12){minrmsd12=nrmsd12;}
}
pol2->score = -5*minrmsd12;
//printf("Score: %8.3f: %d\n",pol2->score,pol2->randtot);
if(pol2->score <= score){
id=1;
}else{
id=0;
double ratio=exp(-1*(pol2->score))/exp(-1*score);
float rand1 = rand()/((double)RAND_MAX);
if(rand1 < ratio){
id=1;
}
}
}
}else{
id=0;
}
//Make a copy of pol2 and save
if(id!=0){
if(pol2->score < minscore){
prevbone=backbone;
printf("minscore: %8.3f: %d \n",pol2->score,prevbone);
minscore = pol2->score;
//Copy to outpol
//printf("Min Score: %8.3f\n", minscore);
for(int w=0; w<pol2->num_mols; w++){
outpol->mols[w]-> z     =pol2->mols[w]->z;
outpol->mols[w]-> x     =pol2->mols[w]->x;
outpol->mols[w]-> y     =pol2->mols[w]->y; 
outpol->mols[w]-> group =pol2->mols[w]->group;
outpol->mols[w]-> seq   =pol2->mols[w]->seq;
outpol->mols[w]-> bx     =pol2->mols[w]-> bx  ;
outpol->mols[w]-> by     =pol2->mols[w]-> by  ;
outpol->mols[w]-> bz    =pol2->mols[w]-> bz  ;
outpol->mols[w]-> b2x =pol2->mols[w]-> b2x ;
outpol->mols[w]-> b2y   =pol2->mols[w]-> b2y ;
outpol->mols[w]-> b2z   =pol2->mols[w]-> b2z ;
outpol->mols[w]-> id   =pol2->mols[w]-> id  ;
}
outpol->num_mols=pol2->num_mols;
outpol->randtot=pol2->randtot;
outpol->score = pol2->score;
for(int t=0; t<pol2->full_atoms.size(); t++){outpol->full_atoms[t]=pol2->full_atoms[t];}
}

//Add to last five



for(int w=0; w<pol2->num_mols; w++){
mols[w]-> z     =pol2->mols[w]->z;
mols[w]-> x     =pol2->mols[w]->x;
mols[w]-> y     =pol2->mols[w]->y; 
mols[w]-> group =pol2->mols[w]->group;
mols[w]-> seq   =pol2->mols[w]->seq;
mols[w]-> bx     =pol2->mols[w]-> bx  ;
mols[w]-> by     =pol2->mols[w]-> by  ;
mols[w]-> bz    =pol2->mols[w]-> bz  ;
mols[w]-> b2x =pol2->mols[w]-> b2x ;
mols[w]-> b2y   =pol2->mols[w]-> b2y ;
mols[w]-> b2z   =pol2->mols[w]-> b2z ;
mols[w]-> id   =pol2->mols[w]-> id  ;
}
num_mols=pol2->num_mols;
randtot=pol2->randtot;
score = pol2->score;
for(int t=0; t<pol2->full_atoms.size(); t++){full_atoms[t]=pol2->full_atoms[t];}
bout=0;
if(score < (-5*rmsd_res*1.5)){bout=5;}
break;
}

}//Loop through all matches
}

}//If groups are 3 apart
if(bout==0 || bout==5){break;}
}//Loop through random group pairs
if(bout==1 || bout==5){
//No more moves are possible 
break;
}

//if(minscore < (-5*(rmsd_res*2))){break;} 
}//Iteration

for(int t=0; t< num_mols; t++){
mols[t]-> z     =outpol->mols[t]->z;
mols[t]-> x     =outpol->mols[t]->x;
mols[t]-> y     =outpol->mols[t]->y; 
mols[t]-> seq   =outpol->mols[t]->seq;
mols[t]-> group   =outpol->mols[t]->group;
mols[t]-> bx    =outpol->mols[t]->bx;
mols[t]-> by    =outpol->mols[t]->by;
mols[t]-> bz    =outpol->mols[t]->bz;
mols[t]-> b2x   =outpol->mols[t]->b2x;
mols[t]-> b2y   =outpol->mols[t]->b2y;
mols[t]-> b2z   =outpol->mols[t]->b2z;
mols[t]-> id    =outpol->mols[t]->id;
}
for(int t=0; t<outpol->full_atoms.size(); t++){full_atoms[t]=outpol->full_atoms[t];}
randtot = outpol->randtot;
score = outpol->score;
printf("Final Structure: %8.3f Replaced: %d Backbone: %d\n",outpol->score, outpol->randtot,prevbone);

for(int m=0; m<pol2->mols.size(); m++){
delete pol2->mols[m];
}
pol2->full_atoms.clear();
pol2->mols.clear();
delete pol2;

for(int m=0; m<outpol->mols.size(); m++){
delete outpol->mols[m];
}
outpol->full_atoms.clear();
outpol->mols.clear();
delete outpol;
minscore=minscore/5;
return minscore; 
} 

void Polymere::mixit(int nmax, int movegroup, const std::vector< Fragment *>  frags,Fragclass * dat,const Polymere * pol,std::vector< Fragatoms *> atoms,std::vector<Fragall *> batoms){
using namespace std;
Polymere  * pol2 = new Polymere();
Polymere  * outpol = new Polymere();
Polymere  * backout = new Polymere();
int bout=2;
//interations = Maximum number of iterations to perform
//max_rand    = Maximum amount of randomization 
//frags       = Fragment database
//dat         = Indexed fragments
//group       = Movement group to randomize
//pol         = This is the master polymere that everything is compared to
std::vector<Matches> mat; 
mat.reserve(7000);
double id=0;
double pid=0;
//uniqpol is a vector of structures that results from a single closed moves search for a given structure
std::vector<Polymere *> uniqpol;
uniqpol.reserve(1000);
std::vector<Polymere *> branches;
branches.reserve(20);
std::vector< std::vector<int> > groups;
int maxg=0;
//Read in the movement group from the ct file 
for(int g=0; g<num_mols; g++){
int gc = pol->mols[g]->group;
cout << "group: " << gc << " bp1: " << pol->mols[g]->id << endl;
if(gc>maxg){
maxg=gc;
}
}
maxg++;
//cout << "Max Groups: " << maxg << endl;
for(int i=1; i < maxg; i++){
std::vector<int> group;
for(int g=0; g<num_mols; g++){
int gc = pol->mols[g]->group;
if(gc==i){
group.push_back(g);
}
}
groups.push_back(group);
}

int sig = groups[movegroup].size();
//cout << "Group Size: " << movegroup << ":" << sig << endl;
//cout << "Total Movement Groups: " << groups.size() << endl;
for(int g=0; g<groups.size(); g++){
//cout << "Group: " << (g+1) << endl;
for(int g1=0; g1 < groups[g].size(); g1++){
//cout << "    pos: " << groups[g][g1] << endl;
}
}

//Set counts to 0 for each position
for(int g=0; g<groups[movegroup].size(); g++){
mols[groups[movegroup][g]]->group=0;
}

int tot;

//Vectors containing the indices of valid movement groups
std::vector<int> groupOrder1;
std::vector<int> groupOrder2;
std::vector<int> groupOrder0;

for(int i=1; i < sig-2; i++){
groupOrder1.push_back(i);
groupOrder2.push_back(i);
}
int nt=(sig-3)*(sig-3);
for(int i=0; i <nt; i++){groupOrder0.push_back(i);} 
//Set initial Score and RMSD threshold
//cout << "Starting randomization: " << nmax << ":" << sig <<  endl;
int rand_threshold=0;
//Foreach starting structure score them based on some criteria for each move.  Need to loop through 
//a certain number of cycles or until the simulation converges
//Add loop to perform multiple generations of movements
printf("Nmax: %d\n",nmax);
for(int iteration=0; iteration < nmax; iteration++){
cout << "Iteration: " << iteration << endl;
////////////////////////////////////////////////////////////
// Move current polymere, score, add improvements to uniq //
////////////////////////////////////////////////////////////
//Generate Random order for movements
int n = groupOrder1.size();
int k;
int temp;
while(n>1){
k=(rand() % n);
--n;
temp=groupOrder1[n];
groupOrder1[n]=groupOrder1[k];
groupOrder1[k]=temp;
}

n=groupOrder2.size();
while(n>1){
k=(rand() % n);
--n;
temp=groupOrder2[n];
groupOrder2[n]=groupOrder2[k];
groupOrder2[k]=temp;
}
n=groupOrder1.size()*groupOrder2.size();
while(n>1){
k=(rand() % n);
--n;
temp=groupOrder0[n];
groupOrder0[n]=groupOrder0[k];
groupOrder0[k]=temp;
}
//cout << "Group1: ";
for(int i=0; i< groupOrder1.size(); i++){
//cout << groupOrder1[i] << " ";
}
//cout << endl;
int maxout=(sig-3)*(sig-3);
int breakout = 0;
int uthres=5;
if(iteration < 20){
uthres=1;
}
for(int tout=0; tout < maxout; tout++){
int i=groupOrder0[tout] % (sig-3);
int j= (groupOrder0[tout]-i)/(sig-3);

//Clear Uniqpols From the previous Run
for(int u=0; u<uniqpol.size(); u++){
for(int m=0; m<uniqpol[u]->mols.size(); m++){
delete uniqpol[u]->mols[m];
}
uniqpol[u]->mols.clear();
uniqpol[u]->full_atoms.clear();
delete uniqpol[u];
}
uniqpol.clear();
//Check to make sure the groups are in focus
//Decide if the group should be moved
float prob=rand()/(float)(RAND_MAX);
float totalprob = (float)(-1*(mols[groups[movegroup][groupOrder1[i]]]->group + mols[ groups[movegroup][groupOrder2[j]] ]->group))/((float)(2*iteration+groups[movegroup].size()));
totalprob=totalprob*(float)groups[movegroup].size()/10;
//cout << "Total prob: " << totalprob << " Rand: " << prob << " : " << groups[movegroup][groupOrder1[i]] << " : " << groups[movegroup][groupOrder2[j]] <<  endl;

if((groupOrder1[i]+3) < groupOrder2[j] && totalprob < prob){
//Loop through possible position combinations
//Only use positions that aren't set to -1 
//If position is less than -1 that is the cluster number to use
//As soon as move is identified make the move
//Find center of each group movement and save to branches
mat=chooseMoves(groups[movegroup][groupOrder1[i]],groups[movegroup][groupOrder2[j]],dat,-1);
tot=mat.size();
//printf("Identified Moves: %d\n",tot);
//Identify the center of cluster of moves  
double id=0;
double pid=0;
for(int w=0; w < tot; w++){
//cout << "Creating new structure " << endl;
//Make a copy of the molecule array
//Use Group to track randomization
int nmol=pol->num_mols;
pol2->num_mols=nmol;
pol2->randtot=0;
for(int t=0; t< nmol; t++){
pol2->mols[t]-> z     =mols[t]->z;
pol2->mols[t]-> x     =mols[t]->x;
pol2->mols[t]-> y     =mols[t]->y; 
pol2->mols[t]-> seq   =mols[t]->seq;
//Adjust the probability for the movement group
if(t==groups[movegroup][groupOrder1[i]] || t==groups[movegroup][groupOrder2[j]]){
pol2->mols[t]-> group   = mols[t]->group - 1;
}else{pol2->mols[t]-> group   =mols[t]->group;}
pol2->mols[t]-> bx    =mols[t]->bx;
pol2->mols[t]-> by    =mols[t]->by;
pol2->mols[t]-> bz    =mols[t]->bz;
pol2->mols[t]-> b2x   =mols[t]->b2x;
pol2->mols[t]-> b2y   =mols[t]->b2y;
pol2->mols[t]-> b2z   =mols[t]->b2z;
pol2->mols[t]-> id    =mols[t]->id;
if(pol2->mols[t]->group < 0){pol2->randtot=pol2->randtot++;}
}
for(int t=0; t<full_atoms.size(); t++){pol2->full_atoms[t]=full_atoms[t];}

pol2->changeFrag2(mat[w].pos1,mat[w].frag1,frags);
pol2->changeFrag2(mat[w].pos2,mat[w].frag2,frags);
pol2->coordSysNt(mat[w].pos1,&atoms[mat[w].frag1]->atoms);
pol2->coordSysNt(mat[w].pos2,&atoms[mat[w].frag2]->atoms);
pol2->updateFull();
pol2->coordSysAT(mat[w].pos1,&batoms[mat[w].frag1]->atoms);
pol2->coordSysAT(mat[w].pos2,&batoms[mat[w].frag2]->atoms);

/*
printf("c1 batoms: %8.3f %8.3f %8.3f\n",batoms[mat[w].frag2]->atoms[11].x,batoms[mat[w].frag2]->atoms[23].x,batoms[mat[w].frag2]->atoms[35].x); 
printf("c1 atoms: %8.3f %8.3f %8.3f\n",atoms[mat[w].frag2]->atoms[0].x,atoms[mat[w].frag2]->atoms[3].x,atoms[mat[w].frag2]->atoms[6].x ); 
*/




//Check base-pair differences
int c1dist1=pol2->c1dist();
if(c1dist1!=0){id=pol2->checkhbonds();}else{id=0;
//cout << "Failed C1dist: " << c1dist1 << endl;
}
//Make a copy of pol2 and save
if(id!=0){
//cout << "Matches: " << mat[w].pos1 << ":" << mat[w].pos2 << endl;
uniqpol.push_back(new Polymere());
for(int w=0; w<pol2->num_mols; w++){
uniqpol.back()->mols[w]-> z     =pol2->mols[w]->z;
uniqpol.back()->mols[w]-> x     =pol2->mols[w]->x;
uniqpol.back()->mols[w]-> y     =pol2->mols[w]->y; 
uniqpol.back()->mols[w]-> group =pol2->mols[w]->group;
uniqpol.back()->mols[w]-> seq   =pol2->mols[w]->seq;
uniqpol.back()->mols[w]-> bx     =pol2->mols[w]-> bx  ;
uniqpol.back()->mols[w]-> by     =pol2->mols[w]-> by  ;
uniqpol.back()->mols[w]-> bz    =pol2->mols[w]-> bz  ;
uniqpol.back()->mols[w]-> b2x =pol2->mols[w]-> b2x ;
uniqpol.back()->mols[w]-> b2y   =pol2->mols[w]-> b2y ;
uniqpol.back()->mols[w]-> b2z   =pol2->mols[w]-> b2z ;
uniqpol.back()->mols[w]-> id   =pol2->mols[w]-> id  ;
}
uniqpol.back()->num_mols=pol2->num_mols;
uniqpol.back()->randtot=pol2->randtot;
uniqpol.back()->score = pol2->score;
for(int t=0; t<pol2->full_atoms.size(); t++){uniqpol.back()->full_atoms[t]=pol2->full_atoms[t];}
}

}
if(uniqpol.size()>uthres){
//cout << "Uniqpol Size: " << uniqpol.size();
////////////////////////////////////////////////////////////
//Identify the center of structures in uniqpol
////////////////////////////////////////////////////////////
//Find structure with minimal average rmsd from all others
int center_ind;
float minrmsd=99999;
if(uniqpol.size() > 2){
center_ind=(rand() % (uniqpol.size()-1))+1;
/*
//Cluster and take the center structure
for(int i=0; i<uniqpol.size(); i++){
float avermsd=0;
for(int j=0; j<uniqpol.size(); j++){
float rmsd_diff = uniqpol[i]->rmsAll(uniqpol[j]);
avermsd=avermsd+rmsd_diff;
}
avermsd=avermsd/uniqpol.size();

if(avermsd<minrmsd){
center_ind=i;
minrmsd=avermsd;
}

}
*/

//Pick random structure
//Save structure to outPol
backout->randtot=uniqpol[0]->randtot;
backout->num_mols=num_mols;
for(int t=0; t< num_mols; t++){
backout->mols[t]-> z     =uniqpol[0]->mols[t]->z;
backout->mols[t]-> x     =uniqpol[0]->mols[t]->x;
backout->mols[t]-> y     =uniqpol[0]->mols[t]->y; 
backout->mols[t]-> seq   =uniqpol[0]->mols[t]->seq;
backout->mols[t]-> group   =uniqpol[0]->mols[t]->group;
backout->mols[t]-> bx    =uniqpol[0]->mols[t]->bx;
backout->mols[t]-> by    =uniqpol[0]->mols[t]->by;
backout->mols[t]-> bz    =uniqpol[0]->mols[t]->bz;
backout->mols[t]-> b2x   =uniqpol[0]->mols[t]->b2x;
backout->mols[t]-> b2y   =uniqpol[0]->mols[t]->b2y;
backout->mols[t]-> b2z   =uniqpol[0]->mols[t]->b2z;
backout->mols[t]-> id    =uniqpol[0]->mols[t]->id;
}
backout->randtot = uniqpol[0]->randtot;
backout->score = uniqpol[0]->score;
for(int t=0; t<uniqpol[0]->full_atoms.size(); t++){backout->full_atoms[t]=uniqpol[0]->full_atoms[t];}

}else{
//If only two of them exist take the first one
center_ind=0;
}
minrmsd=uniqpol[center_ind]->rmsAllc(pol);
cout << "Thres: " << uthres << " Moved " << tot << ": " << center_ind << " RMSD: " << minrmsd << " Randomization: " << uniqpol[center_ind]->randtot << endl;
//Pdb * pdbn = new Pdb(uniqpol[center_ind],'j');
//pdbn->print();
//Set mols to uniqpol[center_ind]
if(uniqpol[center_ind]->randtot >= rand_threshold){
//Save structure to outPol
rand_threshold=uniqpol[center_ind]->randtot;
outpol->num_mols=num_mols;
for(int t=0; t< num_mols; t++){
outpol->mols[t]-> z     =uniqpol[center_ind]->mols[t]->z;
outpol->mols[t]-> x     =uniqpol[center_ind]->mols[t]->x;
outpol->mols[t]-> y     =uniqpol[center_ind]->mols[t]->y; 
outpol->mols[t]-> seq   =uniqpol[center_ind]->mols[t]->seq;
outpol->mols[t]-> group   =uniqpol[center_ind]->mols[t]->group;
outpol->mols[t]-> bx    =uniqpol[center_ind]->mols[t]->bx;
outpol->mols[t]-> by    =uniqpol[center_ind]->mols[t]->by;
outpol->mols[t]-> bz    =uniqpol[center_ind]->mols[t]->bz;
outpol->mols[t]-> b2x   =uniqpol[center_ind]->mols[t]->b2x;
outpol->mols[t]-> b2y   =uniqpol[center_ind]->mols[t]->b2y;
outpol->mols[t]-> b2z   =uniqpol[center_ind]->mols[t]->b2z;
outpol->mols[t]-> id    =uniqpol[center_ind]->mols[t]->id;
}
outpol->randtot = uniqpol[center_ind]->randtot;
outpol->score = uniqpol[center_ind]->score;
for(int t=0; t<uniqpol[center_ind]->full_atoms.size(); t++){outpol->full_atoms[t]=uniqpol[center_ind]->full_atoms[t];}
}

for(int t=0; t< num_mols; t++){
mols[t]-> z     =uniqpol[center_ind]->mols[t]->z;
mols[t]-> x     =uniqpol[center_ind]->mols[t]->x;
mols[t]-> y     =uniqpol[center_ind]->mols[t]->y; 
mols[t]-> seq   =uniqpol[center_ind]->mols[t]->seq;
mols[t]-> group   =uniqpol[center_ind]->mols[t]->group;
mols[t]-> bx    =uniqpol[center_ind]->mols[t]->bx;
mols[t]-> by    =uniqpol[center_ind]->mols[t]->by;
mols[t]-> bz    =uniqpol[center_ind]->mols[t]->bz;
mols[t]-> b2x   =uniqpol[center_ind]->mols[t]->b2x;
mols[t]-> b2y   =uniqpol[center_ind]->mols[t]->b2y;
mols[t]-> b2z   =uniqpol[center_ind]->mols[t]->b2z;
mols[t]-> id    =uniqpol[center_ind]->mols[t]->id;
}
randtot = uniqpol[center_ind]->randtot;
score = uniqpol[center_ind]->score;


for(int t=0; t<uniqpol[center_ind]->full_atoms.size(); t++){full_atoms[t]=uniqpol[center_ind]->full_atoms[t];}
//cout << "\n\n" << endl;
//printf("Full Atom Size: %d\n", uniqpol[center_ind]->full_atoms.size());
//cout << "\n\n" << endl;


//Break out of loop
//Pdb * pdb = new Pdb(uniqpol[center_ind],'A');
//pdb->print();
breakout=1;
bout=0;
uthres=5;
break;
}else{
if(tout == (maxout-1) && bout==0){
bout=1;
uthres=3;
cout << "Back Step" << endl;
for(int t=0; t< num_mols; t++){
mols[t]-> z     =backout->mols[t]->z;
mols[t]-> x     =backout->mols[t]->x;
mols[t]-> y     =backout->mols[t]->y; 
mols[t]-> seq   =backout->mols[t]->seq;
mols[t]-> group   =backout->mols[t]->group;
mols[t]-> bx    =backout->mols[t]->bx;
mols[t]-> by    =backout->mols[t]->by;
mols[t]-> bz    =backout->mols[t]->bz;
mols[t]-> b2x   =backout->mols[t]->b2x;
mols[t]-> b2y   =backout->mols[t]->b2y;
mols[t]-> b2z   =backout->mols[t]->b2z;
mols[t]-> id    =backout->mols[t]->id;
}
randtot = backout->randtot;
score = backout->score;
for(int t=0; t<backout->full_atoms.size(); t++){full_atoms[t]=backout->full_atoms[t];}
}else{
if(tout == maxout-1){
bout=2;
}
}


}


}

}
if(bout==2){break;}
}

//Changes Values of current polymere to outpol
/*
for(int t=0; t< num_mols; t++){
mols[t]-> z     =outpol->mols[t]->z;
mols[t]-> x     =outpol->mols[t]->x;
mols[t]-> y     =outpol->mols[t]->y; 
mols[t]-> seq   =outpol->mols[t]->seq;
mols[t]-> group   =outpol->mols[t]->group;
mols[t]-> bx    =outpol->mols[t]->bx;
mols[t]-> by    =outpol->mols[t]->by;
mols[t]-> bz    =outpol->mols[t]->bz;
mols[t]-> b2x   =outpol->mols[t]->b2x;
mols[t]-> b2y   =outpol->mols[t]->b2y;
mols[t]-> b2z   =outpol->mols[t]->b2z;
mols[t]-> id    =outpol->mols[t]->id;
}
randtot = outpol->randtot;
score = outpol->score;
*/
for(int m=0; m<pol2->mols.size(); m++){
delete pol2->mols[m];
}
pol2->full_atoms.clear();
pol2->mols.clear();
//Remove outpol
for(int m=0; m<outpol->mols.size(); m++){
delete outpol->mols[m];
}
outpol->full_atoms.clear();
outpol->mols.clear();
for(int u=0; u<uniqpol.size(); u++){
for(int m=0; m<uniqpol[u]->mols.size(); m++){
delete uniqpol[u]->mols[m];
}
uniqpol[u]->full_atoms.clear();
uniqpol[u]->mols.clear();
}
uniqpol.clear();
delete outpol;
delete pol2;
//Clear uniqpol;
}



void Polymere::copyOver(Polymere * pol2,int into_pos,int from_pos,int totlen){
//Copies the region between from_pos and from_pos+totlen 
//into the current polymere starting at position into_pos

int tot = into_pos+totlen;
while(mols.size() < tot){mols.push_back(new Base());}
int tt=from_pos;
for(int t=into_pos; t<(into_pos+totlen); t++){
mols[t]-> z     =pol2->mols[tt]->z;
mols[t]-> x     =pol2->mols[tt]->x;
mols[t]-> y     =pol2->mols[tt]->y; 
mols[t]-> seq   =pol2->mols[tt]->seq;
mols[t]-> group =pol2->mols[tt]->group;
mols[t]-> bx    =pol2->mols[tt]->bx;
mols[t]-> by    =pol2->mols[tt]->by;
mols[t]-> bz    =pol2->mols[tt]->bz;
mols[t]-> b2x   =pol2->mols[tt]->b2x;
mols[t]-> b2y   =pol2->mols[tt]->b2y;
mols[t]-> b2z   =pol2->mols[tt]->b2z;
mols[t]-> id    =pol2->mols[tt]->id;
mols[t]-> prob    =pol2->mols[tt]->prob;
tt++;
}
while(full_atoms.size() < ((into_pos+totlen)*81)){full_atoms.push_back(0);}
tt=from_pos*81;
for(int t=(into_pos*81); t<(into_pos*81+totlen*81); t++){full_atoms[t]=pol2->full_atoms[tt++];}
}

void Polymere::copy(Polymere * pol2){
score = pol2->score; 
num_mols=pol2->num_mols;
randtot=pol2->randtot;
int tot = pol2->mols.size();
while(mols.size() < tot){mols.push_back(new Base());}
for(int t=0; t<pol2->mols.size(); t++){
mols[t]-> z     =pol2->mols[t]->z;
mols[t]-> x     =pol2->mols[t]->x;
mols[t]-> y     =pol2->mols[t]->y; 
mols[t]-> seq   =pol2->mols[t]->seq;
mols[t]-> group =pol2->mols[t]->group;
mols[t]-> bx    =pol2->mols[t]->bx;
mols[t]-> by    =pol2->mols[t]->by;
mols[t]-> bz    =pol2->mols[t]->bz;
mols[t]-> b2x   =pol2->mols[t]->b2x;
mols[t]-> b2y   =pol2->mols[t]->b2y;
mols[t]-> b2z   =pol2->mols[t]->b2z;
mols[t]-> id    =pol2->mols[t]->id;
mols[t]-> prob    =pol2->mols[t]->prob;
}
while(full_atoms.size() < pol2->full_atoms.size()){full_atoms.push_back(0);}
for(int t=0; t<pol2->full_atoms.size(); t++){full_atoms[t]=pol2->full_atoms[t];}
}
void writePDB(char * file){
//Create a new PDB from the full atoms vector

}

