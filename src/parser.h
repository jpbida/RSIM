/*
Expression parser
lang:   C++ code, project developed with wxDev-C++
author: Jos de Jong, 2007-12-22
site:   www.speqmath.com

Features:
    Operators:
        & | << >>
        = <> < > <= >=
        + -
        * / % ||
        ^
        !

    Functions:
        Abs, Exp, Sign, Sqrt, Log, Log10
        Sin, Cos, Tan, ASin, ACos, ATan
        Factorial

    Variables:
        Pi, e
        you can define your own variables
        
    Other:
        Scientific notation supported
        Error handling supported

*/


#ifndef PARSER_H
#define PARSER_H

// declarations
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <cmath>

#include <iostream>
#include <string>

#include "constants.h"
#include "error.h"
#include "score.h"
#include "functions.h"
#include "variablelist.h"

using namespace std;

class Parser
{
    // public functions
    public:
        Parser(const char * file);
        Parser(string file);
//Reads in an expression and evaluates it 
        char* parse(const char expr[]);
//Uses the score that was read in from the file and evaluates it 
        char* parse();
//Loops through all residues in the score and returns a score for each of them
        void score(Score * score, std::vector<double> * scorevector);
  int  s0;// score-s0[vind]=gorder[vind];
  int  s1;// score-s1[vind]=gtype1[vind];
  int  s2;// score-s2[vind]=gtype2[vind];
  int  s3;// score-s3[vind]=gtype3[vind];
  int  s4;// score-s4[vind]=gtype4[vind];
  int  s5;// score-s5[vind]=rtype[vind];
float  s6;//  score-s6[vind]=volumes[vind];
  int  s7;//  score-s7[vind]=unpacked[vind];
float  s8;//  score-s8[vind]=surface[vind];
float  s9;//  score-s9[vind]=bvolumes[vind];
  int s10;//  score-s10[vind]=bunpacked[vind];
float s11;//  score-s11[vind]=bsurface[vind];
float s12;//  score-s12[vind]=surfG[vind];
float s13;//  score-s13[vind]=surfA[vind];
float s14;//  score-s14[vind]=surfU[vind];
float s15;//  score-s15[vind]=surfC[vind] ;
float s16;//  score-s16[vind]=bsurfG[vind];
float s17;//  score-s17[vind]=bsurfA[vind];
float s18;//  score-s18[vind]=bsurfU[vind];
float s19;//  score-s19[vind]=bsurfC[vind];
   int s20;// score-s20[vind]=orderG[vind];
   int s21;// score-s21[vind]=orderA[vind];
   int s22;// score-s22[vind]=orderU[vind];
   int s23;// score-s23[vind]=orderC[vind];
   int s24;// score-s24[vind]=borderG[vind];
   int s25;// score-s25[vind]=borderA[vind];
   int s26;// score-s26[vind]=borderU[vind];
   int s27;// score-s27[vind]=borderC[vind];
float s28;//  score-s28[vind]=vatom0[vind];
float s29;//  score-s29[vind]=vatom1[vind];
float s30;//  score-s30[vind]=vatom2[vind];
float s31;//  score-s31[vind]=vatom3[vind];
float s32;//  score-s32[vind]=vatom4[vind];
float s33;//  score-s33[vind]=vatom5[vind];
float s34;//  score-s34[vind]=vatom6[vind];
float s35;//  score-s35[vind]=vatom7[vind];
float s36;//  score-s36[vind]=vatom8[vind];
float s37;//  score-s37[vind]=vatom9[vind];
float s38;//  score-s38[vind]=vatom10[vind];
float s39;//  score-s39[vind]=vatom11[vind];
float s40;//  score-s40[vind]=vatom12[vind];
float s41;//  score-s41[vind]=vatom13[vind];
float s42;//  score-s42[vind]=vatom14[vind];
float s43;//  score-s43[vind]=vatom15[vind];
float s44;//  score-s44[vind]=vatom16[vind];
float s45;//  score-s45[vind]=vatom17[vind];
float s46;//  score-s46[vind]=vatom18[vind];
float s47;//  score-s47[vind]=vatom19[vind];
float s48;//  score-s48[vind]=vatom20[vind];
float s49;//  score-s49[vind]=vatom21[vind];
float s50;//  score-s50[vind]=vatom22[vind];
float s51;//  score-s51[vind]=vatom23[vind];
float s52;//  score-s52[vind]=vatom24[vind];
float s53;//  score-s53[vind]=vatom25[vind];
float s54;//  score-s54[vind]=vatom26[vind];
  int s55;//  score-s55[vind]=atom0[vind];
  int s56;//  score-s56[vind]=atom1[vind];
  int s57;//  score-s57[vind]=atom2[vind];
  int s58;//  score-s58[vind]=atom3[vind];
  int s59;//  score-s59[vind]=atom4[vind];
  int s60;//  score-s60[vind]=atom5[vind];
  int s61;//  score-s61[vind]=atom6[vind];
  int s62;//  score-s62[vind]=atom7[vind];
  int s63;//  score-s63[vind]=atom8[vind];
  int s64;//  score-s64[vind]=atom9[vind];
  int s65;//  score-s65[vind]=atom10[vind];
  int s66;//  score-s66[vind]=atom11[vind];
  int s67;//  score-s67[vind]=atom12[vind];
  int s68;//  score-s68[vind]=atom13[vind];
  int s69;//  score-s69[vind]=atom14[vind];
  int s70;//  score-s70[vind]=atom15[vind];
  int s71;//  score-s71[vind]=atom16[vind];
  int s72;//  score-s72[vind]=atom17[vind];
  int s73;//  score-s73[vind]=atom18[vind];
  int s74;//  score-s74[vind]=atom19[vind];
  int s75;//  score-s75[vind]=atom20[vind];
  int s76;//  score-s76[vind]=atom21[vind];
  int s77;//  score-s77[vind]=atom22[vind];
  int s78;//  score-s78[vind]=atom23[vind];
int s79;//   score-s79[vind]=atom24[vind];
int s80;//   score-s80[vind]=atom25[vind];
int s81;//   score-s81[vind]=atom26[vind];
float s82;// score-s82[vind]=hydros[vind];
float s83;// score-s83[vind]=doublets[vind];

float s84;
float s85;
float s86;
float s87;
float s88;
float s89;



    // enumerations
    private:
    
        enum TOKENTYPE {NOTHING = -1, DELIMETER, NUMBER, VARIABLE, FUNCTION, UNKNOWN};
    
        enum OPERATOR_ID {AND, OR, BITSHIFTLEFT, BITSHIFTRIGHT,                 // level 2
                       EQUAL, UNEQUAL, SMALLER, LARGER, SMALLEREQ, LARGEREQ,    // level 3
                       PLUS, MINUS,                     // level 4
                       MULTIPLY, DIVIDE, MODULUS, XOR,  // level 5
                       POW,                             // level 6
                       FACTORIAL};                      // level 7

    // data
    private:
        char expr[EXPR_LEN_MAX+1];    // holds the expression
        char gexpr[EXPR_LEN_MAX+1];    // holds the expression
        char* e;                      // points to a character in expr
        
        char token[NAME_LEN_MAX+1];   // holds the token
        TOKENTYPE token_type;         // type of the token

        double ans;                   // holds the result of the expression
        char ans_str[255];            // holds a string containing the result 
                                      // of the expression

        Variablelist user_var;        // list with variables defined by user

    // private functions
    private:
        void getToken();
        
        double parse_level1();
        double parse_level2();
        double parse_level3();
        double parse_level4();
        double parse_level5();
        double parse_level6();
        double parse_level7();
        double parse_level8();
        double parse_level9();
        double parse_level10();
        double parse_number();

        int get_operator_id(const char op_name[]);
        double eval_operator(const int op_id, const double &lhs, const double &rhs);
        double eval_function(const char fn_name[], const double &value);
        double eval_variable(const char var_name[]);

        int row();
        int col();
};

#endif
