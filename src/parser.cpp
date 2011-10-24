// declarations
#include "parser.h"
#include <fstream>
#include "score.h"
#include <string>

Parser::Parser(string file)
{
using namespace std;
    expr[0] = '\0';
    gexpr[0] = '\0';
    e = NULL;
    token[0] = '\0';
    token_type = NOTHING;
std::ifstream input;
input.open(file.c_str());
if(input.fail()){
printf("Reading score file failed");
}else{
string line;
int index=0;
string prevchain ("");
string prevres ("");
int prevresnum= -1;
char new_expr[EXPR_LEN_MAX];
while(getline(input, line)) {
string front=line.substr(0,2);

if(front.compare("M:")==0){
int len=line.length()-2;
string model=line.substr(2,len);
strcpy(expr,model.c_str());
//printf("Current Residue Model: \n%s\n",expr);
}

if(front.compare("G:")==0){
int len=line.length()-2;
string model=line.substr(2,len);
strcpy(gexpr,model.c_str());
//printf("Current Global Model: \n%s\n",gexpr);
}

}
}
}
/*
 * constructor. 
 * Initializes all data with zeros and empty strings
 */
Parser::Parser(const char * file)
{
using namespace std;
    expr[0] = '\0';
    gexpr[0] = '\0';
    e = NULL;
    token[0] = '\0';
    token_type = NOTHING;
std::ifstream input;
input.open(file);
if(input.fail()){
printf("Reading score file failed");
}else{
string line;
int index=0;
string prevchain ("");
string prevres ("");
int prevresnum= -1;
char new_expr[EXPR_LEN_MAX];
while(getline(input, line)) {
string front=line.substr(0,2);

if(front.compare("M:")==0){
int len=line.length()-2;
string model=line.substr(2,len);
strcpy(expr,model.c_str());
//printf("Current Residue Model: \n%s\n",expr);
}

if(front.compare("G:")==0){
int len=line.length()-2;
string model=line.substr(2,len);
strcpy(gexpr,model.c_str());
//printf("Current Global Model: \n%s\n",gexpr);
}

}
}
}

void Parser::score(Score * score,std::vector<double> * resvector)
{
//Set the parser variable vector to the values in score
for(int s = 0; s<score->s0.size(); s++){
s0 = score->s0[s];     s10 = score->s10[s];     s20 = score->s20[s];   s30 = score->s30[s];   s40 = score->s40[s];   s50 = score->s50[s];   s60 = score->s60[s];   s70 = score->s70[s];   s80 = score->s80[s];
s1 = score->s1[s];     s11 = score->s11[s];     s21 = score->s21[s];   s31 = score->s31[s];   s41 = score->s41[s];   s51 = score->s51[s];   s61 = score->s61[s];   s71 = score->s71[s];   s81 = score->s81[s];
s2 = score->s2[s];     s12 = score->s12[s];     s22 = score->s22[s];   s32 = score->s32[s];   s42 = score->s42[s];   s52 = score->s52[s];   s62 = score->s62[s];   s72 = score->s72[s];   s82 = score->s82[s];
s3 = score->s3[s];     s13 = score->s13[s];     s23 = score->s23[s];   s33 = score->s33[s];   s43 = score->s43[s];   s53 = score->s53[s];   s63 = score->s63[s];   s73 = score->s73[s];   s83 = score->s83[s];
s4 = score->s4[s];     s14 = score->s14[s];     s24 = score->s24[s];   s34 = score->s34[s];   s44 = score->s44[s];   s54 = score->s54[s];   s64 = score->s64[s];   s74 = score->s74[s];   s84 = score->s84[s];
s5 = score->s5[s];     s15 = score->s15[s];     s25 = score->s25[s];   s35 = score->s35[s];   s45 = score->s45[s];   s55 = score->s55[s];   s65 = score->s65[s];   s75 = score->s75[s];   s85 = score->s85[s];
s6 = score->s6[s];     s16 = score->s16[s];     s26 = score->s26[s];   s36 = score->s36[s];   s46 = score->s46[s];   s56 = score->s56[s];   s66 = score->s66[s];   s76 = score->s76[s];   s86 = score->s86[s];
s7 = score->s7[s];     s17 = score->s17[s];     s27 = score->s27[s];   s37 = score->s37[s];   s47 = score->s47[s];   s57 = score->s57[s];   s67 = score->s67[s];   s77 = score->s77[s];   s87 = score->s87[s];
s8 = score->s8[s];     s18 = score->s18[s];     s28 = score->s28[s];   s38 = score->s38[s];   s48 = score->s48[s];   s58 = score->s58[s];   s68 = score->s68[s];   s78 = score->s78[s];   s88 = score->s88[s];
s9 = score->s9[s];     s19 = score->s19[s];     s29 = score->s29[s];   s39 = score->s39[s];   s49 = score->s49[s];   s59 = score->s59[s];   s69 = score->s69[s];   s79 = score->s79[s];   s89 = score->s89[s];
//Loop through each entry in score
    try
    {
        e = expr;                   // let e point to the start of the expression
        ans = 0;
        getToken();
        if (token_type == DELIMETER && *token == '\0')
        {
            throw Error(row(), col(), 4);
        }
        ans = parse_level1();
        // check for garbage at the end of the expression 
        // an expression ends with a character '\0' and token_type = delimeter
        if (token_type != DELIMETER || *token != '\0')
        {
            if (token_type == DELIMETER)
            {
                // user entered a not existing operator like "//"
                throw Error(row(), col(), 101, token);
            }
            else
            {
                throw Error(row(), col(), 5, token);
            }
        }  
        // add the answer to memory as variable "Ans"
        user_var.add("Ans", ans);
resvector->push_back(ans);
    }
    catch (Error err)
    {
        if (err.get_row() == -1)
        {
            snprintf(ans_str, sizeof(ans_str), "Error: %s (col %i)", err.get_msg(), err.get_col());
        }
        else
        {
            snprintf(ans_str, sizeof(ans_str), "Error: %s (ln %i, col %i)", err.get_msg(), err.get_row(), err.get_col());
        }
    }
}
}


char* Parser::parse()
{
    try
    {

        
        e = expr;                   // let e point to the start of the expression
        ans = 0;

        getToken();
        if (token_type == DELIMETER && *token == '\0')
        {
            throw Error(row(), col(), 4);
        }

        ans = parse_level1();

        // check for garbage at the end of the expression 
        // an expression ends with a character '\0' and token_type = delimeter
        if (token_type != DELIMETER || *token != '\0')
        {
            if (token_type == DELIMETER)
            {
                // user entered a not existing operator like "//"
                throw Error(row(), col(), 101, token);
            }
            else
            {
                throw Error(row(), col(), 5, token);
            }
        }  

        // add the answer to memory as variable "Ans"
        user_var.add("Ans", ans);

        snprintf(ans_str, sizeof(ans_str), "Ans = %g", ans);
    }
    catch (Error err)
    {
        if (err.get_row() == -1)
        {
            snprintf(ans_str, sizeof(ans_str), "Error: %s (col %i)", err.get_msg(), err.get_col());
        }
        else
        {
            snprintf(ans_str, sizeof(ans_str), "Error: %s (ln %i, col %i)", err.get_msg(), err.get_row(), err.get_col());
        }
    }

    return ans_str;
}

char* Parser::parse(const char new_expr[])
{
    try
    {
        // check the length of expr
        if (strlen(new_expr) > EXPR_LEN_MAX)
        {
            throw Error(row(), col(), 200);
        }

        // initialize all variables
        
	strcpy(expr, new_expr);     // copy the given expression to expr
        e = expr;                   // let e point to the start of the expression
        ans = 0;

        getToken();
        if (token_type == DELIMETER && *token == '\0')
        {
            throw Error(row(), col(), 4);
        }

        ans = parse_level1();

        // check for garbage at the end of the expression 
        // an expression ends with a character '\0' and token_type = delimeter
        if (token_type != DELIMETER || *token != '\0')
        {
            if (token_type == DELIMETER)
            {
                // user entered a not existing operator like "//"
                throw Error(row(), col(), 101, token);
            }
            else
            {
                throw Error(row(), col(), 5, token);
            }
        }  

        // add the answer to memory as variable "Ans"
        user_var.add("Ans", ans);

        snprintf(ans_str, sizeof(ans_str), "Ans = %g", ans);
    }
    catch (Error err)
    {
        if (err.get_row() == -1)
        {
            snprintf(ans_str, sizeof(ans_str), "Error: %s (col %i)", err.get_msg(), err.get_col());
        }
        else
        {
            snprintf(ans_str, sizeof(ans_str), "Error: %s (ln %i, col %i)", err.get_msg(), err.get_row(), err.get_col());
        }
    }

    return ans_str;
}


/*
 * checks if the given char c is a minus
 */
bool isMinus(const char c)
{
    if (c == 0) return 0;
    return c == '-';
}



/*
 * checks if the given char c is whitespace
 * whitespace when space chr(32) or tab chr(9)
 */
bool isWhiteSpace(const char c)
{
    if (c == 0) return 0;
    return c == 32 || c == 9;  // space or tab
}

/*
 * checks if the given char c is a delimeter
 * minus is checked apart, can be unary minus
 */
bool isDelimeter(const char c)
{
    if (c == 0) return 0;
    return strchr("&|<>=+/*%^!", c) != 0;
}

/*
 * checks if the given char c is NO delimeter
 */
bool isNotDelimeter(const char c)
{
    if (c == 0) return 0;
    return strchr("&|<>=+-/*%^!()", c) != 0;
}

/*
 * checks if the given char c is a letter or undersquare
 */
bool isAlpha(const char c)
{
    if (c == 0) return 0;
    return strchr("ABCDEFGHIJKLMNOPQRSTUVWXYZ_", toupper(c)) != 0;
}

/*
 * checks if the given char c is a digit or dot
 */
bool isDigitDot(const char c)
{
    if (c == 0) return 0;
    return strchr("0123456789.", c) != 0;
}

/*
 * checks if the given char c is a digit
 */
bool isDigit(const char c)
{
    if (c == 0) return 0;
    return strchr("0123456789", c) != 0;
}


/**
 * Get next token in the current string expr.
 * Uses the Parser data expr, e, token, t, token_type and err
 */
void Parser::getToken()
{
    token_type = NOTHING;
    char* t;           // points to a character in token
    t = token;         // let t point to the first character in token
    *t = '\0';         // set token empty

    //printf("\tgetToken e:{%c}, ascii=%i, col=%i\n", *e, *e, e-expr);

    // skip over whitespaces
    while (*e == ' ' || *e == '\t')     // space or tab
    {
        e++;
    }

    // check for end of expression
    if (*e == '\0')
    {
        // token is still empty
        token_type = DELIMETER;
        return;
    }

    // check for minus
    if (*e == '-')
    {
        token_type = DELIMETER;
        *t = *e;
        e++;
        t++;
        *t = '\0';  // add a null character at the end of token
        return;
    }

    // check for parentheses
    if (*e == '(' || *e == ')')
    {
        token_type = DELIMETER;
        *t = *e;
        e++;
        t++;
        *t = '\0';
        return;
    }

    // check for operators (delimeters)
    if (isDelimeter(*e))
    {
        token_type = DELIMETER;
        while (isDelimeter(*e))
        {
            *t = *e;
            e++;
            t++;
        }
        *t = '\0';  // add a null character at the end of token
        return;
    }
    
    // check for a value
    if (isDigitDot(*e))
    {
        token_type = NUMBER;
        while (isDigitDot(*e))
        {
            *t = *e;
            e++;
            t++;
        }
        
        // check for scientific notation like "2.3e-4" or "1.23e50"
        if (toupper(*e) == 'E')
        {
            *t = *e;
            e++;
            t++;
    
            if (*e == '+' || *e == '-')
            {
                *t = *e;
                e++;
                t++;
            }

            while (isDigit(*e))
            {
                *t = *e;
                e++;
                t++;
            }
        }
        
        *t = '\0';
        return;
    }

    // check for variables or functions
    if (isAlpha(*e))
    {
        while (isAlpha(*e) || isDigit(*e))
        //while (isNotDelimeter(*e))
        {
            *t = *e;
            e++;
            t++;
        }
        *t = '\0';  // add a null character at the end of token

        // check if this is a variable or a function.
        // a function has a parentesis '(' open after the name 
        char* e2 = NULL;
        e2 = e;

        // skip whitespaces
        while (*e2 == ' ' || *e2 == '\t')     // space or tab
        {
            e2++;
        }
        
        if (*e2 == '(') 
        {
            token_type = FUNCTION;
        }
        else
        {
            token_type = VARIABLE;
        }
        return;
    }

    // something unknown is found, wrong characters -> a syntax error
    token_type = UNKNOWN;
    while (*e != '\0')
    {
        *t = *e;
        e++;
        t++;
    }
    *t = '\0';
    throw Error(row(), col(), 1, token);

    return;
}


/*
 * assignment of variable or function
 */
double Parser::parse_level1()
{
    if (token_type == VARIABLE)
    {
        // copy current token
        char* e_now = e;
        TOKENTYPE token_type_now = token_type;
        char token_now[NAME_LEN_MAX+1];
        strcpy(token_now, token);
        
        getToken();
        if (strcmp(token, "=") == 0)
        {
            // assignment
            double ans;
            getToken();
            ans = parse_level2();
            if (user_var.add(token_now, ans) == false)
            {
                throw Error(row(), col(), 300);
            }
            return ans;
        }
        else
        {
            // go back to previous token
            e = e_now;
            token_type = token_type_now;
            strcpy(token, token_now);
        }
    }

    return parse_level2();
}


/*
 * conditional operators and bitshift
 */
double Parser::parse_level2()
{
    int op_id;
    double ans;
    ans = parse_level3();

    op_id = get_operator_id(token);
    while (op_id == AND || op_id == OR || op_id == BITSHIFTLEFT || op_id == BITSHIFTRIGHT)
    {
        getToken();
        ans = eval_operator(op_id, ans, parse_level3());
        op_id = get_operator_id(token);
    }

    return ans;
}

/*
 * conditional operators
 */
double Parser::parse_level3()
{
    int op_id;
    double ans;
    ans = parse_level4();

    op_id = get_operator_id(token);
    while (op_id == EQUAL || op_id == UNEQUAL || op_id == SMALLER || op_id == LARGER || op_id == SMALLEREQ || op_id == LARGEREQ)
    {
        getToken();
        ans = eval_operator(op_id, ans, parse_level4());
        op_id = get_operator_id(token);
    }

    return ans;
}

/*
 * add or subtract
 */
double Parser::parse_level4()
{
    int op_id;
    double ans;
    ans = parse_level5();
    
    op_id = get_operator_id(token);
    while (op_id == PLUS || op_id == MINUS)
    {
        getToken();
        ans = eval_operator(op_id, ans, parse_level5());
        op_id = get_operator_id(token);
    }

    return ans;
}


/*
 * multiply, divide, modulus, xor
 */
double Parser::parse_level5()
{
    int op_id;
    double ans;
    ans = parse_level6();

    op_id = get_operator_id(token);
    while (op_id == MULTIPLY || op_id == DIVIDE || op_id == MODULUS || op_id == XOR)
    {
        getToken();
        ans = eval_operator(op_id, ans, parse_level6());
        op_id = get_operator_id(token);
    }

    return ans;
}


/*
 * power
 */
double Parser::parse_level6()
{
    int op_id;
    double ans;
    ans = parse_level7();

    op_id = get_operator_id(token);
    while (op_id == POW)
    {
        getToken();
        ans = eval_operator(op_id, ans, parse_level7());
        op_id = get_operator_id(token);
    }

    return ans;
}

/*
 * Factorial
 */
double Parser::parse_level7()
{
    int op_id;
    double ans;
    ans = parse_level8();

    op_id = get_operator_id(token);
    while (op_id == FACTORIAL)
    {
        getToken();
        // factorial does not need a value right from the 
        // operator, so zero is filled in.
        ans = eval_operator(op_id, ans, 0.0);
        op_id = get_operator_id(token);
    }

    return ans;
}

/*
 * Unary minus
 */
double Parser::parse_level8()
{
    double ans;
    
    int op_id = get_operator_id(token);    
    if (op_id == MINUS)
    {
        getToken();
        ans = parse_level9();
        ans = -ans;
    }
    else
    {
        ans = parse_level9();
    }
    
    return ans;
}


/*
 * functions
 */
double Parser::parse_level9()
{
    char fn_name[NAME_LEN_MAX+1];
    double ans;

    if (token_type == FUNCTION)
    {
        strcpy(fn_name, token);
        getToken();
        ans = eval_function(fn_name, parse_level10());
    }
    else
    {
        ans = parse_level10();
    }
    
    return ans;
}


/*
 * parenthesized expression or value
 */
double Parser::parse_level10()
{
    // check if it is a parenthesized expression
    if (token_type == DELIMETER)
    {
        if (token[0] == '(' & token[1] == '\0')
        {
            getToken();
            double ans = parse_level2();
            if (token_type != DELIMETER || token[0] != ')' || token[1] || '\0')
            {
                throw Error(row(), col(), 3);
            }
            getToken();
            return ans;
        }
    }

    // if not parenthesized then the expression is a value
    return parse_number();
}


double Parser::parse_number()
{
double ans = 0;

    switch (token_type)
    {
        case NUMBER:
            // this is a number
            ans = strtod(token, NULL);
            getToken();
            break;

        case VARIABLE:
            // this is a variable
            ans = eval_variable(token);
            getToken();  
            break;
            
        default:
            // syntax error or unexpected end of expression
            if (token[0] == '\0')
            {
                throw Error(row(), col(), 6);
            }
            else
            {
                throw Error(row(), col(), 7);
            }
            break;
    }

    return ans;
}


/*
 * returns the id of the given operator
 * treturns -1 if the operator is not recognized
 */
int Parser::get_operator_id(const char op_name[])
{
    // level 2
    if (!strcmp(op_name, "&")) {return AND;}
    if (!strcmp(op_name, "|")) {return OR;}
    if (!strcmp(op_name, "<<")) {return BITSHIFTLEFT;}
    if (!strcmp(op_name, ">>")) {return BITSHIFTRIGHT;}

    // level 3
    if (!strcmp(op_name, "=")) {return EQUAL;}
    if (!strcmp(op_name, "<>")) {return UNEQUAL;}
    if (!strcmp(op_name, "<")) {return SMALLER;}
    if (!strcmp(op_name, ">")) {return LARGER;}
    if (!strcmp(op_name, "<=")) {return SMALLEREQ;}
    if (!strcmp(op_name, ">=")) {return LARGEREQ;}

    // level 4
    if (!strcmp(op_name, "+")) {return PLUS;}
    if (!strcmp(op_name, "-")) {return MINUS;}

    // level 5
    if (!strcmp(op_name, "*")) {return MULTIPLY;}
    if (!strcmp(op_name, "/")) {return DIVIDE;}
    if (!strcmp(op_name, "%")) {return MODULUS;}
    if (!strcmp(op_name, "||")) {return XOR;}

    // level 6
    if (!strcmp(op_name, "^")) {return POW;}

    // level 7
    if (!strcmp(op_name, "!")) {return FACTORIAL;}

    return -1;
}


/*
 * evaluate an operator for given valuess
 */
double Parser::eval_operator(const int op_id, const double &lhs, const double &rhs)
{
    switch (op_id)
    {
        // level 2
        case AND:           return static_cast<int>(lhs) & static_cast<int>(rhs);
        case OR:            return static_cast<int>(lhs) | static_cast<int>(rhs);
        case BITSHIFTLEFT:  return static_cast<int>(lhs) << static_cast<int>(rhs);
        case BITSHIFTRIGHT: return static_cast<int>(lhs) >> static_cast<int>(rhs);

        // level 3
        case EQUAL:     return lhs == rhs;
        case UNEQUAL:   return lhs != rhs;
        case SMALLER:   return lhs < rhs;
        case LARGER:    return lhs > rhs;
        case SMALLEREQ: return lhs <= rhs;
        case LARGEREQ:  return lhs >= rhs;
        
        // level 4
        case PLUS:      return lhs + rhs;
        case MINUS:     return lhs - rhs;
        
        // level 5
        case MULTIPLY:  return lhs * rhs;
        case DIVIDE:    return lhs / rhs;
        case MODULUS:   return static_cast<int>(lhs) % static_cast<int>(rhs); // todo: give a warning if the values are not integer?
        case XOR:       return static_cast<int>(lhs) ^ static_cast<int>(rhs);
        
        // level 6
        case POW:       return pow(lhs, rhs);
        
        // level 7
        case FACTORIAL: return factorial(lhs);
    }

    throw Error(row(), col(), 104, op_id);    
    return 0;
}


/*
 * evaluate a function
 */
double Parser::eval_function(const char fn_name[], const double &value)
{
    try
    {
        // first make the function name upper case
        char fnU[NAME_LEN_MAX+1];
        toupper(fnU, fn_name);

        // arithmetic 
        if (!strcmp(fnU, "ABS")) {return abs(value);}
        if (!strcmp(fnU, "EXP")) {return exp(value);}
        if (!strcmp(fnU, "SIGN")) {return sign(value);}
        if (!strcmp(fnU, "SQRT")) {return sqrt(value);}
        if (!strcmp(fnU, "LOG")) {return log(value);}
        if (!strcmp(fnU, "LOG10")) {return log10(value);}
    
        // trigonometric
        if (!strcmp(fnU, "SIN")) {return sin(value);}
        if (!strcmp(fnU, "COS")) {return cos(value);}
        if (!strcmp(fnU, "TAN")) {return tan(value);}
        if (!strcmp(fnU, "ASIN")) {return asin(value);}
        if (!strcmp(fnU, "ACOS")) {return acos(value);}
        if (!strcmp(fnU, "ATAN")) {return atan(value);}
    
        // probability
        if (!strcmp(fnU, "FACTORIAL")) {return factorial(value);}
    }
    catch (Error err)
    {
        // retrow error, add information about column and row of occurance
        // TODO: doesnt work yet
        throw Error(col(), row(), err.get_id(), err.get_msg());
    }
    
    // unknown function
    throw Error(row(), col(), 102, fn_name);
    return 0;
}


/*
 * evaluate a variable
 */
double Parser::eval_variable(const char var_name[])
{
    // first make the variable name uppercase
    char varU[NAME_LEN_MAX+1];
    toupper(varU, var_name);
    // check for built-in variables
    if (!strcmp(varU, "E")) {return 2.7182818284590452353602874713527;}
    if (!strcmp(varU, "PI")) {return 3.1415926535897932384626433832795;}
    if (!strcmp(varU, "RND")) {float rnd=rand()/(float)RAND_MAX; return rnd; }
    if (!strcmp(varU, "S0")) {return s0;}
    if (!strcmp(varU, "S1")) {return s1;}
    if (!strcmp(varU, "S2")) {return s2;}
    if (!strcmp(varU, "S3")) {return s3;}
    if (!strcmp(varU, "S4")) {return s4;}
    if (!strcmp(varU, "S5")) {return s5;}
    if (!strcmp(varU, "S6")) {return s6;}
    if (!strcmp(varU, "S7")) {return s7;}
    if (!strcmp(varU, "S8")) {return s8;}
    if (!strcmp(varU, "S9")) {return s9;}
    
    if (!strcmp(varU, "S10")) {return s10;}
    if (!strcmp(varU, "S11")) {return s11;}
    if (!strcmp(varU, "S12")) {return s12;}
    if (!strcmp(varU, "S13")) {return s13;}
    if (!strcmp(varU, "S14")) {return s14;}
    if (!strcmp(varU, "S15")) {return s15;}
    if (!strcmp(varU, "S16")) {return s16;}
    if (!strcmp(varU, "S17")) {return s17;}
    if (!strcmp(varU, "S18")) {return s18;}
    if (!strcmp(varU, "S19")) {return s19;}
    if (!strcmp(varU, "S20")) {return s20;}
    if (!strcmp(varU, "S21")) {return s21;}
    if (!strcmp(varU, "S22")) {return s22;}
    if (!strcmp(varU, "S23")) {return s23;}
    if (!strcmp(varU, "S24")) {return s24;}
    if (!strcmp(varU, "S25")) {return s25;}
    if (!strcmp(varU, "S26")) {return s26;}
    if (!strcmp(varU, "S27")) {return s27;}
    if (!strcmp(varU, "S28")) {return s28;}
    if (!strcmp(varU, "S29")) {return s29;}
    if (!strcmp(varU, "S30")) {return s30;}
    if (!strcmp(varU, "S31")) {return s31;}
    if (!strcmp(varU, "S32")) {return s32;}
    if (!strcmp(varU, "S33")) {return s33;}
    if (!strcmp(varU, "S34")) {return s34;}
    if (!strcmp(varU, "S35")) {return s35;}
    if (!strcmp(varU, "S36")) {return s36;}
    if (!strcmp(varU, "S37")) {return s37;}
    if (!strcmp(varU, "S38")) {return s38;}
    if (!strcmp(varU, "S39")) {return s39;}
    if (!strcmp(varU, "S40")) {return s40;}
    if (!strcmp(varU, "S41")) {return s41;}
    if (!strcmp(varU, "S42")) {return s42;}
    if (!strcmp(varU, "S43")) {return s43;}
    if (!strcmp(varU, "S44")) {return s44;}
    if (!strcmp(varU, "S45")) {return s45;}
    if (!strcmp(varU, "S46")) {return s46;}
    if (!strcmp(varU, "S47")) {return s47;}
    if (!strcmp(varU, "S48")) {return s48;}
    if (!strcmp(varU, "S49")) {return s49;}
    if (!strcmp(varU, "S50")) {return s50;}
    if (!strcmp(varU, "S51")) {return s51;}
    if (!strcmp(varU, "S52")) {return s52;}
    if (!strcmp(varU, "S53")) {return s53;}
    if (!strcmp(varU, "S54")) {return s54;}
    if (!strcmp(varU, "S55")) {return s55;}
    if (!strcmp(varU, "S56")) {return s56;}
    if (!strcmp(varU, "S57")) {return s57;}
    if (!strcmp(varU, "S58")) {return s58;}
    if (!strcmp(varU, "S59")) {return s59;}
    if (!strcmp(varU, "S60")) {return s60;}
    if (!strcmp(varU, "S61")) {return s61;}
    if (!strcmp(varU, "S62")) {return s62;}
    if (!strcmp(varU, "S63")) {return s63;}
    if (!strcmp(varU, "S64")) {return s64;}
    if (!strcmp(varU, "S65")) {return s65;}
    if (!strcmp(varU, "S66")) {return s66;}
    if (!strcmp(varU, "S67")) {return s67;}
    if (!strcmp(varU, "S68")) {return s68;}
    if (!strcmp(varU, "S69")) {return s69;}
    if (!strcmp(varU, "S70")) {return s70;}
    if (!strcmp(varU, "S71")) {return s71;}
    if (!strcmp(varU, "S72")) {return s72;}
    if (!strcmp(varU, "S73")) {return s73;}
    if (!strcmp(varU, "S74")) {return s74;}
    if (!strcmp(varU, "S75")) {return s75;}
    if (!strcmp(varU, "S76")) {return s76;}
    if (!strcmp(varU, "S77")) {return s77;}
    if (!strcmp(varU, "S78")) {return s78;}
    if (!strcmp(varU, "S79")) {return s79;}
    if (!strcmp(varU, "S80")) {return s80;}
    if (!strcmp(varU, "S81")) {return s81;}
    if (!strcmp(varU, "S82")) {return s82;}
    if (!strcmp(varU, "S83")) {return s83;}
    if (!strcmp(varU, "S84")) {return s84;}
    if (!strcmp(varU, "S85")) {return s85;}
    if (!strcmp(varU, "S86")) {return s86;}
    if (!strcmp(varU, "S87")) {return s87;}
    if (!strcmp(varU, "S88")) {return s88;}
    if (!strcmp(varU, "S89")) {return s89;}
 
    // check for user defined variables
    double ans;
    if (user_var.get_value(var_name, &ans))
    {
        return ans;
    }

    // unknown variable
    throw Error(row(), col(), 103, var_name);
    return 0;
}



/*
 * Shortcut for getting the current row value (one based)
 * Returns the line of the currently handled expression
 */
int Parser::row()
{
    return -1;
}

/*
 * Shortcut for getting the current col value (one based)
 * Returns the column (position) where the last token starts
 */
int Parser::col()
{
    return e-expr-strlen(token)+1;
}
