/*
Additional functions

*/

#include "error.h"

using namespace std;



/*
 * calculate factorial of value
 * for example 5! = 5*4*3*2*1 = 120
 */
double factorial(double value)
{
    double res;
    int v = static_cast<int>(value);
    
    if (value != static_cast<double>(v))
    {
        throw Error(-1, -1, 400, "factorial");
    }
    
    res = v;
    v--;
    while (v > 1)
    {
        res *= v;
        v--;
    }

    if (res == 0) res = 1;        // 0! is per definition 1
    return res;
}

/* 
 * calculate the sign of the given value
 */
double sign(double value)
{
    if (value > 0) return 1;
    if (value < 0) return -1;
    return 0;
}
