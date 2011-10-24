/*
class for holding the list with user defined variables
*/


#ifndef USER_VARIABLES_H
#define USER_VARIABLES_H

#include <cstdio>
#include <cstring>
#include <vector>

#include "constants.h"

using namespace std;



void toupper(char upper[], const char str[]);


class Variablelist {
    public:
        bool exist(const char* name);
        bool add(const char* name, double value);
        bool del(const char* name);

        bool get_value(const char* name, double* value);
        bool get_value(const int id, double* value);
        int  get_id(const char* name);
        bool set_value(const char* name, const double value);

    private:
        struct VAR
        {
            char name[NAME_LEN_MAX+1];
            double value;
        };

        vector<VAR> var;
};


#endif
