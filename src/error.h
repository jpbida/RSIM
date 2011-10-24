/*
error class for the parser

todo:
    add the line and column where the error occured to the error message
*/


#ifndef ERROR_H
#define ERROR_H

#include "stdio.h"
#include "stdarg.h"
#include <string>

using namespace std;



class Error {
    public:
        Error(const int row, const int col, const int id, ...);

        int get_row() {return err_row;} // Returns the row of the error
        int get_col() {return err_col;} // Returns the column of the error
        int get_id() {return err_id;}   // Returns the id of the error
        char* get_msg() {return msg;}   // Returns a pointer to the error msg

    private:
        int err_row;    // row where the error occured
        int err_col;    // column (position) where the error occured
        int err_id;     // id of the error
        char msg[255];
        
        char* msgdesc(const int id);
};


#endif
