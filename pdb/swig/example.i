%include "exception.i"
%exception {

    try {
        $action
        } catch (const std::exception &e) {
        
        SWIG_exception_fail(SWIG_RuntimeError, e.what());
        }
    }

%module example
%{
/* Put headers and other declarations here */
#include "example.h"
%}

%include "std_vector.i"
%include "std_string.i"

%template(IntVector) std::vector<int>;
%template(DoubleVector) std::vector<double>;
    
%include "example.h"
