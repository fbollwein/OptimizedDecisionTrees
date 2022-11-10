%begin %{
    #define SWIG_PYTHON_CAST_MODE
%}

%module DecisionTree


%{
//#define SWIG_FILE_WITH_INIT
#include "src/DecisionTree.h"
%}

%include <std_string.i>
%include "std_vector.i"


namespace std {
    %template(VectorDouble) vector<double>;
	%template(VectorInt) vector<int>;
    %template(VectorDoubleDouble) vector<vector<double>>;
    %template(VectorString) vector<string>;
};

//%include "numpy.i"

//%init %{
//import_array();
//%}
//%apply (int DIM1) {(int x)};


%include "src/DecisionTree.h"
