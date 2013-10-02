%module smectics
%{
  #define SWIG_FILE_WITH_INIT
  #include "cupgame.h"
%}

%include "numpy.i"
%init %{
  import_array();
%}

%apply (int DIM1, double* INPLACE_ARRAY1) {(int len, double *onpy)};
%apply (int DIM1, double* IN_ARRAY1) {(int len, double *inpy)};
%include "smectics.h"
