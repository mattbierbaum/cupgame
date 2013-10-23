%module cupgamelib
%{
  #define SWIG_FILE_WITH_INIT
  #include "cupgamelib.h"
%}

%include "numpy.i"
%init %{
  import_array();
%}

%apply (int DIM1, int* INPLACE_ARRAY1) {(int NS, int *bounces)};
%apply (int DIM1, double* INPLACE_ARRAY1) {(int NP, double *pos)};
%apply (int DIM1, double* INPLACE_ARRAY1) {(int NV, double *vel)};
%apply (int DIM1, double* INPLACE_ARRAY1) {(int NT, double *traj)};
%include "cupgamelib.h"
