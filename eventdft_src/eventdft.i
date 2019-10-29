%module eventdft_c

%{
#define SWIG_FILE_WITH_INIT
%}

%include "numpy.i"
%include "typemaps.i"

%init %{
    import_array();
%}

%{
#include "utils.h"
#include "eventdft.h"
%}

typedef struct FCOMPLEX {
    float r, i;
} fcomplex;

// setup a typemap for fcomplex to numpy complex
%numpy_typemaps(fcomplex, NPY_CFLOAT, long)

%apply (fcomplex** ARGOUTVIEWM_ARRAY1, long* DIM1) {(fcomplex** vect, long *no)}
%apply (double* INPLACE_ARRAY1, long DIM1) {(double *in, long ni)}

%exception
{
    errno = 0;
    $action

    if (errno != 0)
    {
        switch(errno)
        {
            case ENOMEM:
                PyErr_Format(PyExc_MemoryError, "Failed malloc()");
                break;
            default:
                PyErr_Format(PyExc_Exception, "Unknown exception");
        }
        SWIG_fail;
    }
}

%rename (eventdft) wrap_eventdft;
%inline %{
void wrap_eventdft(double *in, long ni,
                   double lof, double df, long numf,
                   fcomplex** vect, long *no)
{
    *vect = eventdft(in, ni, lof, df, numf);
    *no = numf;
}
%}
%clear (fcomplex **vect, long *no);
%clear (double *in, long ni);

%apply (double* INPLACE_ARRAY1, long DIM1) {(double *in, long nn)};
void prep_eventdft(double *in, long nn, int maxnumharmsum,
                   double lof, double df);

void free_eventdft(void);
