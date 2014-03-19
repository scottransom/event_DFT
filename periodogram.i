%module periodogram_c
%include numpy.i

%apply double* IN_1D_DOUBLE { double *xx, double *tt };
%apply int ARRAYLEN { int numf };
float *periodogram(double *xx, double *tt, int nn, 
		   double lof, double df, int numf);

%apply double* IN_1D_DOUBLE { double *tt };
%apply int ARRAYLEN { int numf };
fcomplex *toafft(double *tt, int nn, 
		 double lof, double df, int numf);

%apply double* IN_1D_DOUBLE { double *tt };
void *prep_toafft(double *tt, int nn, double lof, double df);

void *free_toafft(double *tt, int nn, double lof, double df);

%apply double *OUTPUT { double *freq }; 
dcomplex calc_toafft_point(double *freq);
