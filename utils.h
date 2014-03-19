#define SECPERDAY 86400.0

typedef struct FCOMPLEX{
  float r;
  float i;
} fcomplex;

typedef struct DCOMPLEX{
  double r;
  double i;
} dcomplex;

double *gen_dvect(long length);
fcomplex *gen_cvect(long length);
float *gen_fvect(long length);
double **gen_dmatrix(long nrows, long ncols);
int compare_doubles(const void *a, const void *b);
void davg_dvar(double *x, int n, double *mean, double *var);
int get_filelen(FILE *file, size_t size);
double *read_events(FILE *infile, int bin, int days, int *numevents,
                    double MJD0, double Ttot, double startfrac, double endfrac,
                    double offset);
