#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include "utils.h"

double *gen_dvect(long length)
{
    double *v;
    
    v = (double *) malloc((size_t) (sizeof(double) * length));
    if (!v) {
        perror("\nError in gen_dvect()");
        printf("\n");
        exit(-1);
    }
    return v;
}


fcomplex *gen_cvect(long length)
{
    fcomplex *v;
    
    v = (fcomplex *) malloc((size_t) (sizeof(fcomplex) * length));
    if (!v) {
        perror("\nError in gen_cvect()");
        printf("\n");
        exit(-1);
    }
    return v;
}


float *gen_fvect(long length)
{
    float *v;
    
    v = (float *) malloc((size_t) (sizeof(float) * length));
    if (!v) {
        perror("\nError in gen_fvect()");
        printf("\n");
        exit(-1);
    }
    return v;
}


double **gen_dmatrix(long nrows, long ncols)
{
    /* Note:  To free this matrix, assuming you called it with:    */
    /*             x = gen_dmatrix(10,10);                         */
    /*        all you need to do is the following:                 */
    /*             free(x[0]) ; free(x) ;                          */
    /*        The order is important!                              */
    
    long i;
    double **m;
    
    m = (double **) malloc((size_t) (nrows * sizeof(double *)));
    if (!m) {
        perror("\nError in 1st malloc() in gen_dmatrix()");
        printf("\n");
        exit(-1);
    }
    m[0] = (double *) malloc((size_t) ((nrows * ncols) * sizeof(double)));
    if (!m[0]) {
        perror("\nError in 2nd malloc() in gen_dmatrix()");
        printf("\n");
        exit(-1);
    }
    for (i = 1; i < nrows; i++)
        m[i] = m[i - 1] + ncols;
    return m;
}


int compare_doubles(const void *a, const void *b)
/* qsort comparison function for doubles */
{
    const double *da = (const double *) a;
    const double *db = (const double *) b;
    
    return (*da > *db) - (*da < *db);
}


void davg_dvar(double *x, int n, double *mean, double *var)
/* For a double vector, *x, of length n, this routine  */
/* returns the mean and variance of *x.                */
{
    long i;
    double an=0.0, an1=0.0, dx;
    
    /*  Modified (29 June 98) C version of the following:        */
    /*  ALGORITHM AS 52  APPL. STATIST. (1972) VOL.21, P.226     */
    /*  Returned values were checked with Mathematica 3.01       */
    
    if (n < 1) {
        printf("\vVector length must be > 0 in avg_var().  Exiting\n");
        exit(1);
    } else {
        *mean = (double) x[0];
        *var = 0.0;
    }
    for (i = 1 ; i < n ; i++){
        an = (double) (i + 1);
        an1 = (double) (i);
        dx = (x[i] - *mean) / an;
        *var += an * an1 * dx * dx;
        *mean += dx;
    }
    if (n > 1)
        *var /= an1;
    return;
}


int get_filelen(FILE *file, size_t size)
{
    int filenum, rt;
    struct stat buf;
    
    filenum = fileno(file);
    rt = fstat(filenum, &buf);
    if (rt == -1){
        perror("\nError in get_filelen()");
        printf("\n");
        exit(-1);
    }
    return (int) (buf.st_size / size);
}


double *read_events(FILE *infile, int bin, int days, int *numevents,
                    double MJD0, double Ttot, double startfrac, double endfrac,
                    double offset)
/* This routine reads a set of events from the open file 'infile'.     */
/* It returns a double precision vector of events in seconds from the  */
/* first event.  If 'bin' is true the routine treats the data as       */
/* binary double precision (otherwise text).  If 'days' is 1 then the  */
/* data is assumed to be in days since the 'inf' EPOCH (0 is sec from  */
/* EPOCH in 'inf').  If 'days' is 2, the data are assumed to be MJDs.  */
/* The number of events read is placed in 'numevents', and the raw     */
/* event is placed in 'firstevent'.  MJD0 is the time to use for the   */
/* reference time.  Ttot is the duration of the observation. 'start'   */
/* and 'end' are define the fraction of the observation that we are    */
/* interested in.  'offset' is a time offset to apply to the events.   */
{
    int N=0, nn=0, goodN=0;
    double *ts, *goodts, dtmp, lotime, hitime;
    char line[80];
    
    if (bin){
        N = get_filelen(infile, sizeof(double));
    } else {
        /* Read the input file once to count events */
        while (1){
            fgets(line, 80, infile);
            if (!feof(infile)){
                if (line[0]!='#' && sscanf(line, "%lf", &dtmp)==1) N++;
            } else {
                break;
            }
        }
    }
    
    /* Allocate the event arrays */
    
    ts = (double *)malloc(N * sizeof(double));
    
    /* Rewind and read the events for real */
    
    rewind(infile);
    if (bin){
        fread(ts, sizeof(double), N, infile);
    } else {
        while (1){
            fgets(line, 80, infile);
            if (!feof(infile)){
                if (line[0]!='#' && sscanf(line, "%lf", &ts[nn])==1) nn++;
            } else {
                break;
            }
        }
    }
    
    /* Sort the events  */
    
    qsort(ts, N, sizeof(double), compare_doubles);
    
    /* If there is no offset specified and the data are non-MJD */
    /* days or seconds, then set the offset to be the first event */
    
    if (offset==0.0 && days < 2)
        offset = -ts[0];
    
    /* Convert all the events to MJD */
    
    if (days==0){ /* Events are in seconds since MJD0 */
        for (nn=0; nn<N; nn++)
            ts[nn] = MJD0+(ts[nn]+offset)/SECPERDAY;
    } else if (days==1){ /* Events are in days since MJD0 */
        for (nn=0; nn<N; nn++)
            ts[nn] = MJD0+(ts[nn]+offset);
    } else if (days==2 &&
               offset != 0.0){ /* Events are in MJD with an offset */
        for (nn=0; nn<N; nn++)
            ts[nn] += offset;
    }
    
    /* Count how many events are within our range and only keep them */
    
    if (Ttot==0.0)
        Ttot = ts[N-1];
    lotime = MJD0 + startfrac*Ttot;
    hitime = MJD0 + endfrac*Ttot + 1.0e-15;
    for (nn=0; nn<N; nn++)
        if (ts[nn] >= lotime && ts[nn] < hitime) goodN++;
    if (goodN != N){
        goodts = (double *)malloc(goodN * sizeof(double));
        goodN = 0;
        for (nn=0; nn<N; nn++){
            if (ts[nn] >= lotime && ts[nn] < hitime){
                goodts[goodN] = ts[nn];
                goodN++;
            }
        }
        free(ts);
        ts = goodts;
        N = goodN;
    } else {
        goodts = ts;
    }
    *numevents = N;
    
    /* Convert the events to seconds from MJD0 */
    
    for (nn=0; nn<N; nn++)
        goodts[nn] = (goodts[nn]-MJD0)*SECPERDAY;
    
    return goodts;
}
