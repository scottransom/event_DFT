#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "utils.h"
#include "eventdft.h"

static double normconst, totaltime, lofreq, dfreq, maxharmsum;
static double **coss, **sins, **alph, **beta;
static fcomplex *amplitudes;
static int initialized=0, numevents;
static long long numfreqs;

void free_eventdft()
/* Free the static vectors used by the *_eventdft() routines */
{
    if (initialized){
        free(coss[0]);
        free(coss);
        free(sins[0]);
        free(sins);
        free(alph[0]);
        free(alph);
        free(beta[0]);
        free(beta);
        free(amplitudes);
        initialized = 0;
    }
}

void prep_eventdft(double *events, int numevts, int maxnumharmsum, 
                   double lof, double df)
/* Prepare the static variables and vectors for efficient calculation
   of DFT points based on events.  The FT will begin at frequency
   'lof' (Hz) and continue upwards in frquency by stepsize 'df'.
   There are 'numevents' input events with times 'events' (s).  Up to
   'maxharmsum' harmonics will be summed.  */
{
    int ii, jj;
    double delta, theta, temp, midtime, time;
    
    if (initialized){
        printf("Re-initializing eventdft variables and arrays...\n");
        free_eventdft();
    }
    numevents = numevts;
    maxharmsum = maxnumharmsum;
    numfreqs = 0;
    normconst = 1.0 / sqrt((double) numevents);
    lofreq = lof;
    dfreq = df;
    coss = gen_dmatrix(maxharmsum, numevents);
    sins = gen_dmatrix(maxharmsum, numevents);
    alph = gen_dmatrix(maxharmsum, numevents);
    beta = gen_dmatrix(maxharmsum, numevents);
    amplitudes = gen_cvect(maxharmsum);
    totaltime = events[numevents-1] - events[0];
    midtime = 0.5 * totaltime;
    for (ii = 0; ii < maxharmsum; ii++){
        for (jj = 0; jj < numevents; jj++){
            time = events[jj] - midtime;
            temp = -TWOPI * (ii + 1) * time;
            delta = temp * df;
            theta = temp * lofreq;
            temp = sin(0.5 * delta);
            alph[ii][jj] = -2.0 * temp * temp;
            beta[ii][jj] = sin(delta);
            coss[ii][jj] = cos(theta);
            sins[ii][jj] = sin(theta);
        }
    }
    initialized = 1;
}

fcomplex *calc_eventdft_point(double *freq)
/* Calculate the next frequency of the event DFT. A pointer to the
   vector of normalized Fourier amplitudes (of length maxharmsum) is
   returned.  '*freq' contains the frequency (Hz) that it corresponds
   to (dundamental).  prep_eventdft() must have been called to
   initialize the calcs. */
{
    int ii, jj;
    double aa, bb, cc, ss, *aptr, *bptr, *cptr, *sptr;
    dcomplex result={0.0, 0.0};
    
    (*freq) = lofreq + numfreqs * dfreq;
    for (ii = 0; ii < maxharmsum; ii++){
        aptr = alph[ii];
        bptr = beta[ii];
        cptr = coss[ii];
        sptr = sins[ii];
        result.r = result.i = 0.0;
        for (jj = 0; jj < numevents; jj++){
            aa = aptr[jj];
            bb = bptr[jj];
            cc = cptr[jj];
            ss = sptr[jj];
            result.r += cc;
            result.i += ss;
            cptr[jj] = cc * aa - ss * bb + cc;
            sptr[jj] = ss * aa + cc * bb + ss;
        }
        amplitudes[ii].r = result.r * normconst;
        amplitudes[ii].i = result.i * normconst;
    }
    numfreqs++;
    return amplitudes;
}


fcomplex *eventdft(double *events, int numevents, 
                   double lof, double df, int numf)
/* Return a set of normalized Fourier amplitudes/ calculated using an
   event DFT for a set of 'nn' events 'tt'.  The returned array is
   allocated.  */
{
    int ii;
    double freq;
    fcomplex *famp, *amps;
    
    amps = gen_cvect(numf);
    prep_eventdft(events, numevents, 1, lof, df);
    for (ii = 0; ii < numf; ii++){
        famp = calc_eventdft_point(&freq);
        amps[ii].r = famp->r;
        amps[ii].i = famp->i;
    }
    free_eventdft();
    return amps;
}


float *periodogram(double *xx, double *tt, int nn, 
                   double lof, double df, int numf)
/* Return the normalized Lomb-Scargle Periodogram powers of 'numf'
   frequencies (Hz) from the lowest freq 'lof' upwards by stepsize
   'df'.  There are 'nn' input data points with amplitudes 'xx' and
   times 'tt' (s).  The returned power vector is dynamically
   allocated.  */
{
    int ii, jj;
    float *pows;
    double avg, var, ivar, c, cc, cwtau;
    double s, ss, sumc, sumcxx, sums, sumsh, sumsxx,swtau;
    double wtau, ttavg, ttmax, ttmin;
    double arg, wtemp, *xxnorm, *wi, *wpi, *wpr, *wr;
    
    /* Set-up */
    davg_dvar(xx, nn, &avg, &var);
    if (var==0.0)
        ivar = 0.5;
    else
        ivar = 0.5 / var;
    wr = gen_dvect(nn);
    wi = gen_dvect(nn);
    wpr = gen_dvect(nn);
    wpi = gen_dvect(nn);
    xxnorm = gen_dvect(nn);
    pows = gen_fvect(numf);
    
    /* Scale the times around the midpt */
    ttmax = ttmin = tt[0];
    for (ii = 0; ii < nn; ii++){
        if (tt[ii] > ttmax) ttmax = tt[ii];
        if (tt[ii] < ttmin) ttmin = tt[ii];
    }
    ttavg = 0.5 * (ttmax + ttmin);
    
    /* Generate the trig recurrence values */
    c = cos(TWOPI * lof);
    s = sin(TWOPI * lof);
    for (ii = 0; ii < nn; ii++){
        arg = TWOPI * ((tt[ii] - ttavg) * df);
        wtemp = sin(0.5 * arg);
        wpr[ii] = -2.0 * wtemp * wtemp;
        wpi[ii] = sin(arg);
        wtemp = TWOPI * ((tt[ii] - ttavg) * lof);
        wr[ii] = cos(wtemp);
        wi[ii] = sin(wtemp);
        if (var==0.0)
            xxnorm[ii] = 1.0;
        else
            xxnorm[ii] = xx[ii] - avg;
    }
    
    /* Calculate the periodogram */
    for (ii = 0; ii < numf; ii++){
        sumsh = sumc = 0.0;
        for (jj=0; jj<nn; jj++){
            c = wr[jj];
            s = wi[jj];
            sumsh += s * c;
            sumc += (c - s) * (c + s);
        }
        wtau = 0.5 * atan2(2.0 * sumsh, sumc);
        cwtau = cos(wtau);
        swtau = sin(wtau);
        sums = sumc = sumsxx = sumcxx = 0.0;
        
        /* Step through the data points */
        for (jj = 0; jj < nn; jj++){
            c = wr[jj];
            s = wi[jj];
            ss = s * cwtau - c * swtau;
            cc = c * cwtau + s * swtau;
            sums += ss * ss;
            sumc += cc * cc;
            sumsxx += xxnorm[jj] * ss;
            sumcxx += xxnorm[jj] * cc;
            wr[jj] += c * wpr[jj] - s * wpi[jj];
            wi[jj] += s * wpr[jj] + c * wpi[jj];
        }
        
        /* Set the current power */
        pows[ii] = ivar * (sumcxx * sumcxx / sumc + 
                           sumsxx * sumsxx / sums);
    }
    
    /* Free the temp arrays and return */
    free(wr);
    free(wi);
    free(wpr);
    free(wpi);
    free(xxnorm);
    return(pows);
}
