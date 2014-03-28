#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cdflib.h"

static double extended_equiv_gaussian_sigma(double logp);
static double log_asymtotic_incomplete_gamma(double a, double z);
static double log_asymtotic_gamma(double z);

double incoherent_cand_sigma(double power, int numsum, double numtrials)
/* Return the approximate significance in Gaussian       */
/* sigmas of a candidate of numsum summed powers,        */
/* taking into account the number of independent trials. */
{
    double x = 0.0;
    
    if (power <= 0.0) {
        return 0.0;
    }
    
    if (power > 100.0) {
        double logp;
        
        /* Use some asymtotic expansions for the chi^2 distribution */
        logp = log_asymtotic_incomplete_gamma(numsum, power) -
            log_asymtotic_gamma(numsum);
        /* Now adjust for the number of trials */
        logp += log(numtrials);
        /* Convert to a sigma */
        x = extended_equiv_gaussian_sigma(logp);
    } else {
        int which, status;
        double p, q, bound, mean = 0.0, sd = 1.0, shape, scale = 1.0;
        
        which = 1;
        status = 0;
        shape = (double) numsum;
        x = power;
        /* Determine the basic probability */
        cdfgam(&which, &p, &q, &x, &shape, &scale, &status, &bound);
        if (status) {
            printf("\nError in cdfgam() (candidate_sigma()):\n");
            printf("   status = %d, bound = %g\n", status, bound);
            printf("   p = %g, q = %g, x = %g, shape = %g, scale = %g\n\n",
                   p, q, x, shape, scale);
            exit(1);
        }
        /* Adjust it for the number of trials */
        if (p == 1.0)
            q *= numtrials;
        else
            q = 1.0 - pow(p, numtrials);
        p = 1.0 - q;
        which = 2;
        status = 0;
        /* Convert to a sigma */
        cdfnor(&which, &p, &q, &x, &mean, &sd, &status, &bound);
        if (status) {
            if (status == -2) {
                x = 0.0;
            } else if (status == -3) {
                x = 38.5;
            } else {
                printf("\nError in cdfnor() (candidate_sigma()):\n");
                printf("   status = %d, bound = %g\n", status, bound);
                printf("   p = %g, q = %g, x = %g, mean = %g, sd = %g\n\n",
                       p, q, x, mean, sd);
                exit(1);
            }
        }
    }
    if (x < 0.0)
        return 0.0;
    else
        return x;
}


double coherent_cand_sigma(double power, int numsum, double numtrials)
/* Return the approximate significance in Gaussian sigmas */
/* of a candidate of numsum coherently summed powers,     */
/* taking into account the number of independent trials.  */
{
    /* Note:  Coherently summed candidates are distributed     */
    /*        exponentially as long as the original amplitudes */
    /*        were normalized and power/numsum.                */
    return incoherent_cand_sigma(power/numsum, 1, numtrials);
}


double incoherent_power_for_sigma(double sigma, int numsum, double numtrials)
/* Return the approximate incoherently summed power required */
/* to get a Gaussian significance of 'sigma', after taking   */
/* into account the number of independent trials.            */
{
    int which, status;
    double p, q, x, bound, mean=0.0, sd=1.0, df, scale=1.0;
    
    which = 1;
    status = 0;
    x = sigma;
    cdfnor(&which, &p, &q, &x, &mean, &sd, &status, &bound);
    if (status){
        printf("\nError in cdfnor() (incoherent_power_for_sigma()):\n");
        printf("   cdfstatus = %d, bound = %g\n\n", status, bound);
        printf("   p = %g, q = %g, x = %g, mean = %g, sd = %g\n\n", 
               p, q, x, mean, sd);
        exit(1);
    }
    q = q / numtrials;
    p = 1.0 - q;
    which = 2;
    df = 2.0 * numsum;
    status = 0;
    cdfchi(&which, &p, &q, &x, &df, &status, &bound);
    if (status){
        printf("\nError in cdfchi() (incoherent_power_for_sigma()):\n");
        printf("   status = %d, bound = %g\n", status, bound);
        printf("   p = %g, q = %g, x = %g, df = %g, scale = %g\n\n", 
               p, q, x, df, scale);
        exit(1);
    }
    return 0.5 * x;
}


double coherent_power_for_sigma(double sigma, int numsum, double numtrials)
// Return the approximate coherently summed power required
// to get a Gaussian significance of 'sigma', after taking
// into account the number of independent trials.
{
    /* Note:  Coherently summed candidates are distributed     */
    /*        exponentially as long as the original amplitudes */
    /*        were normalized and power/numsum.                */
    return incoherent_power_for_sigma(sigma, 1, numtrials)*numsum;
}


double extended_equiv_gaussian_sigma(double logp)
// Return the equivalent gaussian sigma corresponding to the log of
// the cumulative gaussian probability logp.  In other words, return x,
// such that Q(x) = p, where Q(x) is the cumulative normal distribution.
// This version uses the rational approximation from Abramowitz and
// Stegun, eqn 26.2.23.  Using the log(P) as input gives a much extended
// range.
{ double t, num, denom;
    
    t = sqrt(-2.0 * logp);
    num = 2.515517 + t * (0.802853 + t * 0.010328);
    denom = 1.0 + t * (1.432788 + t * (0.189269 + t * 0.001308));
    return t - num / denom;
}


double log_asymtotic_incomplete_gamma(double a, double z)
// Return the log of the incomplete gamma function in its
// asymtotic limit as z->infty.  This is from Abramowitz
// and Stegun eqn 6.5.32.
{
    double x = 1.0, newxpart = 1.0, term = 1.0;
    int ii = 1;
    
    while (fabs(newxpart) > 1e-15) {
        term *= (a - ii);
        newxpart = term / pow(z, ii);
        x += newxpart;
        ii += 1;
    }
    return (a - 1.0) * log(z) - z + log(x);
}


double log_asymtotic_gamma(double z)
// Return the log of the gamma function in its asymtotic limit
// as z->infty.  This is from Abramowitz and Stegun eqn 6.1.41.
{
    double x, y;
    
    x = (z - 0.5) * log(z) - z + 0.91893853320467267;
    y = 1.0 / (z * z);
    x += (((-5.9523809523809529e-4 * y
            + 7.9365079365079365079365e-4) * y
           - 2.7777777777777777777778e-3) * y + 8.3333333333333333333333e-2) / z;
    return x;
}
