#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cdflib.h"

double incoherent_cand_sigma(double power, int numsum, double numtrials)
/* Return the approximate significance in Gaussian sigmas */
/* of a candidate of numsum incoherently summed powers,   */
/* taking into account the number of independent trials.  */
{
  double x=0.0;
  
  if (numsum==1 && numtrials==1 && power>2.0){
    /* from Abramowitz and Stegun 26.2.23 (p. 933) */
    /* x error is < 4.5e-4 (absolute)              */
    double t, num, den;

    t = sqrt(2.0*power);
    num = 2.515517 + t * (0.802853 + t * 0.010328);
    den = 1.0 + t * (1.432788 + t * (0.189269 + t * 0.001308));
    x = t - num / den;
  } else {
    int which, status;
    double p, q, bound, mean=0.0, sd=1.0, shape, scale=1.0;

    which = 1;
    status = 0;
    shape = (double) numsum;
    x = power;
    cdfgam(&which, &p, &q, &x, &shape, &scale, &status, &bound);
    if (status){
      printf("\nError in cdfgam() (incoherent_cand_sigma()):\n");
      printf("   status = %d, bound = %g\n", status, bound);
      printf("   p = %g, q = %g, x = %g, shape = %g, scale = %g\n\n", 
	     p, q, x, shape, scale);
      exit(1);
    }
    if (p==1.0)
      q *= numtrials;
    else
      q = 1.0 - pow(p, numtrials);
    p = 1.0 - q;
    which = 2;
    status = 0;
    cdfnor(&which, &p, &q, &x, &mean, &sd, &status, &bound);
    if (status){
      if (status == -2){
	x = 0.0;
      } else if (status == -3){
	x = 38.5;
      } else {
	printf("\nError in cdfnor() (incoherent_cand_sigma()):\n");
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
/* Return the approximate coherently summed power required   */
/* to get a Gaussian significance of 'sigma', after taking   */
/* into account the number of independent trials.            */
{
  /* Note:  Coherently summed candidates are distributed     */
  /*        exponentially as long as the original amplitudes */
  /*        were normalized and power/numsum.                */
  return incoherent_power_for_sigma(sigma, 1, numtrials)*numsum;
}
