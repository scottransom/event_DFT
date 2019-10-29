#define TWOPI 6.2831853071795864769252867665590057683943387987502

/* The structure we will use to describe the candidates */
typedef struct EVENTDFTCAND {
  double freq;      /* Frequency in hz                      */
  double ipow;      /* Normalized incoherently summed power */
  double cpow;      /* Normalized coherently summed power   */
  double sigma;     /* Approx sigma of the candidate        */
  int numharm;      /* Number of harmonics summed           */
  int coherentsum;  /* 0 if incoherent, 1 if coherent       */
} eventdftcand;

/* period.c */
void free_eventdft(void);
void prep_eventdft(double *events, long numevts, int maxnumharmsum, 
		   double lof, double df);
fcomplex *calc_eventdft_point(double *freq);
fcomplex *eventdft(double *events, long numevents, double lof, 
		   double df, long numf);

/* stats.c */
double incoherent_cand_sigma(double power, int numsum, double numtrials);
double coherent_cand_sigma(double power, int numsum, double numtrials);
double incoherent_power_for_sigma(double sigma, int numsum, double numtrials);
double coherent_power_for_sigma(double sigma, int numsum, double numtrials);
