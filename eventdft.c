/*
 * This program searches a portion of the Fourier Transform
 * of a list of events.  The candidates above a certain power 
 * level are sorted and placed placed in an '.out' file
 * along with their significance.
 *
 * Copyright 2014, Scott M. Ransom (sransom@nrao.edu)
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "utils.h"
#include "eventdft.h"
#include "eventdft_cmd.h"

/* Functions at the bottom */

static void print_percent_complete(long long current, long long number);
static int compare_eventdftcands(const void *ca, const void *cb);
static double percolate_eventdftcand(eventdftcand* list, int nlist);
static void calc_incoherent_sums(float *isums, fcomplex *amps, int numsum);
static void calc_coherent_sums(float *csums, fcomplex *amps, int numsum);
void check_cands(float *isums, float *ithresholds, int skip_incoherent,
                 float *csums, float *cthresholds, int skip_coherent,
                 int numsum, double freq, eventdftcand *cands, 
                 int numcands, double numtrials);

/****************************************************************/

int main(int argc, char **argv)
// Based on code written by Scott Ransom in March 2003, which was
// based on an earlier code called 'toafft' that didn't do summing
{
    FILE *iofile;
    char *infilenm, *outfilenm, *cptr;
    int jj, nn, numevents;
    long long ii, numfreqs;
    double *events, T, dfreq, freq;
    float *isums, *csums, *ithresholds, *cthresholds;
    fcomplex *amplitudes;
    eventdftcand *cands;
    Cmdline *cmd;
    
    /* Call usage() if we have no command line arguments */
    
    if (argc == 1) {
        Program = argv[0];
        usage();
        exit(1);
    } else {
        fprintf(stderr, "\n     Event DFT Search Routine\n");
        fprintf(stderr, "       With Harmonic Summing\n");
        fprintf(stderr, "        by Scott M. Ransom\n\n");
    }
    
    /* Parse the command line using the excellent program Clig */
    
    cmd = parseCmdline(argc, argv);
    if (cmd->noincoherentP && cmd->nocoherentP){
        /* Just do a normal non-summing search */
        cmd->noincoherentP = 0;
        cmd->numsum = 1;
    }
    if (cmd->numsum==1)
        cmd->nocoherentP = 1;
    
#ifdef DEBUG
    showOptionValues();
#endif
    
    /* Prep and open our input file */
    
    nn = strlen(cmd->argv[0]);
    infilenm = (char *)malloc((unsigned int)(nn));
    sprintf(infilenm, "%s", cmd->argv[0]);	
    if (cmd->doubleP){
        if ((iofile = fopen(infilenm, "rb")) == NULL){
            perror("\nError opening the input file.\n");
            exit(-1);
        }
    } else {
        if ((iofile = fopen(infilenm, "r")) == NULL){
            perror("\nError opening the input file.\n");
            exit(-1);
        }
    }
    outfilenm = (char *)calloc((unsigned int)(nn + 5), 1);
    cptr = strrchr(infilenm, '.');
    if (cptr != NULL){
        strncpy(outfilenm, infilenm, cptr-infilenm);
        strncpy(outfilenm+(cptr-infilenm), ".out", 4);
    } else
        sprintf(outfilenm, "%s.out", cmd->argv[0]);	
    
    /* Allocate our candidate and statistics structures */
    cands = (eventdftcand *)malloc(sizeof(eventdftcand) * cmd->ncands);
    for (ii=0; ii<cmd->ncands; ii++)
        cands[ii].sigma = 0.0;
    isums = gen_fvect(cmd->numsum);
    csums = gen_fvect(cmd->numsum);
    ithresholds = gen_fvect(cmd->numsum);
    cthresholds = gen_fvect(cmd->numsum);
    for (ii=0; ii<cmd->numsum; ii++){
        isums[ii] = 0.0;
        csums[ii] = 0.0;
        ithresholds[ii] = 0.0;
        cthresholds[ii] = 0.0;
    }
    
    /* Read the events */
    {
        int eventtype=0;
        
        if (cmd->daysP)
            eventtype = 1;
        else if (cmd->mjdsP)
            eventtype = 2;
        events = read_events(iofile, cmd->doubleP, eventtype, &numevents,
                             cmd->MJD0, cmd->T, cmd->startT, cmd->endT, 
                             cmd->offset);
        /* Correct the events for the fdot we are searching */
        if (cmd->fdotbyfP){
            for (ii=0; ii<numevents; ii++)
                events[ii] += 0.5 * events[ii] * events[ii] * cmd->fdotbyf;
        }
    }
    T = events[numevents-1];
    if (!cmd->osampP)
        cmd->osamp = 2 * cmd->numsum;
    dfreq = 1.0 / (T * cmd->osamp);
    numfreqs = (cmd->fmax - cmd->fmin) / dfreq + 1;
    
    /* Calculate the approximate number of independent freqs */
    
    if (!cmd->ifsP)
        cmd->ifs = (cmd->fmax - cmd->fmin) * T;
    /* From Horne and Baliunas, ApJ, 302, 757.                  */
    /* This needs to be seriously modified for eventsearches... */
    /* cmd->ifs = -6.362 + 1.193 * numtoas + 0.00098 * numtoas * numtoas */
    
    /* Output some info to the screen */
    
    fprintf(stderr, "Search information:\n");
    fprintf(stderr, "  Total number of events = %-7d (T = %.15g s)\n", 
            numevents, events[numevents-1]);
    fprintf(stderr, "   Frequencies to search = %.15g to %.15g hz\n", 
            cmd->fmin, cmd->fmax);
    fprintf(stderr, " Frequency stepsize (hz) = %.15g\n", dfreq);
    if (cmd->fdotbyfP){
        fprintf(stderr, " Fdot/F applied (sec^-1) = %.15g\n", cmd->fdotbyf);
    }
    fprintf(stderr, "    Max harmonics summed = %d\n", cmd->numsum);
    fprintf(stderr, "     Oversampling factor = %d\n", cmd->osamp);
    fprintf(stderr, "    Candidates to return = %d\n", cmd->ncands);
    fprintf(stderr, "  Approx # of Ind. Freqs = %.0f\n\n", cmd->ifs);
    
    /* Do the search ... */
    
    prep_eventdft(events, numevents, cmd->numsum, cmd->fmin, dfreq);
    for (ii = 0; ii < numfreqs; ii++){
        if (!cmd->nostatusP) print_percent_complete(ii, numfreqs);

        /* Calculate the exact DFT at the fundamental and harmonics */
        amplitudes = calc_eventdft_point(&freq);
        if (!cmd->noincoherentP)  /* incoherent sums */
            calc_incoherent_sums(isums, amplitudes, cmd->numsum);
        if (!cmd->nocoherentP)    /* coherent sums */
            calc_coherent_sums(csums, amplitudes, cmd->numsum);
        /* See if the candidates are any good */
        check_cands(isums, ithresholds, cmd->noincoherentP,
                    csums, cthresholds, cmd->nocoherentP,
                    cmd->numsum, freq, cands, cmd->ncands, cmd->ifs);
        if (cmd->ioutP){
            fprintf(stdout, "%15.10f", freq);
            for (jj = 0; jj < cmd->numsum; jj++)
                fprintf(stdout, "  %10.4f", isums[jj]);
            fprintf(stdout, "\n");
        }
        if (cmd->coutP){
            fprintf(stdout, "%15.10f", freq);
            for (jj = 0; jj < cmd->numsum; jj++)
                fprintf(stdout, "  %10.4f", csums[jj]);
            fprintf(stdout, "\n");
        }
    }
    free_eventdft();
    
    /* Sort our candidate array so that the most significant */
    /* 'detections' are first.                               */
    
    qsort(cands, cmd->ncands, sizeof(eventdftcand), compare_eventdftcands);
    
    /* Write our (text) output file */ 
    
    if ((iofile = fopen(outfilenm, "w")) == NULL){
        perror("\nError opening the output file.\n");
        exit(-1);
    }
    
    fprintf(iofile, "#\n# Search information:\n");
    fprintf(iofile, "#--------------------\n");
    fprintf(iofile, "#  Total number of events = %-6d (T = %.15g s)\n", 
            numevents, events[numevents-1]);
    fprintf(iofile, "#    Frequencies searched = %.15g to %.15g hz\n", 
            cmd->fmin, cmd->fmax);
    fprintf(iofile, "#    Max harmonics summed = %d\n", cmd->numsum);
    fprintf(iofile, "#     Oversampling factor = %d\n", cmd->osamp);
    fprintf(iofile, "# Frequency stepsize (hz) = %.15g\n", dfreq);
    fprintf(iofile, "#  Approx # of Ind. Freqs = %.0f\n", cmd->ifs);
    fprintf(iofile, "#\n# Results:\n");
    fprintf(iofile, "#---------\n");
    fprintf(iofile, "# Cand     Frequency      Sigma    IPower    CPower   Nharm\n");
    ii = 0;
    while (cands[ii].freq > 0.0 && ii < cmd->ncands){
        fprintf(iofile, " %-4lld  %15.10f  %7.2f   %7.2f   %7.2f  %4d %c\n", 
                ii+1, cands[ii].freq, cands[ii].sigma, 
                cands[ii].ipow, cands[ii].cpow, cands[ii].numharm, 
                cands[ii].coherentsum?'C':'I');
        ii++;
    }
    fprintf(iofile,"\n"); 
    fclose(iofile);
    
    /* Clean-up */
    
    free(events);
    free(cands);
    free(infilenm);
    free(outfilenm);
    free(isums);
    free(csums);
    free(ithresholds);
    free(cthresholds);
    fprintf(stderr, "\nDone.\n\n");
    return (0);								
}									


static double percolate_eventdftcand(eventdftcand* list, int nlist)
/*  Pushes a eventdftcand struct as far up a sorted list of eventdftcands */
/*  as it needs to go to keep the list sorted.  Returns the new low   */
/*  power in the list.                                                */
{
    int ct;
    eventdftcand temp;
    
    for (ct = nlist - 2; ct >= 0; ct--) {
        if (list[ct].sigma < list[ct + 1].sigma) {
            temp = list[ct + 1];
            list[ct + 1] = list[ct];
            list[ct] = temp;
        } else {
            break;
        }
    }
    return list[nlist - 1].sigma;
}


static int compare_eventdftcands(const void *ca, const void *cb)
/*  Used as compare function for qsort() */
{
    eventdftcand *a, *b;
    
    a = (eventdftcand *) ca;
    b = (eventdftcand *) cb;
    if ((b->sigma - a->sigma) < 0.0)
        return -1;
    if ((b->sigma - a->sigma) > 0.0)
        return 1;
    return 0;
}


void calc_incoherent_sums(float *isums, fcomplex *amps, int numsum)
/* Calculate the incoherent power sum for *amps and */
/* return it in *isums.                             */
{
    int ii;
    fcomplex *amp;
    
    amp = amps;
    isums[0] = amp->r*amp->r + amp->i*amp->i;
    for (ii=1; ii<numsum; ii++, amp++)
        isums[ii] = isums[ii-1] + amp->r*amp->r + amp->i*amp->i;
}


void calc_coherent_sums(float *csums, fcomplex *amps, int numsum)
/* Calculate the coherent power sum for *amps and   */
/* return it in *csums.  Note:  It modifies *amps!  */
{
    int ii;
    double phs0, phscorr, ar, ai, cc, ss;
    fcomplex *amp, *prevamp;
    
    amp = amps;
    phs0 = atan2(amp->i, amp->r);
    csums[0] = amp->r*amp->r + amp->i*amp->i;
    for (ii=1, prevamp=amps, amp=amps+1; ii<numsum; ii++, amp++){
        phscorr = phs0 - fmod((ii+1.0)*phs0, TWOPI);
        ar = amp->r;
        ai = amp->i;
        cc = cos(phscorr);
        ss = sin(phscorr);
        amp->r = ar*cc - ai*ss + prevamp->r;
        amp->i = ai*cc + ar*ss + prevamp->i;
        csums[ii] = amp->r*amp->r + amp->i*amp->i;
        prevamp = amp;
    }
}


void check_cands(float *isums, float *ithresholds, int skip_incoherent,
                 float *csums, float *cthresholds, int skip_coherent,
                 int numsum, double freq, eventdftcand *cands, 
                 int numcands, double numtrials)
{
    int ii, goodcand=0, bestsum=0;
    double minsigma, maxsigma=0, sigma;
    eventdftcand *lastcand; 
    
    /* Check to see if any incoherent candidates are better */
    /* than the worst candidate that we have so far.        */
    if (!skip_incoherent){
        for (ii=0; ii<numsum; ii++){
            if (isums[ii] > ithresholds[ii]){
                goodcand=1;
                break;
            }
        }
    }
    /* Check to see if any coherent candidates are better   */
    /* than the worst candidate that we have so far.        */
    if (!goodcand && !skip_coherent){
        for (ii=0; ii<numsum; ii++){
            if (csums[ii] > cthresholds[ii]){
                goodcand=1;
                break;
            }
        }
    }
    
    /* Add the good candidate to out candidate list */
    if (goodcand){
        lastcand = cands+(numcands-1);
        /* Find the best candidate */
        if (!skip_incoherent){
            for (ii=0; ii<numsum; ii++){
                sigma = incoherent_cand_sigma(isums[ii], ii+1, numtrials);
                if (sigma > maxsigma){
                    maxsigma = sigma;
                    bestsum = ii+1;
                    lastcand->ipow = isums[ii];
                    lastcand->cpow = csums[ii];
                    lastcand->sigma = sigma;
                    lastcand->numharm = bestsum;
                    lastcand->coherentsum = 0;
                }
            }
        }
        if (!skip_coherent){
            for (ii=0; ii<numsum; ii++){
                sigma = coherent_cand_sigma(csums[ii], ii+1, numtrials);
                if (sigma > maxsigma){
                    maxsigma = sigma;
                    bestsum = ii+1;
                    lastcand->ipow = isums[ii];
                    lastcand->cpow = csums[ii];
                    lastcand->sigma = sigma;
                    lastcand->numharm = bestsum;
                    lastcand->coherentsum = 1;
                }
            }
        }
        lastcand->freq = freq;
        minsigma = percolate_eventdftcand(cands, numcands);
        /* Raise the thresholds appropriately */
        if (!skip_incoherent){
            for (ii=0; ii<numsum; ii++)
                ithresholds[ii] = incoherent_power_for_sigma(minsigma, ii+1, numtrials);
        }
        if (!skip_coherent){
            for (ii=0; ii<numsum; ii++)
                cthresholds[ii] = coherent_power_for_sigma(minsigma, ii+1, numtrials);
        }
    }
}


static void print_percent_complete(long long current, long long number)
{
    static int newper = 0, oldper = -1;
    
    newper = (int) round((double) current / (double) (number) * 100.0);
    if (newper > oldper) {
        fprintf(stderr, "\rAmount complete = %3d%%", newper);
        fflush(stdout);
        oldper = newper;
    }
}
