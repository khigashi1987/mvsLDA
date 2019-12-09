/*
    learn.h
*/
#ifndef LEARN_H
#define LEARN_H
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include "feature.h"

//#define RANDOM ((double)rand()/(double)RAND_MAX)

int sampling_multinomial(gsl_rng *r, double *p, double *cum_sum_p, int len_p);

extern void mvslda_learn(document *data,
        double **resp,
        double *alpha,
        double beta,
        double nu2, double sigma2,
        int nclass, int nlex, int dlenmax, int nresp,
        int maxiter,
        double **phi, double **theta, double **eta,
        int **n_mz, int **n_zw,
        FILE *likp, FILE *hyperp,
        unsigned long int random_seed);

#endif
