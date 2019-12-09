/*
    likelihood.h
*/
#ifndef LIKELIHOOD_H
#define LIKELIHOOD_H

extern double loglikelihood(int **n_mz, int **n_zw, int *n_m, 
    int nclass, int nlex, int ndocs, int nresp, double **resp,
    double *alpha, double beta, double **eta, double **empirical_z, double nu2, double sigma2);

#endif
