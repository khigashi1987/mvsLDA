/*
    hyper.h
*/
#ifndef HYPER_H
#define HYPER_H

extern void update_alpha(double *alpha, int *n_m, int **n_mz, int ndocs, int nclass);
extern double update_beta(double beta, int *n_z, int **n_zw, int nclass, int nlex);

#endif
