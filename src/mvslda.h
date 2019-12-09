/*
    mvslda.h
    header file for Supervised Latent Dirichlet Allocation with Asymmetric Dirichlet prior.
*/
#ifndef MVSLDA_H
#define MVSLDA_H
#include <stdlib.h>
#include "feature.h"

#define MAXITER_DEFAULT 50
#define NCLASS_DEFAULT 10
#define ALPHA_DEFAULT 0.1
#define BETA_DEFAULT 0.1
#define NU2_DEFAULT 10
#define SIGMA2_DEFAULT 1
#define NRESP_DEFAULT 1
#define RANDOM_SEED_DEFAULT 0

void usage(void);

#endif
